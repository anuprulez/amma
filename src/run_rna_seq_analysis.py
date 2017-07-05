from bioblend.galaxy import GalaxyInstance
import pandas as pd


configfile: "config.yaml"


def check_hist(hist_name):
    '''
    Check if an history exist and return its id if it's exist
    '''
    hist = ''
    for history in histories:
        if history["name"] == hist_name:
            hist = history['id']
    return hist


def get_tool_id(tool_name, version):
    '''
    Retrieve the id of a tool
    '''
    tool_id = ''
    for tool in tools:
        if tool["name"] == tool_name and tool["version"] == version:
            tool_id = tool["id"]
    return tool_id


def get_collection_id(collection_name, hist_id):
    '''
    Retrieve the id of a collection
    '''
    coll_id = ''
    for ds in gi.histories.show_history(hist_id, contents=True, visible = True):
        if ds["history_content_type"] != 'dataset_collection':
            continue
        if ds["name"] == collection_name:
            coll_id = ds["id"]
    return coll_id


def fill_multiqc_inputs(coll_id,hist):
    '''
    Extract the list of element in the collection and format it as input for tool
    '''
    inputs = []
    for el in gi.histories.show_dataset_collection(hist, coll_id)["elements"]:
        inputs.append({'src':'hda','id': el["object"]["id"]})
    return inputs


# Connect to Galaxy and retrieve the history
gi = GalaxyInstance(config["galaxy_url"], config["api_key"])
histories = gi.histories.get_histories()
hist = check_hist(config["hist_name"])
if hist == '':
    hist = gi.histories.create_history(config["hist_name"])["id"]

# Extract the sample names
finename_desc_df = pd.read_csv("data/file_description.csv", index_col = 0)
sample_names = list(finename_desc_df.index)

# Get tools in the Galaxy instance
tools = gi.tools.get_tools()
# Get the id for MultiQC tool
multiqc_id = get_tool_id("multiqc")


rule prepare_files:
    '''
    Import the files from the data library, merge the files sequenced on 2 
    different lanes (for Project_S178 and Project_S225) and move the input files
    into collections
    '''
    run:
        # Find the data library
        lib = gi.libraries.get_libraries(name=config["library_name"])
        assert len(lib) > 0, "No library found for Prinz lab"
        lib_id = lib[0]["id"]
        # Parse the data library datasets
        for ds in gi.libraries.show_library(lib_id, contents=True):
            # Eliminate the folder
            if ds['type'] != 'file':
                continue
            # Eliminate the files from other folders
            if ds["name"].find(config["folder_name"]) == -1:
                continue
            # Add the files to the history
            gi.histories.upload_dataset_from_library(
                hist,
                ds["id"])
        # Retrieve the name of samples to merge
        to_merge = {}
        for sample in sample_names:
            project_id = finename_desc_df["Project id"][sample]
            if project_id not in ["Project_S178", "Project_S225"]:
                continue
            to_merge.setdefault(sample, [])
        # Parse the dataset in history to extract the ids of dataset to merge,
        # rename the other files and add them to a collection
        raw_dataset_ids = []
        for dataset in gi.histories.show_matching_datasets(hist):
            name = dataset['name']
            if not name.endswith("fastq"):
                continue
            sample_name = os.path.splitext(name)[0][:-1]
            if sample_name in to_merge:
                to_merge[sample_name].append(dataset["id"])
                # Hide the file
                gi.histories.update_dataset(
                    hist,
                    dataset["id"],
                    visible = False)
            else:
                # Rename the file and hide it
                gi.histories.update_dataset(
                    hist,
                    dataset["id"],
                    name="%s%s" % (config["name_prefix"]["raw_data"], sample_name),
                    visible = False)
                # Add the file to the collection
                raw_dataset_ids.append({
                    'id': dataset["id"],
                    'name': sample_name,
                    'src': 'hda'
                    })
        # Get concatenate tool
        tool_id = get_tool_id("Concatenate datasets")
        # Merge datasets
        unmerged_dataset_ids = []
        for dataset in to_merge:
            if len(to_merge[dataset]) != 2:
                print("Issue with %s" %(dataset))
                continue
            # Create the input datamap
            datamap = dict()
            datamap["inputs"] = [
                {'src':'hda', 'id': to_merge[dataset][0]},
                {'src':'hda', 'id': to_merge[dataset][1]}]
            # Run the tool
            info = gi.tools.run_tool(hist, tool_id, datamap)
            # Rename the dataset and hide it
            gi.histories.update_dataset(
                hist,
                info['outputs'][0]['id'],
                name="%s%s" % (config["name_prefix"]["raw_data"], dataset),
                visible = False)
            # Add the merge file to the collection of raw dataset
            raw_dataset_ids.append({
                'id': info['outputs'][0]['id'],
                'name': dataset,
                'src': 'hda'
                })
            # Add the input files (before merging) to the collection of unmerged
            # datasets
            unmerged_dataset_ids.append({
                'id': to_merge[dataset][0],
                'name': "%s 0" % dataset,
                'src': 'hda'
                })
            unmerged_dataset_ids.append({
                'id': to_merge[dataset][1],
                'name': "%s 1" % dataset,
                'src': 'hda'
                })
        # Prepare and create the collections
        raw_data_collection = {
            'collection_type': 'list',
            'element_identifiers': raw_dataset_ids,
            'name': config["collection_names"]["raw_data"]
        }
        gi.histories.create_dataset_collection(
            hist,
            raw_data_collection)
        unmerged_collection = {
            'collection_type': 'list',
            'element_identifiers': unmerged_dataset_ids,
            'name': "Unmerged files"
        }
        gi.histories.create_dataset_collection(
            hist,
            unmerged_collection)


rule launch_fastqc:
    '''
    Launch FastQC on raw data and MultiQC report on the FastQC outputs
    '''
    run:
        # Get the id for FastQC tool
        tool_id = get_tool_id("FastQC")
        assert tool_id != '', "No FastQC tool"
        # Search for the collection id with the raw data
        raw_data_coll_id = get_collection_id(
            config["collection_names"]["raw_data"],
            hist)
        assert raw_data_coll_id != '', "No collection for Raw data"
        # Create the input datamap for FastQC
        datamap = {"input_file": {
            'batch': True,
            'values': [{'src':'hdca', 'id': raw_data_coll_id}]}}
        # Run FastQC
        try:
            info = gi.tools.run_tool(hist, tool_id, datamap)
        except:
            print("Issue with FastQC launch")
        # Retrieve the "RawData" collection
        fastqc_data_coll_id = ''
        for ds in gi.histories.show_history(hist, contents=True, visible = True, deleted = False):
            if ds["history_content_type"] != 'dataset_collection':
                continue
            if ds["name"].find("FastQC") == -1:
                continue
            if ds["name"].find("RawData") == -1:
                # Rename the web report collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name="FastQC on Raw data: web report")
            else:
                fastqc_data_coll_id = ds["id"]
                # Rename the raw data report collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name="FastQC on Raw data: raw report")
        assert fastqc_data_coll_id != '', "No collection for FastQC Raw Data"
        # Create and fill the input datamap for MultiQC
        datamap = {
            "results_0|software": "fastqc",
            "results_0|input_file": fill_multiqc_inputs(
                fastqc_data_coll_id,
                hist),
            "results_0|saveLog": "False"}
        # Run MultiQC
        info = gi.tools.run_tool(hist, multiqc_id, datamap)
        gi.histories.update_dataset(
            hist,
            info['outputs'][0]['id'],
            name="MultiQC report of FastQC outputs")


rule launch_trim_galore:
    '''
    Launch Trim Galore! on raw data and MultiQC report on the Trim Galore! 
    outputs
    '''
    run:
        # Get the id for Trim Galore! tool
        tool_id = get_tool_id("Trim Galore!")
        assert tool_id != '', "No Trim Galore! tool"
        # Search for the collection id with the raw data
        raw_data_coll_id = get_collection_id(
            config["collection_names"]["raw_data"],
            hist)
        assert raw_data_coll_id != '', "No collection for Raw data"
        # Create the input datamap for Trim Galore!
        datamap = {
            "singlePaired|sPaired": "single",
            "singlePaired|input_singles" : {'batch': True,'values': [
                {'src':'hdca',
                'id': raw_data_coll_id}]},
            "params|settingsType": "custom",
            "params|settingsType": "custom",
            "params|quality": "20",
            "params|error_rate": "0.1",
            "params|min_length": "20",
            "params|report": "true",
        }
        # Run Trim Galore!
        info = gi.tools.run_tool(hist, tool_id, datamap)
        # Retrieve the "RawData" collection
        trim_galore_data_coll_id = ''
        for ds in gi.histories.show_history(hist, contents=True, visible = True):
            if ds["history_content_type"] != 'dataset_collection':
                continue
            if ds["name"].find("Trim Galore!") == -1:
                continue
            if ds["name"].find("trimmed reads") != -1:
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name="Trim Galore! on Raw data: trimmed reads")
            else:
                trim_galore_data_coll_id = ds["id"]
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name="Trim Galore! on Raw data: report")
        assert trim_galore_data_coll_id != '', "No collection for Trim Galore report"
        # Create and fill the input datamap for MultiQC
        datamap = {
            "results_0|software": "fastqc",
            "results_0|input_file": fill_multiqc_inputs(
                trim_galore_data_coll_id,
                hist),
            "results_0|saveLog": "False"}
        # Run MultiQC
        info = gi.tools.run_tool(hist, multiqc_id, datamap)
        gi.histories.update_dataset(
            hist,
            info['outputs'][0]['id'],
            name="MultiQC report of Trim Galore!")
