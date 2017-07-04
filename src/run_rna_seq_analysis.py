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


def get_tool_id(tool_name):
    '''
    Retrieve the id of a tool
    '''
    tool_id = ''
    for tool in tools:
        if tool["name"] == tool_name:
            tool_id = tool["id"]
    return tool_id


# Connect to Galaxy and retrieve the history
gi = GalaxyInstance(config["galaxy_url"], config["api_key"])
histories = gi.histories.get_histories()
hist = check_hist(config["hist_name"]["raw_data"])

# Extract the sample names
finename_desc_df = pd.read_csv("data/file_description.csv", index_col = 0)
sample_names = list(finename_desc_df.index)

# Get tools in the Galaxy instance
tools = gi.tools.get_tools()
# Get the id for MultiQC tool
multiqc_id = get_tool_id("multiqc")


rule merge_files:
    '''
    Merge files sequenced on 2 different lanes (for Project_S178 and Project_S225)
    '''
    run:
        # Retrieve the name of samples to merge
        to_merge = {}
        for sample in sample_names:
            project_id = finename_desc_df["Project id"][sample]
            if project_id not in ["Project_S178", "Project_S225"]:
                continue
            to_merge.setdefault(sample, [])
        # Parse the dataset in history to extract the ids of dataset to merge
        for dataset in gi.histories.show_matching_datasets(hist):
            name = dataset['name']
            if not name.endswith("fastq"):
                continue
            sample_name = os.path.splitext(name)[0][:-1]
            if sample_name in to_merge:
                to_merge[sample_name].append(dataset["id"])
            else:
                # Rename to remove the last letter
                gi.histories.update_dataset(
                    hist,
                    dataset["id"],
                    name="%s%s" % (config["raw_data_name_prefix"], sample_name))
        # Get concatenate tool
        tool_id = get_tool_id("Concatenate datasets")
        # Merge datasets
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
            # Rename the datasets
            gi.histories.update_dataset(
                hist,
                info['outputs'][0]['id'],
                name="%s%s" % (config["name_prefix"]["raw_data"], dataset))


rule launch_fastqc:
    '''
    Launch FastQC on raw data and MultiQC report on the FastQC outputs
    '''
    run:
        # Get (or create) the history for FastQC outputs
        fastqc_hist = check_hist(config["hist_name"]["fastqc"])
        if fastqc_hist == '':
            fastqc_hist = gi.histories.create_history(
                config["hist_name"]["fastqc"])["id"]
        # Get the id for FastQC tool
        tool_id = get_tool_id("FastQC")
        # Launch FastQC tool on each raw datasets
        for dataset in gi.histories.show_matching_datasets(hist):
            name = dataset['name']
            if not name.startswith(config["name_prefix"]["raw_data"]):
                continue
            if dataset['state'] != "ok" or dataset['deleted']:
                continue
            # Extract sample name
            sample_name = name.split(config["name_prefix"]["raw_data"])[-1]
            # Create the input datamap
            datamap = dict()
            datamap["input_file"] = {'src':'hda', 'id': dataset["id"]}
            # Run the tool
            info = gi.tools.run_tool(fastqc_hist, tool_id, datamap)
            # Rename the datasets
            for output in info['outputs']:
                output_type = output["name"].split(": ")[-1]
                gi.histories.update_dataset(
                    fastqc_hist,
                    output['id'],
                    name="%s%s: %s" % (
                        config["name_prefix"]["fastqc"],
                        sample_name,
                        output_type))
        # Create the input datamap for MultiQC
        datamap = {
            "results_0|software": "fastqc",
            "results_0|input_file": []}
        # Conserve the RawData files ids
        fastqc_raw_data_ids = []
        for dataset in gi.histories.show_matching_datasets(fastqc_hist):
            name = dataset['name']
            if dataset['state'] != "ok" or dataset['deleted']:
                continue
            if name.find("RawData") != -1:
                datamap["results_0|input_file"].append({
                    'src':'hda',
                    'id': dataset["id"]})
        # Run MultiQC tool
        info = gi.tools.run_tool(fastqc_hist, multiqc_id, datamap)

