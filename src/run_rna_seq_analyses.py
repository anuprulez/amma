from bioblend.galaxy import GalaxyInstance
import pandas as pd
import time


def check_hist(hist_name):
    '''
    Check if an history exist and return its id if it's exist
    '''
    hist = ''
    for history in histories:
        if history["name"] == hist_name:
            hist = history['id']
    return hist


def get_tool_id(tool_ids):
    '''
    Retrieve the id of a tool
    '''
    tool_id = ''
    for tool in tools:
        if tool["id"] == tool_ids:
            tool_id = tool["id"]
    return tool_id


def get_working_tool_id(tool_name):
    '''
    Retrieve the id of a tool and test if it exists
    '''
    tool_id = get_tool_id(tool_name)
    assert tool_id != '', "No %s tool" % tool_name
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


def get_working_collection_id(coll_name):
    '''
    Retrieve the id of a collection and test if it exists
    '''
    coll_id = get_collection_id(coll_name, hist)
    assert coll_id != '', "No collection for %s" % coll_name
    return coll_id


def fill_multiqc_inputs(coll_id,hist):
    '''
    Extract the list of element in the collection and format it as input for tool
    '''
    inputs = []
    for el in gi.histories.show_dataset_collection(hist, coll_id)["elements"]:
        inputs.append({'src':'hda','id': el["object"]["id"]})
    return inputs


def run_multiqc(coll_id, software, software_name):
    '''
    Run multiqc tool on the element of a collection given a software
    '''
    # Create and fill the input datamap for MultiQC
    datamap = {
        "results_0|software": software,
        "results_0|input_file": fill_multiqc_inputs(
            coll_id,
            hist),
        "results_0|saveLog": "false"}
    # Run MultiQC
    info = gi.tools.run_tool(hist, multiqc_id, datamap)
    # Rename and/or hide the files
    for output in info['outputs']:
        if output["name"].endswith("Log"):
            gi.histories.update_dataset(
                hist,
                output['id'],
                name="MultiQC log of %s" % software_name,
                visible = False)
        elif output["name"].endswith("Webpage"):
            gi.histories.update_dataset(
                hist,
                output['id'],
                name="MultiQC report of %s" % software_name)


def get_annotation_id():
    '''
    Extract the dataset id in the history
    '''
    annotation_id = ''
    for ds in gi.histories.show_history(hist, contents=True, visible=True):
        if ds["name"].find(config["annotation_name"]) != -1:
            annotation_id = ds["id"]
    assert annotation_id != '', "No annotation file for %s in history" % config["annotation_name"]
    return annotation_id


def extract_dataset_metadate(dataset_name):
    '''
    Extract the type of mice, the age and the gender from the name of the dataset
    '''
    split_dataset_name = dataset_name.split("_")
    mice_type = split_dataset_name[0]
    age = split_dataset_name[1]
    gender = split_dataset_name[2]
    return mice_type, age, gender


def prepare_deseq_datamap():
    '''
    Prepare the DESeq datamap
    '''
    datamap = {
            "tximport|tximport_selector": "count",
            "pdf": "true",
            "normCounts": "true",
            "many_contrasts": "true",
            "fit_type": "1",
            "outlier_replace_off": "false",
            "outlier_filter_off": "false",
            "auto_mean_filter_off": "false"}
    return datamap


def check_col_state(col_id, analysis):
    '''
    Check the state of a collection and wait until ok
    '''
    populated_state = gi.histories.show_dataset_collection(hist, col_id)['populated_state']
    print("%s(%s): %s" % (analysis, col_id, populated_state))
    while populated_state != 'ok' and populated_state != 'error':
        print("%s(%s): %s" % (analysis, col_id, populated_state))
        time.sleep(5)
        populated_state = gi.histories.show_dataset_collection(hist, col_id)['populated_state']
    if len(gi.histories.show_dataset_collection(hist,col_id)['elements']) < 1:
        populated_state = 'error'
    print("%s(%s): %s" % (analysis, col_id, populated_state))
    return populated_state


def get_output_collection_ids(info):
    '''
    Retrieve the ids of the generated collections
    '''
    ids = set()
    print(info.keys())
    for out in info['implicit_collections']:
        print(out)
        ids.add(out['id'])
    assert len(ids) != 0
    ids = list(ids)
    if len(ids) == 1:
        return ids[0]
    else:
        return ids


def rename_generated_collection(info, new_name, visible):
    '''
    Rename a generated collection
    '''
    col_id = get_output_collection_ids(info)
    assert type(col_id) == str
    gi.histories.update_dataset_collection(hist, col_id, name=new_name, visible=visible)
    # Rename the generated files too
    for out in info['outputs']:
        out_id = out['id']
        prov = gi.histories.show_dataset_provenance(hist, out_id)
        print(prov)
        if not 'input|__identifier__' in prov['parameters']:
            print("Issue to find a correct name")
            print(prov)
            continue
        new_name = prov['parameters']['input|__identifier__'].replace('"','')
        gi.histories.update_dataset(hist, out_id, name=new_name)
    return col_id


def launch_deseq_analyses(datamaps, analysis_type):
    '''
    Launch DESeq, change the generated names to have the analysis id in it, and
    extract the differentially expressed genes, the up/down-expressed genes and
    prepare file for Venn diagram

    Can not have one loop: need to force to wait the job is done
    '''
    # Launch DESeq
    raw_output_coll_id = {}
    for an in datamaps:
        info = gi.tools.run_tool(hist, deseq_id, datamaps[an])
        # Rename the generated outputs
        for out in info['outputs']:
            ds_id = out['id']
            ds_name = out['name']
            new_name = "%s: %s" % (an, ds_name)
            gi.histories.update_dataset(hist, ds_id, name=new_name)
        # Rename the generated output collection
        for out in info['output_collections']:
            col_name = out['name']
            raw_output_coll_id[an] = out['id']
            new_name = "%s: %s" % (an, col_name)
            gi.histories.update_dataset_collection(hist, out['id'], name=new_name)
        assert an in raw_output_coll_id
    # Prepare the header line for the Venn diagram formatted file
    venn_header_ds = ''
    for ds in gi.histories.show_history(hist, contents=True, visible=True, deleted=False):
        if ds["name"].find("Venn diagram file header") != -1:
            venn_header_ds = ds['id']
    if venn_header_ds == '':
        # Import the header
        info = gi.tools.paste_content("key,Feature,Description,Gene Name,logFC,adj.P.Val", hist)
        pasted_entries_id = ''
        for out in info["outputs"]:
            ds_id = out['id']
            pasted_entries_id = ds_id
            new_name = "(Non formatted) Venn diagram file header"
            gi.histories.update_dataset(hist, ds_id, name=new_name, visible=False)
        assert pasted_entries_id != ''
        # Change the ',' by '\t' to get a table
        datamap = {
            "input": {'src':'hda','id': pasted_entries_id},
            "strip": "True",
            "condense": "True",
            "convert_from": "C"}
        info = gi.tools.run_tool(hist, convert_id, datamap)
        venn_header_ds = ''
        for out in info['outputs']:
            ds_id = out['id']
            venn_header_ds = ds_id
            new_name = "Venn diagram file header"
            gi.histories.update_dataset(hist, ds_id, name=new_name)
        assert venn_header_ds != ''
    # Filter to conserve the differentially expressed genes
    diff_expr_coll_id = {}
    for an in raw_output_coll_id:
        print("filter on %s" % an)
        col_id = raw_output_coll_id[an]
        print("col_id: %s" % (col_id))
        # Check population state before continuing
        populated_state = check_col_state(col_id, "filter %s" % an) 
        if populated_state == 'error':
            continue
        # Launch filtering
        datamap = {
            "input": {'batch': True,'values': [
                {'src':'hdca',
                'id': col_id}]},
            "cond": "c7<0.05",
            "header_lines": "0"}
        print(datamap)
        info = gi.tools.run_tool(hist, filter_id, datamap)
        new_name = "%s: Differentially expressed genes" % (an)
        diff_expr_coll_id[an] = rename_generated_collection(info, new_name, False)
        #diff_expr_coll_id[an] = new_name
        #rename_generated_collection(info, new_name, True)
        print("%s %s" % (col_id, diff_expr_coll_id[an]))
    # Cut to conserve only the column with the fold change
    cut_diff_expr_coll_id = {}
    for an in diff_expr_coll_id:
        print("cut on %s" % an)
        col_id = diff_expr_coll_id[an]
        #print(diff_expr_coll_id[an])
        #col_id = get_working_collection_id(diff_expr_coll_id[an])
        print("col_id: %s" % (col_id))
        # Check population state before continuing
        populated_state = check_col_state(col_id, "cut %s" % an)
        if populated_state == 'error':
            continue
        datamap = {
            "input": {'batch': True,'values': [
                {'src':'hdca',
                'id': col_id}]},
            "complement": "",
            "delimiter": "",
            "cut_type_options|cut_element": "-f",
            "cut_type_options|list": "c1,c3,c7"}
        print(datamap)
        info = gi.tools.run_tool(hist, cut_id, datamap)
        new_name = "%s: Differentially expressed genes" % (an)
        cut_diff_expr_coll_id[an] = rename_generated_collection(info, new_name, True)
        print("%s %s" % (col_id, cut_diff_expr_coll_id[an]))
    # Filter to conserve the upregulated genes in first level
    for an in cut_diff_expr_coll_id:
        col_id = cut_diff_expr_coll_id[an]
        # Check population state before continuing
        populated_state = check_col_state(col_id, "filter %s" % an) 
        if populated_state == 'error':
            continue
        # Launch filtering
        datamap = {
            "input": {'batch': True,'values': [
                {'src':'hdca',
                'id': col_id}]},
            "cond": "c2>0",
            "header_lines": "0"}
        info = gi.tools.run_tool(hist, filter_id, datamap)
        new_name = "%s: Up-expressed genes" % (an)
        rename_generated_collection(info, new_name, True)
    # Filter to conserve the upregulated genes in first level
    for an in cut_diff_expr_coll_id:
        datamap = {
            "input": {'batch': True,'values': [
                {'src':'hdca',
                'id': cut_diff_expr_coll_id[an]}]},
            "cond": "c2<0",
            "header_lines": "0"}
        info = gi.tools.run_tool(hist, filter_id, datamap)
        new_name = "%s: Down-expressed genes" % (an)
        rename_generated_collection(info, new_name, True)
    # Add a column with the name of the file (to prepare Venn diagram)
    venn_diagram = {}
    for an in cut_diff_expr_coll_id:
        col_id = cut_diff_expr_coll_id[an]
        # Check population state before continuing
        populated_state = check_col_state(col_id, "awk %s" % an) 
        if populated_state == 'error':
            continue
        # Launch script to add input name as a new column
        datamap = {
            "input": {'batch': True,'values': [
                {'src':'hdca',
                'id': col_id}]},
            "header|contains_header": "no"}
        info = gi.tools.run_tool(hist, add_input_name_as_column_id, datamap)
        new_name = "%s: Differentially expressed genes (Venn diagram preparation - 1)" % (an)
        venn_diagram[an] = rename_generated_collection(info, new_name, False)
    # Add a column with the name of the file (to prepare Venn diagram)
    for an in venn_diagram:
        col_id = venn_diagram[an]
        # Check population state before continuing
        populated_state = check_col_state(col_id, "awk %s" % an) 
        if populated_state == 'error':
            continue
        # Launch awk script
        datamap = {
            "infile": {'batch': True,'values': [
                {'src':'hdca',
                'id': col_id}]},
            "code": "{print $4,$1,$1,$1,$2,$3}"}
        info = gi.tools.run_tool(hist, awk_id, datamap)
        new_name = "%s: Differentially expressed genes (Venn diagram preparation - 2)" % (an)
        venn_diagram[an] = rename_generated_collection(info, new_name, False)
    # Replace columns input name by the analysis name (for inter analyses Venn diagrams)
    prepared_collections = {}
    for an in venn_diagram:
        col_id = venn_diagram[an]
        # Check population state before continuing
        populated_state = check_col_state(col_id, "replace %s" % an) 
        if populated_state == 'error':
            continue
        # Launch column name replacing
        datamap = {
            "tabular": {'batch': True,'values': [
                {'src':'hdca',
                'id': col_id}]},
            "column": "c1",
            "find_pattern": "(.+)",
            "replace_pattern": "%s_&" % (an)}
        info = gi.tools.run_tool(hist, replace_id, datamap)
        new_name = "%s: Differentially expressed genes (Venn diagram preparation - 3)" % (an)
        prepared_collections[an] = rename_generated_collection(info, new_name, False)
    # Concatenate header and the files (per DESeq analyses) and 
    # prepare for inter DESeq analysis concatenation
    global_an_id = '%s_global' % analysis_type
    deg_datamaps = {global_an_id:{}}
    for an in prepared_collections:
        col_id = prepared_collections[an]
        # Check population state before continuing
        populated_state = check_col_state(col_id, "concatenate %s" % an) 
        if populated_state == 'error':
            continue
        # Retrieve the datasets to concatenate (in a collection)
        datamap = {"inputs": {'src':'hda','id': venn_header_ds}}
        for idx,el in enumerate(gi.histories.show_dataset_collection(hist,col_id)['elements']):
            ds_id = el['object']['id']
            input_id = 'queries_%s|inputs2' % idx
            datamap[input_id] = {'src':'hda','id': ds_id}
            # Prepare concatenation of differentially expressed genes per age comparison
            an_id = el['element_identifier']
            # Per age
            deg_datamaps.setdefault(an_id, {})
            added_file_nb = len(deg_datamaps[an_id].keys())
            input_id = 'queries_%s|inputs2' % added_file_nb
            deg_datamaps[an_id].setdefault(input_id, {'src':'hda','id': el['object']['id']})
            # Global
            added_file_nb = len(deg_datamaps[global_an_id].keys())
            input_id = 'queries_%s|inputs2' % added_file_nb
            deg_datamaps[global_an_id].setdefault(input_id, {'src':'hda','id': el['object']['id']})
        # Concatenate the header and the files
        info = gi.tools.run_tool(hist, cat_id, datamap)
        new_name = "%s: Differentially expressed genes (formatted for Venn diagram)" % (an)
        for out in info['outputs']:
            ds_id = out['id']
            gi.histories.update_dataset(hist, ds_id, name=new_name)
    # Concatenate inter DESeq analyses
    for an in deg_datamaps:
        deg_datamaps[an]["inputs"] = {'src':'hda','id': venn_header_ds}
        info = gi.tools.run_tool(hist, cat_id, deg_datamaps[an])
        new_name = "%s: Differentially expressed genes (formatted for Venn diagram)" % (an)
        for out in info['outputs']:
            ds_id = out['id']
            gi.histories.update_dataset(hist, ds_id, name=new_name)


def delete_datasets(in_name):
    '''
    Delete all datasets and collections for which in_name is found in the name
    '''
    datasets = gi.histories.show_history(hist, contents=True, visible=False)
    datasets += gi.histories.show_history(hist, contents=True, visible=None)
    for ds in datasets:
        if ds["history_content_type"] == 'dataset':
            if ds['purged']:
                continue
            if ds['name'].find(in_name) != -1:
                print('dataset(%s): %s' % (ds['id'],ds['name']))
                gi.histories.delete_dataset(hist, ds['id'], purge=True)
        if ds["history_content_type"] == 'dataset_collection':
            if ds['deleted']:
                continue
            if ds['name'].find(in_name) != -1:
                print('collection(%s): %s' % (ds['id'],ds['name']))
                gi.histories.delete_dataset_collection(hist, ds['id'])


configfile: "config.yaml"

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
# Get the id for tools
multiqc_id = get_working_tool_id(config["tool_ids"]["multiqc"])
deseq_id = get_working_tool_id(config["tool_ids"]["deseq"])
filter_id = get_working_tool_id(config["tool_ids"]["filter"])
cut_id = get_working_tool_id(config["tool_ids"]["cut"])
add_input_name_as_column_id = get_working_tool_id(config["tool_ids"]["add_input_name_as_column"])
convert_id = get_working_tool_id(config["tool_ids"]["convert"])
replace_id = get_working_tool_id(config["tool_ids"]["replace"])
cat_id = get_working_tool_id(config["tool_ids"]["cat"])
add_column_id = get_working_tool_id(config["tool_ids"]["add_column"])
awk_id = get_working_tool_id(config["tool_ids"]["awk"])


rule prepare_files:
    '''
    Import the files from the data library, merge the files sequenced on 2 
    different lanes (for Project_S178 and Project_S225) and move the input files
    into collections
    '''
    run:
        # Find the data library
        lib = gi.libraries.get_libraries(name=config["library_names"]["input_data"])
        assert len(lib) > 0, "No library found for Prinz lab"
        lib_id = lib[0]["id"]
        # Parse the data library datasets
        for ds in gi.libraries.show_library(lib_id, contents=True):
            # Eliminate the folder
            if ds['type'] != 'file':
                continue
            # Eliminate the files from other folders
            if ds["name"].find(config["folder_names"]["input_data"]) == -1:
                continue
            # Add the files to the history
            gi.histories.upload_dataset_from_library(
                hist,
                ds["id"])
        # Retrieve the name of samples to merge
        to_merge = {}
        for sample in sample_names:
            project_id = finename_desc_df["Project id"][sample]
            if project_id not in ["Project_S178", "Project_S198", "Project_S225"]:
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
        tool_id = get_working_tool_id(config["tool_names"]["merging"])
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
        # Find the data library with the annotations
        lib = gi.libraries.get_libraries(name=config["library_names"]["genome_annotations"])
        assert len(lib) > 0, "No library found for %s lab" % config["library_names"]["genome_annotations"]
        lib_id = lib[0]["id"]
        # Add the annotation file to the history
        for ds in gi.libraries.show_library(lib_id, contents=True):
            # Eliminate the folder
            if ds['type'] != 'file':
                continue
            # Eliminate the files from other folders
            if ds["name"].find(config["folder_names"]["annotation"]) == -1:
                continue
            # Find only the "Mus_musculus.GRCm38.87.gtf (mm10)"
            if ds["name"].find(config["annotation_name"]) == -1:
                continue
            # Add the files to the history
            gi.histories.upload_dataset_from_library(
                hist,
                ds["id"])
            annotation_id = ds["id"]
        assert annotation_id != '', "No annotation file for %s" % config["annotation_name"]


rule launch_fastqc:
    '''
    Launch FastQC on raw data and MultiQC report on the FastQC outputs
    '''
    run:
        # Get the id for FastQC tool
        tool_id = get_working_tool_id(config["tool_names"]["fastqc"])
        # Search for the collection id with the raw data
        raw_data_coll_id = get_working_collection_id(
            config["collection_names"]["raw_data"])
        # Create the input datamap for FastQC
        datamap = {"input_file": {
            'batch': True,
            'values': [{'src':'hdca', 'id': raw_data_coll_id}]}}
        # Run FastQC
        try:
            info = gi.tools.run_tool(hist, tool_id, datamap)
        except Exception as e:
            print("Issue with FastQC launch:\n%s " % str(e))
        # Retrieve the "RawData" collection and rename it
        fastqc_data_coll_id = ''
        for ds in gi.histories.show_history(hist, contents=True, visible=True, deleted=False):
            if ds["history_content_type"] != 'dataset_collection':
                continue
            if ds["name"].find(config["tool_names"]["fastqc"]) == -1:
                continue
            if ds["name"].find("RawData") == -1:
                # Rename the web report collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["fastqc"]["web_report"])
            else:
                fastqc_data_coll_id = ds["id"]
                # Rename the raw data report collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["fastqc"]["raw_report"])
        assert fastqc_data_coll_id != '', "No collection for %s" % config["collection_names"]["fastqc"]["raw_report"]
        # Launch MultiQC
        run_multiqc(
            fastqc_data_coll_id,
            "fastqc",
            config["tool_names"]["fastqc"])


rule launch_trim_galore:
    '''
    Launch Trim Galore! on raw data and MultiQC report on the Trim Galore! 
    outputs
    '''
    run:
        # Get the id for Trim Galore! tool
        tool_id = get_working_tool_id(config["tool_names"]["trim_galore"])
        # Search for the collection id with the raw data
        raw_data_coll_id = get_working_collection_id(
            config["collection_names"]["raw_data"])
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
        # Retrieve the trim Galore collections and rename them
        trim_galore_data_coll_id = ''
        for ds in gi.histories.show_history(hist, contents=True, visible=True):
            if ds["history_content_type"] != 'dataset_collection':
                continue
            if ds["name"].find(config["tool_names"]["trim_galore"]) == -1:
                continue
            if ds["name"].find("trimmed reads") != -1:
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["trim_galore"]["trimmed"])
            else:
                trim_galore_data_coll_id = ds["id"]
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["trim_galore"]["report"])
        assert trim_galore_data_coll_id != '', "No collection for %s" % config["collection_names"]["trim_galore"]["report"]
        # Launch MultiQC
        run_multiqc(
            trim_galore_data_coll_id,
            "cutadapt",
            config["tool_names"]["trim_galore"])


rule launch_preliminary_mapping:
    '''
    Run a preliminary mapping to infer the experiment (library type): extract
    200,000 sequences (the first 800,000 lines), run STAR and then Infer 
    Experiment
    '''
    run:
        # Get the id for the tool to downsample the datasets
        tool_id = get_working_tool_id(config["tool_names"]["seq_extraction"])
        # Search for the collection id with the trimmed data
        input_data_coll_id = get_working_collection_id(
            config["collection_names"]["trim_galore"]["trimmed"])
        # Create the input datamap
        datamap = {
            "infile" : {'batch': True,'values': [
                {'src':'hdca',
                'id': input_data_coll_id}]},
            "complement": "",
            "count": "800000",
        }
        # Run "Select first"
        info = gi.tools.run_tool(hist, tool_id, datamap)
        # Retrieve the generated collections and rename it
        downsample_coll_id = ''
        for ds in gi.histories.show_history(hist, contents=True, visible=True):
            if ds["history_content_type"] != 'dataset_collection':
                continue
            if ds["name"].find(config["tool_names"]["seq_extraction"]) == -1:
                continue
            downsample_coll_id = ds["id"]
            # Rename the collection
            gi.histories.update_dataset_collection(
                hist,
                ds["id"],
                name=config["collection_names"]["preliminary_mapping"]["seq_extraction"])
        assert downsample_coll_id != '', "No collection for %s" % config["collection_names"]["preliminary_mapping"]["seq_extraction"]
        # Get the id for STAR
        tool_id = get_working_tool_id(config["tool_names"]["star"])
        # Extract the annotation dataset id
        annotation_id = get_annotation_id()
        # Create the input datamap for STAR
        datamap = {
            "singlePaired|sPaired": "single",
            "singlePaired|input1" : {'batch': True,'values': [
                {'src':'hdca',
                'id': downsample_coll_id}]},
            "refGenomeSource|geneSource": "indexed",
            "refGenomeSource|GTFconditional|GTFselect": "without-gtf",
            "refGenomeSource|GTFconditional|genomeDir": "mm10",
            "refGenomeSource|GTFconditional|sjdbGTFfile": {
                'src':'hda',
                'id': annotation_id},
            "refGenomeSource|GTFconditional|sjdbOverhang": "100"
        }
        # Run STAR
        info = gi.tools.run_tool(hist, tool_id, datamap)
        # Retrieve the STAR collections and rename (and hide) them
        star_data_coll_id = ''
        for ds in gi.histories.show_history(hist, contents=True, visible=True):
            if ds["history_content_type"] != 'dataset_collection':
                continue
            if ds["name"].find(config["tool_names"]["star"]) == -1:
                continue
            if ds["name"].endswith("log"):
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["preliminary_mapping"]["star"]["log"],
                    visible=False)
            elif ds["name"].endswith("mapped.bam"):
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["preliminary_mapping"]["star"]["mapped"])
                star_data_coll_id = ds["id"]
            elif ds["name"].endswith("splice junctions.bed"):
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["preliminary_mapping"]["star"]["splice_junctions"],
                    visible=False)
        assert star_data_coll_id != '', "No collection for preliminary mapped reads"
        # Get the id for "Infer experiment"
        tool_id = get_working_tool_id(config["tool_names"]["infer_experiment"])
        # Create the input datamap for "Infer experiment"
        datamap = {
            "input" : {'batch': True,'values': [
                {'src':'hdca',
                'id': star_data_coll_id}]},
            "refgene": {'src':'hda','id': annotation_id},
            "sample_size": "200000",
            "mapq": "30"
        }
        # Run "Infer experiment"
        info = gi.tools.run_tool(hist, tool_id, datamap)
        # Retrieve the generated collection and rename it
        for ds in gi.histories.show_history(hist, contents=True, visible=True):
            if ds["history_content_type"] != 'dataset_collection':
                continue
            if ds["name"].find(config["tool_names"]["infer_experiment"]) == -1:
                continue
            gi.histories.update_dataset_collection(
                hist,
                ds["id"],
                name=config["collection_names"]["preliminary_mapping"]["infer_experiment"])


rule launch_star:
    '''
    Launch STAR on the trimmed data and MultiQC report on the STAR logs
    '''
    run:
        # Get the id for the tool to downsample the datasets
        tool_id = get_working_tool_id(config["tool_names"]["star"])
        # Search for the collection id with the trimmed data
        input_data_coll_id = get_working_collection_id(
            config["collection_names"]["trim_galore"]["trimmed"])
        # Extract the annotation dataset id
        annotation_id = get_annotation_id()
        # Create the input datamap
        datamap = {
            "singlePaired|sPaired": "single",
            "singlePaired|input1" : {'batch': True,'values': [
                {'src':'hdca',
                'id': input_data_coll_id}]},
            "refGenomeSource|geneSource": "indexed",
            "refGenomeSource|GTFconditional|GTFselect": "without-gtf",
            "refGenomeSource|GTFconditional|genomeDir": "mm10",
            "refGenomeSource|GTFconditional|sjdbGTFfile": {
                'src':'hda',
                'id': annotation_id},
            "refGenomeSource|GTFconditional|sjdbOverhang": "100"
        }
        # Run STAR
        info = gi.tools.run_tool(hist, tool_id, datamap)
        # Retrieve the generated collections and rename them
        log_coll_id = ''
        for ds in gi.histories.show_history(hist, contents=True, visible=True):
            if ds["history_content_type"] != 'dataset_collection':
                continue
            if ds["name"].find(config["tool_names"]["star"]) == -1:
                continue
            if ds["name"].endswith("log"):
                log_coll_id = ds["id"]
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["star"]["log"])
            elif ds["name"].endswith("mapped.bam"):
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["star"]["mapped"])
            elif ds["name"].endswith("splice junctions.bed"):
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["star"]["splice_junctions"])
        assert log_coll_id != '', "No collection for %s" % config["collection_names"]["star"]["log"]
        # Launch MultiQC on the log data
        run_multiqc(log_coll_id, "rnastar_log",config["tool_names"]["star"])


rule launch_feature_counts:
    '''
    Launch featureCounts on the mapped reads and MultiQC report on the outputs
    '''
    run:
        # Get the id for the tool to downsample the datasets
        tool_id = get_working_tool_id(config["tool_names"]["feature_counts"])
        # Search for the collection id with the trimmed data
        input_data_coll_id = get_working_collection_id(
            config["collection_names"]["star"]["mapped"])
        # Extract the dataset id of the annotation in the history
        annotation_id = get_annotation_id()
        # Create the input datamap for 
        datamap = {
            "alignment" : {'batch': True,'values': [
                {'src':'hdca',
                'id': input_data_coll_id}]},
            "gtf_source|ref_source": "history",
            "gtf_source|reference_gene_sets": {'src':'hda','id': annotation_id},
            "format": "tabdel_short",
            "strand_specificity": "0",
            "multimapping_enabled|multimapping_counts": "",
            "mapping_quality": "10",
            "largest_overlap": "",
            "min_overlap": "1",
            "read_extension_5p": "0",
            "read_extension_3p": "0",
            "read_reduction": "",
            "primary": "",
            "ignore_dup": "",
            "count_split_alignments_only": ""}
        # Run featureCounts
        info = gi.tools.run_tool(hist, tool_id, datamap)
        # Retrieve the generated collections and rename them
        summary_coll_id = ''
        for ds in gi.histories.show_history(hist, contents=True, visible=True):
            if ds["history_content_type"] != 'dataset_collection':
                continue
            if ds["name"].find(config["tool_names"]["feature_counts"]) == -1:
                continue
            if ds["name"].endswith("summary"):
                summary_coll_id = ds["id"]
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["feature_counts"]["summary"])
            else:
                # Rename the collection
                gi.histories.update_dataset_collection(
                    hist,
                    ds["id"],
                    name=config["collection_names"]["feature_counts"]["counts"])
        assert summary_coll_id != '', "No collection for %s" % config["collection_names"]["feature_counts"]["summary"]
        # Launch MultiQC on the log data
        run_multiqc(
            summary_coll_id,
            "featurecounts",
            config["tool_names"]["feature_counts"])


rule run_age_dge_deseq:
    '''
    Run the 4 differential gene expression analyses based on the the age
    '''
    run:
        # Search for the collection id with the count data
        input_data_coll_id = get_working_collection_id(
            config["collection_names"]["feature_counts"]['counts'])
        # Create a dataset dictionary to collect dataset info (id in correct section)
        analysis_datasets = {}
        # Parse the content of the collection and assign datasets to the correct analysis
        coll_content = gi.histories.show_dataset_collection(hist,input_data_coll_id)
        for el in coll_content['elements']:
            el_id = el['element_identifier']
            dataset_id = el['object']['id']
            # Extract metadata and build analysis id
            mice_type, age, gender = extract_dataset_metadate(el_id)
            analysis_id = "age_%s_%s" % (mice_type, gender)
            analysis_datasets.setdefault(analysis_id, {'ages': {}, 'project_lane': {}})
            # Extract the project and lane info
            project_id = finename_desc_df.at[el_id, 'Project id']
            lane = finename_desc_df.at[el_id, 'Lane'].replace("&", "_").replace(" ", "")
            project_lane = "%s_%s" % (project_id, lane)
            # Add dataset id in the correct age section
            analysis_datasets[analysis_id]['ages'].setdefault(age, [])
            analysis_datasets[analysis_id]['ages'][age].append(dataset_id)
            # Add dataset id in the correct project_lane section (and age: to check that )
            analysis_datasets[analysis_id]['project_lane'].setdefault(age, {})
            analysis_datasets[analysis_id]['project_lane'][age].setdefault(project_lane, [])
            analysis_datasets[analysis_id]['project_lane'][age][project_lane].append(dataset_id)
        # Prepare datamap for each analysis
        datamaps = {}
        for an in analysis_datasets:
            # Create the datamap
            datamaps[an] = prepare_deseq_datamap()
            # Retrieve the ages and sort them
            ages = list(analysis_datasets[an]['ages'].keys())
            ages = sorted(ages, key = lambda x: int(x.split("w")[0]))
            # Fill the datamap with the age factor and check if lane must be added as second factor
            # (if it has different datasets than the first factor)
            datamaps[an]["rep_factorName_0|factorName"] = "age"
            datamap_lane_preparation = {}
            add_lane_as_second_factor = False
            age_nb = 0
            for idx,ag in enumerate(ages):
                # Fill the first factor with the age as level and related datasets
                datamaps[an]["rep_factorName_0|rep_factorLevel_%s|factorLevel" % (idx)] = ag
                datamaps[an]["rep_factorName_0|rep_factorLevel_%s|countsFile" % (idx)] = []
                for ds in analysis_datasets[an]['ages'][ag]:
                    datamaps[an]["rep_factorName_0|rep_factorLevel_%s|countsFile" % (idx)].append(
                        {'src':'hda', 'id': ds})
                # Check if there is different project_lane for a same age
                if len(analysis_datasets[an]['project_lane'][ag].keys()) > 1:
                    add_lane_as_second_factor = True
                for pl in analysis_datasets[an]['project_lane'][ag]:
                    datamap_lane_preparation.setdefault(pl, [])
                    for ds in analysis_datasets[an]['project_lane'][ag][pl]: 
                        datamap_lane_preparation[pl].append(ds)
            # Add the project_lane as second factor if needed
            if add_lane_as_second_factor:
                datamaps[an]["rep_factorName_1|factorName"] = "project_lane"
                for idx,pl in enumerate(datamap_lane_preparation):
                    datamaps[an]["rep_factorName_1|rep_factorLevel_%s|factorLevel" % (idx)] = pl
                    datamaps[an]["rep_factorName_1|rep_factorLevel_%s|countsFile" % (idx)] = []
                    for ds in datamap_lane_preparation[pl]:
                        datamaps[an]["rep_factorName_1|rep_factorLevel_%s|countsFile" % (idx)].append(
                            {'src':'hda', 'id': ds})
        # Launch DESeq and change the generated names to have the analysis id in it
        launch_deseq_analyses(datamaps, "age")


rule run_gender_dge_deseq:
    '''
    Run the 6 differential gene expression analyses based on the gender
    '''
    run:
        # Search for the collection id with the count data
        input_data_coll_id = get_working_collection_id(
            config["collection_names"]["feature_counts"]['counts'])
        # Create a dataset dictionary to collect dataset info (id in correct section)
        analysis_datasets = {}
        # Parse the content of the collection and assign datasets to the correct analysis
        coll_content = gi.histories.show_dataset_collection(hist,input_data_coll_id)
        for el in coll_content['elements']:
            el_id = el['element_identifier']
            dataset_id = el['object']['id']
            # Extract metadata and build analysis id
            mice_type, age, gender = extract_dataset_metadate(el_id)
            analysis_id = "gender_%s_%s" % (mice_type, age)
            analysis_datasets.setdefault(analysis_id, {'genders': {}, 'project_lane': {}})
            # Extract the project and lane info
            project_id = finename_desc_df.at[el_id, 'Project id']
            lane = finename_desc_df.at[el_id, 'Lane'].replace("&", "_").replace(" ", "")
            project_lane = "%s_%s" % (project_id, lane)
            # Add dataset id in the correct age section
            analysis_datasets[analysis_id]['genders'].setdefault(age, [])
            analysis_datasets[analysis_id]['genders'][gender].append(dataset_id)
            # Add dataset id in the correct project_lane section (and age: to check that )
            analysis_datasets[analysis_id]['project_lane'].setdefault(age, {})
            analysis_datasets[analysis_id]['project_lane'][gender].setdefault(project_lane, [])
            analysis_datasets[analysis_id]['project_lane'][gender][project_lane].append(dataset_id)
        # Prepare datamap for each analysis
        datamaps = {}
        for an in analysis_datasets:
            # Create the datamap
            datamaps[an] = prepare_deseq_datamap()
            # Retrieve the ages and sort them
            genders = list(analysis_datasets[an]['genders'].keys())
            genders.sort()
            # Fill the datamap with the age factor and check if lane must be added as second factor
            # (if it has different datasets than the first factor)
            datamaps[an]["rep_factorName_0|factorName"] = "age"
            datamap_lane_preparation = {}
            add_lane_as_second_factor = False
            gender_nb = 0
            for idx,ge in enumerate(genders):
                # Fill the first factor with the age as level and related datasets
                datamaps[an]["rep_factorName_0|rep_factorLevel_%s|factorLevel" % (idx)] = ge
                datamaps[an]["rep_factorName_0|rep_factorLevel_%s|countsFile" % (idx)] = []
                for ds in analysis_datasets[an]['ages'][ge]:
                    datamaps[an]["rep_factorName_0|rep_factorLevel_%s|countsFile" % (idx)].append(
                        {'src':'hda', 'id': ds})
                # Check if there is different project_lane for a same age
                if len(analysis_datasets[an]['project_lane'][ge].keys()) > 1:
                    add_lane_as_second_factor = True
                for pl in analysis_datasets[an]['project_lane'][ge]:
                    datamap_lane_preparation.setdefault(pl, [])
                    for ds in analysis_datasets[an]['project_lane'][ge][pl]: 
                        datamap_lane_preparation[pl].append(ds)
            # Add the project_lane as second factor if needed
            if add_lane_as_second_factor:
                datamaps[an]["rep_factorName_1|factorName"] = "project_lane"
                for idx,pl in enumerate(datamap_lane_preparation):
                    datamaps[an]["rep_factorName_1|rep_factorLevel_%s|factorLevel" % (idx)] = pl
                    datamaps[an]["rep_factorName_1|rep_factorLevel_%s|countsFile" % (idx)] = []
                    for ds in datamap_lane_preparation[pl]:
                        datamaps[an]["rep_factorName_1|rep_factorLevel_%s|countsFile" % (idx)].append(
                            {'src':'hda', 'id': ds})
        # Launch DESeq and change the generated names to have the analysis id in it
        launch_deseq_analyses(datamaps, "gender")