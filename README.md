Atlas of Microglia-Microbiota in Aging (AMMA)
=============================================

Microglia, the brain resident macrophages, display high plasticity in response to their environment. Aging of the central nervous system (CNS), where microglial physiology is especially disrupted, is a major risk factor for a myriad of neurodegenerative diseases. Therefore, it is crucial to decipher intrinsic and extrinsic factors, like sex and the microbiome, that potentially modulate this process.

Using transcriptomics, we found that **microglia follow sex-dependent dynamics in aging**. This repository 

# Requirements

- [conda](https://docs.conda.io/en/latest/miniconda.html)
- Creation of the `conda` environment with all the requirements

    ```
    $ conda env create -f environment.yml
    ```

- Launch the `conda` environment

    ```
    $ conda activate amma
    ```

# Run the data analyses

## Prepare files from the sequencing facility

1. Rename the files from the sequencing facility to follow a certain naming convention

    ```
    $ python src/copy_rename_raw_files.py \
        --input_dir <path to input directory> \
        --file_name_description <path to csv file with the correspondance between directory structure and sample name (from the Google drive)> \
        --output_dir <path to output directory>\
    ```

2. Upload the data on Galaxy inside a data library.
3. Update the details in [`config.yaml`](config.yaml), specially the API key
4. Prepare the history in Galaxy (import the files from the data library, merge the files sequenced on 2 different lanes (for Project_S178 and Project_S225) and move the input files into collections)

    ```
    $ python src/prepare_data.py
    ```

## From sequences to gene counts (inside Galaxy)

5. Launch Galaxy workflow to extract gene counts

    ```
    $ python src/extract_gene_counts.py
    ```

    The worklow do:
    1. Quality control and trimming using FastQC and Trim Galore!
    2. Preliminary mapping and experiment inference using STAR and RSeQC
    3. Mapping using STAR
    4. Gene counting using FeatureCounts

The workflow is applied on each dataset (organized into data collection). It can take a while.

Once it is finished, please download the generated count table in `data` as well as the gene length file.

## Differentially Expression Analysis (locally using Jupyter Notebooks)

1. Launch Jupyter

    ```
    $ jupyter notebook
    ```

2. Move to `src` in Jupyter
3. Prepare the differential expression analysis
    1. Open `src/prepare_data.ipynb` and execute all cells
    2. Open `src/dge_analysis.ipynb` and execute all cells
    3. Open `src/pre-visualization.ipynb` and execute all cells
4. Analyze the differentially expressed genes given different comparisons
    1.  DEG analysis

# Generate HTML of the Jupyter Notebooks (for `docs` folder)

```
$ jupyter nbconvert --template=nbextensions --to=html src/*.ipynb --output-dir docs/
```

# Generate the website locally

- Install [Jekyll](https://jekyllrb.com/docs/installation/)
- Move to `docs` folder

    ```
    $ cd docs
    ```
    
- Install the plugins for Jekyll (only once)

    ```
    $ bundle install
    ```

- Serve the website locally

    ```
    $ bundle exec jekyll serve
    ```

