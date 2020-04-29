Atlas of Microglia-Microbiota in Aging (AMMA)
=============================================

Microglia, the brain resident macrophages, display high plasticity in response to their environment. Aging of the central nervous system (CNS), where microglial physiology is especially disrupted, is a major risk factor for a myriad of neurodegenerative diseases. Therefore, it is crucial to decipher intrinsic and extrinsic factors, like sex and the microbiome, that potentially modulate this process.

We found that **microglia follow sex-dependent dynamics in aging**. This repository stores the transcriptomics data analyses and the sources for the website explaining the analysis.

# Transcriptomics data analyses

## Requirements

- [conda](https://docs.conda.io/en/latest/miniconda.html)
- Creation of the `conda` environment with all the requirements

    ```
    $ conda env create -f environment.yml
    ```

- Launch the `conda` environment

    ```
    $ conda activate amma
    ```

## Preparation of files from the sequencing facility

1. Rename the files from the sequencing facility to follow the naming convention

    ```
    $ python src/copy_rename_raw_files.py \
        --input_dir <path to input directory> \
        --file_name_description <path to csv file with the correspondance between directory structure and sample name (from the Google drive)> \
        --output_dir <path to output directory>\
    ```

The naming convention for each sample is `microbiota_age_sex_replicate`

## From sequences to gene counts (in Galaxy)

1. Upload the data on Galaxy (e.g. [https://usegalaxy.eu/](https://usegalaxy.eu/)) inside a data library
2. Update the details in [`config.yaml`](config.yaml), specially the API key
3. Prepare the history in Galaxy (import the files from the data library, merge the files sequenced on 2 different lanes (for Project_S178 and Project_S225) and move the input files into collections)

    ```
    $ python src/prepare_data.py
    ```

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

Once it is finished:

1. Download the generated count table and the gene length file
2. Put these files in the `data` folder

## Differentially Expression Analysis (locally using Jupyter Notebooks)

1. Launch Jupyter

    ```
    $ jupyter notebook
    ```

3. Open [http://localhost:8888/tree/](http://localhost:8888/tree/)
2. Move to `src` in Jupyter
3. Prepare the differential expression analysis
    1. Open `prepare_data.ipynb` and execute all cells
    2. Open `dge_analysis.ipynb` and execute all cells
    3. Open `pre-visualization.ipynb` and execute all cells

4. Analyze the differentially expressed genes given different comparisons
    1.  Effect of the microbiota (GF vs SPF)

        Analysis | Notebook
        --- | ---
        Microbiota effect for both sexes, after controlling for age | `microbiota-effect-sex.ipynb`
        Microbiota effect for the 3 ages, after controlling for sex | `microbiota-effect-age.ipynb`
        Microbiota effect for the 3 ages and both sexes | `microbiota-effect-age-sex.ipynb`

    2. Effect of the sex (Male vs Female)

        Analysis | Notebook
        --- | ---
        Sex effect for both microbiotas, after controlling for age | `sex-effect-microbiota.ipynb`
        Sex effect for the 3 ages, after controlling for microbiota | `sex-effect-age.ipynb`
        Sex effect for the 3 ages and both microbiotas | `sex-effect-microbiota-age.ipynb`

    3. Effect of the ages (Middle-aged vs Young, Old vs Young and Old vs Middle-aged)

        Analysis | Notebook
        --- | ---
        Age effect for both microbiotas, after controlling for sex | `age-effect-microbiota.ipynb`
        Age effect for the both sexes, after controlling for microbiota | `age-effect-sex.ipynb`
        Age effect for the both sexes and both microbiotas | `age-effect-microbiota-sex.ipynb`


# Website

This folder stores the sources of the website describing the analyses in `docs` folder.
Reports from the Jupyter notebooks are available there to show the different steps and images. 

## Generate HTML reports from the Jupyter Notebooks

```
$ jupyter nbconvert --to=html src/*.ipynb --output-dir docs/
```

These reports are stored in the `docs` folder and are linked on the website.

## Generate the website locally

1. Install [Jekyll](https://jekyllrb.com/docs/installation/)
2. Move to `docs` folder

    ```
    $ cd docs
    ```
    
3. Install the plugins for Jekyll (only once)

    ```
    $ bundle install
    ```

4. Serve the website locally

    ```
    $ bundle exec jekyll serve
    ```

5. Open [http://127.0.0.1:4000/amma/](http://127.0.0.1:4000/amma/)