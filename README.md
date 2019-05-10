NeuroMac data analysis
======================

Links
- [Notebook](https://monod.lelab.tailordev.fr/b57db1ec-cc35-47d6-8767-9140d0390bdc#THcB4zTs+1L1Hq9S26HTlJ0IZCjjlCzcnakq7mBGpJo=)
- [Raw data structure](https://docs.google.com/spreadsheets/d/1DL8pEVj5cvGflPIiaSPRXy-dMk2S7CxmnIk6Ubta2xs/edit#gid=0)

# Requirements

- conda
- Creation of the `conda` environment with all the requirements

    ```
    $ conda env create -f environment.yml
    ```

- Add the API key to [https://usegalaxy.eu/](https://usegalaxy.eu/) in [`config.yaml`](config.yaml)

# Run the data analyses

- Launch the `conda` environment

    ```
    $ source activate neuromac
    ```

- Copy and rename the files

    ```
    $ python src/copy_rename_raw_files.py \
        --input_dir <path to input directory> \
        --file_name_description <path to csv file with the correspondance between directory structure and sample name (from the Google drive)> \
        --output_dir <path to output directory>\
    ```

- Upload the data on [http://usegalaxy.eu/](http://usegalaxy.eu/) inside an history named "NeuroMac: GF mices - DGE analysis"

    ```
    $ snakemake --snakefile src/prepare_data.py
    ```

    You can change the name of the history in [`config.yaml`](config.yaml)

- Launch workflow to extract gene counts

    ```
    $ snakemake --snakefile src/extract_gene_counts.py
    ```

    The worklow do:
    1. Quality control and trimming using FastQC and Trim Galore!
    2. Preliminary mapping and experiment inference using STAR and RSeQC
    3. Mapping using STAR
    4. Gene counting using FeatureCounts

    The workflow is applied on each dataset (organized into data collection)

- Do the DGE analyses
    - Launch R and install the libraries

    ```
    $ source("src/install_libraries.R") 
    ```

    - Launch Jupyter

    ```
    $ jupyter notebook
    ```

    - Run the different notebooks in `src`
        1. `src/prepare_data.ipynb`
        2. `src/dge_analysis.ipynb`
        3. `src/pre-visualization.ipynb`
        4. DEG analysis

# Generate HTML of the Jupyter Notebooks (for `docs` folder)

```
$ jupyter nbconvert --template=nbextensions --to=html src/*.ipynb --output-dir docs/
$ for i in $(find docs/ -path "*.html"); do awk '/jquery.min.js/{ print; print "<script src=\"https://cdnjs.cloudflare.com/ajax/libs/three.js/r79/three.min.js\"></script>"; next }1' $i > tmp; mv tmp $i; done
```

# Generate the website locally

- Install [Jekyll](https://jekyllrb.com/docs/installation/)
- Move to `docs` folder

    ```
    $ cd docs
    ```

- Serve the website locally

    ```
    $ bundle exec jekyll serve
    ```

