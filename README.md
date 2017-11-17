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

- Add the API key to [https://galaxy.uni-freiburg.de/](https://galaxy.uni-freiburg.de/) in [`config.yaml`](config.yaml)

# Usage for the data analyses

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

- Upload the data on [http://galaxy.uni-freiburg.de/](http://galaxy.uni-freiburg.de/) inside an history named "NeuroMac: GF mices - DGE analysis"

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

# Usage for the report generation

```
$ pandoc doc/report.md --latex-engine=xelatex --filter pandoc-citeproc  --toc -o doc/report.pdf
```
