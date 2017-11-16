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

    ```
    # Launch FastQC
    $ snakemake --snakefile src/run_rna_seq_analyses.py launch_fastqc
    # Launch TrimGalore!
    $ snakemake --snakefile src/run_rna_seq_analyses.py launch_trim_galore
    ```

- Launch mapping

    ```
    # Launch preliminary mapping
    $ snakemake --snakefile src/run_rna_seq_analyses.py launch_preliminary_mapping
    # Launch STAR mapping
    $ snakemake --snakefile src/run_rna_seq_analyses.py launch_star
    ```

- Launch counting
    
    ```
    $ snakemake --snakefile src/run_rna_seq_analyses.py launch_feature_counts
    ```

- Launch the DGE analyses

    ```
    # Run DGE analses for the ages
    $ snakemake --snakefile src/run_rna_seq_analyses.py run_age_dge
    # Run DGE analses for the gender
    $ 
    # Run DGE analses for the type of mices
    $ 
    ```

# Usage for the report generation

```
$ pandoc doc/report.md --latex-engine=xelatex --filter pandoc-citeproc  --toc -o doc/report.pdf
```