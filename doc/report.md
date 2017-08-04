---
title: Neuromac project on the differential gene expression in germ free mices
author: Bérénice Batut
bibliography: doc/ref.bib
geometry: margin=2cm
---

## Presentation of the biological context

#### Questions

- What are the differentially expressed genes of the microglya cells between germ free (GF) and conventional (SPF) mice at different ages?
- Is there a difference between the ages?
- Is there a difference between the genders?
- Which genes? Which pathways?

#### Hypotheses

- Differential gene expression between SPF and GF already known for 8 weeks mices ([article](http://www.nature.com/neuro/journal/v18/n7/abs/nn.4030.html))
- Expected gender bias in 104w SPF

## Data

- RNA-seq data of microglya cells of mices: (Figure \ref{replicates} and [details](https://docs.google.com/spreadsheets/d/1DL8pEVj5cvGflPIiaSPRXy-dMk2S7CxmnIk6Ubta2xs/edit?usp=sharing)
	- 2 types
    	- Conventional (SPF)
        - Germ free (GF)
    - 3 ages
    	- 8 weeks (8w)
        - 52 weeks (52w)
        - 104 weeks (104w)
    - 2 genders
  		- Male (m)
        - Female (f)
        
![Repartition of the replicates in the different groups\label{replicates}](doc/images/input_data.png){ height=256px }
 
- Library preparation and sequencing
	- Library preparation with Illumina Nextera XT Sample Preparation 
	- First strand
    - Single-end data
    - Sequencing with HiSeq 1000

## Analyses

Workflow applied on each dataset: Figure \ref{workflow}

![Workflow applied on each dataset\label{workflow}](doc/images/workflow.png)

### Quality control and trimming

#### Quality control

Objectives

- Check the quality of the raw sequences

Details

- FastQC on every datasets
- MultiQC [@ewels2016multiqc] report to aggregate the FastQC reports

Results: (details in the [Google doc](https://docs.google.com/spreadsheets/d/1DL8pEVj5cvGflPIiaSPRXy-dMk2S7CxmnIk6Ubta2xs/edit?usp=sharing) in "FastQC report")

- Small sequences; 50 bp
- Between 22.4M and 52.2M of reads
- Good average quality score per sequence (Figure \ref{fastqc_per_sequence_quality_scores_plot})

    ![Per Sequence Quality Scores (generated with FastQC and MultiQC)\label{fastqc_per_sequence_quality_scores_plot}](results/fastqc/fastqc_per_sequence_quality_scores_plot.png){ width=500px }

- No big reduction of quality at the end of the sequences (Figure \ref{fastqc_per_base_sequence_quality_plot})
    
    ![Sequence Quality Histograms (generated with FastQC and MultiQC)\label{fastqc_per_base_sequence_quality_plot}](results/fastqc/fastqc_per_base_sequence_quality_plot.png){ width=500px }

- Per Base Sequence Content: 52 samples with a warning and 2 with an error (non homogenous proportion of each base at the beginning of the sequences)
- Per Sequence GC Content (Figure \ref{fastqc_per_sequence_gc_content_plot}): all samples with a warning (small increase of %GC around 10-20 bp, likely related to the overrepresentation of C and G sequences at the begining seen in the previous check)

    ![Per Sequence GC Content (generated with FastQC and MultiQC)\label{fastqc_per_sequence_gc_content_plot}](results/fastqc/fastqc_per_sequence_gc_content_plot.png){ width=500px }

- Per Base N Content (Figure \ref{fastqc_per_base_n_content_plot}): 10 samples with a warning

    ![Per Base N Content (generated with FastQC and MultiQC)\label{fastqc_per_base_n_content_plot}](results/fastqc/fastqc_per_base_n_content_plot.png){ width=500px }

- Numerous read duplications (Figure \ref{fastqc_sequence_duplication_levels_plot}): expected for RNA seq data

    ![Sequence Duplication Levels (generated with FastQC and MultiQC)\label{fastqc_sequence_duplication_levels_plot}](results/fastqc/fastqc_sequence_duplication_levels_plot.png){ width=500px }

- Overrepresented sequences (Figure \ref{fastqc_overrepresented_sequencesi_plot}): 55 samples with a warning and 1 with an error

    ![Overrepresented sequences (generated with FastQC and MultiQC)\label{fastqc_overrepresented_sequencesi_plot}](results/fastqc/fastqc_overrepresented_sequencesi_plot.png){ width=500px }

    - Need to check with BLAST what are the overrepresented sequences

- Few remaining adapters at the end of the sequences (Figure \ref{fastqc_adapter_content_plot})

    ![Adapter Content (generated with FastQC and MultiQC)\label{fastqc_adapter_content_plot}](results/fastqc/fastqc_adapter_content_plot.png){ width=500px }

#### Trimming

Objectives

- Eliminate the bad quality ends, even if globally quality reports are good. Not so much removed

Details

- Trim Galore! on every datasets
- MultiQC [@ewels2016multiqc] report to aggregate the Trim Galore! reports

Results

- Few 3'-ends bases eliminated
- Few sequences eliminated (too small after trimming)
- To do
    - Add table with raw results
    - Fix MultiQC report aggregation

### Mapping

Objectives

- Map the reads on the reference genome of Drosophila

#### Preliminary mapping

Objectives

- Infer if the library is strand specific or not by
    - Extraction of some reads
    - Mapping them on the reference genome
    - Use an annotation file
    - Infer on gene of which strand the mapping reads fit on

Details

- Downsampling of the dataset: Extraction of 200,000 reads with "Select first" tool
- Mapping with
    - STAR [@dobin2013star]
    - mm10 as reference genome
- Infer the strand with "Infer Experiment" of RSeQC [@wang2012rseqc]

Results

- Unstranded library
- To do
    - Add table with raw results
    - Add MultiQC report aggregation

#### Actual mapping

Objectives

- Map the trimmed reads on the reference genome to annotate them

Details

- Mapping with
    - STAR [@dobin2013star], a splice aware mapper
    - mm10 as reference genome
    - mm10_UCSC_07_15_genes as gene model for splice junctions

Results: (Figure \ref{star_alignment_plot} and details in the [Google doc](https://docs.google.com/spreadsheets/d/1DL8pEVj5cvGflPIiaSPRXy-dMk2S7CxmnIk6Ubta2xs/edit?usp=sharing) in "STAR report")

![Alignment Scores (generated with STAR and MultiQC)\label{star_alignment_plot}](results/star/star_alignment_plot.png){ width=500px }

- Good percentage of reads that are uniquely mapped (> 70.9%)
- Good percentage of reads that are mapped (uniquely + multi-mapped, > 85%)
- Relatively low percentage of unmapped reads
- Keep an eye on
    - `SPF_8w_F_2`
    - `SPF_52w_F_6`
    - `SPF_52w_F_1`
    - `SPF_8w_M_1`

### Gene counting

Objectives

- Count the number of reads that are mapped on genes

Details

- Counting with
    - FeatureCounts [@liao2013featurecounts]
    - mm10_UCSC_07_15_genes as Gene annotation file
    - Unstranded protocol

Results (Figure \ref{featureCounts_assignment_plot} and details in the [Google doc](https://docs.google.com/spreadsheets/d/1DL8pEVj5cvGflPIiaSPRXy-dMk2S7CxmnIk6Ubta2xs/edit?usp=sharing) in "FeatureCounts report")

![Reads assigned (generated with FeatureCounts and MultiQC)\label{featureCounts_assignment_plot}](results/featureCounts/featureCounts_assignment_plot.png){ width=500px }

- Good percentage of assigned reads to genes (> 60 % except 2 samples)
- Keep an eye on
    - `SPF_8w_F_2` (already issue for the mapping)
    - `SPF_8w_F_5`

### Differential expression analyses

Progressive complexification of the analyses

#### Understanding the impact of age: for each gender and each mouse type (GF & SPF)

Questions

- What are the differences between the ages?
- Which genes and pathways are differentially expressed? 
- In which proportion?

Objectives

- To get a better idea of the basal differences between the ages

Details

- Design (Figure \ref{age_impact_analysis})

    ![Analyses to understand the impact of age on gene expression\label{age_impact_analysis}](doc/images/age_impact_analysis.png)

- Analyses with
    - DESeq2 [@love2014moderated]
    - Factor 1: age
    - Factor 2: lane with sequencing project (*e.g.* `project_s195_7`) if different from the age factors
	
#### Understanding the impact of gender: for each age and each type

- Question: How complex are the differences between the gender? Can we have the both gender in a global analysis?
- Why? To get a better idea of the impact of gender on the gene expression, if it is different between the ages and to know if the impact of the gender if more important than the type of mices (GF/SPF)

![Analyses to understand the impact of gender on gene expression](https://framadrive.org/s/Se503gvrZ7apGa8/download)
    
- Factor 1: gender
- Factor 2: lane with sequencing project

#### Understanding the impact of type

- Question: What is the impact of being GF on the expression levels?

1. for each gender and age

	- Question: What are the differentially expressed genes of the microglya cells between GF and SPF mice at different ages?
    - Why? To have limited cofounding factors (age/gender) and still get an idea of the impact of SPF/GF on the gene expression
    
    ![Analyses to understand the impact of GF on gene expression, with limitation of the cofounding factors](https://framadrive.org/s/viE59uhFcH1IRpF/download)
    
    - Factor 1: GF or SPF
    - Factor 2: lane with sequencing project
2. for each gender

	- Question: What are the levels of differential expression between GF and SPF between the ages?
    - Why? To be able to compare the gene expression of GF/SPF between the different ages
    
    ![Analyses to understand the impact of GF on gene expression](https://framadrive.org/s/wSMBMJv8WSTRRVu/download)
    
	- Factor 1: GF or SPF
    - Factor 2: age
    - Factor 3: lane with sequencing project
            
### Post differential expression analyses

To keep in mind

- Comparison of changes between the ages after normalization by the values in 8 weeks
- Comparison of the levels of differential expression between gender comparison and GF/SPF to check if the first ones are not higher than the latter ones to have a big analysis


## References


