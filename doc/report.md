Neuromac project on the differential gene expression in germ free mices
========

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

- RNA-seq data of microglya cells of mices: [details](https://docs.google.com/spreadsheets/d/1DL8pEVj5cvGflPIiaSPRXy-dMk2S7CxmnIk6Ubta2xs/edit?usp=sharing)
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
        
![Repartition of the replicates in the different groups](https://framadrive.org/s/6Fyvkkhh2YQ2WIY/download)       
 
- Library preparation and sequencing
	- Library preparation with Illumina Nextera XT Sample Preparation 
	- First strand
    - Single-end data
    - Sequencing with HiSeq 1000

## Analyses

![Workflow applied on each dataset](https://framadrive.org/s/rWFKBV1HZFUyoZq/download)

Analysis comments:

1. FastQC
	- Small sequences
    - Good global quality
	- Need to keep an eye on the GF_8w_M_2
    - Bad "Per Base Sequence Content" and then "Per Sequence GC Content"
    - 10 samples with warning for the "Per Base N Content"
    - Error for the "Sequence Duplication Levels"
    - Warning for the "Overrepresented sequences"
    	- Checking with BLAST what are the overrepresented sequences
2. Trim Galore!
	- Even if the quality report are good
    	- Small sequences of 50 bp
        - Not too much cutting
3. Preliminary mapping
	- Unstranded data
4. 

### Differential expression analyses with DESeq

Progressive complexification of the analyses

1. Understanding the impact of age: for each gender and each mouse type (GF & SPF)

	- Question: What are the differences between the ages? Which genes and pathways are differentially expressed? In which proportion?
    - Why? To get a better idea of the basal differences between the ages
    
    ![Analyses to understand the impact of age on gene expression](https://framadrive.org/s/DysSOzGaG8C6Aqm/download)

	- Factor 1: age
  	- Factor 2: lane with sequencing project (*e.g.* `project_s195_7`)
	
2. Understanding the impact of gender: for each age and each type

	- Question: How complex are the differences between the gender? Can we have the both gender in a global analysis?
    - Why? To get a better idea of the impact of gender on the gene expression, if it is different between the ages and to know if the impact of the gender if more important than the type of mices (GF/SPF)
    
    ![Analyses to understand the impact of gender on gene expression](https://framadrive.org/s/Se503gvrZ7apGa8/download)
        
	- Factor 1: gender
	- Factor 2: lane with sequencing project

3. Understanding the impact of type
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





