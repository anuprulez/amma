Neuromac project on the differential gene expression in germ free mices
=======================================================================

# Biological context

## Questions

- What are the differentially expressed genes of the microglya cells between germ free (GF) and conventional (SPF) mice at different ages?
- Is there a difference between the ages?
- Is there a difference between the genders?
- Which genes? Which pathways?

## Hypotheses

- Differential gene expression between SPF and GF already known for 8 weeks mices ([article](http://www.nature.com/neuro/journal/v18/n7/abs/nn.4030.html))
- Expected gender bias in 104w SPF
- Differential gene expression between SPF and GF already known for embryonic weeks mices

# Data

We have RNA-seq data of microglya cells of mices for:
- both Conventional (SPF) and Germ free (GF)
- 3 ages: 8 weeks (8w), 52 weeks (52w), 104 weeks (104w)
- 2 genders: Male (m) and Female (f)
        
More details are available [here](data)

# Analyses

1. [Extraction of number of reads mapped on each annotated genes for each samples](gene_count_extraction)
2. Differential expression analyses
    1. [Building of the differential expression analysis (normalization of the counts)](dge_analysis)
    2. [Pre-visualization of the normalized count data before any differential analysis](Pre-visualization)
    3. [Effect of the age on the expressed genes](age-effect)
    4. Effect of the gender on the expressed genes
    5. Effect of the type on the expressed genes



