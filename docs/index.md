Differential gene expression in mice (germ free or normal, different ages, 2 genders)
=====================================================================================

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
    2. [Pre-visualization of the normalized count data before any differential analysis](pre-visualization)
    3. Effect of the age on the expressed genes
        1. [Comparison between the ages (after controlling for type and gender)](age-effect-general) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/age-effect/age)
        2. [Comparison between the ages for the types (after controlling for gender)](age-effect-type) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/age-effect/age_type)
        3. [Comparison between the ages for the genders (after controlling for type)](age-effect-gender) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/age-effect/age_gender)
        4. [Comparison between the ages for the genders and types: [report](age-effect-type-gender) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/age-effect/age_type_gender)
    4. Effect of the type on the expressed genes
        1. [Comparison between the types (after controlling for age and gender)](type-effect-general) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/type-effect/type)
        2. [Comparison between the types for the genders (after controlling for age)](type-effect-gender) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/type-effect/type_gender)
        3. [Comparison between the types for the ages (after controlling for gender)](type-effect-age) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/type-effect/type_age)
        4. [Comparison between the types for the genders and ages: [report](type-effect-age-gender) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/type-effect/type_gender_age)
    5. Effect of the gender on the expressed genes
        1. [Comparison between the genders (after controlling for type and age)](gender-effect-general) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/gender-effect/gender)
        2. [Comparison between the genders for the types (after controlling for age)](gender-effect-type) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/gender-effect/gender_type)
        3. [Comparison between the genders for the ages (after controlling for type)](gender-effect-age) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/gender-effect/gender_age)
        4. [Comparison between the genders for the types and ages](gender-effect-type-age) - [data](https://github.com/bebatut/neuromac_GF_mices/tree/master/results/dge/gender-effect/gender_type_age)



