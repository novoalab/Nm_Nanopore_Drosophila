# Epinano-5mer for Nm detection with Nanopore

![image](https://user-images.githubusercontent.com/44866316/207355327-a730c0e6-42cc-4fb9-9e40-62264c431ae7.png)


## Table of contents
- [General Description](#General-description)
- [Input](#Input)
- [Usage](#Usage)
- [Results](#Results)
- [Dependencies and versions](#Dependencies-and-versions)
- [Citation](#Citation) 
- [Contact](#Contact) 

## General Description
This repository contains the scripts for the analysis of EpiNano results produced from direct RNA sequencing of rRNA of Drosophila upon Knock-Down of Fibrillarin (FBL). The script is customised for the analysis of modifications whose signature in Nanopore sequencing is not limited to the modified position and diffuses in the neighboring positions (5 bases interval - here referred to as "kmer"), such as 2'-O-methylation (Nm). 


## Usage
For 5mer information, we use the Epinano (please visit https://github.com/enovoa/EpiNano for this first step) script TSV_to_Variants_Freq.py3 to analyse the samtotsv output, that generates the .tsv.per.site.var.per_site_var.5mer.csv files for each sample and replicate.

If instead you have EpiNano output at per site level, you can extract information about the neighboring positions at the 5mer level, by running: 

```bash
python Slide_Variants.py file.plus_strand.per.site.csv 5
```

Once you have Epinano 5mer input ready, you can proceed to run Epinano_fivemer_analysis.R interactively in RStudio. 

**!!! The script requires to add custom input and output paths, file names and labels at the beginning !!!**

The script is based on 2 replicates of 2 conditions. If you have more replicates you will have to edit the script accordingly. 
Explanation of the steps of the script is included in the script itself.

## Input

- Epinano 5mer table for each sample ( .tsv.per.site.var.per_site_var.5mer.csv)

- a bed file with predicted modified sites if available (Nm_pos.bed)

You can find examples of the input files in the "input" folder. 

## Results
The script produces as output: 

- FBL_kmer_rep1_colors.bed and FBL_kmer_rep2_colors.bed contain a color code for the EpiNano scores of each kmer centered in each position in the reference. This bed file allows to visualise the scores (yellow to red= low to high score) in IGV for comparison. 

- FBL_kmer_rep1_all.txt and FBL_kmer_rep2_all.txt contain the EpiNano scores of each kmer centered in each position in the reference. Columns 4 and 5 contain the sum of errors for each kmer in non treated and FBL KD respectively, column 6 contains the difference between the sum of errors in the two conditions, referred to as EpiNano score. This file can be converted into a bedgraph by just changing extention from .txt to .bedgraph and can be loaded and visualized in IGV as barplot or heatmap.

- FBL_kmer_replicable_positions_18S.txt and FBL_kmer_replicable_positions_28S.txt contain a list of sites in which a position has a score higher than 3*median of the other positions' score in both replicates for each long rRNA transcript.

- a barplot of the scaled difference of the summed errors per 5mer (score) along the transcript per reference per comparison.


## Dependencies and versions


| Software | Version |
| --- | --- |
| RStudio | 4.2 |
| python | 3 |

## Citation
  
If you find this work useful, please cite: Sklias et al. 2023
  
## Contact
If you have any issues running this code, please go first over previous [issues](https://github.com/novoalab/Nm_Nanopore_Drosophila/issues). If you still can't figure it out based on the prior responses/issues raised, please open a new issue. Thanks!   
