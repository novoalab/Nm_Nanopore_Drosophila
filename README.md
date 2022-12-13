# Epinano-5mer for Nm detection with Nanopore

![image](https://user-images.githubusercontent.com/44866316/207355327-a730c0e6-42cc-4fb9-9e40-62264c431ae7.png)


## Table of contents
- [General Description](#General-description)
- [Input](#Input)
- [Running the code](#Running-the-code)
- [Output](#Output)
- [Dependencies and versions](#Dependencies-and-versions)
- [Citation](#Citation) 
- [Contact](#Contact) 

## General Description
This repository contains the scripts for the analysis of EpiNano results produced from direct RNA sequencing of rRNA of Drosophila upon Knock-Down of Fibrillarin (FBL). The script is customised for the analysis of modifications whose signature in Nanopore sequencing is not limited to the modified position and diffuses in the neighboring positions (5 bases interval - here referred to as "kmer"), such as 2'-O-methylation (Nm). 

For 5mer information, we use the Epinano (https://github.com/enovoa/EpiNano) script TSV_to_Variants_Freq.py3 to analyse the samtotsv output, that generates the .tsv.per.site.var.per_site_var.5mer.csv files for each sample and replicate.

## Running the code

Usage: 
First, you'll need to install dependencies: 
```bash
R
install.packages("argparse")
```

For more details on how to use ModPhred, please see the [GitHub](https://github.com/novoalab/modPhred) repository and the [ReadTheDocs](https://modphred.readthedocs.io/en/latest/install.html) manual.

## Output
In "results" folder: 

- .tsv.per.site.var.per_site_var.5mer.csv files contain the output produced by EpiNano (https://github.com/enovoa/EpiNano) 

- FBL_kmer_rep1_colors.bed contains a color code for the EpiNano scores of each kmer centered in each position in the reference. This bed file allows to visualise the scores (yellow to red= low to high score) in IGV for comparison. 

- FBL_kmer_rep1_all.txt contains the EpiNano scores of each kmer centered in each position in the reference. Columns 4 and 5 contain the sum of errors for each kmer in non treated and FBL KD respectively, column 6 contains the difference between the sum of errors in the two conditions, referred to as EpiNano score.



```bash
Rscript Metagene_Plots.R -i Sample1.bed Sample2.bed -o Test -gtf Annotation.gtf -l Sample-1 Sample-2
```

## Dependencies and versions

Software | Version 
R | 4.2
 | 
 | 
 | 
 | 

## Citation
  
If you find this work useful, please cite: XXX
  
## Contact
If you have any issues running this code, please go first over previous [issues](https://github.com/novoalab/Nm_Epinano_5mer/issues). If you still can't figure it out based on the prior responses/issues raised, please open a new issue. Thanks!   
