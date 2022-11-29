# Nm Nanopore Drosophila
This repository contains the scripts for the analysis of EpiNano results produced from direct RNA sequencing of rRNA of Drosophila upon Knock-Down of Fibrillarin (FBL). 
The script is customised for the analysis of modifications whose signature in Nanopore sequencing is not limited to the modified position and diffuses in the neighboring positions (5 bases interval - here referred to as "kmer"), such as 2'-O-methylation (Nm). 

In "results" folder: 

- .tsv.per.site.var.per_site_var.5mer.csv files contain the output produced by EpiNano (https://github.com/enovoa/EpiNano) 

- FBL_kmer_rep1_colors.bed contains a color code for the EpiNano scores of each kmer centered in each position in the reference. This bed file allows to visualise the scores (yellow to red= low to high score) in IGV for comparison. 

- FBL_kmer_rep1_all.txt contains the EpiNano scores of each kmer centered in each position in the reference. Columns 4 and 5 contain the sum of errors for each kmer in non treated and FBL KD respectively, column 6 contains the difference between the sum of errors in the two conditions, referred to as EpiNano score.
