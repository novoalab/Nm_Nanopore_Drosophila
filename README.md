# Nm Nanopore Drosophila
This repository contains the scripts for the analysis of EpiNano results produced from direct RNA sequencing of rRNA of Drosophila upon Knock-Down of Fibrillarin (FBL). 
The script is customised for the analysis of modifications whose signature in Nanopore sequencing is not limited to the modified position and diffuses in the neighboring positions (5 bases interval - here referred to as "kmer"). 

In "results" folder: 
- Nm_pos.bed is a bed file with rRNA modifications predicted by sequence alignment with mouse and human
- .tsv.per.site.var.per_site_var.5mer.csv files contain the output produced by EpiNano (https://github.com/enovoa/EpiNano)
- FBL_kmer_replicable_positions files contain positions that are differentially modified between non treated and FBL knocked-down cells in 2 replicates

In "scores" folder: 
-FBL_kmer_rep1_colors.bed contains a color code for the EpiNano scores of each kmer centered in each position in the reference. This bed file allows to visualise the scores (yellow to red= low to high score) in IGV for comparison. 

-nanoRMS_score_FBL_kmer_rep1_colors.bed contains a color code for the nanoRMS scores (https://github.com/novoalab/nanoRMS#3-rna-modification-stoichiometry-estimation-using-tombo-resquiggling) of position in the reference. This bed file allows to visualise the scores (yellow to red= low to high score) in IGV for comparison. 

-rep1_FBL_kmer_scores.txt contains the EpiNano scores of each kmer centered in each position in the reference. Columns 4 and 5 contain the sum of errors for each kmer in non treated and FBL KD respectively, column 6 contains the difference between the sum of errors in the two conditions, referred to as EpiNano score.
