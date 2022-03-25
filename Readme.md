# SpliceVault Figures  

Code used to create figures in Dawes et al. 2022(?). Developed in R version 3.6

# Source files  

`src/` contains 300K-RNA as well as information on variants analysed (`fg_vars_events.tsv`).  
  
`src/300KRNA_filt.tsv.gz` is a filtered version of 300K-RNA for variants analysed in the paper. To gain access to the full version use [SpliceVault](https://kidsneuro.shinyapps.io/splicevault/) or see data availability statement in the paper for full details.
  
For code to create the full database see `github.com/kidsneuro-lab/300K-RNA` (file omitted from this repo due to its size).  
  
To use 300K-RNA see [SpliceVault](https://kidsneuro.shinyapps.io/splicevault/)  

`ref/` contains extra files required to create figures. 

# Processing
`1.SAI_scores.Rmd` retrieves SpliceVault predictions for variants analysed in the paper.

`2.SAI_scores.Rmd` is the script we used to retrieve SpliceAI deltas +/-5000 nt of each variant using a custom API. An identical result can be easily achieved using code from the [SpliceAI github](https://github.com/Illumina/SpliceAI), but code is included here for completeness.

`3.SAI_scores_processing.Rmd` filters out deltas < 0.001 and then annotates remaining deltas with predicted mis-splicing events. The source file is not in the repository due to its size but can be retrieved using our source variants and code to score custom sequences available on the [SpliceAI github](https://github.com/Illumina/SpliceAI).
  
`sai_preds_complete.csv` contains the manually annotated SpliceAI-predicted mis-splicing events, using interpretive rules described in the paper.  

`4.merged_predictions.Rmd` pulls SpliceVault and SpliceAI predictions into a format amenable to plotting.

`mendeliome.Rmd` retrieves canonical ensembl transcripts for mendelian disease genes

# Creating figures

Notebooks in `figure_scripts/` can be run to recreate figures in the paper.
