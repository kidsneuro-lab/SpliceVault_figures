# SpliceVault Figures  

Code used to create figures 1C and 2A-E in Dawes et al. 2022(?). Developed in R version 3.6

# Source files  

`src/` contains 300K-RNA as well as information on variants analysed (`fg_vars_events.tsv`).  
  
`300KRNA_filt.tsv.gz` is a filtered version of 300K-RNA for variants analysed in the paper.  
  
For code to create the full database see `github.com/kidsneuro-lab/300K-RNA` (file omitted from this repo due to its size).  
  
To use 300K-RNA see [SpliceVault](https://kidsneuro.shinyapps.io/splicevault/)  

`ref/` contains:  
  
  1. spliceAI deltas +/-5000 nt of each variant (`FG_variants_spliceAI_5000nt_plusindels.tsv`), derived using [this repository](https://github.com/kidsneuro-lab/SpliceAILookup).  
  
  2. MMSplice predictions (`mmsplice_predictions_all.csv`), derived using the author's [GitHub repository](https://github.com/gagneurlab/MMSplice_MTSplice ) 

# SpliceAI delta score processing
1. `SAI_scores.Rmd` filters out deltas < 0.05 and then exports SpliceAI deltas for manual annotation with predicted mis-splicing events (`sai_preds.xlsx`).  
  
2. `sai_preds_complete.xlsx` contains the manually annotated SpliceAI-predicted mis-splicing events, using interpretive rules described in the paper.  
  

# Creating figures
1. data preparation: `processing_predictions.Rmd` merges predictions by 300K-RNA, SpliceAI, MMSplice, along with the actual mis-splicing events observed in RNA studies.    
  
2. Figures: `figures.Rmd` contains code used to create figures. Ouputs to `figs/`   

