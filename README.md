# change-scores

Simulations showing inferential biases that can arise from using change scores.

Code was written by Mark Gilthorpe with a few modifications by myself. The code accompanies a paper that is being submitted. Further details will be added as soon is the paper is available as a preprint.

To run the simulations, execute the file `simulations.R`. Note that **this will generate CSV files in the current directory**. 

To generate an HTML version of Table 2 in the manuscript, knit the file `simulations.Rmd` using the command

```Rscript -e 'library(rmarkdown);render("simulations.Rmd")'```

