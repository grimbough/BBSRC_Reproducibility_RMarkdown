You can copy and paste the following code into an R session to install the packages 
required for the R Markdown session.

```r
if(!requireNamespace(BiocManager))
  install.packages('BiocManager');
  
pkgs <- c('BiocWorkflowTools', 'BiocStyle', 'rticles', 'knitr', 'rmarkdown')
pkgs <- pkgs[!pkgs %in% installed.packages()[,'Package']]
if(length(pkgs))
  BiocManager::install(pkgs)
```