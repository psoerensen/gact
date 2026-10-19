Sys.setenv(RSTUDIO_PANDOC='C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools')
.libPaths(c('C:/Users/au223366/AppData/Local/R/win-library/4.4',.libPaths()))
root<-'C:/Users/au223366/Documents/GitHub/gact/Document'
static<-c('Benchmark_overview.Rmd','Tutorial_timings_and_memory.Rmd')
examples<-c('Gene_ranking_bayesian_pops.Rmd','Gene_set_analysis_bayesian_magma.Rmd')
for(name in static)stopifnot(!any(grepl('^```\\{',readLines(file.path(root,name),warn=FALSE))))
for(name in examples){
 lines<-readLines(file.path(root,name),warn=FALSE);headers<-grep('^```\\{r',lines,value=TRUE)
 stopifnot(any(grepl('opts_chunk\\$set.*eval=FALSE',lines))||all(grepl('eval=FALSE|setup',headers)))
 stopifnot(!any(grepl('eval=TRUE',headers)))
}
for(name in c(static,examples))rmarkdown::render(file.path(root,name),quiet=TRUE,envir=new.env(parent=globalenv()))
cat('Rendered four edited pages with analysis chunks disabled.\n')
