.libPaths(c('C:/Users/au223366/AppData/Local/R/win-library/4.4',.libPaths()))
library(jsonlite)
root<-'work/revision-benchmark/scoring-grid'
paths<-list.files(root,pattern='timing.json',recursive=TRUE,full.names=TRUE)
stopifnot(length(paths)==40)
records<-lapply(paths,read_json,simplifyVector=TRUE)
checks<-list()
for(i in which(vapply(records,function(x)x$threads==1,logical(1)))){
 a<-records[[i]];j<-which(vapply(records,function(x) x$threads==4 && x$individuals==a$individuals && x$markers==a$markers && x$sets==a$sets && x$adjacent_overlap==a$adjacent_overlap,logical(1)))
 stopifnot(length(j)==1)
 x<-readRDS(file.path(dirname(paths[i]),'scores.rds'));y<-readRDS(file.path(dirname(paths[j]),'scores.rds'))
 stopifnot(identical(dimnames(x),dimnames(y)),isTRUE(all.equal(x,y,tolerance=1e-10)))
 checks[[length(checks)+1]]<-list(case=basename(dirname(paths[i])),max_absolute_thread_difference=max(abs(x-y)),rows=nrow(x),columns=ncol(x))
}
write_json(checks,'work/revision-benchmark/scoring-thread-validation.json',pretty=TRUE,auto_unbox=TRUE)
cat('Validated 20 pairs of fixed-input score outputs across thread limits.\n')
