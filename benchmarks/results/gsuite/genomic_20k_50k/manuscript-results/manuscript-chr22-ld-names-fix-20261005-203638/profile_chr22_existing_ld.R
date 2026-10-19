.libPaths(c('build/task-packages/library','build/task-packages/ecosystem/library','build/r-library',.libPaths()))
library(gbase)
out<-'build/benchmarks/manuscript-chr22-filtered-ld-20261005-202118'
args<-commandArgs(TRUE)
destination<-if(length(args))args[1]else out
x<-readRDS(file.path(out,'ld_gbase_filtered-4/output.rds'))
prefix<-x$sparseLD$prefix
measure<-function(label,expr){
  timing<-system.time(force(expr))
  row<-data.frame(component=label,user_seconds=unname(timing[1]),system_seconds=unname(timing[2]),elapsed_seconds=unname(timing[3]))
  print(row);flush.console();row
}
rows<-list(
 measure('identity open',gbase:::cpp_gbase_sparse_ld_identity(prefix)),
 measure('full payload inspection',gbase:::cpp_gbase_sparse_ld_open(prefix)),
 measure('LD scores including validated traversal',gbase:::cpp_gbase_ld_scores(prefix,FALSE,NULL,character(),character())),
 measure('MD5 file checksums',gbase:::.gprep_checksums(x$sparseLD$paths)))
write.csv(do.call(rbind,rows),file.path(destination,'existing-resource-component-timings.csv'),row.names=FALSE)
