args<-commandArgs(TRUE);requested<-as.integer(args[1]);out<-args[2]
Sys.setenv(OPENBLAS_NUM_THREADS=requested,OMP_NUM_THREADS=requested,MKL_NUM_THREADS=requested,OMP_DYNAMIC='FALSE')
library(qgg)
invisible(eigen(diag(2)))
dyn.load('build/benchmarks/blas-thread-probe/probe_blas_threads.dll')
probe<-function(label){x<-.C('probe_blas_threads',counts=integer(3),parallel=integer(3));data.frame(requested,label,dll=c('Rblas.dll','Rlapack.dll','libopenblas.dll'),configured_threads=x$counts,parallel_build=x$parallel)}
before<-probe('before_RhpcBLASctl')
library(RhpcBLASctl);blas_set_num_threads(requested)
after<-probe('after_RhpcBLASctl')
write.csv(rbind(before,after),file.path(out,paste0('dll-threads-',requested,'.csv')),row.names=FALSE)
print(rbind(before,after),row.names=FALSE)
set.seed(20261005)
measurements<-list()
for(m in c(1000L,3000L)){
 B<-matrix(rnorm(503*m),503,m);R<-cor(B)
 t<-proc.time();for(i in 1:3)ev<-eigen(R,symmetric=TRUE)$values;d<-proc.time()-t
 measurements[[length(measurements)+1]]<-data.frame(requested,markers=m,repetitions=3,wall_seconds=unname(d[3]),cpu_seconds=unname(sum(d[1:2])),average_cores=unname(sum(d[1:2])/d[3]))
}
result<-do.call(rbind,measurements);write.csv(result,file.path(out,paste0('eigen-cpu-',requested,'.csv')),row.names=FALSE);print(result,row.names=FALSE)
