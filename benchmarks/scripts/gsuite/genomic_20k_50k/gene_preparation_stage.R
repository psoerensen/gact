# Real gene preparation, preserving the website's eigen(cor(B)) computation.
args<-commandArgs(TRUE);stage<-args[1];threads<-as.integer(args[2]);out<-args[3]
target<-as.integer(sub('gene_preparation_','',stage));stopifnot(target%in%c(100,300,500,1000),threads>=1)
writeLines(as.character(Sys.getpid()),file.path(out,'r-process.pid'))
Sys.setenv(OPENBLAS_NUM_THREADS=threads,OMP_NUM_THREADS=threads,MKL_NUM_THREADS=threads,OMP_DYNAMIC='FALSE')
library(qgg)
if(requireNamespace('RhpcBLASctl',quietly=TRUE))RhpcBLASctl::blas_set_num_threads(threads)
writeLines(c(capture.output(sessionInfo()),paste('qgg',packageVersion('qgg'),find.package('qgg'))),file.path(out,'session.txt'))
cat('START input loading; target',target,'markers, 100 genes\n');flush.console()
f<-readRDS('build/benchmarks/manuscript-20261004-174351/qgg_inputs-1/output.rds')
old<-readRDS('build/benchmarks/manuscript-20261004-194031/qgg_vegas-4/output.rds')$genes
selection<-head(order(abs(old$m-target),old$EnsemblID),100)
genes<-old$EnsemblID[selection];expected<-old$m[selection]
set.seed(20261005);order_genes<-sample(seq_along(genes));genes<-genes[order_genes];expected<-expected[order_genes]
alignment_start<-proc.time()
shared<-Reduce(intersect,lapply(f$vegas_stat,function(s)s$rsids))
sets<-f$marker_sets[genes];flat<-unlist(sets,use.names=FALSE)
present<-!is.na(match(flat,shared));group<-rep(seq_along(sets),lengths(sets))
indices<-match(flat[present],f$Glist$rsids[[1]])
stopifnot(!anyNA(indices));indices<-split(indices,group[present])
stopifnot(identical(as.integer(lengths(indices)),as.integer(expected)))
alignment<-proc.time()-alignment_start
write.csv(data.frame(wall_seconds=alignment[3],cpu_seconds=sum(alignment[1:2])),file.path(out,'alignment-timing.csv'),row.names=FALSE)
G<-f$Glist;rm(f,flat,shared,present,group);gc()
rows<-vector('list',length(genes));spectra<-setNames(vector('list',length(genes)),genes)
measure<-function(fun){a<-proc.time();v<-fun();d<-proc.time()-a;list(value=v,wall=unname(d[3]),cpu=unname(sum(d[1:2])))}
started<-Sys.time();clock<-proc.time();error<-NULL
tryCatch({
for(i in seq_along(genes)){
 k<-indices[[i]];stopifnot(length(unique(G$chr[[1]][k]))==1L)
 combined<-proc.time()
 extraction<-measure(function()qgg::getG(Glist=G,chr=1,cls=k,scale=TRUE))
 B<-extraction$value;extraction$value<-NULL
 correlation<-measure(function()cor(B));R<-correlation$value;correlation$value<-NULL
 decomposition<-measure(function()eigen(R,symmetric=TRUE));ev<-decomposition$value$values
 stopifnot(length(ev)==length(k),all(is.finite(ev)),abs(sum(ev)-length(k))<1e-6*length(k))
 spectra[[i]]<-ev
 elapsed<-proc.time()-combined
 rows[[i]]<-data.frame(gene=genes[i],target_markers=target,markers=length(k),individuals=G$n,
 extraction_seconds=extraction$wall,ld_seconds=correlation$wall,eigen_seconds=decomposition$wall,
 combined_seconds=unname(elapsed[3]),extraction_cpu_seconds=extraction$cpu,
 ld_cpu_seconds=correlation$cpu,eigen_cpu_seconds=decomposition$cpu,
 combined_cpu_seconds=unname(sum(elapsed[1:2])),genotype_matrix_bytes=as.numeric(object.size(B)),
 ld_matrix_bytes=as.numeric(object.size(R)),eigensystem_bytes=as.numeric(object.size(decomposition$value)))
 rm(B,R,decomposition,extraction,correlation,ev)
 if(i%%10==0){write.csv(do.call(rbind,rows[seq_len(i)]),file.path(out,'gene-progress.csv'),row.names=FALSE);cat('Completed',i,'of',length(genes),'genes; target',target,'\n');flush.console()}
}
},error=function(e)error<<-e)
duration<-proc.time()-clock
value<-if(is.null(error))list(gene_timings=do.call(rbind,rows),eigenvalues=spectra,
 eigen_policy='Full eigenvectors computed to match qgg website; only eigenvalues retained',
 marker_policy='Full real gene sets after shared CAD/T2D QC; nearest-size groups, no trimming',
 alignment_seconds=unname(alignment[3]))else error
serialization<-unname(system.time(saveRDS(value,file.path(out,'output.rds')))[3])
write.csv(data.frame(stage,threads,started=format(started,tz='Europe/Copenhagen'),analysis_seconds=unname(duration[3]),
 serialization_seconds=serialization,status=if(is.null(error))'completed'else'failed',
 warnings=if(is.null(error))''else conditionMessage(error)),file.path(out,'timing.csv'),row.names=FALSE)
if(!is.null(error))quit(status=1)
