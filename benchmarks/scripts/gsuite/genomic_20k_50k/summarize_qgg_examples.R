args<-commandArgs(TRUE);runs<-strsplit(args[1],";",fixed=TRUE)[[1]];review<-args[2]
dir.create(review,recursive=TRUE,showWarnings=FALSE)
measurements<-list();inventory<-list()
for(run in runs){
 p<-read.csv(file.path(run,"process-measurements.csv"));p$processWallSeconds<-as.numeric(sub(",",".",p$processWallSeconds,fixed=TRUE))
 for(i in seq_len(nrow(p))){
  timing<-file.path(p$path[i],"timing.csv");if(!file.exists(timing))next
  t<-read.csv(timing)
  measurements[[length(measurements)+1]]<-data.frame(run=basename(run),stage=t$stage,threads=t$threads,analysis_seconds=t$analysis_seconds,serialization_seconds=t$serialization_seconds,process_wall_seconds=p$processWallSeconds[i],peak_working_set_bytes=p$peakWorkingSetBytes[i],sampled_peak_private_bytes=p$sampledPeakPrivateBytes[i],output_bytes=file.info(file.path(p$path[i],"output.rds"))$size,stage_bytes=p$stageDirectoryBytes[i],operation_status=t$status,process_exit=p$exitCode[i],warnings=t$warnings)
 }
 files<-list.files(run,recursive=TRUE,full.names=TRUE);files<-files[!dir.exists(files)]
 inventory[[length(inventory)+1]]<-data.frame(file=files,bytes=file.info(files)$size)
}
m<-do.call(rbind,measurements);write.csv(m,file.path(review,"measurements.csv"),row.names=FALSE)
write.csv(do.call(rbind,inventory),file.path(review,"file-inventory.csv"),row.names=FALSE)
load<-function(s){for(run in rev(runs)){p<-file.path(run,paste0(s,"-1"),"output.rds");if(file.exists(p)){v<-readRDS(p);if(!inherits(v,"error"))return(v)}};NULL}
G<-load("qgg_glist");f<-load("qgg_inputs");v<-load("qgg_vegas");j<-load("qgg_joint");map<-load("qgg_mapping")
checks<-character()
if(!is.null(G)){
 writeLines(capture.output(list(n=G$n,markers=sum(G$mchr),chromosomes=table(unlist(G$chr)))),file.path(review,"reference-dimensions.txt"))
 stopifnot(file.info(G$bedfiles)$size==ceiling(G$n/4)*sum(G$mchr)+3)
 checks<-c(checks,"Filtered BED size agrees with sample and marker counts.")
}
if(!is.null(f)){
 write.csv(data.frame(gene=names(f$marker_sets),raw_markers=lengths(f$marker_sets),fine_mapping=names(f$marker_sets) %in% names(f$fine_sets)),file.path(review,"gene-set-sizes.csv"),row.names=FALSE)
 writeLines(capture.output(list(gwas_markers=sapply(f$vegas_stat,nrow),fine_markers=nrow(f$fine_stat),joint_genes=nrow(f$gene_stat),joint_traits=colnames(f$gene_stat),kegg_sets=length(f$kegg),score_sets=length(f$scoring_sets))),file.path(review,"workload-dimensions.txt"))
}
if(!is.null(v)){
 stopifnot(ncol(v$p)==2,all(is.finite(v$p)),all(v$p>=0 & v$p<=1));write.csv(v$genes,file.path(review,"vegas-gene-dimensions.csv"),row.names=FALSE)
 writeLines(capture.output(list(genes=nrow(v$p),traits=ncol(v$p),zero_pvalues=sum(v$p==0),nonfinite_z=sum(!is.finite(v$z)))),file.path(review,"vegas-output-checks.txt"))
 checks<-c(checks,"VEGAS p-values finite and in [0,1]; tail underflow and Z non-finiteness counted separately.")
}
if(!is.null(j))writeLines(capture.output(str(j,max.level=2)),file.path(review,"joint-output-structure.txt"))
if(!is.null(map)){
 write.csv(map$conv,file.path(review,"mapping-geweke-diagnostics.csv"));write.csv(map$post,file.path(review,"mapping-region-results.csv"))
 checks<-c(checks,paste("Mapping",nrow(map$post),"regions; installed checkConvergence default FALSE, no convergence claim."))
}
for(s in c("qgg_bayesc_fit","qgg_bayesr_fit","qgg_total_c","qgg_total_r","qgg_pathway")){
 x<-load(s);if(is.null(x))next
 writeLines(capture.output(str(x,max.level=2)),file.path(review,paste0(s,"-output.txt")))
 if(is.matrix(x)){stopifnot(nrow(x)==G$n,all(is.finite(x)));checks<-c(checks,paste(s,nrow(x),"people",ncol(x),"finite score columns"))}
}
writeLines(checks,file.path(review,"output-checks.txt"))
file.copy(file.path(tail(runs,1),"hardware.json"),file.path(review,"hardware.json"),overwrite=TRUE)
for(pkg in c("qgg","gact")){
 dll<-list.files(file.path(find.package(pkg),"libs"),pattern="[.]dll$",recursive=TRUE,full.names=TRUE)
 writeLines(c(paste(pkg,packageVersion(pkg),find.package(pkg)),capture.output(packageDescription(pkg)),paste(dll,tools::md5sum(dll))),file.path(review,paste0(pkg,"-provenance.txt")))
}
for(fn in c("vegas","magma","gmap","gbayes","getSparseLD","gscore","gprep","getG")){
 writeLines(deparse(getExportedValue("qgg",fn)),file.path(review,paste0("installed-qgg-",fn,".R")))
}
sources<-c("benchmarks/genomic_20k_50k/qgg_examples_stage.R","benchmarks/genomic_20k_50k/qgg_ld_guard.R","benchmarks/genomic_20k_50k/run_manuscript.ps1","benchmarks/genomic_20k_50k/qgg-examples-protocol.md")
write.csv(data.frame(file=sources,md5=unname(tools::md5sum(sources))),file.path(review,"harness-source-hashes.csv"),row.names=FALSE)
print(m[,c("stage","analysis_seconds","peak_working_set_bytes","process_exit")])
