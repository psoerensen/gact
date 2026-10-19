args<-commandArgs(TRUE); run<-args[1]; review<-args[2]
dir.create(review,recursive=TRUE,showWarnings=FALSE)
.libPaths(c("build/task-packages/library","build/task-packages/ecosystem/library","build/r-library",.libPaths()))
p<-read.csv(file.path(run,"process-measurements.csv")); p$processWallSeconds<-as.numeric(sub(",",".",p$processWallSeconds,fixed=TRUE)); records<-list();checks<-list()
if(identical(p$stage,"database")) {
 path<-p$path[1]; t<-read.csv(file.path(path,"timing.csv")); x<-readRDS(file.path(path,"output.rds"))
 downloads<-read.csv(file.path(path,"downloads.csv")); files<-read.csv(file.path(path,"database-files.csv"))
 stopifnot(all(downloads$status=="success"),all(c("annotation","marker","gsea","gstat","gbayes","atc","dgi") %in% downloads$what),nrow(x$studies)>0,all(file.exists(x$studyfiles)),length(x$rsids)>0,!is.null(x$map),all(file.exists(files$file)),all(files$bytes>0))
 write.csv(data.frame(stage="database",started=t$started,download_install_seconds=t$analysis_seconds,descriptor_save_seconds=t$serialization_seconds,process_wall_seconds=p$processWallSeconds,peak_working_set_bytes=p$peakWorkingSetBytes,sampled_peak_private_bytes=p$sampledPeakPrivateBytes,database_file_bytes=sum(files$bytes),retained_zip_bytes=sum(files$bytes[grepl("[.]zip$",files$file)]),descriptor_bytes=file.info(file.path(path,"output.rds"))$size,stage_directory_bytes=p$stageDirectoryBytes,status=t$status,warnings=t$warnings),file.path(review,"database-measurements.csv"),row.names=FALSE)
 write.csv(downloads,file.path(review,"database-downloads.csv"),row.names=FALSE)
 write.csv(files,file.path(review,"database-files.csv"),row.names=FALSE)
 write.csv(x$studies,file.path(review,"database-studies.csv"),row.names=FALSE)
 writeLines("All seven component categories successful; annotation map and marker identities present; study files and nonempty installed files verified.",file.path(review,"database-checks.txt"))
 print(x$studies[1:2,c("id","trait","reference")]);cat("Database bytes",sum(files$bytes),"; peak GiB",p$peakWorkingSetBytes/1024^3,"\n")
 quit(status=0)
}
for(i in seq_len(nrow(p))) {
 path<-p$path[i]; t<-read.csv(file.path(path,"timing.csv"));x<-readRDS(file.path(path,"output.rds"))
 actual<-switch(t$stage,vegas=x$execution$preparation_threads,joint=x$diagnostics$joint$threads_used,mapping=1,total=x$control$threads,pathway=x$control$threads)
 records[[i]]<-data.frame(stage=t$stage,started=t$started,requested_threads=t$threads,operation_thread_setting=actual,thread_evidence=if(t$stage=="vegas") "reported preparation threads" else if(t$stage=="joint") "reported sampler threads" else if(t$stage=="mapping") "serial R region/trait loop; no thread control" else "native control setting; runtime count not reported",analysis_seconds=t$analysis_seconds,serialization_seconds=t$serialization_seconds,process_wall_seconds=p$processWallSeconds[i],peak_working_set_bytes=p$peakWorkingSetBytes[i],sampled_peak_private_bytes=p$sampledPeakPrivateBytes[i],stage_directory_bytes=p$stageDirectoryBytes[i],output_rds_bytes=file.info(file.path(path,"output.rds"))$size,status=t$status,warnings=t$warnings)
}
for(s in unique(p$stage)) {
 a<-readRDS(file.path(run,paste0(s,"-1"),"output.rds"));b<-readRDS(file.path(run,paste0(s,"-4"),"output.rds"))
 keys<-switch(s,vegas=c("stat"),joint=c("pathways","joint"),mapping=c("estimates","parameters","sets"),total="scores",pathway="scores")
 stopifnot(all(keys %in% names(a)),all(keys %in% names(b)))
 compare<-all.equal(a[keys],b[keys],tolerance=1e-10)
 checks[[s]]<-data.frame(stage=s,check="matched substantive outputs at relative tolerance 1e-10",passed=isTRUE(compare),detail=paste(compare,collapse="; "))
 if(s %in% c("total","pathway")) stopifnot(nrow(a$scores)==2000,all(is.finite(a$scores)),all(is.finite(b$scores)))
 if(s=="vegas") stopifnot(length(a$stat)==2,all(vapply(a$stat,function(x) nrow(x)==1000 && all(is.finite(x$score)) && all(is.finite(x$p_value)),logical(1))))
 if(s=="mapping") stopifnot(nrow(a$estimates)==4000,nrow(a$diagnostic_summary)==40,all(a$diagnostic_summary$converged),all(is.finite(a$estimates$mean)),all(is.finite(a$estimates$pip)))
 if(s=="joint") stopifnot(nrow(a$pathways)==6,all(is.finite(a$pathways$pip)))
}
write.csv(do.call(rbind,records),file.path(review,"measurements.csv"),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(review,"output-checks.csv"),row.names=FALSE)
files<-list.files("build/benchmarks/genomic-20k-50k",full.names=TRUE)
files<-files[basename(files) %in% c("data_simulation.rds","genotypes.bed","genotypes.bim","genotypes.fam","sparse_ld_and_scores.rds","gwas_and_checkstat.rds","vegas.rds","bayesian_regression.rds") | grepl("^LD[.]",basename(files))]
write.csv(data.frame(file=normalizePath(files,winslash="/"),bytes=file.info(files)$size,md5=unname(tools::md5sum(files))),file.path(review,"input-manifest.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(review,"session.txt"))
for(pkg in c("gbase","glma","gsea","gmap","gscore","gbayes")) {
 ns<-asNamespace(pkg);fns<-grep("native_info$",getNamespaceExports(pkg),value=TRUE)
 details<-c(paste("package",pkg,"version",packageVersion(pkg),"library",find.package(pkg)),capture.output(packageDescription(pkg)),unlist(lapply(fns,function(f) capture.output(get(f,ns)()))))
 dll<-list.files(file.path(find.package(pkg),"libs"),pattern="[.]dll$",recursive=TRUE,full.names=TRUE)
 details<-c(details,paste(dll,tools::md5sum(dll)))
 writeLines(details,file.path(review,paste0(pkg,"-provenance.txt")))
}
file.copy(file.path(run,"hardware.json"),file.path(review,"hardware.json"),overwrite=FALSE)
print(do.call(rbind,records));print(do.call(rbind,checks))






