args<-commandArgs(TRUE);run<-args[1];reference_run<-args[2];review<-args[3]
cohort_run<-args[4]
fixture_run<-args[5]
evidence_run<-args[6]
mapping_run<-args[7]
dir.create(review,recursive=TRUE,showWarnings=FALSE)
.libPaths(c("build/task-packages/library","build/task-packages/ecosystem/library","build/r-library",.libPaths()))
p<-read.csv(file.path(run,"process-measurements.csv"));p$processWallSeconds<-as.numeric(sub(",",".",p$processWallSeconds,fixed=TRUE))
q<-read.csv(file.path(reference_run,"process-measurements.csv"));q<-q[q$stage=="reference_import",];q$processWallSeconds<-as.numeric(sub(",",".",q$processWallSeconds,fixed=TRUE))
c<-read.csv(file.path(cohort_run,"process-measurements.csv"));c<-c[c$stage %in% c("founders","bed_conversion"),];c$processWallSeconds<-as.numeric(sub(",",".",c$processWallSeconds,fixed=TRUE))
d<-read.csv(file.path(fixture_run,"process-measurements.csv"));d<-d[d$stage=="preparation",];d$processWallSeconds<-as.numeric(sub(",",".",d$processWallSeconds,fixed=TRUE))
e<-read.csv(file.path(evidence_run,"process-measurements.csv"));e<-e[e$stage %in% c("vegas_real","joint_real"),];e$processWallSeconds<-as.numeric(sub(",",".",e$processWallSeconds,fixed=TRUE))
h<-read.csv(file.path(mapping_run,"process-measurements.csv"));h<-h[h$stage=="mapping_real",];h$processWallSeconds<-as.numeric(sub(",",".",h$processWallSeconds,fixed=TRUE));p<-rbind(q,c,d,e,h,p)
stopifnot(all(p$exitCode==0))
rows<-lapply(seq_len(nrow(p)),function(i){t<-read.csv(file.path(p$path[i],"timing.csv"));data.frame(stage=t$stage,primary=t$stage %in% c("vegas_real","joint_real","mapping_real","total_real","pathway_real"),threads=1,started=t$started,analysis_seconds=t$analysis_seconds,serialization_seconds=t$serialization_seconds,process_wall_seconds=p$processWallSeconds[i],peak_working_set_bytes=p$peakWorkingSetBytes[i],sampled_peak_private_bytes=p$sampledPeakPrivateBytes[i],stage_directory_bytes=p$stageDirectoryBytes[i],output_rds_bytes=file.info(file.path(p$path[i],"output.rds"))$size,status=t$status,warnings=t$warnings)})
measurements<-do.call(rbind,rows);write.csv(measurements,file.path(review,"measurements.csv"),row.names=FALSE)
load<-function(s)readRDS(file.path(if(s=="preparation")fixture_run else if(s %in% c("vegas_real","joint_real"))evidence_run else if(s=="mapping_real")mapping_run else run,paste0(s,"-1"),"output.rds"))
f<-load("preparation");v<-load("vegas_real");j<-load("joint_real");m<-load("mapping_real");total<-load("total_real");pathway<-load("pathway_real")
stopifnot(length(f$individuals)==20000,length(f$genes)==1000,all(lengths(f$genes)==274),length(f$stat)==2,
 all(vapply(v$stat,function(x)nrow(x)==1000 && all(is.finite(x$score)) && all(is.finite(x$p_value)),logical(1))),
 nrow(j$pathways)==6,all(is.finite(j$pathways$pip)),nrow(m$estimates)==4000,all(is.finite(m$estimates$mean)),all(is.finite(m$estimates$pip)),
 nrow(total$scores)==2000,ncol(total$scores)==1,all(is.finite(total$scores)),nrow(pathway$scores)==2000,ncol(pathway$scores)==3,all(is.finite(pathway$scores)))
write.csv(m$diagnostic_summary,file.path(review,"mapping-diagnostics.csv"),row.names=FALSE)
write.csv(m$eigen_attempts,file.path(review,"mapping-eigen-attempts.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(names(v$stat),function(trait)data.frame(trait=trait,v$stat[[trait]][,c("set","effective_rank","clipped_eigenvalues","tail_method","p_available","moment_fallback")]))),file.path(review,"vegas-numerical-diagnostics.csv"),row.names=FALSE)
writeLines(load("weights_real")$provenance,file.path(review,"scoring-weight-provenance.txt"))
writeLines(capture.output(str(f$provenance)),file.path(review,"workload.txt"))
write.csv(data.frame(gene=names(f$genes),markers=lengths(f$genes),first_marker=vapply(f$genes,head,character(1),1),last_marker=vapply(f$genes,tail,character(1),1)),file.path(review,"gene-set-manifest.csv"),row.names=FALSE)
writeLines("Dimensions and finite VEGAS, pathway, mapping and score outputs verified. Mapping convergence statuses are reported separately; no calibration or thread-comparison run was performed.",file.path(review,"output-checks.txt"))
input_files<-c("../qgdata/cache/1000g_phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz","../qgdata/cache/1000g_phase3/genetic_map_GRCh37_chr22.txt.gz","../qgdata/cache/1000g_phase3/integrated_call_samples_v3.20130502.ALL.panel")
write.csv(data.frame(file=normalizePath(input_files,winslash="/"),bytes=file.info(input_files)$size,md5=unname(tools::md5sum(input_files))),file.path(review,"input-manifest.csv"),row.names=FALSE)
files<-c(list.files(run,recursive=TRUE,full.names=TRUE),list.files(file.path(reference_run,"reference_import-1"),recursive=TRUE,full.names=TRUE),list.files(file.path(cohort_run,"founders-1"),recursive=TRUE,full.names=TRUE),list.files(file.path(cohort_run,"bed_conversion-1"),recursive=TRUE,full.names=TRUE),list.files(file.path(fixture_run,"preparation-1"),recursive=TRUE,full.names=TRUE));files<-files[file.exists(files)&!dir.exists(files)]
files<-c(files,list.files(file.path(evidence_run,"vegas_real-1"),recursive=TRUE,full.names=TRUE),list.files(file.path(evidence_run,"joint_real-1"),recursive=TRUE,full.names=TRUE))
files<-c(files,list.files(file.path(mapping_run,"mapping_real-1"),recursive=TRUE,full.names=TRUE))
write.csv(data.frame(file=files,bytes=file.info(files)$size),file.path(review,"file-inventory.csv"),row.names=FALSE)
file.copy(file.path(run,"hardware.json"),file.path(review,"hardware.json"))
for(pkg in c("gsim","gbase","glma","gsea","gmap","gscore","gbayes")){
 ns<-asNamespace(pkg);info<-grep("native_info$",getNamespaceExports(pkg),value=TRUE)
 dll<-list.files(file.path(find.package(pkg),"libs"),pattern="[.]dll$",recursive=TRUE,full.names=TRUE)
 writeLines(c(paste(pkg,packageVersion(pkg),find.package(pkg)),capture.output(packageDescription(pkg)),unlist(lapply(info,function(fn)capture.output(get(fn,ns)()))),paste(dll,tools::md5sum(dll))),file.path(review,paste0(pkg,"-provenance.txt")))
}
print(measurements[,c("stage","analysis_seconds","peak_working_set_bytes")]);print(f$provenance);print(table(m$diagnostic_summary$status))
