# Original qgg website examples: models use installed database statistics.
args<-commandArgs(TRUE);stage<-args[1];threads<-as.integer(args[2]);out<-args[3];root<-dirname(out)
writeLines(as.character(Sys.getpid()),file.path(out,'r-process.pid'))
Sys.setenv(OPENBLAS_NUM_THREADS=threads,OMP_NUM_THREADS=threads,MKL_NUM_THREADS=threads,OMP_DYNAMIC='FALSE')
library(qgg);library(gact)
if(requireNamespace('RhpcBLASctl',quietly=TRUE))RhpcBLASctl::blas_set_num_threads(threads)
set.seed(20261005)
writeLines(c(capture.output(sessionInfo()),paste('qgg',packageVersion('qgg'),find.package('qgg')),paste('gact',packageVersion('gact'),find.package('gact'))),file.path(out,'session.txt'))
inputdir<-file.path(root,'website_inputs-1')
if(stage!='website_inputs')input<-readRDS(file.path(inputdir,if(grepl('pops',stage))'pops-input.rds'else'magma-input.rds'))
started<-Sys.time();clock<-proc.time()[3]
value<-tryCatch(switch(stage,
 website_inputs={
  GA<-readRDS('build/benchmarks/manuscript-20261004-143135/database-1/output.rds')
  pstat<-gact::getVEGAS(GAlist=GA,studyID='GWAS1')
  psets<-gact::getFeatureSets(GAlist=GA,feature='Pathways',minsets=50)
  saveRDS(list(stat=pstat,sets=psets),file.path(out,'pops-input.rds'))
  stat<-gact::getMarkerStat(GAlist=GA,studyID='BLR2')
  stopifnot(all(c('rsids','bm','dm','vm')%in%names(stat)))
  msets<-gact::getMarkerSets(GAlist=GA,feature='Genesplus')
  mapped<-qgg::mapSets(sets=msets,rsids=stat$rsids,index=TRUE)
  scores<-setNames(lapply(c('bm','dm','vm'),function(n)sapply(mapped,function(k)sum(abs(stat[[n]][k])))),c('b','d','v'))
  kegg<-gact::getFeatureSets(GAlist=GA,feature='KEGG',minsets=50)
  saveRDS(list(scores=scores,sets=kegg),file.path(out,'magma-input.rds'))
  data.frame(workload=c('PoPS','BLR-derived MAGMA'),observations=c(nrow(as.matrix(pstat)),length(scores$b)),features=c(length(psets),length(kegg)))
 },
 website_pops_rr=qgg::pops(stat=input$stat,sets=input$sets,method='rr',threshold=.05),
 website_pops_c=qgg::pops(stat=input$stat,sets=input$sets,method='bayesC',pi=.001,nit=5000,nburn=1000),
 website_pops_r=qgg::pops(stat=input$stat,sets=input$sets,method='bayesR',pi=.001,nit=5000,nburn=1000),
 {
  if(!grepl('^website_magma_[bdv]_(joint|bayesr)$',stage))stop('Unknown stage')
  tokens<-strsplit(stage,'_')[[1]]
  if(tokens[4]=='joint')qgg::magma(stat=input$scores[[tokens[3]]],sets=input$sets,type='joint')
  else qgg::magma(stat=input$scores[[tokens[3]]],sets=input$sets,method='bayesR')
 }),error=function(e)e)
seconds<-unname(proc.time()[3]-clock);failed<-inherits(value,'error')
serialization<-unname(system.time(saveRDS(value,file.path(out,'output.rds')))[3])
write.csv(data.frame(stage,threads,started=format(started,tz='Europe/Copenhagen'),analysis_seconds=seconds,serialization_seconds=serialization,status=if(failed)'failed'else'completed',warnings=if(failed)conditionMessage(value)else''),file.path(out,'timing.csv'),row.names=FALSE)
writeLines(capture.output(str(value,max.level=2)),file.path(out,'output-structure.txt'))
if(failed)quit(status=1)
