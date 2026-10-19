# Untimed LD audit; run after time_ld_swap.ps1 has completed.
args <- commandArgs(TRUE); stopifnot(length(args)==1L,dir.exists(args[1]))
.libPaths(c('build/task-packages/library','build/task-packages/ecosystem/library',.libPaths()))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,OMP_DYNAMIC='FALSE')
if(requireNamespace('RhpcBLASctl',quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1L)
input <- readRDS('build/benchmarks/manuscript-20261005-142319/website_sim_stats-1/output.rds')
G <- readRDS('build/benchmarks/manuscript-20261005-231349/cmp_suite_prepare-1/output.rds')
paths <- unique(unlist(c(input$Glist[c('bedfiles','bimfiles','famfiles')],
  G[c('bedfiles','bimfiles','famfiles')]),use.names=FALSE))
write.csv(data.frame(path=paths,bytes=file.info(paths)$size,
  md5=unname(tools::md5sum(paths))),file.path(args[1],'genotype-input-hashes.csv'),row.names=FALSE)
rows <- lapply(names(input$sets),function(region) {
  markers <- input$stat$rsids[input$stat$rsids %in% input$sets[[region]]]
  LD <- gbase::getLD(G,markers,nthreads=1L)
  r2 <- LD[upper.tri(LD)]^2
  data.frame(region,markers=length(markers),maximum_off_diagonal_r2=max(r2),
    pairs_at_r2_0_8=sum(r2>=.8))
})
rows <- do.call(rbind,rows)
write.csv(rows,file.path(args[1],'ld-audit-regions.csv'),row.names=FALSE)
# The preserved genome LD is float32 upper-triangle CSR correlation with
# implicit unit diagonal. Stream its off-diagonal values once.
# No LD preparation or genotype writes take place.
path <- 'build/benchmarks/manuscript-20261005-231349/cmp_suite_ld-1/human.values.f32.bin'
metadata <- readLines(sub('values.f32.bin$','meta.txt',path))
stopifnot('triangle=upper' %in% metadata,'diagonal=implicit_1' %in% metadata,
  'values_type=float32' %in% metadata,'value=r' %in% metadata)
con <- file(path,'rb');entries <- strong <- 0
maximum_r2 <- 0
repeat {
  value <- readBin(con,'numeric',n=1000000L,size=4,endian='little')
  if(!length(value)) break
  stopifnot(all(is.finite(value)))
  entries <- entries+length(value);strong <- strong+sum(value^2>=.8)
  maximum_r2 <- max(maximum_r2,value^2)
}
close(con)
genome <- data.frame(markers=nrow(input$stat),stored_csr_entries=entries,
  strong_off_diagonal_entries=strong,maximum_off_diagonal_r2=maximum_r2)
write.csv(genome,file.path(args[1],'ld-audit-genome.csv'),row.names=FALSE)
print(rows,row.names=FALSE);print(genome,row.names=FALSE)
