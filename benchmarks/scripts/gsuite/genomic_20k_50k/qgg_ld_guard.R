run_website_ld<-function(G,out){
  required<-sum(lengths(G$rsids))*2001*4
  free<-as.numeric(system2("powershell.exe",c("-NoProfile","-Command",shQuote("(Get-PSDrive -Name C).Free")),stdout=TRUE))
  write.csv(data.frame(markers=sum(lengths(G$rsids)),ld_window=1000,required_ld_bytes=required,available_disk_bytes=free),file.path(out,"storage-check.csv"),row.names=FALSE)
  if(!is.finite(free)||required>free-1024^3)stop(sprintf("Website LD needs %.2f GiB; %.2f GiB available. Existing results preserved; no smaller workload substituted.",required/1024^3,free/1024^3))
  qgg::gprep(Glist=G,task="sparseld",msize=1000,rsids=unlist(G$rsids),ldfiles=file.path(out,"g1000_eur_filtered.ld"),overwrite=FALSE,ncores=1)
}
