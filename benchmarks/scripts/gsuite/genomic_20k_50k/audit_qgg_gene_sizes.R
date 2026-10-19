sets<-readRDS('build/benchmarks/manuscript-20261004-143135/database-1/database/hsa.0.0.1/gsets/ensg2rsids.rds')
m<-lengths(sets)
d<-data.frame(gene=names(sets),raw_markers=m,matrix_gib=8*as.double(m)^2/1024^3)
d<-d[order(d$raw_markers,decreasing=TRUE),]
out<-'benchmarks/genomic_20k_50k/manuscript-results/20261004-qgg-website'
write.csv(d,file.path(out,'raw-gene-size-audit.csv'),row.names=FALSE)
print(head(d,10),row.names=FALSE)
print(quantile(m,c(0,.25,.5,.75,.95,.99,1)))
cat('Raw memberships before reference/GWAS QC; matrix estimate excludes eigen workspace.\n')
