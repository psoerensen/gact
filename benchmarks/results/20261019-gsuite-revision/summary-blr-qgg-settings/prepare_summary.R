.libPaths(c('C:/Users/au223366/Documents/GitHub/gsuite/build/source-repository/library','C:/Users/au223366/Documents/GitHub/gsuite/build/r-library','C:/Users/au223366/AppData/Local/R/win-library/4.4',.libPaths()))
library(Matrix);library(gbase);library(jsonlite)
out<-'work/summary-blr-benchmark/input';stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
t0<-proc.time()[3];f<-readRDS('work/revision-benchmark/input-preparation/output.rds');loading<-unname(proc.time()[3]-t0)
t0<-proc.time()[3];genes<-rownames(f$stat);sets<-gbase::mapSets(lapply(f$sets,unique),rsids=genes);sets<-sets[lengths(sets)>0]
X<-sparseMatrix(i=unlist(sets,use.names=FALSE),j=rep(seq_along(sets),lengths(sets)),x=1,dims=c(length(genes),length(sets)),dimnames=list(genes,names(sets)))
y<-as.numeric(scale(f$stat[,1],center=TRUE,scale=TRUE));matrix_seconds<-unname(proc.time()[3]-t0)
t0<-proc.time()[3];n<-nrow(X);m<-ncol(X);mu<-Matrix::colMeans(X);sigma<-sqrt((Matrix::colSums(X^2)-n*mu^2)/(n-1))
stopifnot(all(sigma>0));XX<-as.matrix(Matrix::crossprod(X));XX<-(XX-n*tcrossprod(mu))/outer(sigma,sigma)
Xy<-(as.numeric(Matrix::crossprod(X,y))-mu*sum(y))/sigma;yy<-sum(y^2)
preparation<-unname(proc.time()[3]-t0)
check0<-proc.time()[3];q<-get('computeStat',asNamespace('qgg'))(X=X,y=y,scale=TRUE)
stopifnot(isTRUE(all.equal(unname(XX),unname(q$XX),tolerance=1e-10)),isTRUE(all.equal(unname(Xy),unname(q$Xy),tolerance=1e-10)),isTRUE(all.equal(yy,q$yy)),identical(as.integer(n),as.integer(q$n)))
validation<-unname(proc.time()[3]-check0)
pi_start<-.001;vg<-.5;ve<-.5;vb<-vg/(m*pi_start);ssb<-(4-2)/4*vb;sse<-(4-2)/4*ve
set.seed(10);seed<-sample.int(.Machine$integer.max,1)
writeBin(as.double(XX),file.path(out,'XX.bin'),size=8,endian='little');writeBin(as.double(Xy),file.path(out,'Xy.bin'),size=8,endian='little')
writeLines(names(sets),file.path(out,'features.txt'));writeLines(paste(c(as.character(n),as.character(m),format(c(yy,ve,vb,ssb,sse),digits=17,scientific=FALSE,trim=TRUE),as.character(seed)),collapse=' '),file.path(out,'meta.txt'))
write_json(list(input_loading_seconds=loading,feature_matrix_seconds=matrix_seconds,cross_product_preparation_seconds=preparation,qgg_validation_seconds_excluded=validation,
 genes=n,features=m,feature_nonzeros=length(X@x),max_XX_difference=max(abs(XX-q$XX)),max_Xy_difference=max(abs(Xy-q$Xy)),yy=yy,n=n,
 qgg_defaults=list(method='bayesR',nit=5000,nburn=1000,nthin=1,chains=1,pi_argument=pi_start,actual_BayesR_weights=c(.95,.02,.02,.01),gamma=c(0,.01,.1,1),initial_vg=vg,initial_ve=ve,initial_vb=vb,nub=4,nue=4,nug=4,ssb_prior=ssb,sse_prior=sse,updateB=TRUE,updateE=TRUE,updatePi=TRUE,updateG=TRUE,adjustE=TRUE),seed=seed,
 matching_scope='Same summary scale, one chain, sweeps, initial parameters, variance/weight priors and learning flags. With dense XX every LD neighborhood includes all features, so qgg adjustE contributes zero and updateG cannot affect coefficient sampling. Native marker visit order/RNG consumption and retained output differ; not bitwise qgg trajectory equivalence.'),file.path(out,'preparation.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
writeLines(capture.output(sessionInfo()),file.path(out,'session.txt'))
cat('Sparse cross-product preparation:',preparation,'seconds. Verified against installed qgg computeStat.\n')
