# Reuse the existing gact database, BED reference and stored LD scores.
# Requires compatible development installations of gbase and gcorr.
prepare_gcorr_cad_t2d <- function(GAlist, Glist, partition_file,
                                study_ids=c('GWAS1','GWAS2'),
                                common_study_ids=c('GWAS1','GWAS2','GWAS6')) {
  stopifnot(length(study_ids)>=2L,all(study_ids %in% common_study_ids))
  # Keep the three-study common mask used by the existing gact LDSC example.
  # Only study_ids are fitted below; GWAS6 contributes to this mask.
  summary <- gact::getMarkerStat(GAlist=GAlist,studyID=common_study_ids)
  common_ids <- intersect(rownames(summary$z),rownames(summary$n))
  summary$z <- summary$z[common_ids,,drop=FALSE]
  summary$n <- summary$n[common_ids,,drop=FALSE]
  summary <- summary[c('z','n')]
  if(!is.null(Glist$ldscores)) {
    score_ids <- unlist(if(is.null(Glist$rsidsLD)) Glist$rsids else Glist$rsidsLD,
      use.names=FALSE)
    scores <- unlist(Glist$ldscores,use.names=FALSE)
    stopifnot(length(scores)==length(score_ids))
    names(scores) <- score_ids
  } else {
    scores <- gact::getLDscoresDB(GAlist=GAlist,ancestry='EUR',version='1000G')
  }
  # Convert an existing qgg descriptor using the same physical reference files.
  if (!inherits(Glist,'gs_genotypes')) {
    Glist <- gbase::gprep(bedfiles=Glist$bedfiles,bimfiles=Glist$bimfiles,
      famfiles=Glist$famfiles)
  }
  ids <- unlist(Glist$rsids,use.names=FALSE)
  chr <- as.character(unlist(Glist$chr,use.names=FALSE))
  pos <- unlist(Glist$pos,use.names=FALSE)
  ea <- unlist(Glist$a1,use.names=FALSE);nea <- unlist(Glist$a2,use.names=FALSE)
  at <- match(ids,rownames(summary$z))
  # Glist comes from the website's BED-to-GAlist marker selection. Retain all
  # available reference/summary/score rows; do not add another QC mask.
  retain <- !is.na(at) & ids %in% names(scores)
  index <- which(retain)
  markers <- data.frame(rsids=ids[index],chr=chr[index],pos=pos[index],
    ea=ea[index],nea=nea[index],stringsAsFactors=FALSE)
  labels <- data.table::fread(file.path(GAlist$dirs[['marker']],'markers.txt.gz'),
    select=c('rsids','chr','pos','ea','nea'),data.table=FALSE)
  mapped <- match(markers$rsids,labels$rsids)
  stopifnot(!anyNA(mapped),all(as.character(labels$chr[mapped])==markers$chr),
    all(labels$pos[mapped]==markers$pos))
  direct <- labels$ea[mapped]==markers$ea & labels$nea[mapped]==markers$nea
  swap <- labels$ea[mapped]==markers$nea & labels$nea[mapped]==markers$ea
  stopifnot(all(direct | swap))
  z <- summary$z[at[index],study_ids,drop=FALSE] * ifelse(direct,1,-1)
  n <- summary$n[at[index],study_ids,drop=FALSE]
  stopifnot(all(is.finite(z)),all(is.finite(n)),all(n>0),
    identical(rownames(n),markers$rsids))
  # Preserve the extracted field, including any marker-specific sample sizes.
  # GNOVA and regional HESS validate their constant-N requirement themselves.
  sizes <- vapply(seq_along(study_ids),function(j) {
    if(all(n[,j]==n[1,j]))n[1,j] else NA_real_
  },numeric(1))
  stat <- setNames(lapply(seq_along(study_ids),function(j) {
    x <- markers;x$z <- z[,j];x$n <- n[,j];x
  }),study_ids)
  scores <- scores[markers$rsids]
  stopifnot(all(is.finite(scores)),all(scores>0))
  blocks <- setNames(paste0('delete',pmin(200L,
    ceiling(seq_len(nrow(markers))*200/nrow(markers)))),markers$rsids)
  partition <- read.table(partition_file,header=TRUE)
  partition$chr <- sub('chr','',partition$chr,fixed=TRUE)
  region <- rep(NA_character_,nrow(markers))
  for(c in unique(markers$chr)) {
    rows <- which(markers$chr==c);p <- partition[partition$chr==c,,drop=FALSE]
    at <- findInterval(markers$pos[rows],p$start)
    good <- at>0;good[good] <- markers$pos[rows][good]<p$stop[at[good]]
    region[rows[good]] <- paste0('chr',c,'_',p$start[at[good]],'_',p$stop[at[good]])
  }
  covered <- !is.na(region)
  reference_count <- if(is.null(Glist$physical_marker_counts))
    sum(lengths(Glist$rsids)) else sum(Glist$physical_marker_counts)
  list(statistics=stat,ldscores=scores,blocks=blocks,Glist=Glist,
    regions=split(markers$rsids[covered],region[covered]),
    reference_marker_count=nrow(markers),physical_reference_marker_count=reference_count,
    uncovered_markers=sum(!covered),sample_sizes=setNames(sizes,study_ids))
}
