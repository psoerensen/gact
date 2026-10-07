# Prepare disjoint marker-set LDSC inputs from existing ordinary LD scores.
# For whole chromosomes, this is the within-chromosome LD approximation.
# For arbitrary sets with LD across their boundaries, supply true annotation
# LD scores instead of using this masked-score construction.
prepare_gcorr_set_scores <- function(marker_ids, ldscores, sets,
                                    marker_counts=setNames(lengths(sets),names(sets))) {
  stopifnot(is.character(marker_ids),length(marker_ids)>0L,!anyNA(marker_ids),
    !anyDuplicated(marker_ids),is.numeric(ldscores),!is.null(names(ldscores)),
    !anyDuplicated(names(ldscores)),all(marker_ids %in% names(ldscores)),
    is.list(sets),length(sets)>0L,!is.null(names(sets)),
    !anyNA(names(sets)),all(nzchar(names(sets))),!anyDuplicated(names(sets)),
    all(vapply(sets,function(x)is.character(x)&&length(x)>0L&&!anyNA(x)&&
      all(nzchar(x))&&!anyDuplicated(x),logical(1))))
  members <- unlist(sets,use.names=FALSE)
  stopifnot(!anyDuplicated(members),setequal(members,marker_ids),
    is.numeric(marker_counts),!is.null(names(marker_counts)),
    setequal(names(marker_counts),names(sets)),!anyDuplicated(names(marker_counts)))
  marker_counts <- marker_counts[names(sets)]
  stopifnot(all(is.finite(marker_counts)),all(marker_counts>=lengths(sets)),
    all(marker_counts==floor(marker_counts)))
  ordinary <- ldscores[marker_ids]
  stopifnot(all(is.finite(ordinary)),all(ordinary>0))
  scores <- matrix(0,nrow=length(marker_ids),ncol=length(sets),
    dimnames=list(marker_ids,names(sets)))
  at <- match(members,marker_ids)
  group <- rep.int(seq_along(sets),lengths(sets))
  scores[cbind(at,group)] <- unname(ordinary[at])
  overlap <- diag(as.numeric(marker_counts),nrow=length(sets))
  dimnames(overlap) <- list(names(sets),names(sets))
  list(annotation_ld_scores=scores,annotation_sums=marker_counts,
    annotation_overlap=overlap,reference_marker_count=sum(marker_counts),
    sets=sets,assumption='ordinary scores masked by disjoint set; no cross-set LD contribution')
}
