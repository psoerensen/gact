# Prepare disjoint marker-set LDSC inputs from existing ordinary LD scores.
# For whole chromosomes, this is the within-chromosome LD approximation.
# For arbitrary sets with LD across their boundaries, supply true annotation
# LD scores instead of using this masked-score construction.
prepare_gcorr_set_scores <- function(marker_ids, ldscores, sets,
                                    marker_counts=setNames(lengths(sets),names(sets))) {
  scores <- gcorr::gcorr_set_scores(marker_ids, ldscores, sets)
  stopifnot(is.numeric(marker_counts),!is.null(names(marker_counts)),
    setequal(names(marker_counts),names(sets)),!anyDuplicated(names(marker_counts)))
  marker_counts <- marker_counts[names(sets)]
  stopifnot(all(is.finite(marker_counts)),all(marker_counts>=lengths(sets)),
    all(marker_counts==floor(marker_counts)))
  overlap <- diag(as.numeric(marker_counts),nrow=length(sets))
  dimnames(overlap) <- list(names(sets),names(sets))
  list(annotation_ld_scores=scores,annotation_sums=marker_counts,
    annotation_overlap=overlap,reference_marker_count=sum(marker_counts),
    sets=sets,assumption='ordinary scores masked by disjoint set; no cross-set LD contribution')
}
