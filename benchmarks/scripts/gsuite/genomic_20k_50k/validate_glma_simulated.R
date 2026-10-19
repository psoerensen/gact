# Untimed reference checks, in a separate process from the memory benchmark.
args <- commandArgs(TRUE)
fixture <- readRDS(args[1]); linear <- readRDS(args[2]); out <- args[3]
library(qgg)
frozen <- readRDS('build/benchmarks/manuscript-20261005-142319/website_sim_stats-1/output.rds')
stopifnot(identical(unname(fixture$Glist$ids), unname(unlist(frozen$Glist$ids))),
          identical(as.numeric(fixture$Glist$sample_source_indices_0based),
                    as.numeric(seq_along(fixture$y) - 1L)))
markers <- unique(c(fixture$causal,
  frozen$rsids[round(seq(1, length(frozen$rsids), length.out = 12))]))
W <- qgg::getG(Glist = frozen$Glist, chr = 1L, rsids = markers, scale = FALSE)
stopifnot(nrow(W) == length(fixture$y), ncol(W) == length(markers),
          setequal(colnames(W), markers), !anyNA(W))
reference <- t(vapply(seq_len(ncol(W)), function(j) {
  coef(summary(lm(fixture$y ~ W[, j])))[2, c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')]
}, numeric(4)))
observed <- linear$associations[match(colnames(W), linear$associations$marker), ]
stopifnot(!anyNA(observed$marker), all(observed$status == 'ok'))
validation <- data.frame(marker = colnames(W),
  beta_difference = observed$beta - reference[, 1],
  se_difference = observed$se - reference[, 2],
  statistic_difference = observed$statistic - reference[, 3],
  p_difference = observed$p - reference[, 4])
write.csv(validation, file.path(out, 'lm-reference.csv'), row.names = FALSE)
stopifnot(max(abs(as.matrix(validation[-1]))) < 1e-8)
writeLines(c(capture.output(sessionInfo()),
             paste('qgg', packageVersion('qgg'), find.package('qgg'))),
           file.path(out, 'validation-session.txt'))
cat('Independent lm reference passed for', nrow(validation), 'markers\n')
