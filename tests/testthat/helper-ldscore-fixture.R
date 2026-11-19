ldscore_fixture <- function(build = "GRCh38") {
  root <- tempfile("gact-ldscores-")
  dir.create(root)
  prefix <- file.path(root, "g1000_eur")
  # Four samples; columns include identical and monomorphic SNPs,
  # a missing genotype, another chromosome and an off-backbone SNP.
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, 0x38, 0x38, 0xff, 0x39, 0x38, 0x38)),
           paste0(prefix, ".bed"))
  writeLines(paste("f", paste0("i", 1:4), 0, 0, 0, -9), paste0(prefix, ".fam"))
  markers <- data.frame(rsids = paste0("rs", 1:6), chr = c(1, 1, 1, 1, 2, 2),
    pos = c(10L, 2000000L, 3000000L, 4000000L, 10L, 20L), ea = "A", nea = "G")
  writeLines(with(markers, paste(chr, rsids, 0, pos, ea, nea)), paste0(prefix, ".bim"))
  data.table::fwrite(markers, file.path(root, "markers.txt.gz"))
  list(root = root, prefix = prefix, ga = list(metadata = list(build = build),
    dirs = c(marker = root), rsids = markers$rsids[5:1]))
}
