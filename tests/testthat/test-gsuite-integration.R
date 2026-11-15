test_that("gsuite installation refreshes the official index and verifies the target library", {
  gsuite <- "https://psoerensen.github.io/gtools/gsuite/repository"
  cran <- "https://cloud.r-project.org"
  index <- function(versions) {
    cbind(Package = stats::setNames(names(versions), names(versions)), Version = versions)
  }
  versions <- c(gbase = "0.1.7", gcorr = "0.2.0")
  published <- c(gbase = "0.1.12", gcorr = "0.1.9.1", gsea = "0.1.1.1")
  lib <- tempfile("gsuite-library-")
  requested <- NULL
  refreshes <- 0L
  succeeds <- TRUE
  testthat::local_mocked_bindings(
    available.packages = function(repos, type, ignore_repo_cache = FALSE) {
      expect_identical(type, "source")
      if (identical(repos, gsuite)) {
        expect_true(ignore_repo_cache)
        refreshes <<- refreshes + 1L
        return(index(published))
      }
      expect_identical(unname(repos), cran)
      index(c(gbase = "9.0.0", Rcpp = "1.1.1"))
    },
    installed.packages = function(lib.loc, noCache) {
      expect_identical(lib.loc, lib)
      expect_true(noCache)
      index(versions)
    },
    install.packages = function(pkgs, lib, repos, type, available, ...) {
      requested <<- list(packages = pkgs, lib = lib, repos = repos,
                         type = type, available = available)
      if (succeeds) versions[pkgs] <<- available[pkgs, "Version"]
    }, .package = "utils")
  packages <- c("gbase", "gcorr", "gsea")
  expect_identical(install_gsuite(packages, lib, repos = c(CRAN = "@CRAN@")), packages)
  expect_identical(requested$packages, c("gbase", "gsea"))
  expect_identical(requested$lib, lib)
  expect_identical(unname(requested$repos), c(gsuite, cran))
  expect_identical(requested$type, "source")
  expect_identical(requested$available["gbase", "Version"], "0.1.12")
  expect_true("Rcpp" %in% rownames(requested$available))
  expect_identical(unname(versions["gcorr"]), "0.2.0")

  requested <- NULL
  install_gsuite(packages, lib)
  expect_null(requested)
  expect_identical(refreshes, 2L)

  versions[["gbase"]] <- "0.1.7"
  succeeds <- FALSE
  expect_error(install_gsuite("gbase", lib, repos = c(CRAN = cran)),
               "gbase \\(0\\.1\\.12\\)")
  expect_error(install_gsuite("gmap", lib), "No compatible gsuite source version available for: gmap")
  expect_error(install_gsuite("qgg"), "standalone")
})

test_that("database membership uses compact indices and retains repeated members", {
  expect_identical(gbase::mapSets(list(a=c("m2","absent","m2"),empty="absent"),
    rsids=c("m1","m2")), list(a=c(2L,2L),empty=integer()))
})

test_that("database gene sets select available markers without tutorial bookkeeping", {
  directory <- tempfile("gact-sets-")
  dir.create(directory)
  on.exit(unlink(directory, recursive=TRUE))
  saveRDS(list(gene=c("m3","m1","m1"), absent="missing", other="m2"),
    file.path(directory,"ensg2rsids.rds"))
  database <- list(dirs=c(gsets=directory))
  expect_identical(getMarkerSets(database,feature="Genesplus",rsids=c("m1","m3")),
    list(gene=c("m3","m1")))
})

test_that("summary extraction preserves backbone gaps and ancestry for imputation", {
 directory <- tempfile("gact-summary-"); dir.create(directory)
 on.exit(unlink(directory,recursive=TRUE))
 markers <- data.frame(rsids=c("a","b","c"),chr="1",pos=c(100,200,300),ea="A",nea="G")
 write.table(markers,file.path(directory,"markers.txt.gz"),sep="\t",row.names=FALSE,quote=FALSE)
 file1 <- file.path(directory,"one.txt"); file2 <- file.path(directory,"two.txt")
 one <- cbind(markers[1,,drop=FALSE],b=.2,seb=.1,p=.01,n=100)
 two <- cbind(markers[2,,drop=FALSE],b=-.3,seb=.1,p=.001,n=200)
 data.table::fwrite(one,file1); data.table::fwrite(two,file2)
 database <- list(rsids=markers$rsids,study=list(id=c("S1","S2"),neff=c(100,200),ancestry=c("EUR","EAS")),
   studyfiles=c(S1=file1,S2=file2),dirs=c(marker=directory))
 single <- getMarkerStat(database,"S1",rm.na=FALSE)
 expect_identical(single$rsids,markers$rsids)
 expect_equal(single$z,c(2,NA,NA))
 expect_identical(single$ea,rep("A",3))
 expect_identical(attr(single,"ancestry"),c(S1="EUR"))
 multiple <- getMarkerStat(database,c("S1","S2"),rm.na=FALSE)
 expect_identical(rownames(multiple$z),markers$rsids)
 expect_equal(unname(multiple$z),matrix(c(2,NA,NA,NA,-3,NA),3,2))
 expect_identical(attr(multiple,"ancestry"),c(S1="EUR",S2="EAS"))
 expect_identical(attr(multiple,"marker_alleles")$rsids,markers$rsids)
})
