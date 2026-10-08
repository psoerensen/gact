test_that("gsuite installation is explicit and uses both repositories", {
  requested <- NULL
  ready <- FALSE
  testthat::local_mocked_bindings(install.packages=function(pkgs,lib,repos,...) {
    requested <<- list(packages=pkgs,lib=lib,repos=repos)
    ready <<- TRUE
  }, packageVersion=function(pkg,lib.loc) package_version(if(ready) "0.1.3" else "0.0.0"), .package="utils")
  lib <- dirname(system.file(package="gbase"))
  expect_identical(install_gsuite(packages="gbase",lib=lib,repos=c(CRAN="@CRAN@")), "gbase")
  expect_identical(requested$packages, "gbase")
  expect_identical(unname(requested$repos), c(
    "https://psoerensen.github.io/gtools/gsuite/repository", "https://cloud.r-project.org"))
  requested <- NULL
  install_gsuite(packages="gbase",lib=lib)
  expect_null(requested)
  expect_error(install_gsuite(packages="qgg"), "standalone")
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
