test_that("gsuite installation is explicit and uses both repositories", {
  requested <- NULL
  ready <- FALSE
  testthat::local_mocked_bindings(install.packages=function(pkgs,lib,repos,...) {
    requested <<- list(packages=pkgs,lib=lib,repos=repos)
    ready <<- TRUE
  }, packageVersion=function(pkg,lib.loc) package_version(if(ready) "0.1.2" else "0.0.0"), .package="utils")
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
