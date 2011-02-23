ibs.stats <- function (x) {
  .Call('do_ibs', x, PACKAGE="snpMatrix")
}
