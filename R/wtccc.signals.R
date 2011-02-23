read.wtccc.signals <- function(file, snp.list) {
  .Call("read_signals", file, snp.list, PACKAGE="snpMatrix")
}
