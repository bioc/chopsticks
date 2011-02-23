
#use wtccc.sample.list() against the smallest file (22?) to get a sample list. 
read.snps.chiamo <- function(filename, sample.list, threshold) {
  .Call("read_chiamo", filename, sample.list, threshold,
        PACKAGE="snpMatrix")
}
