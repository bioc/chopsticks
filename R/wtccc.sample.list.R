wtccc.sample.list <- function(infile) {
  if(file.access(infile,mode=4))
    stop(infile, " cannot be read")
  
  filehandle <- gzfile(infile)
  first.line <- readLines(filehandle, n=1)
  close(filehandle)
  
  first.list <- unlist(strsplit(first.line, "\t"))
  len.list <- length(first.list)
  
  if (len.list < 6)
    stop("header has fewer than 6 fields")
  
  # len.list should be 5 + 2 * samples, an odd number
  if ((len.list+1) %% 2)
    stop("sample entries do not seem to be in pairs")

  sample.list <- first.list[seq(from=6, to=len.list, by=2)]

  len.sample.list <- length(sample.list)

  isA <- grep("_A$", sample.list)
  if ( length(isA) != len.sample.list)
    stop("Header does not seems to contain *_A/*_B")

  sample.list <- sub("_A$", "", sample.list)
  sample.list
}
