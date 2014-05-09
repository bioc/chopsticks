read.HapMap.data <- function(url, verbose=FALSE, save=NULL, ...) {
  # in the future we might want to construct the URL
  base <- "http://www.hapmap.org/genotypes"
  build <- "latest"
  strand <- "fwd_strand"
  finish <- "redundant-filtered"

  if (length(save)) {
    saved.file <-  save
  } else {
    saved.file <- tempfile()
  }

  # mode = b is needed for windows
  status <- download.file(url, destfile=saved.file, mode="wb")

  # download.file() supposedly should throw error already, or return 0/1
  if ((status != 0) && (status != 1)) { 
    stop("Download has gone wrong\n");
  }

  snp.result <- .Call("read_hapmap_data", saved.file, verbose, PACKAGE="chopsticks")
  flush.console()
  
  if(length(snp.result)) {
    snp.result$snp.support$dbSNPalleles <- as.factor(snp.result$snp.support$dbSNPalleles)
    snp.result$snp.support$Assignment   <- as.factor(snp.result$snp.support$Assignment  )
    snp.result$snp.support$Chromosome   <- as.factor(snp.result$snp.support$Chromosome  )
    snp.result$snp.support$Strand       <- as.factor(snp.result$snp.support$Strand      )
    if (!length(save)) {  
      if(unlink(saved.file)) { # zero is success for unlink()
        cat("removing temp file failed - this is not fatal\n");
      }
    }
    cat(" ...conversion complete...\n");
  } else {
    cat("conversion failed\ntemp file left at ",
        saved.file,
        " for retrying conversion locally with url file://\n");
  }
  
  snp.result
}

#### Generated with
## zcat www.hapmap.org/downloads/samples_individuals/pedinfo2sample_CEU.txt.gz \
## | cut -f 4,7 | grep '^0' | cut -f 2 | perl -pe 's/:1$//;s/^.+://;' |sort |uniq 

.HapMap.CEU.founders.list <- 
  c("NA06985", "NA06993", "NA06994", "NA07000", "NA07022", "NA07034", "NA07055",
    "NA07056", "NA07345", "NA07357", "NA11829", "NA11830", "NA11831", "NA11832",
    "NA11839", "NA11840", "NA11881", "NA11882", "NA11992", "NA11993", "NA11994",
    "NA11995", "NA12003", "NA12004", "NA12005", "NA12006", "NA12043", "NA12044",
    "NA12056", "NA12057", "NA12144", "NA12145", "NA12146", "NA12154", "NA12155",
    "NA12156", "NA12234", "NA12236", "NA12239", "NA12248", "NA12249", "NA12264",
    "NA12716", "NA12717", "NA12750", "NA12751", "NA12760", "NA12761", "NA12762",
    "NA12763", "NA12812", "NA12813", "NA12814", "NA12815", "NA12872", "NA12873",
    "NA12874", "NA12875", "NA12891", "NA12892")

.HapMap.YRI.founders.list <-
  c("NA18501", "NA18502", "NA18504", "NA18505", "NA18507", "NA18508", "NA18516",
    "NA18517", "NA18522", "NA18523", "NA18852", "NA18853", "NA18855", "NA18856",
    "NA18858", "NA18859", "NA18861", "NA18862", "NA18870", "NA18871", "NA18912",
    "NA18913", "NA19092", "NA19093", "NA19098", "NA19099", "NA19101", "NA19102",
    "NA19116", "NA19119", "NA19127", "NA19128", "NA19130", "NA19131", "NA19137",
    "NA19138", "NA19140", "NA19141", "NA19143", "NA19144", "NA19152", "NA19153",
    "NA19159", "NA19160", "NA19171", "NA19172", "NA19192", "NA19193", "NA19200",
    "NA19201", "NA19203", "NA19204", "NA19206", "NA19207", "NA19209", "NA19210",
    "NA19222", "NA19223", "NA19238", "NA19239")

# place-holders

#read.HapMap.support <- function(url, chromosome=...) {
#}
  
#read.HapMap.support <- function(url, population=...) {
#}
