#I have a couple of small suggestions for your tutorial.  It looks
#really nice and I like ld.with().  Only thing I wanted beyond that was 
#to to make it work with several snps.  See ld.withmany() below (I know 
#you can probably improve it lots!).  The only change I made to the basis 
#which was ld.with() was to include an entry for each snp with itself 
#(default rsq=1, dprime=1, lod=NA), so that binding the results for 
#several snps together made sense.

.ld.withmany <- function (ld.data, snps.wanted,do.join=FALSE) {
  if (class(ld.data) != "snp.dprime")
    stop("data to be must be of class snp.dprime")
  snp.names <- attr(ld.data, "snp.names")
  no.of.snps <- length(snp.names)
  max.depth <- dim(ld.data$dprime)[2]
  max.width <- dim(ld.data$dprime)[1]
  if ((max.width + 1) != no.of.snps)
    stop("malformed ld data object")
#  names.components <- names(ld.data)
#  pos <- which(snp.names %in% snps.wanted)
#  start.snp.pos <- pmax(1, (pos - max.depth))
#  end.snps.pos <- pmin((pos + max.depth), no.of.snps)
#  end.snp.offset <- pmin(max.depth, end.snps.pos - pos)
  nwanted <- length(snps.wanted)
  ld.return <- structure(vector("list",nwanted),names=snps.wanted)
  for(i in 1:nwanted) {
    pos <- which(snp.names == snps.wanted[i])
    start.snp.pos <- max(1, (pos - max.depth))
    end.snps.pos <- min((pos + max.depth), no.of.snps)
    end.snp.offset <- min(max.depth, end.snps.pos - pos)
     
    if (pos > 1) {
      vec1 <- start.snp.pos:(pos - 1)
      vec2 <- pos - vec1
      result.names <- c(snp.names[vec1],snps.wanted[i])
      rsq2 <- c(ld.data$rsq2[cbind(vec1, vec2)],1)
      dprime <- c(ld.data$dprime[cbind(vec1, vec2)],1)
      lod <- c(ld.data$lod[cbind(vec1, vec2)],NA)
    }
    else {
      result.names <- snps.wanted[i]
      rsq2 <- dprime <- 1
      lod <- NA
    }
    if (end.snp.offset) {
      rsq2 <- c(rsq2, ld.data$rsq2[pos, 1:end.snp.offset])
      dprime <- c(dprime, ld.data$dprime[pos, 1:end.snp.offset])
      lod <- c(lod, ld.data$lod[pos, 1:end.snp.offset])
      result.names <- c(result.names, snp.names[(pos + 
                                                 1):end.snps.pos])
    }
    ld.return[[i]] <- data.frame(dprime = dprime, rsq2 = rsq2, lod = 
                                 lod, row.names = result.names)
  }
  if(!do.join)
    return(ld.return)
  
  snps.for.each <- lapply(ld.return,rownames)
  #nsnps.for.each <- unlist(lapply(snps.for.each,length))
  all.snps <- unique(unlist(snps.for.each))
  all.snps <- snp.names[snp.names %in% all.snps] # order correctly
  rsq2 <- dprime <- lod <- 
    matrix(as.numeric(NA),length(all.snps),nwanted,
           dimnames=list(all.snps,snps.wanted))
  for(i in 1:nwanted) {
    snps.this <- rownames(ld.return[[i]])
    rsq2[snps.this,snps.wanted[i]] <- ld.return[[i]]$rsq2
    lod[snps.this,snps.wanted[i]] <- ld.return[[i]]$lod
    dprime[snps.this,snps.wanted[i]] <- ld.return[[i]]$dprime
  }
  return(list(rsq2=rsq2,dprime=dprime,lod=lod))
}

