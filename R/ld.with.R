ld.with <- function(data, snps, include.itself=as.logical(length(snps)-1), signed.r=NULL) {
  if (class(data) != "snp.dprime")
    stop("data to be must be of class snp.dprime")
  snp.names <- attr(data, "snp.names")
  nwanted <- length(snps)
  ld.return <- structure(vector("list",nwanted),names=snps)
  names.components <-  names(data)
  
   for(i in 1:nwanted) {
     ld.return[[i]] <- .ld.with.internal(data, snps[i], include.itself=TRUE, signed.r)
   }

  snps.for.each <- lapply(ld.return,rownames)
  nsnps.for.each <- unlist(lapply(snps.for.each,length))
  all.snps <- unique(unlist(snps.for.each))

  r.out <- dprime <- lod <-
    matrix(as.numeric(NA),length(all.snps),nwanted,
           dimnames=list(all.snps,snps))
  for(i in 1:nwanted) {
    snps.this <- rownames(ld.return[[i]])
    if ("rsq2" %in% names.components) {
      r.out[snps.this,snps[i]] <- ld.return[[i]]$rsq2
    } else {
      r.out[snps.this,snps[i]] <- ld.return[[i]]$r
    }
    lod[snps.this,snps[i]] <- ld.return[[i]]$lod
    dprime[snps.this,snps[i]] <- ld.return[[i]]$dprime
  }
  if ("rsq2" %in% names.components) {
    return(list(rsq2=r.out,dprime=dprime,lod=lod))
  } else {
    return(list(r=r.out,dprime=dprime,lod=lod))
  }
}

.ld.with.internal <- function(data, snps, include.itself=as.logical(length(snps)-1), signed.r=NULL) {
  # signed.r is not used here
  if (!is.null(signed.r))
    stop("Specifying signed.r is not supported for the extractor function ld.with()")
  if (class(data)!="snp.dprime")
    stop("data to be must be of class snp.dprime")
  
  snp.names <- attr(data, 'snp.names')
  no.of.snps <- length(snp.names)
  
  max.depth <- dim(data$dprime)[2]
  max.width <- dim(data$dprime)[1]
  
  if ((max.width +1) != no.of.snps)
    stop ("malformed ld data object")
  
  # the middle one can be rsq2 or r
  names.components <-  names(data)
  if ("rsq2" %in% names.components) {
    r.in <- data$rsq2
  } else {
    r.in <- data$r
  }

  pos <- which(snp.names == snps)
  start.snp.pos <- max(1, (pos - max.depth) )
  end.snps.pos <- min( (pos + max.depth), no.of.snps)  

  if (pos > 1) {
    #because 1:0 is valid
    vec1 <- start.snp.pos:(pos-1)
    vec2 <- pos - vec1
    
    result.names <- snp.names[start.snp.pos:(pos-1)]
    
    r.out <-   r.in[cbind(vec1, vec2)]
    dprime <- data$dprime[cbind(vec1, vec2)]
    lod <-    data$lod[cbind(vec1, vec2)]
  } else {
    # for pos=1, there is nothing before it
    result.names <- character(0)
    r.out <- numeric(0)
    dprime <- numeric(0)
    lod <- numeric(0)
  }

  if (include.itself) {
    result.names <- c(result.names, snps)
    r.out <- c(r.out, 1)
    dprime <- c(dprime, 1)
    lod <- c(lod, NA)          
  }
    
  # because 1:0 is a valid R expression...
  end.snp.offset <- min(max.depth, end.snps.pos - pos)
  if(end.snp.offset) {
  r.out <-   c(r.out,   r.in[pos,1:end.snp.offset]  )
  dprime <- c(dprime, data$dprime[pos,1:end.snp.offset])
  lod <-    c(lod,    data$lod[pos,1:end.snp.offset]   )
  result.names <- c(result.names, snp.names[(pos+1):end.snps.pos]) 
}
  if ("rsq2" %in% names.components) {
    data.frame(dprime=dprime, rsq2=r.out, lod=lod, row.names=result.names)
  } else {
    data.frame(dprime=dprime, r=r.out, lod=lod, row.names=result.names)
  }
}
