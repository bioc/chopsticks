setOldClass(c("haplotype", "genotype"), test=TRUE)

#content should be raw - TODO: put restriction 
setClass("snp.matrix", contains="matrix")

setClass("X.snp.matrix", representation("snp.matrix", Female="logical"),
         contains="snp.matrix")

setClass("snp", contains="raw")
setClass("X.snp", representation(Female="logical"), contains="snp")

setClass("snp.dprime")
         
# constructor for new()
setMethod("initialize", 
          "snp.matrix", 
          function(.Object, ...) {
            if(is.null(.Object)) {
              stop("object is null to constructor");
            }
            .Object <- callNextMethod()
            # Please don't remove the tests - they are to do with promises and not redundant.
            # The accessor methods transparently passes to .Data, the assignment methods are different...
            if(!is.null(dim(.Object)) && (dim(.Object)[1] > 0) && (dim(.Object)[2] > 0)) {
              if (mode(.Object) != "raw") {
                cat("coercing object of mode ", mode(.Object), " to snp.matrix\n")
                mode(.Object@.Data) <- "raw"
              }
              if (is.null(dimnames(.Object))) {
                cat("object has no names - using numeric order for row/column names\n")
                dimnames(.Object@.Data) <- c(list(as.character(1:(dim(.Object)[1]))), list(as.character(1:(dim(.Object)[2]))))
              }
            }
            .Object
          })

setMethod("initialize", 
          "X.snp.matrix", 
          function(.Object, ...) {
            .Object <- callNextMethod()
            # The accessor methods transparently passes to .Data, the assignment methods are different...
            if (length(.Object@Female)==0){
              warning("Sex guessed from heterozygosity")
              guess <- .guessSex(.Object)
              .Object@Female <- guess$Female
              if (any(!guess$Female & (guess$Heterozygosity>0), na.rm=T)) {
                warning("heterozygous calls for males set to NA")
                .Object@.Data <- .forceHom(.Object@.Data, guess$Female)
              }
            }
            else {
              if (length(.Object@Female)!=nrow(.Object))
                stop("female argument incorrect length")
              if (any(is.na(.Object@Female)))
                warning("female argument contains NAs")
              het <- row.summary(.Object)$Heterozygosity
              if (any(!.Object@Female & (!is.na(het)&(het>0)))){
                warning("heterozygous calls for males set to NA")
                .Object@.Data <- .forceHom(.Object@.Data, .Object@Female)
              }                
            }
          .Object
        })

setMethod("[", signature(x="snp.matrix",i="ANY",j="ANY",drop="ANY"),
          function(x, i, j, drop) {
            if (missing(drop))
              drop <- TRUE
            if (missing(i)) {
              if (missing(j))
                return(x)
              else
                x <-  x@.Data[,j,drop=drop]
            }
            else {
              if (missing(j))
                x <- x@.Data[i,,drop=drop]
              else 
                x <- x@.Data[i,j,drop=drop]
            }
            class(x) <-
              if (is.matrix(x)) {
                "snp.matrix"
              } else { 
                "snp"
              }
            # setting an S4 class doesn't not automatically
            # set the object's internal tag; do it manually
            x <- asS4(x)
            x
          }
)

setMethod("[", signature(x="X.snp.matrix",i="ANY",j="ANY",drop="ANY"),
          function(x, i, j, drop) {
            if (missing(drop))
              drop <- TRUE
            female <- x@Female
            names(female) <- rownames(x)
            if (missing(i)) {
              if (missing(j))
                return(x)
              else 
                x <-  x@.Data[,j,drop=drop]
            }
            else {
              female <- female[i]
              if (missing(j))
                x <- x@.Data[i,,drop=drop]
              else 
                x <- x@.Data[i,j,drop=drop]
            }
            if (is.matrix(x)) {
              new("X.snp.matrix", x, Female=female)
            } else {
              x.names <- names(x)
              res<- new("X.snp", x, Female=female)
              # new() seems to destroy names in this case,
              # so we'll copy it back afterwards
              names(res) <- x.names
              res
            }
          })

# The sub-assignment methods

setMethod("[<-", signature(x="X.snp.matrix",i="ANY",j="ANY",value="X.snp.matrix"),
          function(x, i, j, value) {
            # The raw sub assignment should *just work*:
            if (!missing(i) & !missing(j)) {
              x@.Data[i,j] <- value
            } else if (missing(i) & missing(j)) {
              x@.Data[,] <- value
            } else if (missing(i)) {
              x@.Data[,j] <- value
            } else if (missing(j)) {
              x@.Data[i,] <- value
            }
            # All we want to do is to shuffle the rows
            # if a row index is passed
            if(!missing(i) & missing(j)) {
              slot(x, "Female")[i] <- slot(value, "Female")
            }
            x
          })

setAs("snp.matrix", "numeric",
      function(from) {
        to <- unclass(from)
        mode(to) <- "integer"
        if (any(to>3))
          warning("Illegal SNP code(s)")
        to[to==0 | to>3] <- NA
        to-1
      })

setAs("snp", "character",
      function(from) {
        result <- c("", "A/A", "A/B", "B/B")[1+as.integer(from)]
        names(result)<- names(from)
        result
      }
      )

setAs("X.snp", "character",
      function(from) {
        fem <- from@Female
        if (length(fem)==1)
          fem <- rep(fem, length(from))
        ifelse(fem,
               c("", "A/A", "A/B", "B/B")[1+as.integer(from)],
               c("", "A/Y", "Error", "B/Y")[1+as.integer(from)]
               )
      }
      )

setAs("snp.matrix", "character",
      function(from) {
        df <- dim(from)
        to <- c("", "A/A", "A/B", "B/B")[1+as.integer(from)]
        dim(to) <- df
        to
      })

setAs("snp.matrix", "X.snp.matrix",
      function(from) {
        new("X.snp.matrix", from)
      })

setAs("X.snp.matrix", "character",
      function(from) {
        df <- dim(from)
        to <- c("", "A/Y", "Error", "B/Y", "", "A/A", "A/B", "B/B")[
            1+as.integer(from)+ 4*rep(from@Female, ncol(from))]
        dim(to) <- df
        to
      })

setAs("snp", "numeric",
      function(from) ifelse(from==00, NA, as.integer(from)-1))

setAs("X.snp", "numeric",
      function(from) ifelse(from==00, NA, as.integer(from)-1))

#2007-02-27 does this really work with the genetics package loaded?
# Error in callNextMethod() : a call to callNextMethod() appears in
# a call to "asMethod", but the call does not seem to come
# from either a generic function or another 'callNextMethod'

setAs("X.snp", "genotype",
      function(from) {
        warning("X loci are not explicitly dealt with in the genetics package")
        callNextMethod()
      })


setAs("snp", "genotype",
      function(from) {
        to <- as(from,"numeric")+1
        attr(to, "levels") <- c( "A/A", "A/B", "B/B")
        attr(to ,"allele.names") <- c("B", "A")
        attr(to,"allele.map") <- matrix(c("A", "A", "A", "B", "B", "B"),
                                      byrow=TRUE, ncol=2)
        attr(to, "class") <- c("genotype", "factor")
        to
    })
      
setAs("matrix", "snp.matrix",
      function(from) {
        if (is.numeric(from)) {
          valid  <- (from==0 | from==1 | from==2)
          if(!all(is.na(from) | valid)) 
            warning("values other than 0, 1 or 2 set to NA")
          to <- from+1
          to[is.na(from)|!valid] <- 0
          mode(to) <- "raw"
        }
        else if (is.raw(from)) {
          not.valid  <- (from>3)
          if (any(not.valid))
            warning("values other than 01, 02 or 03 set to NA (00)")
          to <- from
          to[not.valid] <- as.raw(0)
        }
        else 
          stop("Can only convert `numeric' or `raw' matrices") 
        new("snp.matrix", to)
      })
           
# Summary of a snp.matrix

setGeneric("summary")

setMethod("summary", "snp.matrix",
   function(object) {
     .Call("snp_summary", object, PACKAGE="snpMatrix")
   })

setMethod("summary", "X.snp.matrix",
   function(object) {
     .Call("X_snp_summary", object, PACKAGE="snpMatrix")
   })

# Imputation

setClass("snp.reg.imputation")


setMethod("summary", "snp.reg.imputation",
   function(object) {
     r2 <- .Call("r2_impute", object, PACKAGE="snpMatrix")
     table(cut(r2[,1], c((0:9)/10, 0.95, 0.99, 1.0)), r2[,2],
           dnn=c("R-squared", "SNPs used"))
   })

# Show and plot methods

setMethod("show", "snp.reg.imputation",
   function(object) {
     to <- names(object)
     for (i in 1:length(object)) {
       cat(to[i], " ~ ", paste(object[[i]]$snps, collapse="+"),
           " (R-squared = ",object[[i]]$r.squared, ")\n", sep="")
     }
   })


setMethod("plot", signature(x="snp.reg.imputation", y="missing"),
   function(x, y, ...) {
     mat <- summary(x)
     if (colnames(mat)[1]=="0") 
       mat <- mat[,-1]
     n <- nrow(mat)
     m <- ncol(mat)
     val <- barplot(t(mat[n:1,]), legend.text=paste(1:m, "tag SNPs"),
                    beside=F, col=heat.colors(m),
                    xlab="r-squared", ylab="Number of SNPs",
                    names.arg=rep("", n), ...)
     mtext(rownames(mat)[n:1], at=val, side=1, las=2, line=0.5)
   })

setMethod("show", "X.snp.matrix",
   function(object) {
     nr <- nrow(object)
     nc <- ncol(object) 
     cat("An X.snp.matrix with ", nr, "rows and ", nc,
         "columns, \nholding data for ")
     fem <- object@Female
     cat(sum(!fem), " males and ", sum(fem), " females\n")
     if (nr>1)
       cat("Row names: ", rownames(object)[1],"...", rownames(object)[nr],"\n")
     else
       cat("Row name: ", rownames(object)[1], "\n")
     if (nc>1)
       cat("Col names: ", colnames(object)[1],"...", colnames(object)[nc],"\n")
     else
       cat("Col name: ", colnames(object)[1],"\n")
       
   })

setMethod("show", "snp.matrix",
   function(object) {
     nr <- nrow(object)
     nc <- ncol(object) 
     cat("A snp.matrix with ", nr, "rows and ", nc,
         "columns\n")
     if (nr>1)
       cat("Row names: ", rownames(object)[1],"...", rownames(object)[nr],"\n")
     else
       cat("Row name: ", rownames(object)[1],"\n")
     if (nc>1)
       cat("Col names: ", colnames(object)[1],"...", colnames(object)[nc],"\n")
     else
       cat("Col name: ", colnames(object)[1],"\n")
   })

setMethod("show", "X.snp",
          function(object) {
            fem <- object@Female
            if (length(fem)==1) {
              cat("Snp(s) on the X chromosome ")
              if (fem)
                cat(" (female)\n")
              else
                cat(" (male)\n")
              print(as(object, "character"))
            }
            else {
              cat("A snp on the X chromosome in",
                  sum(!fem), " males and ", sum(fem), " females\n")
              print(as(object, "character"))
            }
          })

# setMethod("names", signature(x="X.snp"), function(x) { names(as(x,"character"))})

setMethod("show", "snp",
          function(object) {
            cat("Autosomal snp(s):\n")
            print(as(object, "character"))
          })

setMethod("is.na", "snp.matrix", function(x){ x==0})

setGeneric("ld.with")

setMethod("ld.with", signature(data="snp.matrix",snps="character", include.itself="ANY", signed.r="ANY"),
          function(data, snps,include.itself=TRUE, signed.r=FALSE) {
            if(!include.itself)
              stop("ld.with(snp.matrix, snps) always includes the input snps")
            snp.names <- colnames(data)            
            # work out the zero-offset idices, and check every input is found
            snp.idx <- match(snps, snp.names) - 1
            snp.idx <- as.integer(snp.idx) # default real
            if(any(is.na(snp.idx)))
               stop("some of snps is not found in the snp.matrix object. Please check.")
            .Call("ld_with", data, snp.idx, signed.r, PACKAGE="snpMatrix")
          })

.rbind2 <- function(x,y){
  .External("snp_rbind",x, y, PACKAGE="snpMatrix")
}
snp.rbind <- function(...){
  .External("snp_rbind", ..., PACKAGE="snpMatrix")
}
.cbind2 <- function(x,y){
  .External("snp_cbind",x, y, PACKAGE="snpMatrix")
}
snp.cbind <- function(...){
  .External("snp_cbind", ..., PACKAGE="snpMatrix")
}

setMethod("rbind2", signature(x="snp.matrix", y="snp.matrix"), .rbind2)
setMethod("cbind2", signature(x="snp.matrix", y="snp.matrix"), .cbind2)



# Tests

setClass("snp.tests.single", 
         representation(snp.names="character", chisq="matrix",
                        N="integer", N.r2="numeric"))
setClass("snp.tests.single.score", 
         representation("snp.tests.single", U="matrix", V="matrix"),
         contains="snp.tests.single")

setMethod("[",
          signature(x="snp.tests.single", i="ANY",
                    j="missing", drop="missing"),
          function(x, i, j, drop) {
            if (is.character(i)) {
              i <- match(i, x@snp.names)
              if (any(is.na(i))) {
                warning(sum(is.na(i)),
                        " SNPs couldn't be found in snp.tests.single object\n")
                i <- i[!is.na(i)]
              }
            } else if (is.logical(i)) {
              if (length(i)!=length(x@snp.names))
                stop("logical selection array has incorrect length")
              i <- (1:length(i))[i]
            } else if (is.numeric(i)) {
              if (min(i)<0 || max(i)>length(x@snp.names))
                stop("selection index out of range")
            } else {
              stop("illegal selection array")
            }
            if (length(x@N.r2)>0)
              N.r2 <- x@N.r2[i]
            else
              N.r2 <- numeric(0)
            new("snp.tests.single",
              snp.names=x@snp.names[i], chisq=x@chisq[i,, drop=FALSE],
                N=x@N[i], N.r2=N.r2)
          })

setMethod("[",
          signature(x="snp.tests.single.score", i="ANY",
                    j="missing", drop="missing"),
          function(x, i, j, drop) {
            if (is.character(i)) {
              i <- match(i, x@snp.names)
              if (any(is.na(i))) {
                warning(sum(is.na(i)),
                        " SNPs couldn't be found in snp.tests.single object\n")
                i <- i[!is.na(i)]
              }
            } else if (is.logical(i)) {
              if (length(i)!=length(x@snp.names))
                stop("logical selection array has incorrect length")
              i <- (1:length(i))[i]
            } else if (is.numeric(i)) {
              if (min(i)<0 || max(i)>length(x@snp.names))
                stop("selection index out of range")
            } else {
              stop("illegal selection array")
            }
            if (length(x@N.r2)>0)
              N.r2 <- x@N.r2[i]
            else
              N.r2 <- numeric(0)
            new("snp.tests.single.score",
              snp.names=x@snp.names[i], chisq=x@chisq[i,, drop=FALSE],
                N=x@N[i], N.r2=N.r2,
                U=x@U[i,,drop=FALSE], V=x@V[i,,drop=FALSE])
          })
                  
setMethod("summary", "snp.tests.single",
          function(object) {
            if (length(object@N.r2)>0)
              summary(data.frame(N=object@N, N.r2=object@N.r2,
                        Chi2=object@chisq,
                        P.1df=pchisq(object@chisq[,1], df=1, lower.tail=FALSE),
                        P.2df=pchisq(object@chisq[,2], df=2, lower.tail=FALSE)
                        ))
            else
               summary(data.frame(N=object@N, Chi2=object@chisq,
                        P.1df=pchisq(object@chisq[,1], df=1, lower.tail=FALSE),
                        P.2df=pchisq(object@chisq[,2], df=2, lower.tail=FALSE)
                        ))
          })

            
setMethod("show", "snp.tests.single",
          function(object) {
            if (length(object@N.r2)>0)
              print(data.frame(N=object@N, N.r2=object@N.r2, Chi2=object@chisq,
                        P.1df=pchisq(object@chisq[,1], df=1, lower.tail=FALSE),
                        P.2df=pchisq(object@chisq[,2], df=2, lower.tail=FALSE),
                        row.names=object@snp.names))
            else
              print(data.frame(N=object@N, Chi2=object@chisq,
                        P.1df=pchisq(object@chisq[,1], df=1, lower.tail=FALSE),
                        P.2df=pchisq(object@chisq[,2], df=2, lower.tail=FALSE),
                        row.names=object@snp.names))
            })

# There are no standard generics for these, but the call to standardGeneric
# avoids a warning message 
  
setGeneric("p.value", function(x, df) standardGeneric("p.value"),
           useAsDefault=FALSE)

setGeneric("chi.squared", function(x, df) standardGeneric("chi.squared"),
           useAsDefault=FALSE)

setGeneric("pool2", function(x, y, score) standardGeneric("pool2"),
           useAsDefault=FALSE)

setMethod("p.value", signature(x="snp.tests.single", df="numeric"),
          function(x, df) {
            if (df>2)
              stop("df must be 1 or 2")
            p <- pchisq(x@chisq[,df], df=df, lower.tail=FALSE)
            names(p) <- x@snp.names
            p
          })

setMethod("chi.squared", signature(x="snp.tests.single", df="numeric"),
          function(x, df) {
            if (df>2)
              stop("df must be 1 or 2")
            chi2 <- x@chisq[,df]
            names(chi2) <- x@snp.names
            chi2
          })

setMethod("p.value", signature(x="snp.tests.single.score", df="numeric"),
          function(x, df) {
            if (df>2)
              stop("df must be 1 or 2")
            p <- pchisq(x@chisq[,df], df=df, lower.tail=FALSE)
            names(p) <- x@snp.names
            p
          })

setMethod("chi.squared", signature(x="snp.tests.single.score", df="numeric"),
          function(x, df) {
            if (df>2)
              stop("df must be 1 or 2")
            chi2 <- x@chisq[,df]
            names(chi2) <- x@snp.names
            chi2
          })

setMethod("pool2",
          signature(x="snp.tests.single.score",y="snp.tests.single.score",
                    score="logical"),
          function(x, y, score) {
            all.snps <- union(x@snp.names, y@snp.names)
            can.pool <- intersect(x@snp.names, y@snp.names)
            x.only <- setdiff(x@snp.names, can.pool)
            y.only <- setdiff(y@snp.names, can.pool)
            if (length(can.pool)>0) {
              xsel <- match(can.pool, x@snp.names)
              ysel <- match(can.pool, y@snp.names)
              N <- x@N[xsel] + y@N[ysel]
              U <- x@U[xsel,] + y@U[ysel,]
              V <- x@V[xsel,] + y@V[ysel,]
              if (length(x@N.r2)>0) {
                if (length(y@N.r2)>0)
                  Nr2 <- x@N.r2[xsel] + y@N.r2[ysel]
                else
                  Nr2 <- x@N.r2[xsel] + y@N[ysel]
              } else {
                if (length(y@N.r2)>0)
                  Nr2 <- x@N[xsel]+  y@N.r2[ysel]
                else
                  Nr2 <- numeric(0)
              }
            } else {
              N <- NULL
              U <- NULL
              V <- NULL
              Nr2 <- numeric(0)
            }
            if (length(x.only)>0) {
              xsel <- match(x.only, x@snp.names)
              N <- c(N, x@N[xsel])
              U <- rbind(U, x@U[xsel,])
              V <- rbind(V, x@V[xsel,])
              if (length(Nr2>0)) {
                if (length(x@N.r2)>0) 
                  Nr2 <- c(Nr2, x@N.r2[xsel])
                else
                  Nr2 <- c(Nr2, x@N[xsel])
              } else {
                if (length(x@N.r2)>0)
                  Nr2 <- c(N, x@N.r2[xsel])
              }
            }
            if (length(y.only)>0) {
              ysel <- match(y.only, y@snp.names)
              N <- c(N, y@N[ysel])
              U <- rbind(U, y@U[ysel,])
              V <- rbind(V, y@V[ysel,])
              if (length(Nr2>0)) {
                if (length(y@N.r2)>0) 
                  Nr2 <- c(Nr2, y@N.r2[ysel])
                else
                  Nr2 <- c(Nr2, y@N[ysel])
              } else {
                if (length(y@N.r2)>0)
                  Nr2 <- c(N, y@N.r2[ysel])
              }
            }
            chisq <- .Call("chisq_single", list(U=U, V=V, N=N),
                           PACKAGE="snpMatrix")
            if (score)
              res <- new("snp.tests.single.score",
                         snp.names=c(can.pool, x.only, y.only),
                         chisq=chisq, N=N, U=U, V=V, N.r2=Nr2)
            else
              res <- new("snp.tests.single",
                         snp.names=c(can.pool, x.only, y.only),
                         chisq=chisq, N=N, N.r2=Nr2)
            res
          })
