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
              hetm <- !.Object@Female & (!is.na(het)&(het>0))
              if (any(hetm)){
                warning(sum(hetm, na.rm=TRUE),
                        " heterozygous calls for males set to NA")
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
            cl <-
              if (is.matrix(x)) {
                "snp.matrix"
              } else { 
                "snp"
              }
            attr(cl, "package") <- "snpMatrix"
            class(x) <- cl
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
        dimnames(to) <- dimnames(from)
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

setClass("snp.reg.imputation", contains="list")


setMethod("summary", "snp.reg.imputation",
   function(object) {
     r2 <- .Call("r2_impute", object, PACKAGE="snpMatrix")
     table(cut(r2[,1], c((0:9)/10, 0.95, 0.99, 1.0)), r2[,2],
           dnn=c("R-squared", "SNPs used"))
   })

setMethod("[",
       signature(x="snp.reg.imputation", i="ANY", j="missing", drop="missing"),
       function(x, i){
         if (is.character(i))
           i <- match(i, names(x), nomatch=0)
         res <- new("snp.reg.imputation", x@.Data[i])
         names(res) <- names(x)[i]
         res
       })


# Show and plot methods

setMethod("show", "snp.reg.imputation",
   function(object) {
     to <- names(object)
     for (i in 1:length(object)) {
       cat(to[i], " ~ ", paste(object[[i]]$snps, collapse="+"),
           " (MAF = ", object[[i]]$maf, 
           ", R-squared = ", object[[i]]$r.squared,
           ")\n", sep="")
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
                    xlab=expression(r^2), ylab="Number of SNPs",
                    names.arg=rep("", n), cex.lab=1.3, ...)
     mtext(rownames(mat)[n:1], at=val, side=1, las=2, line=0.5, cex=0.8)
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
setClass("snp.tests.glm",
         representation(test.names="character", chisq="numeric", df="integer",
                        N="integer"))
setClass("snp.tests.glm.score", representation("snp.tests.glm", score="list"),
         contains="snp.tests.glm")

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
              snp.names=x@snp.names[i], chisq=x@chisq[i,,drop=FALSE],
                N=x@N[i], N.r2=N.r2,
                U=x@U[i,,drop=FALSE], V=x@V[i,,drop=FALSE])
          })

setMethod("[",
          signature(x="snp.tests.glm", i="ANY",
                    j="missing", drop="missing"),
           function(x, i, j, drop) {
            if (is.character(i)) 
              i <- match(i, x@test.names, nomatch=0)
            new("snp.tests.glm",
                test.names = x@test.names[i],
                chisq = x@chisq[i],
                df = x@df[i],
                N = x@N[i])
          })

setMethod("[",
          signature(x="snp.tests.glm.score", i="ANY",
                    j="missing", drop="missing"),
          function(x, i, j, drop) {
            if (is.character(i)) 
              i <- match(i, x@test.names, nomatch=0)
            new("snp.tests.glm.score",
                test.names = x@test.names[i],
                chisq = x@chisq[i],
                df = x@df[i],
                N = x@N[i],
                score = x@score[i])
          })
                   
setMethod("summary", "snp.tests.single",
          function(object) {
            if (length(object@N.r2)>0)
              summary(data.frame(N=object@N, N.r2=object@N.r2,
                        Chi.squared=object@chisq,
                        P.1df=pchisq(object@chisq[,1], df=1, lower.tail=FALSE),
                        P.2df=pchisq(object@chisq[,2], df=2, lower.tail=FALSE)
                        ))
            else
               summary(data.frame(N=object@N,
                        Chi.squared=object@chisq,
                        P.1df=pchisq(object@chisq[,1], df=1, lower.tail=FALSE),
                        P.2df=pchisq(object@chisq[,2], df=2, lower.tail=FALSE)
                        ))
          })

setMethod("summary", "snp.tests.glm",
          function(object) {
            chi2 <- object@chisq
            df <- object@df
            p <- pchisq(chi2, df=df, lower.tail=FALSE)
            summary(data.frame(row.names=object@test.names,
                               Chi.squared=chi2, Df=df, p.value=p))
          })

setMethod("show", "snp.tests.single",
          function(object) {
            if (length(object@N.r2)>0)
              print(data.frame(N=object@N, N.r2=object@N.r2,
                        Chi.squared=object@chisq,
                        P.1df=pchisq(object@chisq[,1], df=1, lower.tail=FALSE),
                        P.2df=pchisq(object@chisq[,2], df=2, lower.tail=FALSE),
                        row.names=object@snp.names))
            else
              print(data.frame(N=object@N,
                        Chi.squared=object@chisq,
                        P.1df=pchisq(object@chisq[,1], df=1, lower.tail=FALSE),
                        P.2df=pchisq(object@chisq[,2], df=2, lower.tail=FALSE),
                        row.names=object@snp.names))
            })

setMethod("show", "snp.tests.glm",
          function(object) {
            chi2 <- object@chisq
            df <- object@df
            p <- pchisq(chi2, df=df, lower.tail=FALSE)
            print(data.frame(row.names=object@test.names,
                             Chi.squared=chi2, Df=df, p.value=p))
          })

# There are no standard generics for these, but the call to standardGeneric
# avoids a warning message 
  
setGeneric("p.value", function(x, df) standardGeneric("p.value"),
           useAsDefault=FALSE)

setGeneric("chi.squared", function(x, df) standardGeneric("chi.squared"),
           useAsDefault=FALSE)

setGeneric("deg.freedom", function(x) standardGeneric("deg.freedom"),
           useAsDefault=FALSE)

setGeneric("effect.sign", function(x, simplify) standardGeneric("effect.sign"),
           useAsDefault=FALSE)

setGeneric("sample.size", function(x) standardGeneric("sample.size"),
           useAsDefault=FALSE)

setGeneric("pool2", function(x, y, score) standardGeneric("pool2"),
           useAsDefault=FALSE)

setGeneric("names")

setMethod("p.value", signature(x="snp.tests.single", df="numeric"),
          function(x, df) {
            if (df>2)
              stop("df must be 1 or 2")
            p <- pchisq(x@chisq[,df], df=df, lower.tail=FALSE)
            names(p) <- x@snp.names
            p
          })

setMethod("p.value", signature(x="snp.tests.glm", df="missing"),
          function(x, df) {
            p <- pchisq(q=x@chisq, x@df, lower.tail=FALSE)
            names(p) <- x@test.names
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

setMethod("chi.squared", signature(x="snp.tests.glm", df="missing"),
          function(x, df) {
            chi2 <- x@chisq
            names(chi2) <- x@test.names
            chi2
          })


setMethod("deg.freedom", signature(x="snp.tests.glm"),
          function(x) {
            df <- x@df
            names(df) <- x@test.names
            df
          })

setMethod("effect.sign", signature(x="snp.tests.glm", simplify="logical"),
          function(x, simplify=TRUE) {
            res <- sapply(x@score, function(x) sign(x$U))
            names(res) <- x@test.names
            res
         })

setMethod("effect.sign",
          signature(x="snp.tests.single.score", simplify="missing"),
          function(x, simplify) {
            res <- sign(x@U[,1])
            names(res) <- x@snp.names
            res
          })

setMethod("sample.size", signature(x="snp.tests.glm"),
          function(x) {
            res <- x@N
            names(res) <- x@test.names
            res
          })

setMethod("sample.size", signature(x="snp.tests.single"),
          function(x) {
            res <- x@N
            names(res) <- x@snp.tests
            res
          })

setMethod("names", signature(x="snp.tests.single"), function(x) x@snp.names)

setMethod("names", signature(x="snp.tests.glm"), function(x) x@test.names)

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
              U <- x@U[xsel,,drop=FALSE] + y@U[ysel,,drop=FALSE]
              V <- x@V[xsel,,drop=FALSE] + y@V[ysel,,drop=FALSE]
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
              U <- rbind(U, x@U[xsel,,drop=FALSE])
              V <- rbind(V, x@V[xsel,,drop=FALSE])
              if (length(Nr2>0)) {
                if (length(x@N.r2)>0) 
                  Nr2 <- c(Nr2, x@N.r2[xsel])
                else
                  Nr2 <- c(Nr2, x@N[xsel])
              } else {
                if (length(x@N.r2)>0)
                  Nr2 <- c(N, x@N.r2[xsel])
              }
              N <- c(N, x@N[xsel])
            }
            if (length(y.only)>0) {
              ysel <- match(y.only, y@snp.names)
              U <- rbind(U, y@U[ysel,,drop=FALSE])
              V <- rbind(V, y@V[ysel,,drop=FALSE])
              if (length(Nr2>0)) {
                if (length(y@N.r2)>0) 
                  Nr2 <- c(Nr2, y@N.r2[ysel])
                else
                  Nr2 <- c(Nr2, y@N[ysel])
              } else {
                if (length(y@N.r2)>0)
                  Nr2 <- c(N, y@N.r2[ysel])
              }
              N <- c(N, y@N[ysel])
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

# need to edit this

setMethod("pool2",
          signature(x="snp.tests.glm.score",y="snp.tests.glm.score",
                    score="logical"),
          function(x, y, score) {
            nm.x <- x@test.names
            nm.y <- y@test.names
            if (is.null(nm.x) || is.null(nm.y)) {
              if (length(x)!=length(y))
                stop("Cannot pool unnamed snp.test.glm objects of different lengths")
              res <- .Call("pool2_glm", x, y, score, PACKAGE="snpMatrix")
              return(res)
            }
            to.pool <- intersect(nm.x, nm.y)
            if (length(to.pool)>0) {
              res <- .Call("pool2_glm",
                           x[to.pool],
                           y[to.pool], score, 
                           PACKAGE="snpMatrix")
 
            }
            else {
              res <- NULL
            }
            x.only <- setdiff(nm.x, to.pool)
            y.only <- setdiff(nm.y, to.pool)
            if (length(x.only)==0 && length(y.only)==0)
              return(res);
            ix <- match(x.only, x@test.names, nomatch=0)
            iy <- match(y.only, y@test.names, nomatch=0)
            if (score)
              res <- new("snp.tests.glm.score",
               test.names=c(res@test.names, x@test.names[ix], y@test.names[iy]),
               chisq=c(res@chisq, x@chisq[ix], y@chisq[iy]),        
               df=c(res@df, x@df[ix], y@df[iy]),        
               N=c(res@N, x@N[ix], y@N[iy]),        
               score=append(res@score, append(x@score[ix], y@score[iy])))
            else
               res <- new("snp.tests.glm",
               test.names=c(res@test.names, x@test.names[ix], y@test.names[iy]),
               chisq=c(res@chisq, x@chisq[ix], y@chisq[iy]),        
               df=c(res@df, x@df[ix], y@df[iy]),        
               N=c(res@N, x@N[ix], y@N[iy]))       
            res
          })

# Switch allele methods 

setGeneric("switch.alleles",
           function(x, snps) standardGeneric("switch.alleles"),
           useAsDefault=FALSE)

setMethod("switch.alleles", signature(x="snp.matrix", snps="ANY"),
          function(x, snps) {
            if (is.character(snps)) 
              snps <- match(snps, colnames(x))
            else if (is.logical(snps)) {
              if (length(snps)!=ncol(x))
                stop("logical snp selection vector is wrong length")
              snps <- (1:ncol(x))[snps]
            }
            else if (!is.integer(snps))
              stop("snp selection must be character, logical or integer")
            if (any(is.na(snps) | snps>ncol(x) | snps<1))
              stop("illegal snp selection")
            .Call("smat_switch", x, snps, PACKAGE="snpMatrix")
          })

setMethod("switch.alleles", signature(x="snp.tests.single.score", snps="ANY"),
          function(x, snps) {
            if (is.character(snps)) {
              snps <- match(snps, x@snp.names)
              if (any(is.na(snps)))
                warning(sum(is.na(snps)),
                        " SNP names were not found in tests object")
              snps <- snps[!is.na(snps)]
            } 
            ntest <- length(x@snp.names)
            if (is.logical(snps)) {
              if (length(snps)!=ntest)
                stop("incompatible arguments")
              if (sum(snps)==0)
                return(x)
            } else if (is.numeric(snps)) {
              if (length(snps)==0)
                return(x)
              if (max(snps)>ntest || min(snps)<1)
                stop("incompatible arguments")
            } else {
              stop("SNPs to be switched must be indicated by name, position, or by a logical vector")
            }
            res <- x
            res@U[snps,1] <- -x@U[snps,1]
            if (ncol(x@U)==3) {
              res@U[snps,2] <- -x@U[snps,2]
              res@U[snps,3] <- x@U[snps,3] - x@U[snps,2]
              res@V[snps,4] <- x@V[snps,4] - 2*x@V[snps,3] +
                x@V[snps,2]
              res@V[snps,3] <- x@V[snps,3] - x@V[snps,2]
            } else {
              res@U[snps,2] <- x@U[snps,2] - x@U[snps,1]
              res@V[snps,3] <- x@V[snps,3] - 2*x@V[snps,2] +
                x@V[snps,1]
              res@V[snps,2] <- x@V[snps,2] - x@V[snps,1]
            }
            res
          })

setMethod("switch.alleles", signature(x="snp.tests.glm.score",
                                      snps="character"),
          function(x, snps) {
            new("snp.tests.glm.score", test.names=x@test.names,
                chisq=x@chisq, df=x@df, N=x@N,
                score=lapply(x@score, sw1.glm, snps))
          })

sw1.glm <- function(x, snps) {
  to.switch <- names(x$U) %in% snps
  if (!any(to.switch))
    return(x)
  res <- x
  u <- ifelse(to.switch, -1, 1)
  v <- u %*% t(u)
  v <- v[!lower.tri(v)]
  res$U <- u*res$U
  res$V <- v*res$V
  res
}


pool <- function(..., score=FALSE) {
  argl <- list(...)
  na <- length(argl)
  while (na > 0 && is.null(argl[[na]])) {
    argl <- argl[-na]
    na <- na - 1
  }
  if (na<2)
    stop("need at least two test sets to pool")
  if (na>2) {
    p2 <- pool2(..1, ..2, score=TRUE)
    r <- do.call(pool, c(p2, argl[3:na], score=score))
  }
  else
    r <- pool2(..1, ..2, score=score)
  r
}

# To do

# names method for snp and X.snp

