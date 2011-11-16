setOldClass(c("haplotype", "genotype"), test=TRUE)

#content should be raw - TODO: put restriction 
setClass("snp.matrix", contains="matrix")

setClass("X.snp.matrix", representation(Female="logical"),
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
              if (any(!guess$Female&(guess$Heterozygosity>0))) {
                warning("heterozygous calls for males set to NA")
                .Object@.Data <- .forceHom(.Object@.Data, guess$Female)
              }
            }
            else {
              if (length(.Object@Female)!=nrow(.Object))
                stop("female argument incorrect length")
              if (any(is.na(.Object@Female)))
                stop("female argument contains NAs")
              if (any(!.Object@Female&(row.summary(.Object)$Heterozygosity>0))){
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
              names(res@.Data) <- x.names
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
        if(!all(is.na(from) | from==0 | from==1 | from==2))
          stop("values must be NA, 0, 1 or 2")  # FIXME: is this really right? invalid entries should auto-convert to NA
        to <- from+1
        to[is.na(from)] <- 0
        mode(to) <- "raw"
        new("snp.matrix", to)
      })
           
# Summary of a snp.matrix

setGeneric("summary")

setMethod("summary", "snp.matrix",
   function(object) {
     .Call("snp_summary", object, PACKAGE="chopsticks")
   })

setMethod("summary", "X.snp.matrix",
   function(object) {
     .Call("X_snp_summary", object, PACKAGE="chopsticks")
   })

# Show methods

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
            cat("An autosomal snp:\n")
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
            .Call("ld_with", data, snp.idx, signed.r, PACKAGE="chopsticks")
          })

.rbind2 <- function(x,y){
  .External("snp_rbind",x, y, PACKAGE="chopsticks")
}
snp.rbind <- function(...){
  .External("snp_rbind", ..., PACKAGE="chopsticks")
}
.cbind2 <- function(x,y){
  .External("snp_cbind",x, y, PACKAGE="chopsticks")
}
snp.cbind <- function(...){
  .External("snp_cbind", ..., PACKAGE="chopsticks")
}

setMethod("rbind2", signature(x="snp.matrix", y="snp.matrix"), .rbind2)
setMethod("cbind2", signature(x="snp.matrix", y="snp.matrix"), .cbind2)

