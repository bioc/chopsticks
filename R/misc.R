row.summary <- function(object) {
   if (inherits(object, "snp.matrix"))
     .Call("row_summary", object, PACKAGE="snpMatrix")
   else
     stop("not a snp.matrix object")
}

col.summary <- function(object) {
   if (inherits(object, "snp.matrix")) {
     if (inherits(object, "X.snp.matrix")) 
       .Call("X_snp_summary", object, PACKAGE="snpMatrix")
     else  
       .Call("snp_summary", object, PACKAGE="snpMatrix")
   }
   else
     stop("not a snp.matrix object")
 }

