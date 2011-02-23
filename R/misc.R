row.summary <- function(object) {
     .Call("row_summary", object, PACKAGE="snpMatrix")
   }
