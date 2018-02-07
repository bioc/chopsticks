.onLoad <- function(libname, pkgname) {
  library.dynam("chopsticks", pkgname, libname)
  ## methods:::bind_activation(TRUE)
}

## .Last.lib <- function(libname, package) {
##   methods:::bind_activation(FALSE)
## }

.onUnload <- function(libpath) {
  library.dynam.unload("chopsticks", libpath)
}
