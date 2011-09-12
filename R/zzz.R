.onLoad <- function(lib, pkg) {
  library.dynam("chopsticks", package=pkg)
  methods:::bind_activation(TRUE)
}

.Last.lib <- function(libname, package) {
  methods:::bind_activation(FALSE)
}
