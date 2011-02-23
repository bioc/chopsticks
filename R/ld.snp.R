#
# Copyright (C) 2006  Hin-Tak Leung
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 51 Franklin Street
# Fifth Floor, Boston, MA 02110-1301  USA.
#

#
# The 4 C routines used here are all described in sdfpw.h
# and implemented in sdfpw.c

# SEXP snp_pairwise(SEXP v, SEXP i, SEXP j);
# SEXP snp_pair_graphics(SEXP v, SEXP fileoutput, SEXP i, SEXP j, SEXP depth, SEXP do_notes);
# SEXP snp_dprime_draw(SEXP list_in, SEXP fileoutput, SEXP scheme, SEXP do_notes);
# SEXP snp_pair_range(SEXP v, SEXP i, SEXP j, SEXP depth, SEXP signed_r);


# The order of the arguments here is different from the C code,
# because depth is the one which most people would *not* leave as default
#  - 100 choosen arbitrarily as a "good" value.

ld.snp <- function(snpdata, depth=100, start=1, end=dim(snpdata)[2], signed.r=FALSE) {
  # sanity check of all the arguments
  if ((class(snpdata)!="snp.matrix") && (class(snpdata)!="X.snp.matrix") && (class(snpdata)!="snp.data.frame"))
    stop("snps argument must be of class snp.matrix/X.snp.matrix")

  # type check: we check/coerce to int in the C code as well, this is
  # just being paranoid 
  start <- as.integer(start)
  depth <- as.integer(depth)
  end  <- as.integer(end)
  signed.r <-  as.logical(signed.r)

  # range check: need a few more here
  if (start < 1)
    start <- 1

  if (end > dim(snpdata)[2])
    end <- dim(snpdata)[2]
  
  if ((depth > (end - start)) || depth == 0)
      depth <- end - start

  ans <- .Call("snp_pair_range", snpdata, start, end, depth, signed.r,
               PACKAGE="snpMatrix")
  # ans returns fully formed as a "snp.dprime" object
  # if we need further work we'll insert it here.
  ans
}

plot.snp.dprime <- function(x, filename, scheme="standard", do.notes=FALSE, metric=NULL, ...) {
  # sanity check arguments
  if (class(x)!="snp.dprime")
    stop("data to be drawn must be of class snp.dprime")
  
  filename <- as.character(filename)
  scheme <- as.character(scheme)
  do.notes <- as.logical(do.notes)

  if (is.null(filename))
    stop("you must specify an output file name!")
  
  #converting exotic types to types that the C code expects
  
  int.scheme <- 0
  if (scheme == "rsq") {
  # currently only two schemes are implemented, the "standard" and the "rsq",
  # see inside of the C code; assume everything non-rsq is standard.
    int.scheme <- 1
  }

  int.notes <- 0
  if (do.notes) {
    int.notes <- 1
  } # everything else is assumed FALSE
  
  #if (length(x$snp.names) > 1200) {
  #  cat("Acrobat reader has an implemention limit of 200 inches\n")
  #  cat("You will need a different pdf reader, e.g. xpdf/kpdf/evince/gsview\n")
  #}
  
  .Call("snp_dprime_draw", x, filename, int.scheme, int.notes, metric,
        PACKAGE="snpMatrix") 
  # plotting operates entirely on side effects, emit a reminder before returning
  # if (do.notes)
  #   cat("Don't forget to run ps2pdf -dEPSCrop", filename, ".\n")  
  invisible()
}

# This is the all-in-one version, it doesn't keep anything in memory,
# so is useful where one want to draw big diagrams
epsout.ld.snp <- function (snpdata, filename, start, end, depth, do.notes=FALSE) {
  # anybody calling this should know what he is doing,
  # so we are not setting any defaults here

  # sanity check of all the arguments
  if ((class(snpdata)!="snp.matrix") && (class(snpdata)!="X.snp.matrix") && (class(snpdata)!="snp.data.frame"))
    stop("snps argument must be of class snp.matrix/X.snp.matrix")

  # type check: we check/coerce to int in the C code as well, this is
  # just being paranoid 
  start <- as.integer(start)
  depth <- as.integer(depth)
  end  <- as.integer(end)
  
  .Call("snp_pair_graphics", snpdata, filename, start, end, depth, do.notes,
        PACKAGE="snpMatrix")
  invisible()
}

# this is mainly a debug routime for single pair result
pair.result.ld.snp <- function(snpdata, loc.snpA, loc.snpB) {
  if ((class(snpdata)!="snp.matrix") && (class(snpdata)!="X.snp.matrix") && (class(snpdata)!="snp.data.frame"))
    stop("snps argument must be of class snp.matrix/X.snp.matrix")
  
  .Call("snp_pairwise", snpdata, loc.snpA, loc.snpB,
        PACKAGE="snpMatrix")
  invisible()
}

#S3 generic
niceprint <- function(x, ...) NextMethod("print", x, ...)

print.snp.dprime <- function(x, ..., output="") {
  if (class(x)!="snp.dprime")
    stop("data to be must be of class snp.dprime")
  
  do.rsq2 <- ("rsq2" %in% names(x))
  
  snp.names <- attr(x, 'snp.names')
  
  if (do.rsq2) {
    cat (file=output, "#M1", "#M2", "rsq2", "Dprime", "lod","\n", sep="\t")
    r.maybe <- x$rsq2
  } else {
    cat (file=output, "#M1", "#M2", "r", "Dprime", "lod","\n", sep="\t")
    r.maybe <- x$r
  }
  
  max.depth <- dim(x$dprime)[2]
  
  for (i.snp in c( 1:(length(snp.names)-1) )) {
    for (j.snp in c( (i.snp+1):length(snp.names) )) {
      #cat(i.snp, j.snp, "\n")
      step <- j.snp - i.snp
      if (step > max.depth) {
        # the snp.dprime object has no data beyond max.depth
        break
      }
      cat(file=output,
          snp.names[i.snp], snp.names[j.snp],
          r.maybe[i.snp, step],
          x$dprime[i.snp, step],
          x$lod[i.snp, step],
          "\n", sep="\t", append=TRUE)
    }
  }
}
