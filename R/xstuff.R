.guessSex <- function(snps) {
  het <- row.summary(snps)$Heterozygosity
  fem <- (het>0.1)
  repeat {
    female <- fem
    fm <- mean(het[female], na.rm=TRUE)
    mm <- mean(het[!female], na.rm=TRUE)
    if (is.na(fm) || is.na(mm))
      break
    fd <- (het-fm)^2/(fm*(1-fm))
    md <- (het-mm)^2/((mm+0.001)*(1-mm))
    fem <- (fd<md)
    if (all(fem==female, na.rm=T))
      break
  }
  list(Female=female, Heterozygosity=het)
}

.forceHom <- function(xsnps, female) {
  .Call("force_hom", xsnps, female, PACKAGE="snpMatrix")
}
