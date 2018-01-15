source('R/dist.R')

sim <- function(extremes, dist = 'gpa', trim = c(0, 0)) { # stationary index flood method ;)
  
  l.atsite <- t(apply(extremes, 2, samlmu, trim = trim))
  l.atsite[,2] <- l.atsite[,2]/l.atsite[,1]
  w <- apply(extremes, 2, function(x) length(x[!is.na(x)]))
  
  w.l <- apply(l.atsite[,-1], 2, weighted.mean, w = w)

  list(data = as.data.frame(extremes), 
       scaling_factor = l.atsite[,1], 
       REG = do.call(paste0('pel', dist), args = list(c(1, w.l))), 
       dist = dist)
  
}


# zero - rgev, nezavis.
# one - NA jsou nepar. nasamplovany, radky rmvnorm, hodnoty jsou po sloupcich rgev, orankovani ...
# two - cela matice rmvnorm, NA jsou vlozeny do matice, zbytek jako metoda one
# three - cela matice rmvnorm, NA jsou vlozeny do matice, hodnoty jsou pres residua prevedeny na gev

smp.gpa <- function(sim, length = 1, type = 'nonpar') { # testovaci nonpar sample (pro test fitovani)
  
  dta <- as.data.table(sim[[1]])
  
  i <- 1:length
  
  if(type == 'nonpar') {
    out <- mapply(function(i) {
      m <- dta[sample(1:nrow(dta), nrow(dta), replace = TRUE), ]
      return(m)
    }, i, SIMPLIFY = FALSE)
  }
  
  if(type == 'zero') {
    out <- mapply(function(i) {
      m <- apply(dta, 2, function(x) {x[which(!is.na(x))] <- rgpa(which(!is.na(x)), sim$REG); x})
      return(m)
    }, i, SIMPLIFY = FALSE)
  }
  
  names(out) <- paste0('sample_',i)
  out
}