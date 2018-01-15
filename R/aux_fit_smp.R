source('R/dist.R')

sim <- function(extremes, dist = 'gpa', trim = c(0, 0)) { # stationary index flood method ;)
  
  l.atsite <- t(apply(extremes, 2, samlmu, trim = trim))
  l.atsite[,2] <- l.atsite[,2]/l.atsite[,1]
  w <- apply(extremes, 2, function(x) length(x[!is.na(x)]))
  
  w.l <- apply(l.atsite[,-1], 2, weighted.mean, w = w)

  structure(.Data = list(data = as.data.frame(extremes), 
                         scaling_factor = l.atsite[,1], 
                         REG = do.call(paste0('pel', dist), args = list(c(1, w.l))), 
                         dist = dist),
            sim.call = match.call(),
            class = 'sim')
  
}

# zero - rgev, nezavis.
# one - NA jsou nepar. nasamplovany, radky rmvnorm, hodnoty jsou po sloupcich rgev, orankovani ...
# two - cela matice rmvnorm, NA jsou vlozeny do matice, zbytek jako metoda one
# three - cela matice rmvnorm, NA jsou vlozeny do matice, hodnoty jsou pres residua prevedeny na gev

sample <- function(...) {UseMethod('sample')}

sample.default <- base::sample

sample.sim <- function(sim, length = 1, type = 'nonpar') { # testovaci nonpar sample (pro test fitovani)
  
  dta <- as.data.table(sim$data)
  
  i <- 1:length
  
  if(type == 'nonpar') {
    out <- mapply(function(i) {
      m <- dta[sample(1:nrow(dta), nrow(dta), replace = T), ]
      return(m)
    }, i, SIMPLIFY = FALSE)
  }
  
  if(type == 'zero') {
    out <- mapply(function(i) {
      m <- apply(dta, 2, function(x) {x[which(!is.na(x))] <- rgpa(which(!is.na(x)), sim$REG); x})
      return(m)
    }, i, SIMPLIFY = FALSE)
  }
  
  structure(.Data = out,
            names = paste0('b_sample_',i),
            class = 'simsample')

}

fit <- function(...) {UseMethod('fit')}

fit.simsample <- function(smp, sim) {
  
  cl <- attr(sim, 'sim.call')
  lapply(smp, function(x) {sim(x, dist = cl$dist, trim = as.numeric(as.character(cl$trim)[2:3]))})
}