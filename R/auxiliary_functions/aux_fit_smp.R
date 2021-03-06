sim <- function(extremes, dist = 'gpa', trim = c(0, 0)) { # stationary index flood method ;)
  
  l.atsite <- t(apply(extremes, 2, samlmu, trim = trim))
  
  if(dim(l.atsite)[1] > 2)  {
    
    l.atsite[,2] <- l.atsite[,2]/l.atsite[,1]
    
    w <- apply(extremes, 2, function(x) length(x[!is.na(x)]))
    w.l <- apply(l.atsite[,-1], 2, weighted.mean, w = w)
    
    structure(.Data = list(data = as.data.frame(extremes), 
                           scaling_factor = l.atsite[,1], 
                           REG = do.call(paste0('pel', dist), args = list(c(1, w.l))), 
                           dist = dist),
              sim.call = match.call(),
              class = 'sim')
  } else {
    
    structure(.Data = list(data = as.data.frame(extremes), 
                           scaling_factor = 1, 
                           REG = do.call(paste0('pel', dist), args = list(l.atsite)), 
                           dist = dist),
              sim.call = match.call(),
              class = 'sim')
  }
}

# zero - rgev, nezavis.
# one - NA jsou nepar. nasamplovany, radky rmvnorm, hodnoty jsou po sloupcich rgev, orankovani ...
# two - cela matice rmvnorm, NA jsou vlozeny do matice, zbytek jako metoda one
# three - cela matice rmvnorm, NA jsou vlozeny do matice, hodnoty jsou pres residua prevedeny na gev

sample <- function(...) {UseMethod('sample')}

sample.default <- base::sample

sample.sim <- function(model_object, length = 1, type = 'nonpar', na = T) { # testovaci nonpar sample (pro test fitovani)
  
  dta <- as.data.table(model_object$data)
  para <- model_object$REG
  
  i <- 1:length
  
  if(type == 'nonpar') {
    
    out <- mapply(function(i) {
      
      m <- resid.sim(model_object)
      m <- m[sample(1:dim(m)[1], dim(m)[1], replace = T),]
      m <- backtodata(m, model_object)
      return(m)
    }, i, SIMPLIFY = FALSE)
  }
  
  if(type == 'zero') {
    
    out <- mapply(function(i) {
      
      m <- apply(dta, 2, function(x) {x[which(!is.na(x))] <- rgpa(which(!is.na(x)), para); x})
      return(m)
    }, i, SIMPLIFY = FALSE)
  }
  
  if(type == 'para_cor') {
    
    out <- mapply(function(i) {
      cvr <- cov(resid.sim(model_object), use = 'pairwise.complete.obs')
      
      crr <- cov2cor(cvr)
      crr[] <- mean(crr[upper.tri(crr)], na.rm = TRUE)
      diag(crr) <- 1
      sdMat <- diag(sqrt(diag(cvr)))
      cvr <- sdMat %*% crr %*% t(sdMat)
      
      norm.res <- data.table(mvtnorm::rmvnorm(dim(res)[1], sigma = cvr, method = 'chol'))
      
      if(na) {
        
        norm.res <- norm.res*(dta/dta) 
      }
      
      gpa.val <- sapply(norm.res, function(x) {qgpa(pnorm(x), para)}) ###################
      
      return(gpa.val)
    }, i, SIMPLIFY = FALSE)
  }
  
  
  structure(.Data = out,
            names = paste0('b_sample_',i),
            class = 'simsample')

}

fit <- function(...) {UseMethod('fit')}

fit.simsample <- function(smp, model_object) {
  
  cl <- attr(model_object, 'sim.call')
  
  if(all(!is.na(as.numeric(as.character(cl$trim)[2:3])))) {
    
    trim <- as.numeric(as.character(cl$trim)[2:3])
  } else {
    
    trim <- as.numeric(as.character(formals(sim)$trim)[2:3])
  }
  
  structure(.Data = lapply(smp, function(x) {sim(x, dist = cl$dist, trim = trim)}),
            class = 'simsample')
}

