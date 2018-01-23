# library(nim)
# 
# data('precip_max')
# 
# dta <- precip_max[,c(1,12)]
# extremes(dta) = -1
# 
# n <- nim( ~1, data = dta)
# 
# smp <- sample(n, length = 50)
# f <- fit(n, smp)

quantile.nim <- function(nim, p = NULL, Tm = c(2, 5, 10, 50), at_site = TRUE){
  
  cl = match.call()
  qO = Vectorize(function(p){
    r = data.frame(XI*(1 - G/K*(1 -(-log(p))^(-K))))
    if (model_info(nim)$cvrt == 'I') {
      r = data.frame(I = 1, r[1, ])
      if (is.null(names(XI))) {
        names(r) = c(names(r)[1], 'value')
      } else {
        names(r) = c(names(r)[1], names(XI))#if (!at_site) {names(r)[2] = 'regional'}
      }
    } else {
      r = data.frame(nim$REG[[1]], r, check.names = FALSE)
      names(r)[1] = model_info(nim)$cvrt
      names(r)[2:ncol(r)] = if (at_site) (names(nim$XI)) else ('regional')
    }
    return(r)
  }, SIMPLIFY = FALSE)
  
  if (is.null(p)) p = 1 - 1 / Tm
  
  if(at_site) {
    XI <- outer(regional(nim)$XI, atsite(nim))
    G <- matrix(data = rep(x = exp(nim$REG$G), times = dim(XI)[2]),
                nrow = dim(XI)[1])
    K <- matrix(data = rep(x = nim$REG$K, times = dim(XI)[2]),
                nrow = dim(XI)[1])
  } else {
    XI <- regional(nim)$XI
    G <- exp(nim$REG$G)
    K <- nim$REG$K
  }
  
  res = qO(p)
  names(res) = p
  res = rbindlist(res, idcol = 'p')
  res
}

AD.test.multidist <- function(val, location = 0, scale = 1, shape = 0, dist = 'gev') {
  
  val <- unlist(val)
  val <- val[!is.na(val)]
  
  location <- unlist(location)
  scale <- unlist(scale)
  shape <- unlist(shape)
  
  p <- do.call(what = paste0('p', dist), args = list(q = val, para = c(location, scale, shape)))
  
  u <- sort(p[(p != 0) & (p != 1)])
  nr <- length(u)
  -nr - 1/nr*sum((2*1:nr - 1)*log(u) + (2*nr - 2*1:nr + 1)*log(1 - u))
  
}

n.max <- function(x, n = 2){
  l <- length(x)
  if (n > l){
    warning('vole !!!')
    n <- length(x)
  }
  if (n == 0) {
    x[NULL]
  } else {
    sort(x)[(l - n + 1):l]
  }
}

n.min <- function(x, n = 2){
  l <- length(x)
  if(n > l){
    warning('vole !!!')
    n <- length(x)
  }  
  if (n == 0) {
    x[NULL]
  } else {
    sort(x)[1:n]
  }
}

# pelgpa(samlmu(x))
# gpa.para(moje.lm(x))

################################### resid

resid.sim <- function(model_object){
  
  dta <- model_object$data
  para <- model_object$REG
  sf <- model_object$scaling_factor
  
  resi <- suppressMessages(melt(dta))
  resi <- data.table(variable = names(dta), sf = sf, t(para))[resi, on = c('variable')]
  resi <- resi[, `:=`(resi = 1/k*log(1 + k*((value/sf)/alpha)),
                      aux = seq_len(.N)), by = 'variable'] # coles - sigma = alpha, xi = kappa #######
  
  dcast(resi, aux ~ variable, value.var = 'resi')[,-1]
}

backtodata <- function(res, model_object){
  
  para <- model_object$REG
  sf <- model_object$scaling_factor
  
  res <- suppressWarnings(melt(res))
  res <- data.table(variable = names(model_object$data), sf = sf, t(para))[res, on = c('variable')]
  dta <- res[, `:=`(val = ((alpha*(exp(value*k) - 1))/k)*sf,
                    aux = seq_len(.N)), by = 'variable']
  
  dcast(dta, aux ~ variable, value.var = 'val')[,-1]
}