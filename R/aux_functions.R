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

gp.new <- function(n) {
  
  res.gp <- data.table(attributes(n)$data)
  
  res.gp <- melt(res.gp, id.var = 1)
  res.gp <- data.table(variable = names(n$XI), XI = n$XI)[res.gp, on = c('variable')]
  res.gp <- res.gp[!is.na(value), p := (rank(value) - .3) / (length(value) + .4), by = variable]
  res.gp <- res.gp[, p := -log(-log(as.numeric(p)))]
  res.gp <- res.gp[, val.xi := value/XI]
  
  p <- exp(-exp(-seq(min(res.gp$p, na.rm = T), max(res.gp$p, na.rm = T), .01)))
  
  ifelse(mean(n$REG$K) == 0,
         assign('val', mean(n$XI) - exp(mean(n$REG$G))*mean(n$XI)*log(-log(p))), 
         assign('val', mean(n$XI) + exp(mean(n$REG$G))*mean(n$XI)*((-log(p))^(-mean(n$REG$K)) - 1)/mean(n$REG$K)))
  
  gev <- data.table(x = -log(-log(p)), y = val/mean(n$XI))
  
  return(ggplot(res.gp) +
           geom_line(aes(x = p, 
                         y = val.xi, 
                         group = variable), 
                     colour = 'steelblue1', 
                     alpha = .25) +
           geom_point(aes(x = p, 
                          y = val.xi), 
                      colour = 'steelblue4',
                      fill = 'steelblue4',
                      alpha = .5) +
           geom_line(data = gev, 
                     aes(x = x, 
                         y = y), 
                     col = 'red4', 
                     lwd = .75) + 
           theme_bw() +
           labs(x = '-log(-log(p))', 
                y = 'Value'))
  
}

gc.new <- function (n, f, ribbon.1 = c(0.05, 0.95), ribbon.2 = c(0.25, 0.75)) {
  
  prbs <- sort(c(ribbon.1, ribbon.2))
  
  qntl <- data.table(do.call(rbind, lapply(f, FUN = function(x)quantile.nim(x, p = seq(.01,.99,.01), at_site = T))))
  
  names(qntl)[3:dim(qntl)[2]] <- paste0('station.',1:(dim(qntl)[2] - 2))
  
  suppressWarnings(qntl.m <- melt(qntl[,-2]))
  qntl.m[, p := as.numeric(p)]

  res.gc <- qntl.m[, .(quantile = quantile(value, probs = prbs)), by = p]
  res.gc <- data.table(cbind(res.gc, c(paste0('q.',prbs[1]),paste0('q.',prbs[2]),paste0('q.',prbs[3]),paste0('q.',prbs[4]))))
  names(res.gc) <- c('p', 'value', 'q')
  res.gc <- res.gc[, value := value/mean(n$XI)]
  res.gc <- res.gc[, p := -log(-log(as.numeric(p)))]
  res.gc <- dcast(res.gc, p ~ q, value.var = 'value')
  
  p <- exp(-exp(-seq(min(res.gc$p, na.rm = T), max(res.gc$p, na.rm = T), .01)))
  
  ifelse(mean(n$REG$K) == 0,
         assign('val', mean(n$XI) - exp(mean(n$REG$G))*mean(n$XI)*log(-log(p))), 
         assign('val', mean(n$XI) + exp(mean(n$REG$G))*mean(n$XI)*((-log(p))^(-mean(n$REG$K)) - 1)/mean(n$REG$K)))
  
  gev <- data.table(x = -log(-log(p)), y = val/mean(n$XI))
  # gev <- data.table(x = -log(-log(p)), y = val)
  
  
  return(ggplot(res.gc) +
           geom_ribbon(aes_string(x = colnames(res.gc)[1], 
                                  ymin = colnames(res.gc)[2], 
                                  ymax = colnames(res.gc)[5]), 
                       fill = 'steelblue1', 
                       alpha = .5) +
           geom_ribbon(aes_string(x = colnames(res.gc)[1], 
                                  ymin = colnames(res.gc)[3], 
                                  ymax = colnames(res.gc)[4]), 
                       fill = 'steelblue4', 
                       alpha = .5) +
           geom_line(data = gev, 
                     aes(x = x, 
                         y = y), 
                     col = 'red4', 
                     lwd = .75) + 
           theme_bw() +
           labs(x = '-log(-log(p))', 
                y = 'Value'))

}

# gc_new(n, f)
# 
# gp.new(n)

gp.gpa <- function(dta, para, scaling.factor) {
  
  gpa.gp <- melt(dta, id.var = 1)
  gpa.gp <- data.table(variable = names(dta), l1 = scaling.factor)[gpa.gp, on = c('variable')]
  gpa.gp <- gpa.gp[!is.na(value), p := (rank(value) - .3) / (length(value) + .4), by = variable]
  gpa.gp <- gpa.gp[, val.scaled := value/l1]
  
  p <- seq(min(gpa.gp$p, na.rm = T), max(gpa.gp$p, na.rm = T), .001)
  
  gpa.q <- data.table(x = p,
                      y = if(para[3] == 0) {
                        para[1] + para[2]*(-log(1 - p))
                      } else {
                        para[1] + para[2]/para[3]*(1 - (1 - p)^(para[3]))
                      })
  
  gp <- ggplot(gpa.gp, aes(x = -log(-log(p)), y = val.scaled, group = variable)) +
    geom_line(colour = 'steelblue1', alpha = .25) +
    geom_point(colour = 'royalblue4', fill = 'royalblue4', alpha = .5) +
    geom_line(data = gpa.q, aes(x = -log(-log(x)), y = y, group = 1), col = 'red4', lwd = .75) +
    theme_bw() +
    labs(x = '-log(-log(p))', y = 'Value')

  return(gp)
}


AD.test.multidist <- function(val, location = 0, scale = 1, shape = 0, distr = 'gev') {
  
  val <- unlist(val)
  val <- val[!is.na(val)]
  
  location <- unlist(location)
  scale <- unlist(scale)
  shape <- unlist(shape)
  
  if (distr == 'gev') {
    # if (shape == 0) {
    #   p <- exp(-exp(-(val - location)/scale))
    # } else {
    #   p <- exp(-pmax(1 + shape*(val - location)/scale, 0)^(-1/shape))
    # }
    p <- lmom::cdfgev(x = val, c(location, scale, shape))
  }
  
  if (distr == 'gpa') {
    # if (shape == 0) {
    #   p <- 1 - exp(-pmax(val - location, 0)/scale)
    # } else {
    #   p <- pmax(1 + shape*pmax(val - location, 0)/scale, 0)
    #   p <- 1 - p^(-1/shape)
    # }
    p <- lmom::cdfgpa(x = val, c(location, scale, shape))
  }
  
  u <- sort(p[(p != 0) & (p != 1)]) ################################# ?
  # plot(u)
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

moje.lm <- function (x) {
  
  x <- sort(x)
  
  n <- length(x)
  nn <- rep(n - 1, n)
  
  pp <- seq(0, n - 1)
  p1 <- pp/nn
  p2 <- p1*(pp - 1)/(nn - 1)
  p3 <- p2*(pp - 2)/(nn - 2)
  
  b0 <- sum(x)/n
  b1 <- sum(p1*x)/n
  b2 <- sum(p2*x)/n
  b3 <- sum(p3*x)/n
  
  l1 <- b0
  l2 <- 2*b1 - b0
  l3 <- 2*(3*b2 - b0)/(2*b1 - b0) - 3
  l4 <- 5*(2*(2*b3 - 3*b2) + b0)/(2*b1 - b0) + 6
  
  unlist(mget(paste0('l', 1:4)))
}

gpa.para <- function(l) {
  
  k <- (1 - 3*l[3])/(1 + l[3])
  alpha <- (1 + k)*(2 + k)*l[2]
  xi <- l[1] - (2 + k)*l[2]
  
  setNames(c(xi, alpha, k), c('xi', 'alpha', 'k'))
}

# pelgpa(samlmu(x))
# gpa.para(moje.lm(x))

gg.MRD <- function(taus) {
  
  num <- seq(min(taus[,1])*.5, max(taus[,1])*1.2, .01)
  mr <- data.table(t3 = num, t4 = num*(1 + 5*num)/(5 + num))
  
  names(taus) <- names(mr)
  
  lmrd <- ggplot(data = NULL, aes(x = t3, y = t4)) +
    geom_line(data = mr, colour = 'red4') +
    geom_point(data = taus, colour = 'grey15', fill = 'steelblue4', shape = 21) +
    theme_classic() +
    labs(x = 'L-skewness', y = 'L-kurtosis', title = 'GPA L-moment ratio diagram')
  
  return(lmrd)
}

gg.homo <- function(lmoms) {
  
  names(lmoms) <- c('l1', 'l2', 't3', 't4')
  
  homo <- ggplot(lmoms) +
    geom_point(aes(x = t3, y = l2/l1), colour = 'grey15', fill = 'steelblue4', shape = 21) +
    theme_classic() +
    labs(x = 'L-skweness', y = 'L-CV')
  
  return(homo)
}