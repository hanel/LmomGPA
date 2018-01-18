######## GEV ########

rgev <- function(n, para = c(xi, alpha, k)) {
  
  qgev(runif(n), para)
}

qgev <- function(p, para = c(xi, alpha, k)) {
  
  if (para[3] == 0) {
    para[1] - para[2]*log(-log(p))
  } else {
    para[1] + para[2]/para[3]*(1 - (-log(p))^para[3])
  }
}

pgev <- function(q, para = c(xi, alpha, k)) {
  
  if (para[3] == 0) {
    out <- (q - para[1])/para[2]
  } else {
    out <- -1/para[3]*log(pmax(0, 1 - para[3]*(q - para[1])/para[2]))
  }
  
  return(exp(-exp(-out)))
}

######## GPA ########

rgpa <- function(n, para = c(xi, alpha, k)) {
  
  qgpa(runif(n), para)
}

qgpa <- function(p, para = c(xi = 0, alpha = 1, k = 0)) {
  
  if (para[3] == 0) {
    para[1] + para[2]*(-log(1 - p))
  } else {
    para[1] + para[2]/para[3]*(1 - (1 - p)^(para[3]))
  }
}

pgpa <- function(q, para = c(xi, alpha, k)) {
  
  if (para[3] == 0) {
    out <- (q - para[1])/para[2]
  } else {
    out <- -1/para[3]*log(pmax(0, 1 - para[3]*(q - para[1])/para[2]))
  }
  
  return(1 - exp(-pmax(out, 0)))
}

######## L-mom ########

lmanual <- function (x) {
  
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
