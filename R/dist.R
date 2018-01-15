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
    out <- -1/para[3] * log(pmax(0, 1 - para[3] * (q - para[1])/para[2]))
  }
  
  return(1 - exp(-pmax(out, 0)))
}
