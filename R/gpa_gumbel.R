library(data.table)
library(ggplot2)
library(lmomRFA)

source('R/aux_functions.R')
source('R/aux_fit_smp.R')

trim <- c(0,0)
cluster.number <- 1

{

mx <- data.table(readRDS('data/DV.rds'))

mx <- mx[which(cluster %in% cluster.number),]

mx.trim <- mx[mx[, .I[which(!(dV %in% n.min(dV, trim[1])) & !(dV %in% n.max(dV, trim[2])))], by = SP_ID]$V1]
mx.trim[, p :=  (rank(dV)-.3)/(length(dV) + .4), by = SP_ID]

MX <- data.table(dcast(mx.trim[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
MX <- MX[, year := NULL]

lmom.atsite.t <- as.data.frame(t(apply(MX, 2, function(x) samlmu(x, trim = trim))))
at.site.para <- t(apply(lmom.atsite.t, 1, pelgpa))

head(lmom.atsite.t)
head(at.site.para)

# gg.MRD(lmom.atsite.t[,3:4]) + ggtitle(paste('Cluster', cluster.number))
# gg.homo(lmom.atsite.t)

# temp <- as.regdata(data.frame(name = names(MX),
#                               n = apply(MX, 2, function(x) length(which(!is.na(x)))),
#                               mean = lmom.atsite.t[,1],
#                               t = lmom.atsite.t[,2]/lmom.atsite.t[,1],
#                               t_3 = lmom.atsite.t[,3],
#                               t_4 = lmom.atsite.t[,4]))
# 
# para.reg <- regfit(temp, dist = 'gpa')

dta.fit <- sim(MX)
gp.gpa(dta.fit$data, dta.fit$REG, dta.fit$scaling_factor) + ggtitle(paste('Cluster', cluster.number))

}

s <- smp.gpa(dta.fit, length = 100,  type = 'zero')
f <- lapply(s, sim)

gc.gpa(dta.fit, f)
