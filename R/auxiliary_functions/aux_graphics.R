gumbelplot <- function(model_object, dist = NULL, method = 'base') {
  
  # res (init res will need to be rewritten for nim)
  
  dta <- as.data.table(model_object$data)
  para <- model_object$REG
  scaling.factor <- model_object$scaling_factor
  
  if(is.null(dist)) {dist <- attr(model_object, 'sim.call')$dist}
  
  # dta <- data.table(attributes(n)$data)
  
  res.gp <- suppressWarnings(melt(dta))
  res.gp <- data.table(variable = names(dta), sf = scaling.factor)[res.gp, on = c('variable')]
  res.gp <- res.gp[!is.na(value), p := (rank(value) - .3) / (length(value) + .4), by = variable]
  res.gp <- res.gp[, gumbel.variate := -log(-log(as.numeric(p)))]
  res.gp <- res.gp[, scaled.value := value/sf]
  
  p <- seq(min(res.gp$p, na.rm = T), max(res.gp$p, na.rm = T), .001)
  
  regional <- data.table(x = -log(-log(p)), y = do.call(paste0('q', dist), list(p, para)))
  
  # graphics
  
  if(method == 'base') {
    
    sres.gp <- split(res.gp[!is.na(res.gp$scaled.value),], res.gp$variable[!is.na(res.gp$scaled.value)])
    
    plot(NULL,
         xlim = c(min(regional$x), max(regional$x)*1.15),
         ylim = c(min(res.gp$scaled.value, na.rm = T), max(res.gp$scaled.value, na.rm = T)),
         bty = 'n',
         xlab = expression(-log(-log(p))),
         ylab = 'Value',
         main = 'Gumbel plot')
    
    grid()
    
    lapply(sres.gp, function(x) {
      points(sort(x$gumbel.variate),
             sort(x$scaled.value),
             pch = 21,
             col = 'grey15', 
             bg = '#36648b50',
             cex = .75)
      lines(sort(x$gumbel.variate),
            sort(x$scaled.value),
            col = '#36648b50')
    })
    
    lines(regional,
          type = 'l',
          col = 'red4',
          lwd = .75)
  }
  
  if(method %in% c('ggplot', 'plotly')) {
    
    if(method == 'ggplot') {
      
      gp <- ggplot2::ggplot(res.gp) +
        ggplot2::geom_line(ggplot2::aes(x = gumbel.variate, y = scaled.value, group = variable), colour = 'steelblue4', alpha = .5) +
        ggplot2::geom_point(ggplot2::aes(x = gumbel.variate, y = scaled.value, group = variable), colour = 'grey15', fill = 'steelblue4', alpha = .5, shape = 21) +
        ggplot2::geom_line(data = regional, ggplot2::aes(x = x, y = y), col = 'red4', lwd = .75) + 
        ggplot2::theme_bw() +
        ggplot2::labs(x = '-log(-log(p))', y = 'Value', title = 'Gumbel plot') +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5))
    }
    
    if(method == 'plotly') {
      
      gp <- plotly::ggplotly(gp)
    }
    
    return(gp)
  }
}

growthcurve <- function (model_object, fitted_bootstrap, dist = NULL, method = 'base', ribbon.1 = c(0.05, 0.95), ribbon.2 = c(0.25, 0.75)) {
  
  # res (init res will need to be rewritten for nim)
  
  prbs <- sort(c(ribbon.1, ribbon.2))
  para <- model_object$REG
  
  if(is.null(dist)) {dist <- attr(model_object, 'sim.call')$dist}
  
  qs <- seq(.01, .995, .001)
  qaux <- data.table(rbindlist(lapply(fitted_bootstrap, function(x) {data.frame(q = do.call(paste0('q',dist), list(qs, x$REG)))}),
                               idcol = 'sample'), 
                     probs = seq_along(qs))
  q <- qaux[, .(val = quantile(q, c(.05, .25, .75, .95)),
                q = c('rib_1_min', 'rib_2_min', 'rib_2_max', 'rib_1_max')), 
            by = probs]
  res.gc <- cbind(dcast(q, probs ~ q, value.var = 'val'),
                  data.table(gumbel.variate = -log(-log(qs)),
                             scaled.value = qgpa(qs, para)))
  
  # graphics
  
  if(method == 'base') {
    
    plot(NULL,
         xlim = c(min(res.gc$gumbel.variate), max(res.gc$gumbel.variate)),
         ylim = c(min(res.gc[,c('rib_1_min', 'rib_2_min', 'rib_2_max', 'rib_1_max')]), 
                  max(res.gc[,c('rib_1_min', 'rib_2_min', 'rib_2_max', 'rib_1_max')])),
         bty = 'n',
         xlab = expression(-log(-log(p))),
         ylab = 'Value',
         main = 'Growth curve')
    
    grid()
    
    polygon(c(res.gc$gumbel.variate, rev(res.gc$gumbel.variate)), 
            c(res.gc$rib_1_max, rev(res.gc$rib_1_min)),
            col = '#36648b40', border = NA)
    
    polygon(c(res.gc$gumbel.variate, rev(res.gc$gumbel.variate)), 
            c(res.gc$rib_2_max, rev(res.gc$rib_2_min)),
            col = '#36648b80', border = NA)
    
    lines(res.gc$gumbel.variate,
          res.gc$scaled.value,
          type = 'l',
          col = 'red4',
          lwd = .75)
  }
    
  if(method %in% c('ggplot', 'plotly')) {
    
    if(method == 'ggplot') {
      
      gc <- ggplot2::ggplot(res.gc) +
        ggplot2::geom_ribbon(ggplot2::aes(x = gumbel.variate, ymin = rib_1_min, ymax = rib_1_max), fill = 'steelblue4', alpha = .4) +
        ggplot2::geom_ribbon(ggplot2::aes(x = gumbel.variate, ymin = rib_2_min, ymax = rib_2_max), fill = 'steelblue4', alpha = .8) +
        ggplot2::geom_line(ggplot2::aes(x = gumbel.variate, y = scaled.value), col = 'red4', lwd = .75) + 
        ggplot2::theme_bw() +
        ggplot2::labs(x = '-log(-log(p))', y = 'Value', title = 'Growth curve') +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5))
    }
    
    if(method == 'plotly') {
      
      gc <- plotly::ggplotly(gc)
    }
    
    return(gc)
  }
}

qq <- function(model_object, dist = NULL, method = 'base') {
  
  dta <- as.data.table(model_object$data)
  para <- model_object$REG
  scaling.factor <- model_object$scaling_factor
  
  if(is.null(dist)) {dist <- attr(dta.fit, 'sim.call')$dist}
  
  res.qq <- suppressWarnings(melt(dta))
  res.qq <- data.table(variable = names(dta), sf = scaling.factor)[res.qq, on = c('variable')]
  res.qq <- res.qq[, scaled.value := value/sf]
  
  
  if(method == 'base') {
    
    inipar <- par()
    
    par(pty = 's')
    
    sres.qq <- split(res.qq[!is.na(res.qq$scaled.value),], res.qq$variable[!is.na(res.qq$scaled.value)])
    
    plot(NULL, 
         xlim = c(0, max(res.qq$scaled.value, na.rm = T)*1.15),
         ylim = c(0, max(res.qq$scaled.value, na.rm = T)*1.15),
         pch = 21,
         col = 'grey15', 
         bg = '#36648b90',
         bty = 'n',
         xlab = 'theoretical',
         ylab = 'sample',
         main = 'qqplot')
    
    grid()
    
    lapply(sres.qq, function(x) {
      points(sort(x$scaled.value),
             sort(rgpa(length(x$scaled.value), dta.fit$REG)),       
             pch = 21,
             col = 'grey15', 
             bg = '#36648b90')
    })
    
    abline(0,1, col = 'red4')
    
    suppressWarnings(par(inipar))
  }
  
  if(method %in% c('ggplot', 'plotly')) {
    
    if(method == 'ggplot') {
      
      qq <- ggplot(res.qq) +
        geom_qq(aes(sample = scaled.value, group = variable), geom = 'point', distribution = noquote(paste0('q', dist)), dparams = list(para), colour = 'grey15', fill = 'steelblue4', shape = 21) +
        geom_abline(colour = ('red4')) +
        coord_fixed() +
        lims(x = c(0, max(gpa.qq$value/gpa.qq$l1, na.rm = T)),
             y = c(0, max(gpa.qq$value/gpa.qq$l1, na.rm = T))) +
        theme_bw()
    }
    
    if(method == 'plotly') {
      
      qq <- plotly::ggplotly(qq)
    }
    
    return(qq)
  }
}

ratiodiagram <- function(taus) {
  
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
