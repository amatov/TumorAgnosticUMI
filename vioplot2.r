#This is a custom version of pcks vioplot to allow list of vectors and multiple
#colors.
df <- data.frame(id=sample(c('a','b'), 100, replace=T), val = rnorm(100, 5, 20))
x_list <- split(df$val,df$id)
aux1 <- pdata[,1]
aux1<-aux1[!is.na(aux1)]
aux2 <- pdata[,2]
aux2<-aux2[!is.na(aux2)]
aux3 <- pdata[,3]
aux3<-aux3[!is.na(aux3)]
aux4 <- pdata[,4]
aux4<-aux4[!is.na(aux4)]
aux5 <- pdata[,5]
aux5<-aux5[!is.na(aux5)]
aux6 <- pdata[,6]
aux6<-aux6[!is.na(aux6)]
aux7 <- pdata[,7]
aux7<-aux7[!is.na(aux7)]
aux8 <- pdata[,8]
aux8<-aux8[!is.na(aux8)]
aux9 <- pdata[,9]
aux9<-aux9[!is.na(aux9)]
aux10 <- pdata[,10]
aux10<-aux10[!is.na(aux10)]
aux11 <- pdata[,11]
aux11<-aux11[!is.na(aux11)]
aux12 <- pdata[,12]
aux12<-aux12[!is.na(aux12)]
aux13 <- pdata[,13]
aux13<-aux13[!is.na(aux13)]
aux14 <- pdata[,14]
aux14<-aux14[!is.na(aux14)]
aux15 <- pdata[,15]
aux15<-aux15[!is.na(aux15)]
aux16 <- pdata[,16]
aux16<-aux16[!is.na(aux16)]
x_list$a <- aux1
x_list$b <- aux2
x_list$c <- aux3
x_list$d <- aux4
x_list$e <- aux5
x_list$f <- aux6
x_list$g <- aux7
x_list$h <- aux8
x_list$p <- aux9
x_list$q <- aux10
x_list$k <- aux11
x_list$l <- aux12
x_list$m <- aux13
x_list$n <- aux14
x_list$r <- aux15
x_list$s <- aux16
vioplot2(x = x_list, col=c("red", "blue", "green","red", "blue", "green","red", "blue", "green","red", "blue", "green"))

vioplot2 <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                      horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                      lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                      at, add = FALSE, wex = 1, drawRect = TRUE){
  
  # EXAMPLE:
  # df <- data.frame(id=sample(c('a','b'), 100, replace=T), val = rnorm(100, 5, 20))
  # x_list <- split(df$val,df$id)
  # vioplot2(x = x_list, col=c("red", "blue"))
  
  #To load dependencies easy just to load vioplot
  require(vioplot)
  
  #Change from vioplot so we can accept lists of vectors
  if(!is.list(x)){
    datas <- list(x, ...)
  } else{
    datas<-x
  }
  
  
  n <- length(datas)
  
  
  #Change from vioplot so we can plot multiple colors per violin. In polygons
  #below we cycle through color by col[i]
  if(length(col)==1) col <- rep(col, n)
  
  
  
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}
