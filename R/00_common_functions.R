
Derivative <- function(y, x){
  # function estimates derivative of y with respect to x
  # will only return a meaningful result if y is a function only of x
  
  c(NA, diff(y)) / c(NA, diff(x))
  
}

CalcMfx <- function(object, X, pred_fun = predict, predictors = colnames(X), 
                    max_pts = 100, ...){
  
  ### Check consistency of inputs ----
  
  ### Prepare some global variables ----

  ### Do assuming 1 predictor and expand later ----
  if (length(predictors) > 1) {
    
    result <- textmineR::TmParallelApply(predictors, function(p){
      CalcMfx(object = object, X = X, pred_fun = pred_fun, predictors = p,
              max_pts = max_pts)
    }, export = c("object", "X", "pred_fun", "max_pts"), ...)
    
    names(result) <- predictors
    
    class(result) <- "Mfx_list"
    
  } else {
    
    p <- predictors
    
    # get the sequence to iterate over
    if (is.na(max_pts)) {
      
      pts <- sort(unique(X[[ p ]]))
      
    } else {
      
      pts <- seq(min(X[[ p ]]), max(X[[ p ]]), length.out = max_pts)
      
    }
    
    # get predictions for each point in the sequence
    yhat <- lapply(pts, function(point){
      
      X_new <- cbind(X[ , setdiff(names(X), p) ], rep(point, nrow(X)))
      
      colnames(X_new)[ ncol(X_new) ] <- p
      
      yhat <- pred_fun(object, X_new)
      
      yhat
    })
    
    yhat <- do.call(rbind, yhat)
    
    # store a variable used for plot method 
    yh0 <- yhat[ 1 , ]
    
    # get the derivative
    dy <- apply(yhat, 2, function(y) c(NA, diff(y)))
    
    dx <- c(NA, diff(pts))
    
    # get the mfx
    dydx <- dy / dx
    
    mfx <- mean(dydx, na.rm = TRUE)
    
    stdd <- sd(dydx, na.rm = TRUE)
    
    # confidence interval (totally overconfident, does not account for non-linearity)
    conf <- c(mfx - 1.96 * stdd / sqrt(length(dydx)),
              mfx + 1.96 * stdd / sqrt(length(dydx)))
    
    
    # return the result
    result <- list(mfx = mfx, conf = conf,
                   dy = dy,
                   dx = dx,
                   x = pts,
                   yh0 = yh0,
                   varname = p)
    
    class(result) <- "Mfx"
    
  }
  
  return(result)
  
}

plot.Mfx_list <- function(mfx_list, ...) {
  
  # get current par setting
  opar <- par()
  
  # run plots, pausing between
  for (j in seq_along(mfx_list)) {
    plot(mfx_list[[ j ]], ...)
    par(ask = TRUE)
  }
  
  # return par to previous setting
  par(ask = opar$ask)
  
}


plot.Mfx <- function(mfx, type = c("derivative", "level"), centered = TRUE, ...){

  ### Below still only works with one variable
  
  ### check consistency of inputs ----
  # TODO
  
  ### Do calculations to get the chosen plot ----
  if (type[ 1 ] == "level") {
    plotmat <- rbind(mfx$yh0, mfx$dy[ -1 , ])
    
    plotmat <- apply(plotmat, 2, cumsum)
    
    mfx_pred <- mfx$x * mfx$mfx + 
      mean(plotmat[ which.min(abs(mfx$x)) , ], na.rm = TRUE) # I'm not sure about this constant
    
    ylab <- "partial y-hat"
    
    if (centered) {
      
      plotmat <- t(t(plotmat) - plotmat[ 1 , ])
      
      mfx_pred <- mfx_pred - mfx_pred[ 1 ]
      
      ylab <- paste(ylab, "(centered)")
      
    }
    
  } else {
    
    plotmat <- mfx$dy / mfx$dx
    
    mfx_pred <- rep(mfx$mfx, length(mfx$x))
    
    ylab <- "derivative of partial y-hat"
  }
  
  # get our plot options straight
  dots <- list(...)
  
  if ( ! "ylim" %in% names(dots)) {
    dots$ylim <- range(plotmat, na.rm = TRUE)
  }
  
  if (! "col" %in% names(dots)) {
    dots$col <- rgb(0,0,0,0.3)
  }
  
  if (! "ylab" %in% names(dots)) {
    dots$ylab <- ylab
  }
  
  if (! "xlab" %in% names(dots)) {
    dots$xlab <- mfx$varname
  }
  
  if ("type" %in% names(dots)) {
    warning("reverting type to type = 'l'")
  }
  
  dots$type <- "l"
  
  ### do the plotting ----
  
  # plot the ICE curves
  dots$x <- mfx$x
  dots$y <- plotmat[ , 1 ]
  
  do.call(plot, dots)
  
  for (j in 2:ncol(plotmat)) {
    dots$y <- plotmat[ , j ]
    do.call(lines, dots)
  }
  
  # plot the curve predicted by the mfx
  dots$y <- mfx_pred
  dots$col <- "red"
  dots$lty = 1
  dots$lwd = 3
  
  do.call(lines, dots)
  
  
  # plot the average ICE curve
  dots$y <- rowMeans(plotmat)
  dots$col <- "blue"
  dots$lty = 2
  dots$lwd = 3
  
  do.call(lines, dots)
  
}

summary.Mfx_list <- function(mfx_list){
  tab <- data.frame(variable = names(mfx_list),
                    effect = sapply(mfx_list, function(x) x$mfx))
  
  conf <- do.call(rbind, lapply(mfx_list, function(x) x$conf))
  
  colnames(conf) <- c("95_conf_low", "95_conf_high")
  
  tab <- cbind(tab, conf)
  
  print(tab)
}

print.Mfx_list <- function(mfx_list){
  summary(mfx_list)
}
