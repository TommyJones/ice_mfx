#' Calculate marginal effects and ICE curves
#' @description calculate marginal effects of features on predicted outcomes 
#' from an arbitrary prediction model. Also calculates the necessary components 
#' to construct individual conditional expectation (ICE) curves.
#'
#' @param object a model object 
#' @param X an object of class \code{data.frame} containing data from which 
#' marginal effects should be calculated
#' @param pred_fun a function that accepts two arguments corresponding to 
#' \code{object} and \code{X}, above. The function should return a vector of 
#' predicted responses. If this argument is not passed, the generic 
#' \code{predict} function will be tried.
#' @param predictors a character vector of column names of \code{X} for which
#' marginal effects are to be calculated
#' @param max_pts integer of the maximum number of points to be used for 
#' calculating ICE curves. If a predictor has fewer than \code{max_pts} unique
#' values, then only the predictor's unique values will be used.
#' @param dydx_mean logical indicating whether to calculate the marginal effect
#' based on the mean of the derivative of ICE curves. If \code{FALSE}, then
#' the marginal effect is calculated at each actual observation then averaged.
#' @param ... arguments to be passed on to \code{textmineR}'s function
#' code{\link[textmineR]{TmParallelApply}}.
#'
#' @return
#' Returns a \code{list} of class \code{Mfx} or \code{Mfx_list} depending on
#' whether the \code{predictors} argument is of length 1 or greater. If length 1,
#' then the result is of class \code{Mfx}. If greater than 1, the result is of
#' class \code{Mfx_list}, each element of which is of class \code{Mfx}. An object
#' of class \code{Mfx} has the following slots:
#' 
#' 
#' @examples
#' data(mtcars)
#'
#' folds <- CreateKFolds(data = mtcars, k = 10)
#' @export
CalcMfx <- function(object, X, pred_fun = predict, predictors = colnames(X), 
                    max_pts = 100, dydx_mean = FALSE, ...){
  
  ### Check consistency of inputs ----
  # TODO
  
  
  ### Gets to calculating! ----
  if (length(predictors) > 1) {
    
    result <- textmineR::TmParallelApply(predictors, function(p){
      CalcMfx(object = object, X = X, pred_fun = pred_fun, predictors = p,
              max_pts = max_pts, dydx_mean = dydx_mean)
    }, export = c("object", "X", "pred_fun", "max_pts", "dydx_mean"), ...)
    
    names(result) <- predictors
    
    class(result) <- "Mfx_list"
    
  } else {
    
    p <- predictors
    
    # get the sequence to iterate over
    if (is.factor(X[[ p ]])) {
      
      pts <- levels(X[[ p ]])
      
    } else if (is.na(max_pts)) {
      
      pts <- sort(unique(X[[ p ]]))
      
    } else {
      
      pts <- seq(min(X[[ p ]]), max(X[[ p ]]), length.out = max_pts)
      
      # add a lower bound so we can get non-infinate values for
      # delta x and delta y at the bottom
      pts <- c(pts[ 1 ] - mean(diff(pts)), pts)
    }
    
    # get predictions for each point in the sequence
    yhat <- lapply(pts, function(point){
      
      X_new <- X
      
      X_new[[ p ]] <- point
      
      yhat <- pred_fun(object, X_new)
      
      yhat
    })
    
    yhat <- do.call(rbind, yhat)
    
    # store a variable used for plot method 
    yh0 <- yhat[ 1 , ]
    
    # get the derivative
    if (is.factor(X[[ p ]])) {
      
      warning("Factor detected. CalcMfx for factors is in prototype. 
              I have reason to believe it doesn't calculate the right thing.
              Caveat emptor!")
      
      ### More research is needed to get factors right...
      # since each level of p represents a binary variable, I think I can just
      # take rowMeans of yhat to get the marginal effect. But I'm probably wrong
      
      # dy <- yhat
      
      dy <- apply(yhat, 2, function(y){
        ans <- numeric(length(y))
        
        for (j in seq_along(ans)) {
          ans[ j ] <- y[ j ] - mean(y[ -j ])
        }
        ans
      })
      
      dx <- rep(1, length(pts))
      
      dydx <- dy
      
      mfx <- rowMeans(dydx, na.rm = TRUE)
      
      stdd <- apply(dydx, 1, function(y) sd(y, na.rm = TRUE))
      
      # confidence interval (totally overconfident, does not account for non-linearity)
      conf <- cbind(mfx - 1.96 * stdd / sqrt(nrow(dydx)),
                    mfx + 1.96 * stdd / sqrt(nrow(dydx)))
      
      names(mfx) <- levels(X[[ p ]])
      rownames(conf) <- names(mfx)
      
      is_factor <- TRUE
      
    } else {
      # dy/dx of curves
      dy <- apply(yhat, 2, function(y) c(NA, diff(y)))
      
      dx <- c(NA, diff(pts))
      
      dydx <- dy / dx
      
      # dy/dx of true values
      small_x <- sapply(X[[ p ]], function(x) max(pts[ pts < x ], na.rm = T))
      
      dx_true <- X[[ p ]] - small_x
      
      X_new <- X
      X_new[[ p ]] <- small_x
      
      yhat_true <- pred_fun(object, X)
      
      dy_true <- yhat_true - pred_fun(object, X_new)
      
      dydx_true <- dy_true / dx_true
      
      # get the mfx
      if (dydx_mean) {
        # If you want to use mean values of the curves
        mfx <- mean(dydx, na.rm = TRUE)
        
        stdd <- sd(dydx, na.rm = TRUE)
        
        se <- stdd / sqrt(length(dydx))
        
        # confidence interval (totally overconfident, does not account for non-linearity)
        conf <- c(mfx - 1.96 * se,
                  mfx + 1.96 * se)
      } else {
        # If you only want to consider values at actual data points
        # dy/dx of true values
        
        mfx <- mean(dydx_true, na.rm = TRUE)
        
        stdd <- sd(dydx_true, na.rm = TRUE)
        
        se <- stdd / sqrt(length(dydx_true))
        
        # confidence interval (totally overconfident, does not account for non-linearity)
        conf <- c(mfx - 1.96 * se,
                  mfx + 1.96 * se)
      }
      
      
      
      is_factor <- FALSE
    }
    
    # return the result
    result <- list(mfx = mfx, 
                   se = se,
                   conf = conf,
                   dy = dy,
                   dx = dx,
                   x = pts,
                   yh0 = yh0,
                   true_values = data.frame(x = X[[ p ]],
                                            dx = dx_true,
                                            y = yhat_true,
                                            dy = dy_true,
                                            stringsAsFactors = FALSE),
                   varname = p,
                   is_factor = is_factor)
    
    class(result) <- "Mfx"
    
  }
  
  return(result)
  
}