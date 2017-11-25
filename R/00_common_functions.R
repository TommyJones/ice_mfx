
source("R/Mfx_functions.R")

# Declare a function to calculate AUC
Trapezoid <- function(x,y){
  len = length(x)
  w <- x[2:len] - x[1:(len-1)]
  h <- (y[2:len] + y[1:(len-1)])/2
  sum(h*w)
}

CalcClassificationStats <- function(predicted_probabilities, true_values){
  
  # error checking
  if( length(unique(true_values)) != 2){
    stop("true_values must be a vector with two levels. Currently, 0/1 and TRUE/FALSE are supported")
  }
  
  # put true_values into a 0/1 vector
  true_values <- as.numeric(true_values)
  
  # get a list of unique thresholds to the thousandths place
  thresholds <- seq(0, 1, 0.01)
  
  # get TP, TN, FP, FN & calc other stats
  result <- lapply(thresholds, function(T){
    preds <- as.numeric(predicted_probabilities >= T)
    
    tp <- sum(preds == 1 & true_values == 1)
    fp <- sum(preds == 1 & true_values == 0)
    tn <- sum(preds == 0 & true_values == 0)
    fn <- sum(preds == 0 & true_values == 1)
    
    sens_recall <- tp / (tp + fn)
    spec_tnr <- tn / (fp + tn)
    prec_ppv <- tp / (tp + fp)
    npv <- tn / (tn + fn)
    fpr <- fp / (fp + tn)
    fdr <- fp / (fp + tp)
    fnr <- fn / (fn + tp)
    acc <- (tp + tn) / (tp + tn + fp + fn)
    f1 <- 2 * tp / (2 * tp + fp + fn)
    # mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    
    data.frame(threshold=T, TP=tp, FP=fp, TN=tn, FN=fn,
               Sens_Recall=sens_recall,
               Spec_TNR=spec_tnr,
               Prec_PPV=prec_ppv,
               NPV=npv,
               FPR=fpr,
               FDR=fdr,
               FNR=fnr,
               Accuracy=acc,
               F1=f1,
               # MCC=mcc,
               stringsAsFactors=FALSE)
  })
  
  result <- do.call(rbind, result)
  
  roc <- data.frame(xOneMinusSpec=1 - result$Spec_TNR,
                    ySensitivity=result$Sens_Recall,
                    stringsAsFactors=FALSE)
  auc <- Trapezoid(x=1 - roc$xOneMinusSpec, y=roc$ySensitivity)
  
  
  result <- list(stats=result, roc=roc, auc=auc)
  
  return(result)
}


#' Partition a dataset into K folds for cross-validation
#' @description This function takes a dataset as input and returns a list with
#'  indices for observations in each fold and the dataset itself.
#'
#' @param data A matrix or data frame of the data to be partitioned
#' @param k An iteger indicating the number of folds
#' @param strat_var A character or numeric indexing the column of 
#'  \code{data} on which to stratify. Note that \code{data[[ strat_var ]]} must
#'  be categorical. i.e. It must be a character or factor.
#' @param seed An integer to set the random seed if desired. Defaults to \code{NULL}.
#'
#' @return
#' Returns a \code{list} of length K. Each element of this list has
#' indices for training and testing at each fold. For example, \code{training}
#' contains indices for K - 1 folds for training a model and \code{test} contains
#' indices for the K-th fold for testing.
#' @export
#' @examples
#' data(mtcars)
#'
#' folds <- CreateKFolds(data = mtcars, k = 10)
SampleKFolds <- function(data, k, strat_var = NULL, seed = NULL) {
  
  # set seed if desired / not already set
  if( ! is.null(seed) ) set.seed(seed)
  
  # if data doesn't have rownames, give it to them
  if (is.null(rownames(data))) {
    rownames(data) <- 1:nrow(data)
  }
  
  # declare a cv function
  MakeCv <- function(data, k) {
    # set up basics
    rows <- rownames(data)
    
    size <- floor(nrow(data)/k)
    
    remainder <- nrow(data) - size * k # make each fold as balanced as possible
    
    # get the indexes of each fold into a list
    folds <- vector(mode="list", length=k)
    
    for(j in 1:k){
      folds[[ j ]] <- sample(rows, size=size, replace=FALSE)
      rows <- rows[ ! rows %in% folds[[ j ]] ]
    }
    
    if( remainder > 0){
      for(j in 1:length(rows)){
        folds[[ j ]][ length(folds[[ j ]]) + 1 ] <- sample(rows, size=1)
        rows <- rows[ ! rows %in% folds[[ j ]] ]
      }
    }
    
    # combinations of each fold to use for training sets / test sets
    combos <- combn(x=1:k, m=(k-1))
    
    # pull indices into format that's easier to use
    cv.folds <- vector(mode="list", length=k)
    
    for( j in 1:k ){
      cv.folds[[ j ]]$training <- unlist( folds[ combos[ , j ] ] )
      cv.folds[[ j ]]$test <- folds[[ (1:k)[ ! (1:k) %in% combos[ , j ] ]  ]]
    }
    
    result <- cv.folds
    return(result)
  }
  
  # Check whether or not we're stratifying. If yes, do it. If no, just do CV on whole.
  if ( ! is.null(strat_var) ) {
    
    # is the input correctly formatted?
    if ( ! class(data[[ strat_var ]]) %in% c("factor", "character"))
      stop("class(data[[ strat_var ]] must be one of factor or character to stratify")
    
    # stratify then sample and re-combine
    # Note: this could be WAY more efficient, I'm sure
    result <- by(data, INDICES = data[[ strat_var ]], function(x){
      MakeCv(data = x, k = k)
    })
    
    result <- lapply(seq_len(k), function(fold){
      
      training <- unlist(lapply(result, function(strat){
        strat[[ fold ]]$training
      }))
      
      test <- unlist(lapply(result, function(strat){
        strat[[ fold ]]$test
      }))
      
      list(training = training, test = test)
    })
  } else {
    result <- MakeCv(data = data, k = k)
  }
  
  return(result)
  
}

Factor2Binary <- function(y){
  ynames <- levels(y)
  
  Y <- data.frame(y = y)
  
  Y <- model.matrix(~y, Y)
  
  colnames(Y) <- ynames
  
  # model.matrix always makes the first column all "1", fix this
  Y[ , 1 ] <- as.numeric(rowSums(Y[ , -1 ]) == 0)
  
  
  Y
}

Standardize <- function(x, sd_x = NULL, mean_x = NULL){
  
  if(is.null(sd_x) | is.null(mean_x)){
    # remove any columns for which there is no variation
    sd_x <- apply(x, 2, function(y) sd(y, na.rm = TRUE))
    
    x <- x[ , sd_x > 0 ]
    
    sd_x <- sd_x[ sd_x > 0 ]
    
    # get means to store
    mean_x <- colMeans(x, na.rm = TRUE)
    
    # standardize 
    result <- apply(x, 2, function(y){
      (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
    })
    
    # return result
    return(list(data = result, 
                sd_x = sd_x,
                mean_x = mean_x))
  } else {
    
    result <- mapply(function(d, m, s){
      (d - m) / s
    }, d = as.list(x[ , names(mean_x) ]), m = mean_x, s = sd_x,
    SIMPLIFY = FALSE)
    
    return(as.data.frame(result, stringsAsFactors = FALSE))
    
  }

}

CalcR2 <- function(y, yhat) {
  
  sse <- sum((y - yhat) ^ 2)
  
  sst <- sum((y - mean(y)) ^ 2)
  
  1 - sse / sst
  
}


### Recyclable functions for fitting and predicting things ---------------------
# lasso
FitLasso <- function(y, x, family = "binomial", ...){
  glmnet::glmnet(y = y, x = as.matrix(x), family = family)
}

PredictLasso <- function(object, newdata, type = "response"){
  out <- predict(object, as.matrix(newdata), type = type)
  out[ , ncol(out) ]
}

# random forest
FitRf <- function(y, x){
  # will automatically do regression or classification depending on format of y
  randomForest::randomForest(y = y, x = x)
}

PredictRf <- function(object, newdata) {
  
  # will automatically predict linear or logit probabilities depending on format
  # of object
  result <- try(predict(object, newdata, type = "vote")[ , "1" ],
                silent = TRUE)
  
  if(class(result) == "try-error")
    result <- predict(object, newdata)
  
  result
}
