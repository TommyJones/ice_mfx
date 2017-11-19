
rm(list = ls())

source("R/00_common_functions.R")

### Load wine quality data ----

url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/"

red <- read.delim(paste0(url, "winequality-red.csv"), stringsAsFactors = FALSE, sep = ";")
white <- read.delim(paste0(url, "winequality-white.csv"), stringsAsFactors = FALSE, sep = ";")

wine <- rbind(red, white)

wine <- wine[ sample(1:nrow(wine), nrow(wine)) , ] # shuffle rows


### Declare training and predicting functions -------------

# lasso
FitLasso <- function(y, x, family = "binomial"){
  glmnet::glmnet(y = y, x = as.matrix(x), family = family)
}

PredictLasso <- function(object, newdata, type = "response"){
  out <- predict(object, as.matrix(newdata), type = type)
  out[ , ncol(out)]
}

# random forest
FitRf <- function(y, x){
  randomForest::randomForest(y = y, x = x)
}

PredictRf <- function(object, newdata) {

  # new model
  result <- try(predict(object, newdata, type = "vote")[ , "1" ],
                silent = TRUE)
  
  if(class(result) == "try-error")
    result <- predict(object, newdata)
  
  result
}

# deep neural net
library(magrittr)
library("keras")

FitNn <- function(y, x, activation = "relu", epochs = 50) {
  
  x <- as.matrix(x)
  y <- matrix(y, ncol = 1)
  
  net <- keras_model_sequential()
  
  net %>% 
    layer_dense(units = 36, activation = activation, input_shape = ncol(x)) %>% 
    layer_dropout(rate = 0.4) %>% 
    # layer_dense(units = 27, activation = activation) %>% 
    # layer_dropout(rate = 0.4) %>% 
    layer_dense(units = 20, activation = activation) %>%
    layer_dropout(rate = 0.4) %>%
    layer_dense(units = 15, activation = activation) %>% 
    layer_dropout(rate = 0.4) %>% 
    # layer_dense(units = 11, activation = activation) %>% 
    # layer_dropout(rate = 0.4) %>% 
    layer_dense(units = 8, activation = activation) %>%
    layer_dropout(rate = 0.4) %>%
    # layer_dense(units = 6, activation = activation) %>% 
    # layer_dropout(rate = 0.4) %>% 
    # layer_dense(units = 4, activation = activation) %>% 
    # layer_dropout(rate = 0.4) %>% 
    layer_dense(units = ncol(y), activation = "sigmoid")
  
  # summary(net)
  
  net %>% compile(
    loss = 'binary_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- net %>% fit(
    x = x, y = y, 
    epochs = epochs, batch_size = 128, 
    validation_split = 0.2
  )
  
  net
}

PredictNn <- function(object, newdata) {
  # re-calculating interactions for ICE plots
  # newdata <- RePrep(newdata)
  
  as.numeric(predict(object, as.matrix(newdata)))
}

### training and test set -------------------
wine$quality <- as.character(wine$quality)

cv <- SampleKFolds(data = wine, k = 5, strat_var = "quality")

wine$quality <- as.numeric(wine$quality)

training <- wine[ cv[[ 1 ]]$training, ] 

test <- wine[ cv[[ 1 ]]$test , ]


### Fit models and get mfx ---------

lasso <- FitLasso(y = training$quality, 
                  x = training[ , setdiff(names(training), "quality") ],
                  family = "gaussian")

lasso_mfx <- CalcMfx(object = lasso, 
                     X = test[ , setdiff(names(test), "quality") ],
                     pred_fun = PredictLasso,
                     predictors = setdiff(names(test), "quality"))

lasso$mfx <- summary(lasso_mfx, print = F)

rf <- FitRf(y = training$quality, 
            x = training[ , setdiff(names(training), "quality") ])

rf_mfx <- CalcMfx(object = rf, 
                  X = test[ , setdiff(names(test), "quality") ],
                  pred_fun = PredictRf,
                  predictors = setdiff(names(test), "quality"))

rf$mfx <- summary(rf_mfx, print = F)