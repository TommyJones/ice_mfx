################################################################################
# This script runs the MFX expiriment for simple 3-variable models with one
# interaction term and one quadratic term. The data generating processes are 
# linear and logistic binomial. The models tested are OLS (baseline), logit 
# (baseline), random forest, and a 1-layer neural network.
################################################################################


rm(list = ls())

source("R/00_common_functions.R")

### Generate random data -------------------------------------------------------

set.seed(666)

# Some continuous variables
x1 <- rnorm(n = 1000, mean = 1, sd = 2)
x2 <- rnorm(n = 1000, mean = 1.5, sd = 2.5)
x3 <- rnorm(n = 1000, mean = .6, sd = 1.5)

# Round them to three decimal places to ease computation downstream
x1 <- round(x1, 3)
x2 <- round(x2, 3)
x3 <- round(x3, 3)

z <- 1 + 1.2 * x1 + 4.5 * x2 + 2.8 * x3 + 4.2 * x1 * x3 - 2 * x2 ^ 2

# linear combination
y_linear <- z + rnorm(n = 1000, mean = 0, sd = 1)

# binomial combination
y_logit <- rbinom(10000,1, prob = 1 / (1 + exp(-z)))


# set training and validation data sets
X <- as.data.frame(cbind(x1, x2, x3, y_linear, y_logit))

X$x1sq <- X$x1 ^ 2

X$x2sq <- X$x2 ^ 2

X$x3sq <- X$x3 ^ 2

X$x1x2 <- X$x1 * X$x2

X$x1x3 <- X$x1 * X$x3

X$x2x3 <- X$x2 * X$x3

X_test <- X[ 1:1000 , ]

X <- X[ -1 * 1:1000 , ]

### Get a function to re-prepare interactions for ICE --------------------------
RePrep <- function(newdata) {
  newdata$x1x2 <- newdata$x1 * newdata$x2
  newdata$x1x3 <- newdata$x1 * newdata$x3
  newdata$x2x3 <- newdata$x2 * newdata$x3
  
  newdata$x1sq<- newdata$x1 ^ 2
  newdata$x2sq<- newdata$x2 ^ 2
  newdata$x3sq<- newdata$x3 ^ 2
  
  newdata
}

### Create a random linear model -----------------------------------------------

# correct model
model_linear <- lm(y_linear ~ x1 + x2 + x3 + x1x3 + x2sq, data = X) 

# misspecified model
model_linear_all <- lm(y_linear ~ x1 + x2 + x3 + x1x2 + x1x3 + x2x3 + 
                         x1sq + x2sq + x3sq, data = X)

# predict function
PredictLinear <- function(object, newdata) {
  
  # re-calculating interactions for ICE plots
  # newdata <- RePrep(newdata)
  
  # new model
  predict(object, newdata)
  
}

### Create a random logistic model ---------------------------------------------

# correct model
model_logit <- glm(y_logit ~ x1 + x2 + x3 + x1x3 + x2sq, data = X,
                    family = binomial("logit"), x = TRUE) 
model_logit$mfx <- erer::maBina(w = model_logit, 
                                x.mean = FALSE) # x.mean = FALSE averages mfx across data points

# misspecified model
model_logit_all <- glm(y_logit ~ x1 + x2 + x3 + x1x2 + x1x3 + x2x3 + 
                          x1sq + x2sq + x3sq, data = X,
                        family = binomial("logit"), x = TRUE)

model_logit_all$mfx <- erer::maBina(w = model_logit_all, 
                                    x.mean = FALSE) # x.mean = FALSE averages mfx across data points

# predict function
PredictLogit <- function(object, newdata) {
  
  # re-calculating interactions for ICE plots
  # newdata <- RePrep(newdata)
  
  # new model
  predict(object, newdata, type = "response")
  
}

### Fit random forest models ---------------------------------------------------

rf_linear <- randomForest::randomForest(y_linear ~ x1 + x2 + x3 + x1x2 + x1x3 + x2x3 + 
                                          x1sq + x2sq + x3sq,
                                        data = X)

rf_logit <- randomForest::randomForest(y = factor(X$y_logit),
                                       x = X[ , grep("^x", colnames(X), value = TRUE) ])

PredictRf <- function(object, newdata) {
  
  # re-calculating interactions for ICE plots
  # newdata <- RePrep(newdata)
  
  # new model
  result <- try(predict(object, newdata, type = "vote")[ , 1 ])
  
  if(class(result) == "try-error")
    result <- predict(object, newdata)
  
  result
}


### Fit neural network model ---------------------------------------------------

nn_logit <- nnet::nnet(y_logit ~ x1 + x2 + x3 + x1x2 + x1x3 + x2x3 + 
                         x1sq + x2sq + x3sq,
                       data = X, size = 5,
                       entropy = TRUE)

PredictNn <- function(object, newdata) {
  # re-calculating interactions for ICE plots
  # newdata <- RePrep(newdata)
  
  predict(object, newdata, type = "raw")[ , 1 ]
}


### Get mfx from ICE curves ----------------------------------------------------
mfx_linear <- CalcMfx(object = model_linear, X = X_test, pred_fun = PredictLinear,
                      predictors = c("x1", "x2", "x3", "x1x3", "x2sq"))

mfx_linear_all <- CalcMfx(object = model_linear_all, X = X_test, pred_fun = PredictLinear,
                          predictors = grep("^x", colnames(X_test), value = T))

mfx_logit <- CalcMfx(object = model_logit, X = X_test, pred_fun = PredictLogit,
                     predictors = intersect(names(X_test), colnames(model_logit$x)))

mfx_logit_all <- CalcMfx(object = model_logit_all, X = X_test, pred_fun = PredictLogit,
                          predictors = grep("^x", colnames(X_test), value = T))

mfx_rf_linear <- CalcMfx(object = rf_linear, X = X_test, pred_fun = PredictRf,
                         predictors = grep("^x", colnames(X_test), value = T))

mfx_rf_logit <- CalcMfx(object = rf_logit, X = X_test, pred_fun = PredictRf,
                         predictors = grep("^x", colnames(X_test), value = T))

mfx_nn_logit <- CalcMfx(object = nn_logit, X = X_test, pred_fun = PredictNn,
                        predictors = grep("^x", colnames(X_test), value = T))

save.image("data_derived/one_interaction.RData")

