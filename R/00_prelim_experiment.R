################################################################################
# This script runs some preliminary expiriments getting mfx from ICE curves
#
# If (reasonably) successful, I'll copy it to Rmarkdown for the paper/real stuff
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

z <- 1 + 1.2 * x1 + 4.5 * x2 + 2.8 * x3 # + 4.2 * x1 * x3 

# linear combination
y_linear <- z + rnorm(n = 1000, mean = 0, sd = 1)

# binomial combination
y_logit <- rbinom(10000,1, prob = 1 / (1 + exp(-z)))

X <- as.data.frame(cbind(x1, x2, x3, y_linear, y_logit))

X_test <- X[ 1:1000 , ]

X <- X[ -1 * 1:1000 , ]

### Create a random linear model -----------------------------------------------

y_hat_linear <- lm(y_linear ~ x1 + x2 + x3 , data = X) # + I(x1 * x3)

### Create a random logistic model ---------------------------------------------

y_hat_logit <- glm(y_logit ~ x1 + x2 + x3, data = X, # + I(x1 * x3) 
                   family = binomial("logit"), x = TRUE)


### Get true marginal effects for logistic from erer ---------------------------
y_hat_logit$mfx <- erer::maBina(w = y_hat_logit, 
                                x.mean = FALSE) # x.mean = FALSE averages mfx across data points


### Get mfx from ICE curves ----------------------------------------------------
PredictLogistic <- function(object, newdata){
  predict(object, newdata, type = "response")
}

mfx_logit <- CalcMfx(object = y_hat_logit, 
                     X = X_test,
                     pred_fun = PredictLogistic, 
                     predictors = c("x1", "x2", "x3"))

mfx_linear <- CalcMfx(object = y_hat_linear, 
                      X = X_test,
                      pred_fun = predict, 
                      predictors = c("x1", "x2", "x3"))


# Note: May fail in the face of interactions, or the I(x1 * x3) fails
# Either way, more work is necessary

### Create Random Forest logistic and linear models ----------------------------
rf_linear <- randomForest::randomForest(y_linear ~ x1 + x2 + x3,
                                        data = X)

rf_logit <- randomForest::randomForest(y = factor(X$y_logit),
                                       x = X[ , c("x1", "x2", "x3") ])

### Random Forest ICE curves ---------------------------------------------------
PredictRfBinom <- function(object, newdata){
  predict(object, newdata, type = "prob")[ , "1" ]
}

mfx_rf_logit <- CalcMfx(object = rf_logit, 
                        X = X_test,
                        pred_fun = PredictRfBinom,
                        predictors = c("x1", "x2", "x3"))

mfx_rf_linear <- CalcMfx(object = rf_linear, 
                         X = X_test,
                         pred_fun = predict,
                         predictors = c("x1", "x2", "x3"))

save.image("data_derived/prelim_experiment.RData")

# Notes on closing: 
# * Can use standard prediction from logistic to recover coefficients
# * Not sure how to get that with random forest though...
# * This doesn't cover interactions: you should handle those explicitly
# * There are many different kinds of regression you should include (e.g. poisson)
# * You'll also need a wider range of simulations to cover the effects/properties
# * You should be getting ICE curves on held-out data for all models
