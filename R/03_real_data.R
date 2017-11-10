################################################################################
# This script calculates various models for crash data and gets MFX on those
# models. The models are binomial logistic, random forest, and deep neural net
################################################################################

rm(list = ls())

source("R/00_common_functions.R")

load("data_raw/csdata_formatted_v2.RData")

### Do some pre-curation of the data -------------------------------------------

crash_characteristics <- c("bmode", "delta_v", "maxc1", "intrus_p", 
                           "multi", "pdo1", "shl1", "dvl", "dvd")

vehicle_characteristics <- c("curbwgt_kg", "btype")

occupant_characteristics <- c("belt_use", "age", "height", "weight", "bmi", 
                              "seatrack", "seatpos")

# quadratic terms
csdata$delta_v_sq <- csdata$delta_v ^ 2 

csdata$maxc1_sq <- csdata$maxc1 ^ 2

csdata$intrus_p_sq <- csdata$intrus_p ^ 2

csdata$pdo1_sq <- csdata$pdo1 ^ 2

csdata$dvl_sq <- csdata$dvl ^ 2

csdata$dvd_sq <- csdata$dvd ^ 2

sq_terms <- grep("_sq$", names(csdata), value = TRUE)

# interaction terms
csdata$i_dv_maxc1 <- csdata$delta_v * csdata$maxc1

csdata$i_dv_intrusp <- csdata$delta_v * csdata$intrus_p

csdata$i_dv_pdo1 <- csdata$delta_v * csdata$pdo1

csdata$i_dv_dvl <- csdata$delta_v * csdata$dvl

csdata$i_dv_dvd <- csdata$delta_v * csdata$dvd

csdata$i_dvl_dvd <- as.numeric(csdata$dvl * csdata$dvd)

i_terms <- grep("^i_", names(csdata), value = TRUE)

csdata$mais3pl <- as.numeric(csdata$mais3pl)

# keep only the columns we need moving forward
csdata <- csdata[ , c("mais3pl", "ratwgt2", 
                      crash_characteristics, 
                      vehicle_characteristics,
                      occupant_characteristics,
                      i_terms, sq_terms) ]

### Get cross-validation indices -----------------------------------------------

# stratifying by crash mode
cv <- SampleKFolds(data = csdata, k = 10, strat_var = "bmode", seed = 666)

### Re-code all factor variables as binary -------------------------------------

f_vars <- sapply(csdata, class)

f_vars <- names(f_vars)[ f_vars == "factor" ]

names(f_vars) <- f_vars

f_vars <- lapply(f_vars, function(x) Factor2Binary(y = csdata[[ x ]]))

for (j in seq_along(f_vars)) 
  colnames(f_vars[[ j ]]) <- paste(names(f_vars)[ j ], colnames(f_vars[[ j ]]), sep = "_")

# take out those factor variables
csdata <- csdata[ , setdiff(names(csdata), names(f_vars)) ]

# append binary variables
csdata <- cbind(csdata, do.call(cbind, f_vars))

### Do k-fold cross validation to get baseline model accuracies ----------------

# resampling procedure to make case weights like that of population
Resample <- function(vec, weights){
  sample(vec, length(vec), 
         prob = weights, 
         replace = TRUE)
}

# logistic regression
PredictLogit <- function(object, newdata) {
  
  # re-calculating interactions for ICE plots
  # newdata <- RePrep(newdata)
  
  # new model
  predict(object, newdata, type = "response")
  
}

model_logit <- parallel::mclapply(cv, function(x){
  samp <- Resample(x$training, csdata[ x$training, "ratwgt2" ])
  
  X <- csdata[ samp , setdiff(names(csdata), c("mais3pl", "ratwgt2")) ]
  y <- csdata[ samp, "mais3pl" ]
  
  f <- as.formula(paste("mais3pl ~ -1 +", paste(colnames(X), collapse = " + ")))
  
  model <- glm(f, data = cbind(mais3pl = y, X), x = TRUE, family = binomial("logit"))
  
  samp <- Resample(x$test, csdata[ x$test, "ratwgt2" ])
  
  p <- predict(model, newdata = csdata[ samp , colnames(X) ], type = "response")
  
  result <- CalcClassificationStats(predicted_probabilities = p, 
                                    true_values = csdata[ samp , "mais3pl" ])
  
}, mc.cores = 4)




# random forest
PredictRf <- function(object, newdata) {
  
  # re-calculating interactions for ICE plots
  # newdata <- RePrep(newdata)
  
  # new model
  result <- try(predict(object, newdata, type = "vote")[ , "1" ])
  
  if(class(result) == "try-error")
    result <- predict(object, newdata)
  
  result
}

rf_logit <- parallel::mclapply(cv, function(x){
  samp <- Resample(x$training, csdata[ x$training, "ratwgt2" ])
  
  X <- csdata[ samp , setdiff(names(csdata), c("mais3pl", "ratwgt2")) ]
  y <- csdata[ samp, "mais3pl" ]
  
  model <- randomForest::randomForest(y = factor(y), x = X)
  
  samp <- Resample(x$test, csdata[ x$test, "ratwgt2" ])
  
  p <- PredictRf(model, csdata[ samp , colnames(X) ])
  
  result <- CalcClassificationStats(predicted_probabilities = p, 
                                    true_values = csdata[ samp , "mais3pl" ])
  
}, mc.cores = 4)

# deep neural net
library(magrittr)
library("keras")

FitNn <- function(y, x, activation = "sigmoid", epochs = 50) {
  
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
    layer_dense(units = 4, activation = activation) %>% 
    layer_dropout(rate = 0.4) %>% 
    layer_dense(units = ncol(y), activation = "sigmoid")
  
  # summary(net)
  
  net %>% compile(
    loss = 'mean_squared_error',
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

nn_logit <- lapply(cv, function(x){
  samp <- Resample(x$training, csdata[ x$training, "ratwgt2" ])
  
  X <- csdata[ samp , setdiff(names(csdata), c("mais3pl", "ratwgt2")) ]
  y <- csdata[ samp, "mais3pl" ]
  
  model <- FitNn(y = y, x = X, activation = "sigmoid", epochs = 30)
  
  samp <- Resample(x$test, csdata[ x$test, "ratwgt2" ])

    p <- PredictNn(model, csdata[ samp , colnames(X) ])
  
  result <- CalcClassificationStats(predicted_probabilities = p,
                                    true_values = csdata[ samp , "mais3pl" ])
  
})

# aggregate results into a list

### Get final models for each type ---------------------------------------------

# get training rows

