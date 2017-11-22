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

# offset from centerline variables measured in absolute value
csdata$dvd <- abs(csdata$dvd)

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

csdata$i_dv_beltuse <- csdata$delta_v * csdata$belt_use

csdata$i_beltuse_pdo1 <- csdata$belt_use * csdata$pdo1

csdata$i_beltuse_dvd <- csdata$belt_use * csdata$dvd

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

FitLasso <- function(y, x){
  glmnet::glmnet(y = y, x = as.matrix(x), family = "binomial")
}

PredictLasso <- function(object, newdata){
  out <- predict(object, as.matrix(newdata), type = "response")
  out[ , ncol(out)]
}

model_lasso <- parallel::mclapply(cv, function(x){
  samp <- Resample(x$training, csdata[ x$training, "ratwgt2" ])
  
  X <- csdata[ samp , setdiff(names(csdata), c("mais3pl", "ratwgt2")) ]
  y <- csdata[ samp, "mais3pl" ]
  
  model <- FitLasso(y = y, x = X)
  
  samp <- Resample(x$test, csdata[ x$test, "ratwgt2" ])
  
  p <- PredictLasso(object = model, newdata = csdata[ samp, colnames(X) ])
  
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

model_rf <- parallel::mclapply(cv, function(x){
  samp <- Resample(x$training, csdata[ x$training, "ratwgt2" ])
  
  X <- csdata[ samp , setdiff(names(csdata), c("mais3pl", "ratwgt2")) ]
  y <- csdata[ samp, "mais3pl" ]
  
  model <- randomForest::randomForest(y = factor(y), x = X)
  
  samp <- Resample(x$test, csdata[ x$test, "ratwgt2" ])
  
  p <- PredictRf(model, csdata[ samp , colnames(X) ])
  
  result <- CalcClassificationStats(predicted_probabilities = p, 
                                    true_values = csdata[ samp , "mais3pl" ])
  
}, mc.cores = 4)


# aggregate results into a list
model_accuracy <- list(lasso = model_lasso,
                       rf = model_rf)

### Get final models for each type ---------------------------------------------

# get training and testing sets
set.seed(666)

test_rows <- sample(1:nrow(csdata), round(.2 * nrow(csdata)))

X_training <- csdata[ -test_rows , ]

X_test <- csdata[ test_rows , ]

# resample
samp <- Resample(rownames(X_training), weights = X_training$ratwgt2)

X_training <- X_training[ samp , ]

samp <- Resample(rownames(X_test), weights = X_test$ratwgt2)

X_test <- X_test[ samp , ]

# logistic 
model_lasso <- FitLasso(y = X_training$mais3pl, 
                        x = X_training[ , setdiff(names(X_training), c("mais3pl", "ratwgt2"))])

mfx_lasso <- CalcMfx(object = model_lasso, 
                     X = X_test[ , setdiff(names(X_test), c("mais3pl", "ratwgt2")) ],
                     pred_fun = PredictLasso,
                     predictors = setdiff(names(X_test), c("mais3pl", "ratwgt2")))

# random forest
model_rf <- randomForest::randomForest(y = as.factor(X_training$mais3pl), 
                        x = X_training[ , setdiff(names(X_training), c("mais3pl", "ratwgt2"))])

mfx_rf <- CalcMfx(object = model_rf, 
                  X = X_test[ , setdiff(names(X_test), c("mais3pl", "ratwgt2")) ],
                  pred_fun = PredictRf,
                  predictors = setdiff(names(X_test), c("mais3pl", "ratwgt2")))


save(model_accuracy, 
     model_lasso, model_rf, 
     mfx_lasso, mfx_rf,
     file = "data_derived/real_data.RData")
