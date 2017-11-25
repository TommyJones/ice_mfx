rm(list = ls())

source("R/00_common_functions.R")

load("data_derived/real_data.RData")

eval <- lapply(model_accuracy, function(x){
  lapply(x, function(y) y$stats$F1)
})

plot(0, ylim = c(0,1), xlim = c(0,1), col = "white",
     ylab = "F1", xlab = "Threshold")

for(j in 1:10){
  lines(seq(0,1, length = length(eval[[ 1 ]][[ j ]])),
        eval[[ 1 ]][[ j ]], col = rgb(1,0,0,0.3))
  lines(seq(0,1, length = length(eval[[ 2 ]][[ j ]])),
        eval[[ 2 ]][[ j ]], col = rgb(0,0,1,0.3))
}

lines(seq(0,1, length = length(eval[[ 1 ]][[ j ]])),
      Reduce("+", eval[[ 1 ]]) / 10, col = rgb(1,0,0,1), lwd = 2)
lines(seq(0,1, length = length(eval[[ 1 ]][[ j ]])),
      Reduce("+", eval[[ 2 ]]) / 10, col = rgb(0,0,1,1), lwd = 2)

legend("topright", legend = names(eval), col = c("red", "blue"), lty = 1)
