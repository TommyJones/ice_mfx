################################################################################
# This script runs some specific plots and saves them in output_charts_tables
################################################################################

rm(list = ls())

source("R/00_common_functions.R")

### plots from 02_one_interaction.R ----------------
load("data_derived/one_interaction.RData")

# compare x2 on RF model

png("output_charts_tables/x2sq.png",
    width = 8, height = 6, units = "in", res = 300)
plot(mfx_rf_logit[[ "x2sq" ]])
dev.off()

png("output_charts_tables/x2sq_deriv.png",
    width = 8, height = 6, units = "in", res = 300)
plot(mfx_rf_logit[[ "x2sq" ]], type = "derivative")
dev.off()

png("output_charts_tables/x2sq_deriv_full_range.png",
    width = 8, height = 6, units = "in", res = 300)
plot(mfx_rf_logit[[ "x2sq" ]], type = "derivative", 
     ylim = range(mfx_rf_logit$x2sq$true_values$dy / mfx_rf_logit$x2sq$true_values$dx))
dev.off()

png("output_charts_tables/x2sq_deriv_overlay.png",
    width = 8, height = 6, units = "in", res = 300)
plot(mfx_rf_logit[[ "x2sq" ]], type = "derivative")
abline(h = mean(mfx_rf_logit$x2sq$dy / mfx_rf_logit$x2sq$dx), 
       col = "red", lty = 3, lwd = 5)
dev.off()


mfx_new <- mean(mfx_rf_logit$x2sq$dy / mfx_rf_logit$x2sq$dx)

int <- 0.565

mfx_line_new <- mfx_rf_logit$x2sq$x * mfx_new + int

png("output_charts_tables/x2sq_overlay.png",
    width = 8, height = 6, units = "in", res = 300)
plot(mfx_rf_logit[[ "x2sq" ]])
lines(mfx_rf_logit$x2sq$x,
      mfx_line_new, lty = 3, lwd = 5, col = "red")
dev.off()

### plots from 03_real_data.R ----------------
load("data_derived/real_data.RData")

eval <- lapply(model_accuracy, function(x){
  lapply(x, function(y) y$stats$F1)
})

png("output_charts_tables/real_data_f1.png", 
    width = 8, height = 6, units = "in", res = 300)

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

dev.off()

plot(mfx_lasso$delta_v, centered = T, main = "Lasso")

plot(mfx_rf$delta_v, centered = T, main = "Random Forest")
