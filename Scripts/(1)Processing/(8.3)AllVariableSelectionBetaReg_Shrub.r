## Beta Regression
# Shrub

library(betareg)
library(furrr)
library(future)
library(caret)
library(corrplot)
library(dplyr)
library(future.apply)

plan(multisession)
options(future.globals.maxSize = +Inf)

#------------------------------------------------------------------------------
# LOAD AND CLEAN DATA
#------------------------------------------------------------------------------

csv <- read.csv("G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\Plot_Level_Datasets\\responsematrix_06-19-25.csv")

regen <- subset(csv, shrub.percentages > 0)
train <- regen[, 7:ncol(regen)]
train <- train[, !grepl("propor", names(train))]

colnames(train)[1] <- "y"
train$success <- as.numeric(as.character(train$y))  # beta response
train <- train[, !(colnames(train) %in% c("y"))]

hist(train$success)

#------------------------------------------------------------------------------
# REMOVE CORRELATED VARIABLES
#------------------------------------------------------------------------------

var_df <- train[, !(colnames(train) %in% c("success", "plot_id"))]
cor_matrix <- cor(var_df)
drop_idx <- findCorrelation(cor_matrix, cutoff = 0.7)

if (length(drop_idx) > 0) {
  noncor_vars <- var_df[, -drop_idx]
} else {
  noncor_vars <- var_df
}

# Visual check (optional)
corrplot(cor(noncor_vars),
         method = "color", type = "upper",
         order = "hclust", addCoef.col = "black",
         tl.col = "black", tl.cex = 0.8,
         number.cex = 0.7, diag = FALSE)

#------------------------------------------------------------------------------
# GENERATE COMBINATIONS
#------------------------------------------------------------------------------

variables <- colnames(noncor_vars)

get_combinations <- function(variables, n) {
  unlist(lapply(1:n, function(x) combn(variables, x, simplify = FALSE)), recursive = FALSE)
}

combinations <- get_combinations(variables, 5)

#------------------------------------------------------------------------------
# FIT MODELS AND CALCULATE AIC
#------------------------------------------------------------------------------

pb <- txtProgressBar(min = 0, max = length(combinations), style = 3)

for (i in seq_along(combinations)) {
  setTxtProgressBar(pb, i)
  vars <- combinations[[i]]
  form <- as.formula(paste("success ~", paste(vars, collapse = "+")))
  
  model <- tryCatch(betareg(form, data = train), error = function(e) NULL)
  aic <- if (!is.null(model)) AIC(model) else 999999
  
  aic_list[[i]] <- data.table(
    index = i,
    variables = paste(vars, collapse = ", "),
    aic = aic
  )
}
close(pb)

aic_scores <- rbindlist(aic_list)

write.csv(aic_scores, "G:\\YOSE_Regen_Analysis\\Runs\\AICscoresOrdReg_conifer_085525.csv", row.names = FALSE)

#------------------------------------------------------------------------------
# k-Fold Cross Validation
#------------------------------------------------------------------------------

# Get the best row
best_row <- aic_scores[which.min(aic_scores$aic), ]
# View it
print(best_row)

best_vars <- unlist(strsplit(best_row$variables, ",\\s*"))
form <- as.formula(paste("success ~", paste(best_vars, collapse = " + ")))


# Ensure response is in (0,1) for beta regression
train$success <- as.numeric(train$success)
epsilon <- 1e-4
train$success <- pmin(pmax(train$success, epsilon), 1 - epsilon)

# Define variables
k <- 5  # number of folds
set.seed(42)
folds <- createFolds(train$success, k = k, list = TRUE)

cv_results <- data.frame(fold = integer(), RMSE = numeric(), MAE = numeric())

for (i in seq_along(folds)) {
  test_idx <- folds[[i]]
  train_fold <- train[-test_idx, ]
  test_fold <- train[test_idx, ]
  
  model <- tryCatch(betareg(form, data = train_fold), error = function(e) NULL)
  
  if (!is.null(model)) {
    preds <- predict(model, newdata = test_fold, type = "response")
    actual <- test_fold$success
    
    rmse <- sqrt(mean((preds - actual)^2))
    mae <- mean(abs(preds - actual))
    
    cv_results <- rbind(cv_results, data.frame(fold = i, RMSE = rmse, MAE = mae))
  }
}

mean_rmse <- mean(cv_results$RMSE, na.rm = TRUE)
mean_mae <- mean(cv_results$MAE, na.rm = TRUE)

print(cv_results)
cat("Mean RMSE:", mean_rmse, "\n")
cat("Mean MAE:", mean_mae, "\n")

final_model <- betareg(form, data = train)
summary(final_model)

#------------------------------------------------------------------------------
# Residuals
#------------------------------------------------------------------------------

# Extract residuals
res <- residuals(final_model, type = "pearson")   # or type = "deviance", "response"

# Basic histogram
hist(res, breaks = 20, main = "Histogram of Residuals", xlab = "Residuals")

# Residuals vs Fitted
plot(predict(final_model, type = "response"), res,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)

# Q-Q plot for residuals
qqnorm(res)
qqline(res, col = "blue")

par(mfrow = c(3, 2))
suppressWarnings(RNGversion("3.5.0"))
set.seed(123)
plot(final_model, which = 1:4, type = "pearson")
plot(final_model, which = 5, type = "deviance", sub.caption = "")
plot(final_model, which = 1, type = "deviance", sub.caption = "")
par(mfrow = c(1, 1))


