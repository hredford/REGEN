################################################################################
#
#
#
#
#                     Residual Plotting
#
#
#
#
#
#
################################################################################

library(dplyr)
library(ggplot2)
library(ISLR)
library(gridExtra)
library(arm)
library(tidyr)
library(corrplot)
library(car)
library(pROC)

csv = read.csv("G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\responsematrix_04-25-25.csv")

conifer_matrix = dplyr::select(csv, c(regencount, linear_75to95, iqr_90to99, b11w, ndviw))# regen$y[regen$y > 0] = 1
colnames(conifer_matrix)[1] = "y"
conifer_matrix$y[conifer_matrix$y > 0] = 1
conifer_matrix$y[conifer_matrix$y < 1] = 0
conifer_matrix$success = conifer_matrix$y
conifer_matrix = conifer_matrix[,-1]
conifer_matrix$success = as.factor(conifer_matrix$success)

oak_matrix = dplyr::select(csv, c(oakcount, anisotr_85to95, blue, b9, b3s))
colnames(oak_matrix)[1] = "y"
oak_matrix$y[oak_matrix$y > 0] = 1
oak_matrix$y[oak_matrix$y < 1] = 0
oak_matrix$success = oak_matrix$y
oak_matrix = oak_matrix[,-1]
oak_matrix$success = as.factor(oak_matrix$success)

shrub_matrix = dplyr::select(csv, c(shrub.percentages, rumple_85to95, blue, b10s, ndirs))
colnames(shrub_matrix)[1] = "y"
shrub_matrix$y[shrub_matrix$y == 0] = 0
shrub_matrix$y[shrub_matrix$y > 0] = 1
shrub_matrix$success = shrub_matrix$y
shrub_matrix = shrub_matrix[,-1]
shrub_matrix$success = as.factor(shrub_matrix$success)

############################# PLOT MAIN VARS ###################################

# Conifer Metrics
grid.arrange(
ggplot(conifer_matrix, aes(success, linear_75to95)) + 
  geom_boxplot() +
  ggtitle("Linearity") + 
  ylab("Linearity in the 75th to 95th Percentile"),
ggplot(conifer_matrix, aes(success, iqr_90to99)) + 
  geom_boxplot() +
  ggtitle("Interquartile Range") + 
  ylab("IQR in the 90th to 99th Percentile"),
ggplot(conifer_matrix, aes(success, b11w)) + 
  geom_boxplot() +
  ggtitle("Band 11 Reflectance") + 
  ylab("Band 11 (October Imagery)"),
ggplot(conifer_matrix, aes(success, ndviw)) + 
  geom_boxplot() +
  ggtitle("NDVI") + 
  ylab("NDVI (October Imagery)"),
ncol = 4)

numeric_vars <- conifer_matrix[, sapply(conifer_matrix, is.numeric)]
cor_matrix <- cor(numeric_vars, use = "complete.obs")
corrplot(cor_matrix, method = "circle", tl.cex = 0.8, addCoef.col = "black")


# Oak Metrics
grid.arrange(
  ggplot(oak_matrix, aes(success, anisotr_85to95)) + 
    geom_boxplot() +
    ggtitle("Anisotropy") + 
    ylab("Anisotropy in the 75th to 95th Percentile"),
  ggplot(oak_matrix, aes(success, blue)) + 
    geom_boxplot() +
    ggtitle("Blue Reflectance") + 
    ylab("NAIP Blue Band"),
  ggplot(oak_matrix, aes(success, b9)) + 
    geom_boxplot() +
    ggtitle("Band 9 Reflectance") + 
    ylab("Band 9 (Seasonal Difference)"),
  ggplot(oak_matrix, aes(success, b3s)) + 
    geom_boxplot() +
    ggtitle("B3 Imagery") + 
    ylab("B3 (June Imagery)"),
  ncol = 4)

numeric_vars <- oak_matrix[, sapply(oak_matrix, is.numeric)]
cor_matrix <- cor(numeric_vars, use = "complete.obs")
corrplot(cor_matrix, method = "circle", tl.cex = 0.8, addCoef.col = "black")

# Shrub Metrics
grid.arrange(
  ggplot(shrub_matrix, aes(success, rumple_85to95)) + 
    geom_boxplot() +
    ggtitle("Rumple") + 
    ylab("Rumple in the 85th to 95th Percentile"),
  ggplot(shrub_matrix, aes(success, blue)) + 
    geom_boxplot() +
    ggtitle("Blue Reflectance") + 
    ylab("NAIP Blue Band"),
  ggplot(shrub_matrix, aes(success, b10s)) + 
    geom_boxplot() +
    ggtitle("Band 10 Reflectance") + 
    ylab("Band 10 (June Imagery)"),
  ggplot(shrub_matrix, aes(success, ndirs)) + 
    geom_boxplot() +
    ggtitle("Normalized Difference Infrared Index") + 
    ylab("NDII (June Imagery)"),
  ncol = 4)

numeric_vars <- shrub_matrix[, sapply(shrub_matrix, is.numeric)]
cor_matrix <- cor(numeric_vars, use = "complete.obs")
corrplot(cor_matrix, method = "circle", tl.cex = 0.8, addCoef.col = "black")

#################################################################################




############################## SPLIT DATASETS ##################################

split_data <- function(matrix_data, train_frac = 0.7) {
  idx <- sample(nrow(matrix_data), size = floor(train_frac * nrow(matrix_data)))
  train <- matrix_data[idx, ]
  test  <- matrix_data[-idx, ]
  return(list(train = train, test = test))
}

conifer_split <- split_data(conifer_matrix)
oak_split     <- split_data(oak_matrix)
shrub_split   <- split_data(shrub_matrix)

conifer_train <- conifer_split$train
conifer_test  <- conifer_split$test

oak_train     <- oak_split$train
oak_test      <- oak_split$test

shrub_train   <- shrub_split$train
shrub_test    <- shrub_split$test


################################## MODEL #######################################

conifer_model <- glm(success ~ linear_75to95 + iqr_90to99 + b11w + ndviw, data = conifer_train, family = binomial)
oak_model     <- glm(success ~ ., data = oak_train,     family = binomial)
shrub_model   <- glm(success ~ ., data = shrub_train,   family = binomial)

summary(conifer_model)
summary(oak_model)
summary(shrub_model)


################################# PLOT RESIDUALS ######################################

par(cex = 1.2,        # overall text size multiplier
    cex.lab = 1.75,    # axis label size
    cex.axis = 1.3,   # axis number size
    cex.main = 1.5)   # title size

residualPlots(conifer_model)
residualPlots(oak_model)
residualPlots(shrub_model)

##################################### CHI Square ###############################
drop1(conifer_model, test = "Chisq")


#################################### ROC CURVE #################################

conifer_test$prob <- predict(conifer_model, newdata = conifer_test, type = "response")
roc_obj <- roc(conifer_test$success, conifer_test$prob)
plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve - Conifer Model")
abline(a = 0, b = 1, lty = 2, col = "gray")  # Diagonal line
auc(roc_obj)

oak_test$prob <- predict(oak_model, newdata = oak_test, type = "response")
roc_obj <- roc(oak_test$success, oak_test$prob)
plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve - Conifer Model")
abline(a = 0, b = 1, lty = 2, col = "gray")  # Diagonal line
auc(roc_obj)

shrub_test$prob <- predict(shrub_model, newdata = shrub_test, type = "response")
roc_obj <- roc(shrub_test$success, shrub_test$prob)
plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve - Conifer Model")
abline(a = 0, b = 1, lty = 2, col = "gray")  # Diagonal line
auc(roc_obj)

################################################################################

# List of predictors to plot
#predictors <- c("linear_75to95", "iqr_90to99", "b11w", "ndviw")
#predictors <- c("anisotr_85to95", "blue", "b9", "b3s")
predictors <- c("rumple_85to95", "blue", "b10s", "ndirs")

train = shrub_train
model = shrub_model

# Get means of all variables for baseline
means <- colMeans(train[sapply(train, is.numeric)], na.rm = TRUE)

# Build prediction dataframe for each variable
plot_data <- bind_rows(lapply(predictors, function(var) {
  seq_vals <- seq(min(train[[var]], na.rm = TRUE), 
                  max(train[[var]], na.rm = TRUE), 
                  length.out = 100)
  
  new_data <- as.data.frame(t(rep(means, 100)))
  names(new_data) <- names(means)
  new_data <- new_data[rep(1, 100), ]
  new_data[[var]] <- seq_vals
  new_data$variable <- var
  new_data$xval <- seq_vals
  
  preds <- predict(model, newdata = new_data, type = "response", se.fit = TRUE)
  new_data$fit <- preds$fit
  new_data$upper <- preds$fit + 1.96 * preds$se.fit
  new_data$lower <- preds$fit - 1.96 * preds$se.fit
  
  return(new_data[, c("xval", "fit", "upper", "lower", "variable")])
}))


ggplot(plot_data, aes(x = xval, y = fit)) +
  geom_line(size = 1.2, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.6) +
  facet_wrap(~variable, scales = "free_x", ncol = 2) +
  labs(x = NULL, y = "Predicted Probability of Conifer Success") +
  theme_minimal(base_size = 13)










