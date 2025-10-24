library(dplyr)
library(tidyr)
library(purrr)
library(glmnet)
library(MuMIn)
library(Rcpp)
library(furrr)
library(dplyr)
library(parallel)
library(future)
library(yardstick)
library(caret)
library(ordinalNet)
library(corrplot)

csv = read.csv("G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\Plot_Level_Datasets\\responsematrix_06-19-25.csv")



################# CONIFER VARIABLES ############################################
train = csv[,c(6, 8:115)]
train = train[, !grepl("propor", names(train))]


colnames(train)[1] = "success"
train = train[train$success > 0.1,]
train$plot_id = 1:length(train[,1])
n_plots = length(train[,1])

train$regenCL = ifelse(train$success < 5, "Low", 
                       ifelse(train$success >= 5 & train$success < 25, "Medium",
                              ifelse(train$success >= 25, "High", "NA")))


counts <- table(train$regenCL)
barplot(counts)

train = train[-1]

colnames(train)[colnames(train) == "regenCL"] <- "success"
train$success <- factor(train$success, levels = c("Low", "Medium", "High"), ordered = TRUE)
barplot(table(train$success))

str(train$success)

#------------------------------------------------------------------------------#

# --- Step 2: Prepare X (predictors) and y (response) for ordinalNet ---
y <- train$success
X <- as.matrix(train[, !(colnames(train) %in% c("success", "plot_id"))])

#------------------------------------------------------------------------------

kruskal_results <- apply(X, 2, function(x) kruskal.test(x ~ y)$p.value)
top_vars <- names(sort(kruskal_results))[1:25]  # Select top 25
X_reduced <- X[, top_vars]


# --- Step 3: Fit Penalized Ordinal Logistic Regression with Cross-Validation ---
set.seed(123)
cvfit <- ordinalNetCV(
  x = scale(X_reduced),
  y = y,
  family = "cumulative",
  link = "logit",
  parallelTerms = TRUE,
  nonparallelTerms = TRUE,  # <- enables semi-parallel model
  alpha = 0.5,
  tuneMethod = "cvLoglik",
  maxiterOut = 500,
  maxiterIn = 500
)

# --- Step 4: Get Selected Variables at Best Lambda ---
best_model <- cvfit$fit
summary(best_model)
best_lambda_index <- which.min(cvfit$misclass)
coefs <- coef(best_model, whichLambda = best_lambda_index, matrix = TRUE)

# Display non-zero coefficients
nonzero_vars <- rownames(coefs)[apply(coefs, 1, function(x) any(x != 0))]
nonzero_vars

#------------------------------------------------------------------------------
# Remove correlated variables
nonzero_vars_clean <- setdiff(nonzero_vars, c("(Intercept)"))
var_df <- dplyr::select(train, all_of(nonzero_vars_clean))

cor_matrix <- cor(var_df)  # X = predictor matrix (no response var)
corrplot(cor_matrix,
         method = "color",       # use color shading
         type = "upper",         # show only upper triangle
         order = "hclust",       # hierarchical clustering order
         addCoef.col = "black",  # add correlation coefficients
         tl.col = "black",       # text label color
         tl.cex = 0.8,           # text label size
         number.cex = 0.7,       # correlation number size
         diag = FALSE)  
drop_idx <- findCorrelation(cor_matrix, cutoff = 0.7)
noncor_vars <- var_df[, -drop_idx]  # keep only uncorrelated predictors
colnames(noncor_vars)
#------------------------------------------------------------------------------
X_nocor <- as.matrix(noncor_vars)
y <- train$success

set.seed(123)
cvfit <- ordinalNetCV(
  x = X_nocor,
  y = y,
  family = "cumulative",
  link = "logit",
  parallelTerms = TRUE,
  nonparallelTerms = TRUE,
  alpha = 0.5,
  tuneMethod = "cvLoglik",
  maxiterOut = 500
)

best_lambda_index <- which.min(cvfit$misclass)
coefs <- coef(cvfit$fit, whichLambda = best_lambda_index, matrix = TRUE)
nonzero_vars <- rownames(coefs)[apply(coefs, 1, function(x) any(x != 0))]
print(nonzero_vars)

#------------------------------------------------------------------------------#
# Predict
# Predict class labels on training data
predicted_class <- predict(cvfit$fit, 
                           newx = X_nocor, 
                           whichLambda = which.min(cvfit$misclass), 
                           type = "class")

# Confusion matrix
table(Observed = y, Predicted = predicted_class)

level_labels <- c("Low", "Medium", "High")

y <- factor(y, levels = level_labels, ordered = TRUE)
predicted_class <- factor(predicted_class, 
                          levels = 1:3, 
                          labels = level_labels, 
                          ordered = TRUE)
confusionMatrix(data = predicted_class, reference = y)

