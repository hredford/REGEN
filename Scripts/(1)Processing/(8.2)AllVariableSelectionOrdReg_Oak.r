## Ordinal Logistic Regression
# Oak

library(dplyr)
library(tidyr)
library(purrr)
library(ordinalNet)
library(corrplot)


csv = read.csv("G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\Plot_Level_Datasets\\responsematrix_06-19-25.csv")

train = csv[,c(5, 8:115)]
train = train[, !grepl("propor", names(train))]
train = subset(train, oakcount >= 1)

train$success = ifelse(train$oakcount < 4, "Low", 
                       ifelse(train$oakcount >= 4 & train$oakcount < 8, "Medium",
                              ifelse(train$oakcount >= 8, "High", "NA")))

# Check that the levels are equal
counts <- table(train$success)
barplot(counts)

##################################################################################
train = train[-1] # remove the original count data

# Make sure the levels are ordered correctly
train$success <- factor(train$success, levels = c("Low","Medium", "High"), ordered = TRUE)
barplot(table(train$success))
str(train$success)

#------------------------------------------------------------------------------#
# --- Prepare X (predictors) and y (response) for ordinalNet ---

var_df = train[, !(colnames(train) %in% c("success", "plot_id"))]

y <- train$success
X <- as.matrix(var_df)

#------------------------------------------------------------------------------
# --- Remove Correlated Variables --
cor_matrix <- cor(var_df)  # X = predictor matrix (no response var)
drop_idx <- findCorrelation(cor_matrix, cutoff = 0.7) # this will group the variables and keep the one with the highest influence

if (length(drop_idx) > 0) {         # keep only uncorrelated predictors
  noncor_vars <- var_df[, -drop_idx]
} else {
  noncor_vars <- var_df  # keep all variables
} 
colnames(noncor_vars)
cor_matrix <- cor(noncor_vars) 

corrplot(cor_matrix,
         method = "color",       # use color shading
         type = "upper",         # show only upper triangle
         order = "hclust",       # hierarchical clustering order
         addCoef.col = "black",  # add correlation coefficients
         tl.col = "black",       # text label color
         tl.cex = 0.8,           # text label size
         number.cex = 0.7,       # correlation number size
         diag = FALSE)  

#------------------------------------------------------------------------------
# --- Kruskal-Wallis Test ---
# Non-parametric test identifying variables with significant distribution differences across oak classes 

X = noncor_vars # make new X df

kruskal_results <- apply(X, 2, function(x) kruskal.test(x ~ y)$p.value)
top_vars <- names(sort(kruskal_results))[c(1:9)]  # Select top 10
top_vars = c("emid_55to85",   "emin_90to99",   "linear_90to99", "emin__TAO",     "iqr95_90to99",  "planar_90to99") #ONE THROUGH NINE ABOVE
X_reduced <- as.matrix(X[, top_vars]) 

#------------------------------------------------------------------------------
# --- Fit Penalized Ordinal Logistic Regression with Cross-Validation ---
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

# --- Get Selected Variables at Best Lambda ---
best_model <- cvfit$fit
summary(best_model)
best_lambda_index <- which.min(cvfit$misclass)
coefs <- coef(best_model, whichLambda = best_lambda_index, matrix = TRUE)

# Display non-zero coefficients
nonzero_vars <- rownames(coefs)[apply(coefs, 1, function(x) any(x != 0))]
nonzero_vars

selected_vars <- setdiff(nonzero_vars, "(Intercept)")
X_final <- scale(noncor_vars[, selected_vars])

#------------------------------------------------------------------------------#
# Predict

# scale the variables for prediction
X_reduced_scaled <- scale(X_reduced)

# Predict class labels on training data
predicted_class <- predict(
  cvfit$fit,
  newx = X_reduced_scaled,
  whichLambda = best_lambda_index,
  type = "class"
)

# Confusion matrix
table(Observed = y, Predicted = predicted_class)

level_labels <- c("Low", "Medium", "High")

y <- factor(y, levels = level_labels, ordered = TRUE)
predicted_class <- factor(predicted_class, 
                          levels = 1:3, 
                          labels = level_labels, 
                          ordered = TRUE)
confusionMatrix(data = predicted_class, reference = y)



