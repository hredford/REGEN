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

csv = read.csv("G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\Plot_Level_Datasets\\responsematrix_06-19-25.csv")



################# SHRUB VARIABLES ############################################
train = csv[,c(7, 8:115)]
train = train[, !grepl("propor", names(train))]


colnames(train)[1] = "y"
train$y[train$y > 0.001] = 1
train$y[train$y < 0.001] = 0

train$success = train$y
train = train[,-1]
train$success = as.numeric(train$success)
train$plot_id = 1:length(train[,1])
n_plots = length(train[,1])

get_combinations <- function(variables, n) {
  unlist(lapply(1:n, function(x) combn(variables, x, simplify = FALSE)), recursive = FALSE)
}


plan(multisession, workers = 15)  # or use multicore if you're on Unix

# Setup
variables <- colnames(train)[!colnames(train) %in% c("success", "plot_id")]
combo3 = get_combinations(variables, 2)
combinations <- get_combinations(variables, 4)
combinations = combinations[c(length(combo3):length(combinations))]
rm(combo3)


# Pre-sample a fixed set of combination indices for each bootstrap
n_boots <- 100  # or more
n_combos <- 5000  # combos per bootstrap
set.seed(123)
combo_indices <- replicate(n_boots, sample(seq_along(combinations), n_combos), simplify = FALSE)


options(future.globals.maxSize = +Inf)


run_bootstrap_static <- function(seed, combo_indices_subset, combo_indices_id, n_plots = 113) {
  
  set.seed(seed)
  selected_plots <- sample(unique(train$plot_id), n_plots, replace = TRUE)
  train_sub <- train[train$plot_id %in% selected_plots, ]
  
  results <- vector("list", length(combo_indices_subset))
  
  for (j in seq_along(combo_indices_subset)) {
    idx <- combo_indices_subset[j]
    vars <- combinations[[idx]]
    form <- as.formula(paste("success ~", paste(vars, collapse = "+")))
    
    model <- tryCatch(glm(form, data = train_sub, family = binomial), error = function(e) NULL)
    
    if (!is.null(model)) {
      aic <- AIC(model)
      
      # Predict and classify at 0.5 threshold
      probs <- predict(model, type = "response")
      pred_class <- ifelse(probs >= 0.5, 1, 0)
      
      # AUC
      auc <- tryCatch({
        if (length(unique(train_sub$success)) > 1) {
          as.numeric(pROC::auc(train_sub$success, probs))
        } else NA
      }, error = function(e) NA)
      
      sig_vars <- character(0)
      if (!is.na(auc) && auc >= 0.7) {
        sig_vars <- tryCatch({
          summ <- summary(model)
          coef_pvals <- coef(summ)[, "Pr(>|z|)"]
          coef_pvals <- coef_pvals[names(coef_pvals) %in% vars]  # remove intercept
          names(coef_pvals)[coef_pvals < 0.05]
        }, error = function(e) character(0))
      }
      
      # Cohen's Kappa
      kappa_val <- tryCatch({
        cm <- confusionMatrix(factor(pred_class), factor(train_sub$success), positive = "1")
        as.numeric(cm$overall["Kappa"])
      }, error = function(e) NA)
      
      # Max correlation
      predictor_data <- train_sub[, vars, drop = FALSE]
      cor_mat <- tryCatch(cor(predictor_data, use = "pairwise.complete.obs"), error = function(e) NA)
      max_corr <- if (is.matrix(cor_mat)) {
        max(abs(cor_mat[lower.tri(cor_mat)]), na.rm = TRUE)
      } else NA
      
      results[[j]] <- data.frame(
        seed = seed,
        combo_index = idx,
        combo_set = combo_indices_id,
        variables = paste(vars, collapse = ", "),
        vars = I(list(vars)),
        sig_vars = I(list(sig_vars)),
        aic = aic,
        max_corr = max_corr,
        auc = auc,
        kappa = kappa_val,
        stringsAsFactors = FALSE
      )
    }
  }
  
  dplyr::bind_rows(results)
}


# Ensure future can export large objects
options(future.globals.maxSize = +Inf)

#----------------------------------------------------------------------------
# RUN THE BOOTSTRAPPING FUNCTION — furrr will now export combinations and train
#----------------------------------------------------------------------------
all_results <- future_map_dfr(1:n_boots, function(i) {
  run_bootstrap_static(
    seed = i,
    combo_indices_subset = combo_indices[[i]],
    combo_indices_id = i
  )
}, .options = furrr_options(seed = TRUE))

#------------------------------------------------------------------------------
# Get the relevant statitistics
top_sig <- all_results %>%
  filter(!is.na(kappa)) %>%  # Only keep successful model runs
  unnest_longer(sig_vars) %>%  # Expand variables to one per row
  group_by(sig_vars) %>%
  summarize(
    significant_count = n(),
    mean_kappa = mean(kappa, na.rm = TRUE),
    max_kappa = max(kappa, na.rm = TRUE),
    pct_above_0_4 = mean(kappa > 0.4, na.rm = TRUE),
    mean_auc = mean(auc, na.rm = TRUE),
    mean_aic = mean(aic, na.rm = TRUE),
    .groups = "drop"
  )

sig_var_counts <- all_results %>%
  filter(!is.na(auc) & auc >= 0.7) %>%
  unnest(sig_vars) %>%
  group_by(sig_vars) %>%
  summarise(significant_count = n(), .groups = "drop")

variable_votes <- all_results %>%
  unnest(vars) %>%
  group_by(vars) %>%
  summarise(total_count = n(), .groups = "drop")

variable_importance <- full_join(variable_votes, sig_var_counts,
                                 by = c("vars" = "sig_vars")) %>%
  mutate(significant_count = replace_na(significant_count, 0)) %>%
  arrange(desc(significant_count))

write.csv(top_sig, "G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\Plot_Level_Datasets\\VariableVotes_Shrub.csv", row.names = FALSE)


################################## Find best variables #########################

library(ggplot2)
library(reshape2)

sub_variable = variable_importance #[-c(2,3,5,6,7,10:16),]

# 1. Select top N variables
top_vars <- sub_variable %>%
  arrange(desc(significant_count)) %>%
  slice(1:10) %>%
  pull(vars)

# 2. Subset your training data to these variables
predictors_top <- train[, top_vars]

# 3. Compute the correlation matrix
cor_mat <- cor(predictors_top, use = "pairwise.complete.obs")

# 4. Reshape for plotting
cor_df <- melt(cor_mat)

ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), size = 3, color = "black") +  # <-- Add this line
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limit = c(-1, 1), name = "Correlation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  coord_fixed()


####################### PLOT METRICS ###################################



df <- as_tibble(all_results)
df_filtered <- df[!is.na(df$auc) & df$auc >= 0.7, ]
df_selected <- df_filtered[, c("sig_vars", "aic", "auc", "kappa")]
df_unnested <- unnest(df_selected, cols = c(sig_vars))

# Order variables
top_sig_summary <- summarise(
  group_by(df_unnested, sig_vars),
  significant_count = n(),
  mean_aic = mean(aic),
  mean_auc = mean(auc),
  mean_kappa = mean(kappa)
)

top_sig_summary <- top_sig_summary[order(-top_sig_summary$significant_count), ]
top_sig_summary$sig_vars <- factor(top_sig_summary$sig_vars, levels = rev(top_sig_summary$sig_vars))


top_sig_summary$auc_scaled <- top_sig_summary$mean_auc * 100
top_sig_summary$kappa_scaled <- top_sig_summary$mean_kappa * 100

top_sig_summary <- top_sig_summary[c(1:min(10, nrow(top_sig_summary))), ]

# Select and reshape to long format
plot_df <- top_sig_summary[, c("sig_vars", "mean_aic", "auc_scaled", "kappa_scaled")]

# Pivot longer
plot_df <- pivot_longer(
  data = plot_df,
  cols = c("mean_aic", "auc_scaled", "kappa_scaled"),
  names_to = "metric",
  values_to = "value"
)

# Rename metric levels
plot_df$metric <- recode(plot_df$metric,
                         mean_aic = "AIC",
                         auc_scaled = "AUC",
                         kappa_scaled = "Kappa")

# Define color and shape
metric_colors <- c("AUC" = "darkgreen", "Kappa" = "purple", "AIC" = "red")
metric_shapes <- c("AUC" = 16, "Kappa" = 17, "AIC" = 15)

# Final plot
ggplot(plot_df, aes(x = value, y = sig_vars, color = metric, shape = metric)) +
  geom_point(size = 3) +
  scale_color_manual(values = metric_colors) +
  scale_shape_manual(values = metric_shapes) +
  scale_x_continuous(
    name = "AIC",
    sec.axis = sec_axis(~./100, name = "AUC / Kappa")
  ) +
  labs(
    y = "Variable",
    color = " ",
    shape = " ",
    title = "Top Significant Shrub Variables: AIC, AUC & Kappa"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.title.x.top = element_text(color = "black"),
    axis.title.x.bottom = element_text(color = "red"),
    legend.position = "top"
  )


########################## Get the Statistics ##################################

library(bestglm)


regen = dplyr::select(train, c(success, rumple_90to99, blue, b10, ndvis))
#regen = regen[c(2:9, 11:30),]
hist(regen$regencount)
regen = na.omit(regen)

colnames(regen)[1] = "y"

regen$y[regen$y > 0] = 1
regen$y[regen$y < 1] = 0

regen$success = regen$y
#oaks$failure = 1-oaks$y
regen = regen[,-1]
regen$success = as.numeric(regen$success)


model = bestglm(regen, family = binomial, IC = "AIC")
summary(model)

model$BestModels

truepos = 0
trueneg = 0
falspos = 0
falsneg = 0



for (i in c(1:length(regen[,1]))) {
  
  df = regen[-c(i),]
  
  #
  glmmodel = glm(success~ rumple_90to99 + blue + b10 + ndvis , data = df, family = binomial)
  
  pred = predict(glmmodel, regen[i,], type = 'response')
  print(pred)
  
  if (pred > 0.5) {
    if (regen$success[i] == 1) {
      truepos = truepos +1
      
    } else{
      falspos = falspos +1
    }
    
    
  } 
  
  if (pred <=  0.5) {
    if (regen$success[i] == 1) {
      falsneg = falsneg +1
      
    } else {
      trueneg = trueneg + 1
    }
  }
  
  
  
}


precision = truepos/(truepos+falspos)
recall = truepos/(truepos+falsneg)

fscore = 2*(precision*recall)/ (precision+recall)

summary(glmmodel)

##################################################################################
##################### K-Fold Cross Validaiton #############################

regen = dplyr::select(train, c(success, rumple_90to99, blue, b10, ndvis))
regen$success <- factor(regen$success, levels = c(0, 1), labels = c("no", "yes"))

true_labels <- regen$success

set.seed(42)

customSummary <- function(data, lev = NULL, model = NULL) {
  require(pROC)
  require(caret)
  roc_val <- tryCatch({
    roc_obj <- roc(response = data$obs, predictor = data$yes)
    auc(roc_obj)
  }, error = function(e) NA)
  
  kappa_val <- confusionMatrix(data = data$pred, reference = data$obs)$overall["Kappa"]
  
  out <- c(ROC = as.numeric(roc_val),
           Kappa = as.numeric(kappa_val),
           Sens = sensitivity(data$pred, data$obs, positive = "yes"),
           Spec = specificity(data$pred, data$obs, positive = "yes"))
  return(out)
}

ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = customSummary,
  sampling = "up"
)

final_model <- train(
  success ~ rumple_90to99 + blue + b10 + ndvis,
  data = regen,
  method = "glm",
  family = binomial,
  metric = "ROC",
  trControl = ctrl
)

final_model

final_model$results
c = confusionMatrix(final_model)

# Extract probabilities for "yes"
pred_probs <- predict(final_model, type = "prob")[, "yes"]

# Build ROC and threshold table
roc_obj <- roc(true_labels, pred_probs)
threshold_df <- as.data.frame(coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity")))

# Calculate Kappa at each threshold
kappa_vals <- sapply(threshold_df$threshold, function(t) {
  pred_class <- factor(ifelse(pred_probs > t, "yes", "no"), levels = c("no", "yes"))
  
  if (length(unique(pred_class)) == 2) {
    cm <- confusionMatrix(pred_class, true_labels, positive = "yes")
    as.numeric(cm$overall["Kappa"])
  } else {
    NA
  }
})

# Find best threshold and kappa
valid_indices <- which(!is.na(kappa_vals))
best_index <- valid_indices[which.max(kappa_vals[valid_indices])]
best_threshold <- threshold_df$threshold[best_index]
best_kappa <- kappa_vals[best_index]

# Output
cat("🔍 Best threshold for Kappa:", round(best_threshold, 3), "\n")
cat("📈 Max Kappa at that threshold:", round(best_kappa, 3), "\n")

final_pred_class <- factor(ifelse(pred_probs > best_threshold, "yes", "no"), levels = c("no", "yes"))
confusionMatrix(final_pred_class, true_labels, positive = "yes")


#----------------Get the Prediction Table ------------------------

label <- character(length(true_labels))

for (i in seq_along(true_labels)) {
  if (true_labels[i] == "yes" && final_pred_class[i] == "yes") {
    label[i] <- "True Positive"
  } else if (true_labels[i] == "no" && final_pred_class[i] == "no") {
    label[i] <- "True Negative"
  } else if (true_labels[i] == "no" && final_pred_class[i] == "yes") {
    label[i] <- "False Positive"
  } else if (true_labels[i] == "yes" && final_pred_class[i] == "no") {
    label[i] <- "False Negative"
  }
}

# Step 2: Combine into final prediction table
prediction_table <- data.frame(
  plot = csv$plot,
  observed = true_labels,
  predicted = final_pred_class,
  predicted_prob = round(pred_probs, 3),
  classification = label,
  oakcount = csv$oakcount
)

# Preview
head(prediction_table)
write.csv(prediction_table, "G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\Plot_Level_Datasets\\PredTable_Shrub.csv", row.names = FALSE)
