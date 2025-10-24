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

csv = read.csv("G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\responsematrix_04-25-25.csv")



################# CONIFER VARIABLES ############################################
train = csv[,c(6, 8:115)]
train = train[, !grepl("propor", names(train))]


colnames(train)[1] = "y"
train$y[train$y > 0] = 1
train$y[train$y < 1] = 0

train$success = train$y
train = train[,-1]
train$success = as.numeric(train$success)
train$plot_id = 1:length(train[,1])
n_plots = length(train[,1])

get_combinations <- function(variables, n) {
  unlist(lapply(1:n, function(x) combn(variables, x, simplify = FALSE)), recursive = FALSE)
}


plan(multisession)  # or use multicore if you're on Unix

# Setup
variables <- colnames(train)[!colnames(train) %in% c("success", "plot_id")]
combo3 = get_combinations(variables, 2)
combinations <- get_combinations(variables, 4)
combinations = combinations[c(length(combo3):length(combinations))]
rm(combo3)


# Pre-sample a fixed set of combination indices for each bootstrap
n_boots <- 100  # or more
n_combos <- 50000  # combos per bootstrap
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
      probs <- predict(model, type = "response")
      auc <- tryCatch({
        if (length(unique(train_sub$success)) > 1) {
          as.numeric(pROC::auc(train_sub$success, probs))
        } else NA
      }, error = function(e) NA)
      
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
        aic = aic,
        max_corr = max_corr,
        auc = auc,
        stringsAsFactors = FALSE
      )
    }
  }
  
  dplyr::bind_rows(results)
}


# Ensure future can export large objects
options(future.globals.maxSize = +Inf)

# Run safely — furrr will now export combinations and train
all_results <- future_map_dfr(1:n_boots, function(i) {
  run_bootstrap_static(
    seed = i,
    combo_indices_subset = combo_indices[[i]],
    combo_indices_id = i
  )
}, .options = furrr_options(seed = TRUE))


variable_votes <- all_results %>%
  separate_rows(vars, sep = ", ") %>%
  group_by(vars) %>%
  summarize(count = n()) %>%
  arrange(desc(count))


write.csv(variable_votes, "G:\\YOSE_Regen_Analysis\\Runs\\April18_run_CURRENT\\VariableVotes_Conifer.csv", row.names = FALSE)

################################## Find best variables #########################

# Set thresholds for defining "top" models
top_n_per_boot <- 50  # number of top models per bootstrap
auc_thresh <- 0.75     # define high AUC
corr_thresh <- 0.6     # acceptable collinearity

# Unnest all variables into rows
vars_long <- all_results %>%
  dplyr::select(seed, aic, auc, max_corr, vars) %>%
  unnest_longer(vars) %>%
  rename(variable = vars)

# Count appearances overall
var_counts <- vars_long %>%
  count(variable, name = "n_total")


top_aic_vars <- all_results %>%
  group_by(seed) %>%
  slice_min(aic, n = top_n_per_boot, with_ties = FALSE) %>%
  unnest_longer(vars) %>%
  count(vars, name = "n_top_aic") %>%
  rename(variable = vars)
n_distinct(top_aic_vars$variable) 

top_auc_vars <- all_results %>%
  filter(auc >= auc_thresh) %>%
  unnest_longer(vars) %>%
  count(vars, name = "n_high_auc") %>%
  rename(variable = vars)

low_corr_vars <- all_results %>%
  filter(max_corr <= corr_thresh | is.na(max_corr)) %>%
  unnest_longer(vars) %>%
  count(vars, name = "n_low_corr") %>%
  rename(variable = vars)


var_summary <- var_counts %>%
  full_join(top_auc_vars, by = "variable") %>%
  full_join(low_corr_vars, by = "variable") %>%
  replace(is.na(.), 0) %>%
  arrange(desc(n_total))

############################ Plot the Vars ####################################
library(ggplot2)

top_vars <- var_summary %>%
  arrange(desc(n_total)) %>%
  slice_head(n = 20) %>%
  pivot_longer(cols = c(n_total, n_high_auc, n_low_corr),
               names_to = "category", values_to = "count")

ggplot(top_vars, aes(x = reorder(variable, -count), y = count, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Variable Importance Across Bootstrapped Models",
       x = "LiDAR Variable",
       y = "Count of Appearances",
       fill = "Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

########################## Get the Statistics ##################################

library(bestglm)


# first strata 1 to 4
regen = dplyr::select(train, c(linear_75to95, iqr_90to99, b11w, ndviw, success))
#regen = regen[c(2:9, 11:30),]
# hist(regen$regencount)
# regen = na.omit(regen)
# 
# colnames(regen)[1] = "y"
# 
# regen$y[regen$y > 0] = 1
# regen$y[regen$y < 1] = 0
# 
# regen$success = regen$y
# #oaks$failure = 1-oaks$y
# regen = regen[,-1]
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
  glmmodel = glm(success~ linear_75to95 + iqr_90to99 + b11w + ndviw, data = df, family = binomial)
  
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






























