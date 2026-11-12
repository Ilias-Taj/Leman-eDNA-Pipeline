################################# PACKAGES USED ################################
library(readr)
library(Biostrings)
library(BiocManager)
library(tidyverse)
library(glmmTMB)
################################################################################

################################################################################
######################## ANALYSIS ON WATER 18S TAXONOMY ########################
################################################################################

setwd("C:/Users/info/OneDrive/Desktop/Genorobotics")
water_taxonomy_18S <- read_csv("Water_eDNA_18S_COI_14_01_26_20261004_161409/taxonomy_summary/18S/pr2/comprehensive_taxonomy_18S.csv")
View(water_taxonomy_18S)

# 1. Check the dependence structure between confidence levels and other variables (most importantly, reads per otu)
read_per_otu <- readDNAStringSet("Water_eDNA_18S_COI_14_01_26_20261004_161409/temp_clustering/consensus_18S_clean.fasta")
read_per_otu
names(read_per_otu)

read_per_otu_df <- data.frame(
  centroid = sub("centroid=", "", sub(";.*", "", names(read_per_otu))),
  n_reads = as.numeric(sub(".*seqs=", "", names(read_per_otu)))
)
head(read_per_otu_df)
# Now, we have to merge the two datasets together. 
water_taxonomy_18S <- water_taxonomy_18S %>%
  left_join(read_per_otu_df, by = c("OTU_ID" = "centroid"))
# Now, we have all the information we need.

# We can do some exploratory analysis first.
# Let us start with the dependent variable, PR2_Domain_Conf.

mean(water_taxonomy_18S$PR2_Domain_Conf == 1, na.rm = TRUE)
mean(water_taxonomy_18S$PR2_Domain_Conf == 0, na.rm = TRUE)
mean(water_taxonomy_18S$PR2_Domain_Conf >= 0.8, na.rm = TRUE)

summary(water_taxonomy_18S$PR2_Domain_Conf)
sort(table(water_taxonomy_18S$PR2_Domain_Conf), decreasing = TRUE)[1:20]

# We have a lot of really large values, with 62% of observations exactly equal to 1 and
# 99.3% of observations that have a confidence level larger than the threshold (0.8).
# Therefore, it doesn't even make sense to analyse it.

# What about the Supergroup category?
mean(water_taxonomy_18S$PR2_Supergroup_Conf == 1, na.rm = TRUE)
mean(water_taxonomy_18S$PR2_Supergroup_Conf == 0, na.rm = TRUE)
mean(water_taxonomy_18S$PR2_Supergroup_Conf >= 0.8, na.rm = TRUE)
# We have 26.76% of observations equal to 1, which is still quite a lot, but we only
# have 48.59% of observations below the target confidence level.

summary(water_taxonomy_18S$PR2_Supergroup_Conf)
sort(table(water_taxonomy_18S$PR2_Supergroup_Conf), decreasing = TRUE)[1:20]
# Many of the most common values lie below 0.80.
# We should consider a mixture of two distributions (some 1-inflated model, for instance).
# As a first thing, we can look at the corresponding histogram:
ggplot(water_taxonomy_18S, aes(x = PR2_Supergroup_Conf,
                               fill = PR2_Supergroup)) +
  geom_histogram(bins = 10, colour = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0.80, linetype = "dashed",
             colour = "red", linewidth = 1) +
  labs(
    x = "PR2 domain confidence",
    y = "Number of OTUs",
    title = "Distribution of PR2 domain confidence",
    subtitle = "Using a threshold of 80% as reference",
    fill = "PR2 Domain"
  ) +
  theme_minimal()
# Does the distribution change for different domains? Let us have a look:
ggplot(water_taxonomy_18S[!is.na(water_taxonomy_18S$PR2_Supergroup_Conf), ], 
       aes(x = PR2_Supergroup_Conf, fill = PR2_Supergroup)) +
  geom_histogram(bins = 30, colour = "black", linewidth = 0.5) +
  geom_vline(
    xintercept = 0.80,
    linetype = "dashed",
    colour = "red"
  ) +
  facet_wrap(~ PR2_Supergroup) +
  labs(
    x = "PR2 supergroup confidence",
    y = "Number of OTUs",
    title = "Distribution of PR2 supergroup confidence by domain",
    fill = "PR2 Supergroup"
  ) +
  theme_minimal()
# The supergroups with a sample size suifficiently large to analyse them are Amoebozoa,
# excavata, obazoa and TSAR. They all show slightly different distributions, so it
# might be worth it to add the supergroup factor in R.
# Due to the general skewness of the distributions, it might be useful to introduce a ù
# log-transformation of the dependent variable:
water_taxonomy_18S$log_Supergroup_Conf <- log(water_taxonomy_18S$PR2_Supergroup_Conf)
# What is the histogram going to look like now?
ggplot(water_taxonomy_18S, aes(x = log_Supergroup_Conf,
                               fill = PR2_Supergroup)) +
  geom_histogram(bins = 10, colour = "black", linewidth = 0.5) +
  geom_vline(xintercept = log(0.80), linetype = "dashed",
             colour = "red", linewidth = 1) +
  labs(
    x = "PR2 Supergroup confidence",
    y = "Number of OTUs",
    title = "Distribution of the logarithm of PR2 Supergroup confidence",
    subtitle = "Using a threshold of log(0.80) as reference",
    fill = "PR2 Supergroup"
  ) +
  theme_minimal()
# A strong suspect from experts (Ilias) is that the most plausible explanation for the 
# small confidence probability is the number of reads for each otu. We can 
# check whether this is plausible at the Supergroup level or not, using a scatterplot:
ggplot(water_taxonomy_18S,
       aes(x = n_reads, y = PR2_Supergroup_Conf)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_minimal()
# We see that there are two units that have a really large value for x and that 
# influence heavily the fitted smooth.
# Luckily for us, the number of reads per OTU are already in decreasing order, so 
# the first two units are the ones we should remove.
ggplot(water_taxonomy_18S[-c(1,2),],
       aes(x = n_reads, y = PR2_Supergroup_Conf)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_minimal()

# Since we have very large values for the number of reads, we can use the logarithm:
water_taxonomy_18S$log10_nreads <- log10(water_taxonomy_18S$n_reads)
# Now, we can start with the modelling procedure:

# point a: linear model is not good, strongly skewed distribution;
# point b: the log-transformation of the dependent variable doesn't work out well,
#          the skewness stays there;
# point c: since the response variable is a probability ranging from 0 to 1, we could use
#          a Beta distribution, but it doesn't handle observations at the boundary;
# point d: to solve this issue, we could try to fit a mixture of models:
#          a Binomial GLM to analyse the factors that determine why the probability goes to 1,
#          a Beta GLM to analyse what factors influence the behaviour of units in the interior of the support.

# For our analysis, we should keep observations with valid confidence and read counts values:
model_18S_data_Supergroup <- water_taxonomy_18S %>%
  filter(
    !is.na(PR2_Supergroup_Conf),
    !is.na(n_reads),
    n_reads > 0,
    PR2_Supergroup_Conf >= 0,
    PR2_Supergroup_Conf <= 1
  ) %>%
  mutate(
    # Create the binary response: 1 if the confidence is exactly 1
    at_max = as.integer(PR2_Supergroup_Conf == 1),
    
    # Log-transform the predictor, not the response
    log_reads = log10(n_reads)
  )

# We can quickly check the data
summary(model_18S_data_Supergroup$PR2_Supergroup_Conf)
table(model_18S_data_Supergroup$at_max)
mean(model_18S_data_Supergroup$at_max)
# Now, we are ready to model the data, starting from the binomial distribution.
mod_supergroup_binomial <- glm(
  at_max ~ log_reads + factor(PR2_Supergroup) + Total_Abundance,
  data = model_18S_data_Supergroup,
  family = binomial(link = "logit")
)
# Let's look at the summary for the model
summary(mod_supergroup_binomial)
# We can see that:
# - the total abundance is not significant;
# - the supegroup factor is (very slightly) significant only for two categories,
#   Archaeplastida and Excavata, that have a very small sample size, so we can consider 
#   removing it as well.
# This leads to these 2 new models:
mod_supergroup_bin_2 <- glm(
  at_max ~ log_reads + factor(PR2_Supergroup),
  data = model_18S_data_Supergroup,
  family = binomial(link = "logit")
)
# Let's look at the summary for the model
summary(mod_supergroup_bin_2)

mod_supergroup_bin_3 <- glm(
  at_max ~ log_reads,
  data = model_18S_data_Supergroup,
  family = binomial(link = "logit")
)
# Let's look at the summary for the model
summary(mod_supergroup_bin_3)
# In terms of AIC, model 2 seems to be the best, so we can define it as:
mod_supergroup_bin_opt <- mod_supergroup_bin_2
# Hpw do interpret this?
exp(coef(mod_supergroup_bin_opt)["log_reads"])
# Since the odds ratio is larger than 1, it means that an increase in the number 
# of reads is associated with an expected larger possibility of reaching maximum
# confidence.

# Now, let us take care of the "Beta" part of the model:
model_18S_beta_data_Supergroup <- model_18S_data_Supergroup %>%
  filter(
    at_max != 1
  )

model_supergroup_beta <- glmmTMB(
  PR2_Supergroup_Conf ~ log_reads + factor(PR2_Supergroup) + Total_Abundance,
  data = model_18S_beta_data_Supergroup,
  family = beta_family(link = "logit")
)
summary(model_supergroup_beta)
# Once again, the Total_Abundance variable is not significant, so we can fit a model
# without that variable:
model_supergroup_beta_2 <- glmmTMB(
  PR2_Supergroup_Conf ~ log_reads + factor(PR2_Supergroup),
  data = model_18S_beta_data_Supergroup,
  family = beta_family(link = "logit")
)
summary(model_supergroup_beta_2)
# This has a lower AIC and only contains variables with significant effects, so 
# we can choose this as an optimal plot

# We can develop a function that is able to analyse the variable of interest exists,
# the percentage of observations equal to 0 and to 1 and, if it is lower than 0.10,
# transforms those observations to 0.01 and to 0.99 respectively:

preprocess_boundary <- function(var_interest, data = water_taxonomy_18S,
                                threshold = 0.10) {
  
  # Check that the variable of interest exists
  if (!var_interest %in% names(data)) {
    stop("Variable '", var_interest, "' was not found in the dataset.")
  }
  
  x <- data[[var_interest]]
  
  # Exclude missing observations from the denominator
  valid <- !is.na(x)
  n <- sum(valid)
  
  if (n == 0) {
    stop("The variable of interest contains no non-missing observations.")
  }
  
  # Calculate the percentages at each boundary
  percentage_zero <- mean(x[valid] == 0)
  percentage_one  <- mean(x[valid] == 1)
  percentage_boundary <- percentage_zero + percentage_one
  
  # Report the percentages separately
  cat("Variable of interest:", var_interest, "\n")
  cat("Percentage equal to 0:",
      round(100 * percentage_zero, 2), "%\n")
  cat("Percentage equal to 1:",
      round(100 * percentage_one, 2), "%\n")
  cat("Combined boundary percentage:",
      round(100 * percentage_boundary, 2), "%\n\n")
  
  # Apply the appropriate rule
  if (percentage_boundary < threshold) {
    
    data[[var_interest]][valid & x == 0.00] <- 0.01
    data[[var_interest]][valid & x == 1.00] <- 0.99
    
    message(
      "Low percentage of observations on the support boundary: ",
      "you may proceed with modelling"
    )
    
  } else if (percentage_boundary > threshold) {
    
    message(
      "Large percentage of observations on the support boundary: ",
      "consider a mixture of two models"
    )
    
  } else {
    
    message(
      "The percentage of observations on the support boundary ",
      "is exactly at the threshold: consider your modelling approach"
    )
  }

}
preprocess_boundary("PR2_Supergroup_Conf")
# This works. We can do it for every classification level:
preprocess_boundary("PR2_Division_Conf")
preprocess_boundary("PR2_Class_Conf")
preprocess_boundary("PR2_Order_Conf")
preprocess_boundary("PR2_Family_Conf")
preprocess_boundary("PR2_Genus_Conf")
preprocess_boundary("PR2_Species_Conf")

# ow, we can proceed by fitting a two-part model for each variable, except from the Species confidence level:

fit_2part_model <- function(data = water_taxonomy_18S, response, covariates) {
  
  # 1. Check whether the inputs are valid
  if (!is.data.frame(data))
    stop("'data' must be a data frame.")
  
  if (length(response) != 1 || !response %in% names(data))
    stop("Response variable not in dataset.")
  
  if (!all(covariates %in% names(data)))
    stop("Covariates not found: ",
         paste(setdiff(covariates, names(data)), collapse = ", "))
  
  if (response %in% covariates)
    stop("The response cannot also be a covariate.")
  
  if (!is.numeric(data[[response]]))
    stop("Response must be numeric and bounded between 0 and 1.")
  
  if (any(data[[response]] < 0 | data[[response]] > 1,
          na.rm = TRUE))
    stop("Response must lie between 0 and 1.")
  
  if (!"log10_nreads" %in% covariates)
    stop("The logarithm of the number of reads must be included in covariates.")
  
  if (any(data$n_reads <= 0, na.rm = TRUE))
    stop("'n_reads' must be positive to calculate log10(n_reads).")

  # 2. Create transformed read-count predictor
  model_covariates <- covariates
  
  # 3. Keep complete cases
  required_vars <- unique(c(response, model_covariates))
  
  model_data <- data[
    complete.cases(data[, required_vars, drop = FALSE]),
    ,
    drop = FALSE
  ]
  
  if (nrow(model_data) == 0)
    stop("No complete observations remain.")
  
  # 4. Create the binary indicator
  model_data$at_max <- as.integer(model_data[[response]] == 1)
  
  y <- model_data[[response]]
  
  cat("Response:", response, "\n")
  cat("Complete observations:", nrow(model_data), "\n")
  cat("Equal to 0:", sum(y == 0), "\n")
  cat("Equal to 1:", sum(y == 1), "\n")
  cat("Strictly between 0 and 1:", sum(y > 0 & y < 1), "\n")
  
  # 5. Construct formulas
  formula_binomial <- reformulate(
    model_covariates, response = "at_max"
  )
  
  formula_beta <- reformulate(
    model_covariates, response = response
  )
  
  # 6. Fit binomial model
  if (length(unique(model_data$at_max)) < 2)
    stop("Binomial model requires both 0 and 1 outcomes.")
  
  model_binomial <- glm(
    formula_binomial,
    data = model_data,
    family = binomial(link = "logit")
  )
  
  # 7. Fit beta regression on interior observations
  beta_data <- model_data[
    model_data[[response]] > 0 &
      model_data[[response]] < 1,
    ,
    drop = FALSE
  ]
  
  if (nrow(beta_data) < 2)
    stop("Too few observations strictly between 0 and 1.")
  
  model_beta <- glmmTMB::glmmTMB(
    formula_beta,
    data = beta_data,
    family = glmmTMB::beta_family(link = "logit")
  )
  
  # 8. Report n_reads effects
  cat("\n========== EFFECT OF READ COUNT ==========\n")
  
  # Binomial model
  binom_coef <- summary(model_binomial)$coefficients
  
  b_name <- "log10_n_reads"
  
  if (b_name %in% rownames(binom_coef)) {
    b1 <- binom_coef[b_name, "Estimate"]
    se1 <- binom_coef[b_name, "Std. Error"]
    p1 <- binom_coef[b_name, "Pr(>|z|)"]
    
    cat("\nBINOMIAL MODEL: probability of confidence = 1\n")
    cat("Estimate (log-odds):", b1, "\n")
    cat("Odds ratio for a 10-fold increase in reads:",
        exp(b1), "\n")
    cat("95% CI for odds ratio:",
        exp(b1 + c(-1, 1) * 1.96 * se1), "\n")
    cat("P-value:", p1, "\n")
  }
  
  # Beta model
  beta_coef <- summary(model_beta)$coefficients$cond
  
  if (b_name %in% rownames(beta_coef)) {
    b2 <- beta_coef[b_name, "Estimate"]
    se2 <- beta_coef[b_name, "Std. Error"]
    p2 <- beta_coef[b_name, "Pr(>|z|)"]
    
    cat("\nBETA MODEL: interior confidence values\n")
    cat("Estimate (logit scale):", b2, "\n")
    cat("P-value:", p2, "\n")
    
    # Average predictions before and after a 10-fold increase
    reference <- beta_data
    increased <- reference
    increased$log10_n_reads <- increased$log10_n_reads + 1
    
    pred_before <- predict(
      model_beta,
      newdata = reference,
      type = "response",
      re.form = NA
    )
    
    pred_after <- predict(
      model_beta,
      newdata = increased,
      type = "response",
      re.form = NA
    )
    
    delta <- pred_after - pred_before
    
    cat("Mean predicted confidence before:", mean(pred_before), "\n")
    cat("Mean predicted confidence after 10-fold increase:",
        mean(pred_after), "\n")
    cat("Average change in expected confidence:", mean(delta), "\n")
    cat("95% interval of individual prediction changes:",
        quantile(delta, c(0.025, 0.975)), "\n")
  }
  
  return(list(
    model_binomial = model_binomial,
    model_beta = model_beta,
    data_binomial = model_data,
    data_beta = beta_data,
    formula_binomial = formula_binomial,
    formula_beta = formula_beta
  ))
}

# Let's see if this works:
results_supergroup <- 
  fit_2part_model(response = "PR2_Supergroup_Conf", covariates = c("log10_nreads", "PR2_Supergroup"))
results_supergroup
# For the supergroup confidence, we have two models:
# - Binomial model, for which the logarithm of the number of reads has significance:
supergroup_bin_sig <- summary(results_supergroup$model_binomial)$coefficients[
  "log10_nreads", "Pr(>|z|)"
  ]
supergroup_bin_sig
#   and for which we can compute a 95% Wald CI as follows:
#   1) First, we extract the estimate of the coefficient and the corresponding standard error
binom_supergroup_table <- summary(results_supergroup$model_binomial)$coefficients
est_sg <- binom_supergroup_table["log10_nreads", "Estimate"]
se_sg <- binom_supergroup_table["log10_nreads", "Std. Error"]
#   2) Then, we compute the 95% CI on the log-odds scale
supergroup_bin_CI <- c("LB" = est_sg-qnorm(0.975)*se_sg,
                       "Est" = est_sg, 
                       "UB" = est_sg+qnorm(0.975)*se_sg)
#   3) Taking the exponent, we can convert it to odds-ratio scale
supergroup_bin_OR_ci <- exp(supergroup_bin_CI)
supergroup_bin_OR_ci

# - Beta regression, for which we can do the same thing as before
supergroup_beta_sig <- summary(results_supergroup$model_beta)$coefficients$cond[
  "log10_nreads", "Pr(>|z|)"
  ]
supergroup_beta_sig

beta_supergroup_table <- summary(results_supergroup$model_beta)$coefficients$cond
est2_sg <- beta_supergroup_table["log10_nreads", "Estimate"]
se2_sg <- beta_supergroup_table["log10_nreads", "Std. Error"]
# Notice that in this case the estimate represents the change in the log-odds of the
# expected confidence level following a unit increase in log10_nreads.
supergroup_beta_CI <- c("LB" = est2_sg-qnorm(0.975)*se2_sg,
                        "Est" = est2_sg, 
                        "UB" = est2_sg+qnorm(0.975)*se2_sg)
supergroup_beta_CI

# We can do the same for every other classification level, except for Species, for which 
# the number of observations equal to 1 is low enough to treat that directly with a beta regression.
# DIVISION CONFIDENCE
results_division <- 
  fit_2part_model(response = "PR2_Division_Conf", covariates = c("log10_nreads", "PR2_Division"))
results_division
# For the division confidence, we have two models:
# - Binomial model, for which the logarithm of the number of reads has significance:
division_bin_sig <- summary(results_division$model_binomial)$coefficients[
  "log10_nreads", "Pr(>|z|)"
]
division_bin_sig
#   and for which we can compute a 95% Wald CI as follows:
#   1) First, we extract the estimate of the coefficient and the corresponding standard error
binom_division_table <- summary(results_division$model_binomial)$coefficients
est_div <- binom_division_table["log10_nreads", "Estimate"]
se_div <- binom_division_table["log10_nreads", "Std. Error"]
#   2) Then, we compute the 95% CI on the log-odds scale
division_bin_CI <- c("LB" = est_div-qnorm(0.975)*se_div,
                       "Est" = est_div, 
                       "UB" = est_div+qnorm(0.975)*se_div)
#   3) Taking the exponent, we can convert it to odds-ratio scale
division_bin_OR_ci <- exp(division_bin_CI)
division_bin_OR_ci
# Beta regression
division_beta_sig <- summary(results_division$model_beta)$coefficients$cond[
  "log10_nreads", "Pr(>|z|)"
]
division_beta_sig

beta_division_table <- summary(results_division$model_beta)$coefficients$cond
est2_div <- beta_division_table["log10_nreads", "Estimate"]
se2_div <- beta_division_table["log10_nreads", "Std. Error"]
# Notice that in this case the estimate represents the change in the log-odds of the
# expected confidence level following a unit increase in log10_nreads.
division_beta_CI <- c("LB" = est2_div-qnorm(0.975)*se2_div,
                      "Est" = est2_div, 
                      "UB" = est2_div+qnorm(0.975)*se2_div)
division_beta_CI


# CLASS CONFIDENCE (in this case, we have 42 unique levels for the factor PR2_Class,
#                   so an average of around 3.5 observation per level, hence we can remove
#                   the class factor from the analysis).
# (Total_abundance has been considered, but it is highly non-significant).
results_class <- 
  fit_2part_model(response = "PR2_Class_Conf", covariates = c("log10_nreads"))
results_class
# For the class confidence, we have two models:
# - Binomial regression
class_bin_sig <- summary(results_class$model_binomial)$coefficients[
  "log10_nreads", "Pr(>|z|)"
]
class_bin_sig

binom_class_table <- summary(results_class$model_binomial)$coefficients
est_class <- binom_class_table["log10_nreads", "Estimate"]
se_class <- binom_class_table["log10_nreads", "Std. Error"]
class_bin_CI <- c("LB" = est_class-qnorm(0.975)*se_class,
                     "Est" = est_class, 
                     "UB" = est_class+qnorm(0.975)*se_class)
class_bin_OR_ci <- exp(class_bin_CI)
class_bin_OR_ci
# - Beta regression
class_beta_sig <- summary(results_class$model_beta)$coefficients$cond[
  "log10_nreads", "Pr(>|z|)"
]
class_beta_sig

beta_class_table <- summary(results_class$model_beta)$coefficients$cond
est2_class <- beta_class_table["log10_nreads", "Estimate"]
se2_class <- beta_class_table["log10_nreads", "Std. Error"]
# Notice that in this case the estimate represents the change in the log-odds of the
# expected confidence level following a unit increase in log10_nreads.
class_beta_CI <- c("LB" = est2_class-qnorm(0.975)*se2_class,
                   "Est" = est2_class, 
                   "UB" = est2_class+qnorm(0.975)*se2_class)
class_beta_CI

# ORDER CONFIDENCE
results_order <- 
  fit_2part_model(response = "PR2_Order_Conf", covariates = c("log10_nreads"))
results_order
# For the order confidence, we have two models:
# - Binomial model
order_bin_sig <- summary(results_order$model_binomial)$coefficients[
  "log10_nreads", "Pr(>|z|)"
]
order_bin_sig

binom_order_table <- summary(results_order$model_binomial)$coefficients
est_ord <- binom_order_table["log10_nreads", "Estimate"]
se_ord <- binom_order_table["log10_nreads", "Std. Error"]
order_bin_CI <- c("LB" = est_ord-qnorm(0.975)*se_ord,
                  "Est" = est_ord, 
                  "UB" = est_ord+qnorm(0.975)*se_ord)
order_bin_OR_ci <- exp(order_bin_CI)
order_bin_OR_ci
# Beta regression
order_beta_sig <- summary(results_order$model_beta)$coefficients$cond[
  "log10_nreads", "Pr(>|z|)"
]
order_beta_sig

beta_order_table <- summary(results_order$model_beta)$coefficients$cond
est2_ord <- beta_order_table["log10_nreads", "Estimate"]
se2_ord <- beta_order_table["log10_nreads", "Std. Error"]
# Notice that in this case the estimate represents the change in the log-odds of the
# expected confidence level following a unit increase in log10_nreads.
order_beta_CI <- c("LB" = est2_ord-qnorm(0.975)*se2_ord,
                   "Est" = est2_ord, 
                   "UB" = est2_ord+qnorm(0.975)*se2_ord)
order_beta_CI

# FAMILY CONFIDENCE
results_family <- 
  fit_2part_model(response = "PR2_Family_Conf", covariates = c("log10_nreads"))
results_family
# For the family confidence, we have two models:
# - Binomial model
family_bin_sig <- summary(results_family$model_binomial)$coefficients[
  "log10_nreads", "Pr(>|z|)"
]
family_bin_sig

binom_family_table <- summary(results_family$model_binomial)$coefficients
est_fam <- binom_family_table["log10_nreads", "Estimate"]
se_fam <- binom_family_table["log10_nreads", "Std. Error"]
family_bin_CI <- c("LB" = est_fam-qnorm(0.975)*se_fam,
                  "Est" = est_fam, 
                  "UB" = est_fam+qnorm(0.975)*se_fam)
family_bin_OR_ci <- exp(family_bin_CI)
family_bin_OR_ci
# Beta regression
family_beta_sig <- summary(results_family$model_beta)$coefficients$cond[
  "log10_nreads", "Pr(>|z|)"
]
family_beta_sig

beta_family_table <- summary(results_family$model_beta)$coefficients$cond
est2_fam <- beta_family_table["log10_nreads", "Estimate"]
se2_fam <- beta_family_table["log10_nreads", "Std. Error"]
# Notice that in this case the estimate represents the change in the log-odds of the
# expected confidence level following a unit increase in log10_nreads.
family_beta_CI <- c("LB" = est2_fam-qnorm(0.975)*se2_fam,
                   "Est" = est2_fam, 
                   "UB" = est2_fam+qnorm(0.975)*se2_fam)
family_beta_CI

# GENUS CONFIDENCE
results_genus <- 
  fit_2part_model(response = "PR2_Genus_Conf", covariates = c("log10_nreads"))
results_genus
# For the genus confidence, we have two models:
# - Binomial model
genus_bin_sig <- summary(results_genus$model_binomial)$coefficients[
  "log10_nreads", "Pr(>|z|)"
]
genus_bin_sig

binom_genus_table <- summary(results_genus$model_binomial)$coefficients
est_gen <- binom_genus_table["log10_nreads", "Estimate"]
se_gen <- binom_genus_table["log10_nreads", "Std. Error"]
genus_bin_CI <- c("LB" = est_gen-qnorm(0.975)*se_gen,
                   "Est" = est_gen, 
                   "UB" = est_gen+qnorm(0.975)*se_gen)
genus_bin_OR_ci <- exp(genus_bin_CI)
genus_bin_OR_ci
# Beta regression
genus_beta_sig <- summary(results_genus$model_beta)$coefficients$cond[
  "log10_nreads", "Pr(>|z|)"
]
genus_beta_sig

beta_genus_table <- summary(results_genus$model_beta)$coefficients$cond
est2_gen <- beta_genus_table["log10_nreads", "Estimate"]
se2_gen <- beta_genus_table["log10_nreads", "Std. Error"]
# Notice that in this case the estimate represents the change in the log-odds of the
# expected confidence level following a unit increase in log10_nreads.
genus_beta_CI <- c("LB" = est2_gen-qnorm(0.975)*se2_gen,
                    "Est" = est2_gen, 
                    "UB" = est2_gen+qnorm(0.975)*se2_gen)
genus_beta_CI

# SPECIES CONFIDENCE
# As we said before, the species confidence can be treated directly with a Beta 
# regression model. 

# 1. Prepare the data
species_model_data <- water_taxonomy_18S %>%
  filter(
    !is.na(PR2_Species_Conf),
    !is.na(log10_nreads)
  ) %>%
  mutate(
    species_conf_adj = if_else(
      PR2_Species_Conf == 1,
      0.99,
      PR2_Species_Conf
    )
  ) %>%
  # Beta regression requires values strictly between 0 and 1
  filter(
    species_conf_adj > 0,
    species_conf_adj < 1
  )

# Report the number of observations adjusted or excluded
cat("Observations used:", nrow(species_model_data), "\n")
cat("Original values equal to 1:",
    sum(species_model_data$PR2_Species_Conf == 0.99), "\n")
cat("Original values equal to 0:",
    sum(water_taxonomy_18S$PR2_Species_Conf == 0, na.rm = TRUE), "\n")
# No values are outside the support of the beta model
# 2. Fit beta regression
model_species_beta <- glmmTMB(
  species_conf_adj ~ log10_nreads,
  data = species_model_data,
  family = beta_family(link = "logit")
)

# 3. Extract coefficient estimate, standard error and p-value
coef_table <- summary(model_species_beta)$coefficients$cond

estimate <- coef_table["log10_nreads", "Estimate"]
std_error <- coef_table["log10_nreads", "Std. Error"]
species_beta_sig <- coef_table["log10_nreads", "Pr(>|z|)"]
species_beta_sig
# 4. Calculate 95% Wald confidence interval
wald_ci <- estimate + c(-1, 1) * qnorm(0.975) * std_error

species_beta_CI <- c("LB" = wald_ci[1],
                     "Est" = estimate,
                     "UB" = wald_ci[2])
species_beta_CI

# 5. Report results
cat("\nBeta regression: PR2_Species_Conf\n")
cat("Coefficient estimate:", estimate, "\n")
cat("Standard error:", std_error, "\n")
cat("95% Wald CI:", wald_ci[1], "to", wald_ci[2], "\n")
cat("P-value:", format.pval(species_beta_sig, digits = 4), "\n")

if (format.pval(species_beta_sig, digits = 4) < 0.05) {
  cat("The coefficient is statistically significant at the 5% level.\n")
} else {
  cat("The coefficient is not statistically significant at the 5% level.\n")
}

# Now, we can summarize our results into a single dataframe.
nreads_result <- data.frame("Classification_Level" = c("Domain", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"),
                            "Binomial_OR_LB" = c(NA, supergroup_bin_OR_ci[1], division_bin_OR_ci[1], class_bin_OR_ci[1],
                                                 order_bin_OR_ci[1], family_bin_OR_ci[1], genus_bin_OR_ci[1], NA),
                            "Binomial_OR_Est" = c(NA, supergroup_bin_OR_ci[2], division_bin_OR_ci[2], class_bin_OR_ci[2],
                                                 order_bin_OR_ci[2], family_bin_OR_ci[2], genus_bin_OR_ci[2], NA),
                            "Binomial_OR_UB" = c(NA, supergroup_bin_OR_ci[3], division_bin_OR_ci[3], class_bin_OR_ci[3],
                                                 order_bin_OR_ci[3], family_bin_OR_ci[3], genus_bin_OR_ci[3], NA),
                            "Binomial_pval" = c(NA, supergroup_bin_sig, division_bin_sig, class_bin_sig,
                                                order_bin_sig, family_bin_sig, genus_bin_sig, NA),
                            "Beta_LB" = c(NA, supergroup_beta_CI[1], division_beta_CI[1], class_beta_CI[1],
                                          order_beta_CI[1], family_beta_CI[1], genus_beta_CI[1], species_beta_CI[1]),
                            "Beta_Est" = c(NA, supergroup_beta_CI[2], division_beta_CI[2], class_beta_CI[2],
                                          order_beta_CI[2], family_beta_CI[2], genus_beta_CI[2], species_beta_CI[2]),
                            "Beta_UB" = c(NA, supergroup_beta_CI[3], division_beta_CI[3], class_beta_CI[3],
                                          order_beta_CI[3], family_beta_CI[3], genus_beta_CI[3], species_beta_CI[3]),
                            "Beta_pval" = c(NA, supergroup_beta_sig, division_beta_sig, class_beta_sig,
                                            order_beta_sig, family_beta_sig, genus_beta_sig, family_beta_sig))

# Let us use this dataset to produce some interesting plots.
# How does the p-value changes across classification levels?
nreads_result <- nreads_result %>%
  mutate(
    bin_significance = case_when(
      is.na(Binomial_pval) ~ NA_character_,
      Binomial_pval < 0.05 ~ "Significant (p < 0.05)",
      TRUE ~ "Not significant (p >= 0.05)"
    ),
    beta_significance = case_when(
      is.na(Beta_pval) ~ NA_character_,
      Beta_pval < 0.05 ~ "Significant (p < 0.05)",
      TRUE ~ "Not significant (p >= 0.05)"
    ),
    Classification_Level = factor(
      Classification_Level,
      levels = unique(as.character(Classification_Level))
    )
  )
# For the binomial model:
ggplot(
  nreads_result,
  aes(x = Classification_Level, y = Binomial_pval, fill = bin_significance)
) +
  geom_col(width = 0.7, na.rm = TRUE) +
  geom_text(
    aes(label = ifelse(
      is.na(Binomial_pval),
      "",
      format.pval(Binomial_pval, digits = 3, eps = 0.001)
    )),
    vjust = -0.5,
    na.rm = TRUE,
    size = 3.5
  ) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "Significant (p < 0.05)" = "seagreen3",
      "Not significant (p >= 0.05)" = "indianred2"
    ),
    na.value = "grey80",
    name = "Significance"
  ) +
  scale_y_continuous(
    name = "P-value",
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Binomial model significance by classification level",
    x = "Classification level"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )
# We notice that there is a general increasing trend in terms of significance: i.e.,
# if we have a finer classification the number of reads has a less and less significant
# effect on the probability that the confidence level is equal to 100%. 
# This regression was not the main interest, but we can conclude that this loss of 
# significance may be due to the fact that the proportion of observations equal to 1
# becomes smaller and smaller.

# In fact, the main interest of this small research is the following result:
ggplot(
  nreads_result,
  aes(x = Classification_Level, y = Beta_pval, fill = beta_significance)
) +
  geom_col(width = 0.7, na.rm = TRUE) +
  geom_text(
    aes(label = ifelse(
      is.na(Beta_pval),
      "",
      format.pval(Beta_pval, digits = 3, eps = 0.001)
    )),
    vjust = -0.5,
    na.rm = TRUE,
    size = 3.5
  ) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "Significant (p < 0.05)" = "seagreen3",
      "Not significant (p >= 0.05)" = "indianred2"
    ),
    na.value = "grey80",
    name = "Significance"
  ) +
  scale_y_continuous(
    name = "P-value",
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Beta model significance by classification level",
    x = "Classification level"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )
# Here, we see that the number of reads is strongly significant for every classification 
# level, with a notable decreasing trend. So, the greater the detail, the more important
# the number of reads becomes.

# What about the estimates of the coefficients?
ggplot(
  nreads_result,
  aes(
    x = Classification_Level,
    y = Beta_Est,
    group = 1
  )
) +
  # Connect coefficient estimates
  geom_line(
    color = "steelblue",
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  
  # Dashed vertical confidence intervals with horizontal caps
  geom_errorbar(
    aes(
      ymin = Beta_LB,
      ymax = Beta_UB
    ),
    width = 0.15,
    linetype = "dashed",
    color = "steelblue",
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  
  # Point estimates
  geom_point(
    color = "steelblue",
    size = 2.8,
    na.rm = TRUE
  ) +
  
  # Reference line at zero
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  labs(
    title = "Beta regression coefficients across taxonomic levels",
    subtitle = "Points indicate estimates; dashed intervals indicate 95% Wald confidence intervals",
    x = "Classification level",
    y = "Beta coefficient estimate (logit scale)"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# As we can see, the value of the coefficient estimate (at the logit scale) remains
# almost constant and around the interval (0.50,0.75).


