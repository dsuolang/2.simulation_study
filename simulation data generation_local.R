#-------simulation data generation--------#



#--------required packages--------#

library(haven)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(simstudy)
library(Matrix)
library(ltm)
library(psych)
library(mice)
library(MASS)
library(norm)
library(VIM)
library(ggplot2)
library(naniar)
library(devtools)

options(scipen = 999)
setwd('/Users/dsuolang/Desktop/Study2')
set.seed(2024)



#--------prepare existing data--------#

create_binary_variables <- function(real_data) {
  binary_real_data <- data.frame(matrix(NA, nrow = nrow(real_data), ncol = 0))
  for (col in names(real_data)) {
    real_data[[col]] <- factor(real_data[[col]])
    levels(real_data[[col]]) <- gsub("[ /</>/-]", "_", levels(real_data[[col]]))
    binary_vars <-
      model.matrix( ~ . - 1, data = real_data[, col, drop = FALSE])
    colnames(binary_vars) <- paste0(col, "_", levels(real_data[[col]]))
    binary_real_data <- cbind(binary_real_data, binary_vars)
  }
  real_data <- cbind(real_data, binary_real_data)
  return(real_data)
}

#real_data<-read_dta('/nfs/turbo/isr-bwest1/dsuolang/iveware3.imputed_all.dta')
real_data <- read_dta('nhanes_clean2.dta')
real_data$mvpa_total_acc<-sqrt(real_data$modpa_total_acc + real_data$vigpa_total_acc)
real_data$modpa_total<-sqrt(real_data$modpa_total)
real_data$vigpa_total<-sqrt(real_data$vigpa_total)
real_data$modpa_zero<-ifelse(real_data$modpa_total==0,0,1)
real_data$vigpa_zero<-ifelse(real_data$vigpa_total==0,0,1)

cont_vars_mixed <-  # semi-continuous
  real_data[, c("modpa_total", "vigpa_total")]
cont_vars_norm <-  # continuous
  real_data[, c("mvpa_total_acc", "fPCA1", "bmi")]
cont_vars_unif <-  # uniform
  real_data[, c("age")]
cat_vars_binary <-  # binary
  real_data[, c("modpa_zero", "vigpa_zero", 
                "gender", "work", "insurance", "srvy_yr", "hypertension", "heartdiseases", "stroke", "cancers", "diabetes", "obesity")]
cat_vars_binary$srvy_yr<-ifelse(cat_vars_binary$srvy_yr=='201112',0,1)
cat_vars_nonbinary <-  # transform non-binary to binary
  real_data[, c("race", "edu", "marital", "poverty", "self_reported_health", "smoker", "alcohol_cat")]
cat_vars_nonbinary <- create_binary_variables(cat_vars_nonbinary)
cat_vars_binary <- cbind (cat_vars_binary, cat_vars_nonbinary[, 8:ncol(cat_vars_nonbinary)])
cat_vars_nonbinary<-NULL

real_data<-cbind(cont_vars_mixed, cont_vars_norm, cont_vars_unif, cat_vars_binary)



#--------simulate data --------#

compute_correlation_matrix <- function(df, continuous_indices, binary_indices) {
    cor_matrix <- matrix(NA, ncol = ncol(real_data), nrow = ncol(real_data),
                         dimnames = list(names(real_data), names(real_data)))
    for (i in 1:ncol(real_data)) {
      for (j in i:ncol(real_data)) {
        if (i %in% continuous_indices && j %in% continuous_indices) {
          # Both variables are continuous
          cor_matrix[i, j] <-
            cor_matrix[j, i] <- cor(real_data[[i]], real_data[[j]])
        } else if (i %in% binary_indices && j %in% binary_indices) {
          # Both variables are binary
          cor_matrix[i, j] <-
            cor_matrix[j, i] <-
            cor(real_data[[i]], real_data[[j]], method = "pearson")
        } else {
          # At least one variable is binary; place the binary one in the y position
          continuous_var_index <-
            ifelse(i %in% continuous_indices, i, j)
          binary_var_index <- ifelse(i %in% binary_indices, i, j)
          # Calculate biserial correlation using psych::biserial
          cor_matrix[i, j] <-
            cor_matrix[j, i] <-
            biserial(real_data[[continuous_var_index]], real_data[[binary_var_index]])
        }
      }
    }
    return(cor_matrix)
  }

generateSimDefinition <- function(cont_vars_mixed, cont_vars_norm, cont_vars_unif, cat_vars_binary) {
  def <- defData(varname = "placeholder", formula = 0)
  
  for (varname in colnames(cont_vars_norm)) {
    mu_real <- mean(cont_vars_norm[[varname]])
    sd_real <- sd(cont_vars_norm[[varname]])
    def <- defDataAdd(def, varname = varname, 
                      formula = mu_real, 
                      variance = sd_real^2,
                      dist = "normal")
  }
  
  for (varname in colnames(cont_vars_unif)) {
    min_val <- min(cont_vars_unif[[varname]])
    max_val <- max(cont_vars_unif[[varname]])
    range_str <- paste(min_val, max_val, sep = ";")
    def <- defDataAdd(def, varname = varname, 
                      formula = range_str,
                      dist = "uniform")
  }
  
  for (varname in colnames(cat_vars_binary)) {
    prob <- mean(cat_vars_binary[[varname]])
    def <- defDataAdd(def, varname = varname, 
                      formula = prob,
                      variance = prob * (1 - prob),
                      dist = "binary")
  }
  for (varname in colnames(cont_vars_mixed)) {
    mu_real <- mean(cont_vars_mixed[[varname]][cont_vars_mixed[[varname]] != 0])
    sd_real <- sd(cont_vars_mixed[[varname]][cont_vars_mixed[[varname]] != 0])
    print(cont_vars_mixed[[varname]][cont_vars_mixed[[varname]] != 0])
    def <- defDataAdd(def, varname = varname, 
                      formula = mu_real, 
                      variance = sd_real^2,
                      dist = "normal")
  }
  
  def <- def[!def$varname == "placeholder",]
  return(def)
}

generateSimData<-function(n, def, cor_matrix) {
  sim_data <- genCorFlex(n, def, cor_matrix)
  sim_data$modpa_total <- ifelse(sim_data$modpa_zero == 0, 0, sim_data$modpa_total)
  sim_data$vigpa_total <- ifelse(sim_data$vigpa_zero == 0, 0, sim_data$vigpa_total)
  return(sim_data)
}

continuous_indices <- 1:6
binary_indices <- 7:47
cor_matrix <- compute_correlation_matrix(real_data, continuous_indices, binary_indices) # multivariate correlation matrix
cor_matrix <- as.matrix(nearPD(cor_matrix, corr = TRUE, keepDiag = TRUE)$mat)

def<-generateSimDefinition(cont_vars_mixed, cont_vars_norm, cont_vars_unif, cat_vars_binary)
sim_data<-generateSimData(10000, def, cor_matrix)


sim_data$mvpa_total_acc<-abs(sim_data$mvpa_total_acc)
cor(sim_data$mvpa_total_acc, sim_data$gender)
hist(real_data$mvpa_total_acc)
hist(sim_data$mvpa_total_acc)
save.image()

#--------R-squared--------#

# level 1 information
model1 <- lm(mvpa_total_acc ~ . - bmi - age - modpa_zero - vigpa_zero - gender - work - insurance  - hypertension - heartdiseases - stroke - cancers - diabetes - obesity - race_Hispanic - race_NHblack - race_NHwhite - race_Other - edu_college_and_above - edu_highschool_ged - edu_less_than_highschool - edu_some_college - marital_living_with_partner - marital_married - marital_never_married - marital_separated_divorced_widowed - poverty_1_to_2 - poverty_2_to_4 - poverty_4_or_higher - poverty_less_than_1 - poverty_missing - self_reported_health_1 - self_reported_health_2 - self_reported_health_3 - self_reported_health_4 - self_reported_health_5 - smoker_current - smoker_former - smoker_not_smoker - alcohol_cat__1_per_week - alcohol_cat__4_per_week - alcohol_cat_1_4_per_week - alcohol_cat_missing, 
             data = real_data)
summary(model1)$r.squared

# level 2 information
model2 <- lm(mvpa_total_acc ~ . - bmi  - modpa_zero - vigpa_zero - hypertension - heartdiseases - stroke - cancers - diabetes - obesity - self_reported_health_1 - self_reported_health_2 - self_reported_health_3 - self_reported_health_4 - self_reported_health_5 - smoker_current - smoker_former - smoker_not_smoker - alcohol_cat__1_per_week - alcohol_cat__4_per_week - alcohol_cat_1_4_per_week - alcohol_cat_missing, 
             data = real_data)
summary(model2)$r.squared

# level 3 information
model3 <- lm(mvpa_total_acc ~ . - modpa_zero - vigpa_zero - hypertension - heartdiseases - stroke - cancers - diabetes - obesity,
             data = real_data)
summary(model3)$r.squared

# level 4 information
model4 <- lm(mvpa_total_acc ~ . - modpa_zero - vigpa_zero, data = real_data)
summary(model4)$r.squared

colnames(real_data)
#--------introducing missingness (MAR)--------#

mu.X <- c(1, 1)
Sigma.X <- matrix(c(1, 1, 1, 4), nrow = 2)
n <- 100
X.complete.cont <- mvrnorm(n, mu.X, Sigma.X)

lambda <- 0.5
X.complete.discr <- rpois(n, lambda)

n.cat <- 5
X.complete.cat <- rbinom(n, size=5, prob = 0.5)

X.complete <- data.frame(cbind(X.complete.cont, X.complete.discr, X.complete.cat))
X.complete[,4] <- as.factor(X.complete[,4])
levels(X.complete[,4]) <- c("F", "E", "D", "C", "B", "A")


source_url('https://raw.githubusercontent.com/R-miss-tastic/website/master/static/how-to/generate/amputation.R')
patterns <- matrix(0, nrow = ncol(sim_data), ncol = 2)
colnames(patterns) <- c("influenceMissingness", "subjectToMissingness")
variable_names <- colnames(sim_data)
patterns[colnames(sim_data) %in% c("age", "gender", "race_Hispanic", "race_NHwhite", 
                               "race_NHblack", "race_Other"), "influenceMissingness"] <- 1
patterns[colnames(sim_data) == "mvpa_total_acc", "subjectToMissingness"] <- 1

X.miss <-  produce_NA(sim_data, 
                      mechanism="MAR", 
                      perc.missing = 0.2, 
                      idx.incomplete = patterns[,2], 
                      idx.covariates = patterns[,1])
miss <- produce_NA(X.complete, mechanism="MAR", perc.missing = 0.2, idx.incomplete = c(1, 0, 0, 0), idx.covariates = c(0,1,0,0), weights.covariates = c(0,1,0,0))
X.complete
rownames(sim_data)<-NULL
rownames()colSums(is.na(X.miss$data.incomp))
cor_matrix["mvpa_total_acc", "age"]
sim_miss <- function(data, miss_perc) {
  alpha0 <- switch(
    miss_perc,
    "50" = -3.3610,
    # 5000
    "70" = -9.1990,
    # 7000
    "80" = -12.3995,
    # 8000
    "90" = -16.3890,
    # 9000
    "95" = -19.4500
  ) # 9500
  lambda <-
    alpha0 + cor_matrix["mvpa_total_acc", "age"] * data[[age]] + cor_matrix["mvpa_total_acc", "gender"] 
  e_lambda <- exp(lambda)
  lgs_lambda <- e_lambda / (1 + e_lambda)
  print(summary(lgs_lambda))
  # missingness indicator R
  set.seed(2024)
  data[[paste0("miss_mar_", miss_perc)]] <-
    as.numeric(lgs_lambda >= runif(nrow(data)))
  data <-
    data[,!(names(data) %in% c("lambda", "e_lambda", "lgs_lambda"))]
  return(data)
}
simulated_data <- sim_miss(sim_data, miss_perc = "50")
simulated_data <- sim_miss(simulated_data, miss_perc = "70")
simulated_data <- sim_miss(simulated_data, miss_perc = "80")
simulated_data <- sim_miss(simulated_data, miss_perc = "90")
simulated_data <- sim_miss(simulated_data, miss_perc = "95")



#--------inducing missingness (MCAR)--------#

generate_missing <- function(data, column_name, percentage) {
  n <- nrow(data)
  set.seed(2024)
  miss_indices <- sample(seq_len(n), size = percentage * n)
  data[[column_name]] <- 1
  data[[column_name]][miss_indices] <- 0
  return(data)
}

simulated_data <-
  generate_missing(simulated_data, "miss_mcar_70", 0.7)
simulated_data <-
  generate_missing(simulated_data, "miss_mcar_80", 0.8)
simulated_data <-
  generate_missing(simulated_data, "miss_mcar_90", 0.9)
simulated_data <-
  generate_missing(simulated_data, "miss_mcar_95", 0.95)

write_csv(simulated_data, 'simulated_data.csv')
