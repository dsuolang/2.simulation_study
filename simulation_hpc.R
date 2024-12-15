# Via simulation, this study will examine the impacts of four key factors on imputation performance: 
#  Missing data mechanism
#  Missing rate
#  Amount of shared auxiliary variables
#  Model misspecification. 
#  Modeling approach 
#  Total scenarios = 2 * 3 * 3 * 2 * 2 = 72
#  Last updated 10/23/2024

############################################
#--------SRMI--------#
#work_dir <- "/Users/dsuolang/Downloads/study2/srmi" # CHANGE TEST
work_dir <- "/nfs/turbo/isr-bwest1/dsuolang/study2/srmi"
############################################

############################################
#--------PMM--------#
#work_dir <- "/Users/dsuolang/Downloads/study2/pmm" # CHANGE TEST
#work_dir <- "/nfs/turbo/isr-bwest1/dsuolang/study2/pmm"
############################################
setwd(work_dir)



#--------required packages--------#
############################################
library(haven)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(devtools)
library(foreign)
library(nnet)
library(pscl)
library(doParallel)
library(readr)
library(stringr)
library(fastDummies)
library(vcd)
library(magrittr)
library(mice)
library(caret) 
library(pscl)
library(pROC)
options(scipen = 999)

#source("/Users/dsuolang/Downloads/study2/functions.R")
source("/nfs/turbo/isr-bwest1/dsuolang/study2/functions.R") # CHANGE TEST



#--------Draw samples from a synthetic population--------#
############################################
synthpop <- read_dta("/nfs/turbo/isr-bwest1/dsuolang/study2/synthpop.dta") # CHANGE TEST
#synthpop <- read_dta("testdata.dta")
summary(sqrt(synthpop$mvpa_total_acc))
synthpop$mvpa_total_acc <- synthpop$modpa_total_acc + synthpop$vigpa_total_acc
synthpop[c("mvpa_total_acc_sqrt", "modpa_total_sqrt", "vigpa_total_sqrt")] <-
  sqrt(synthpop[c("mvpa_total_acc", "modpa_total", "vigpa_total")])

synthpop <- synthpop %>%
  mutate(self_reported_health = as.character(self_reported_health)) %>%
  mutate(self_reported_health = recode(self_reported_health,
                                       `1` = "excellent",
                                       `2` = "very good",
                                       `3` = "good",
                                       `4` = "fair",
                                       `5` = "poor"))

synthpop <- synthpop %>% mutate(
  srvy_yr = as.factor(srvy_yr),activity_pattern = as.factor(activity_pattern)
)

summary(synthpop$mvpa_total_acc) # Variable of interest 1: MVPA duration
table(synthpop$activity_pattern)/nrow(synthpop)# Variable of interest 2: activity pattern

# creating pseudo variable
synthpop$activity_pattern <- as.factor(synthpop$activity_pattern)

synthpop <- synthpop %>%
  dummy_cols(select_columns = c("activity_pattern","race", "gender", "edu", "poverty", "self_reported_health"))  
synthpop$mvpa_total_acc_scaled <- scale(synthpop$mvpa_total_acc)
synthpop$age_scaled<- scale(synthpop$age)

z1 <- 0 +
  0.15 * synthpop$mvpa_total_acc_scaled -       
  0.3 * synthpop$activity_pattern_2 -         
  0.2 * synthpop$activity_pattern_3 -          
  0.1 * synthpop$activity_pattern_4 +          
  0.1 * synthpop$age_scaled +                   
  0.1 * synthpop$gender_1 -                    
  0.1 * synthpop$race_Hispanic -               
  0.15 * synthpop$race_NHblack -                
  0.2 * synthpop$race_Other +                   
  0.3 * synthpop$`poverty_1 to 2` +            
  0.2 * synthpop$`poverty_2 to 4` +           
  0.1 * synthpop$`poverty_4 or higher` -
  0.1 * synthpop$poverty_missing  +
  rlnorm(nrow(synthpop), meanlog = 0, sdlog = 0.1) 

probabilities <- 1 / (1 + exp(-z1))
#hist(probabilities)
synthpop$fitness_access<- cut(probabilities, breaks = quantile(probabilities, probs = seq(0, 1, length.out = 3)), labels = c("yes", "no"), include.lowest = TRUE)
#table(synthpop$fitness_access)


z2 <- 0 +
  0.3 * synthpop$mvpa_total_acc_scaled -       
  0.6 * synthpop$activity_pattern_2 -         
  0.4 * synthpop$activity_pattern_3 -          
  0.2 * synthpop$activity_pattern_4 +          
  0.1 * synthpop$`edu_highschool/ged` +                 
  0.2 * synthpop$`edu_some college` +                    
  0.3 * synthpop$`edu_college and above` +  
  0.3 * synthpop$self_reported_health_excellent +
  0.2 * synthpop$`self_reported_health_very good` + 
  0.2 * synthpop$self_reported_health_good +
  0.1 * synthpop$self_reported_health_fair +
  rlnorm(nrow(synthpop), meanlog = 0, sdlog = 0.1) 

probabilities <- 1 / (1 + exp(z2))
#hist(probabilities)
synthpop$health_literacy<- cut(probabilities, breaks = quantile(probabilities, probs = seq(0, 1, length.out = 4)), labels = c("basic", "intermediate", "advanced"), include.lowest = TRUE)
#table(synthpop$health_literacy)

#assocstats(table(synthpop$fitness_access, synthpop$health_literacy))$cramer

variables <- c("cluster", "seqn",
               "modpa_total", "vigpa_total", "mvpa_total_acc", "modpa_total_sqrt", "vigpa_total_sqrt", "mvpa_total_acc_sqrt",
               "activity_pattern", "source", "srvy_yr",
               "age", "race","gender", "marital", "edu","poverty","work","insurance",
               "self_reported_health","bmi", "smoker", "alcohol_cat",
               "hypertension", "diabetes", "heartdiseases", "cancers","stroke", 
               "fitness_access", "health_literacy")
synthpop<- synthpop %>%
  select(all_of(variables))

set.seed(2024)
n_samples <- 10000
n_datasets <- 100
#n_datasets <- 10 # CHANGE TEST

sampled_datasets <- vector("list", n_datasets)

for (i in 1:n_datasets) {
  sampled_datasets[[i]] <- synthpop[sample(nrow(synthpop), n_samples, replace = FALSE), ]
}
names(sampled_datasets) <- paste0("dataset_", 1:n_datasets)
length(sampled_datasets)

#save.image(file = file.path(work_dir, ".RData"))

#--------Evaluate R-squared values with varying predictors--------#
############################################
# model1 <- list(
#   lm(mvpa_total_acc_sqrt ~ modpa_total_sqrt + vigpa_total_sqrt + activity_pattern + source + srvy_yr +
#        age + race + gender + marital, data = synthpop),
#   multinom(activity_pattern ~ modpa_total_sqrt + vigpa_total_sqrt + mvpa_total_acc_sqrt + source + srvy_yr +
#              age + race + gender + marital, data = synthpop))
# print(paste("Adjusted R-squared for linear model1:",summary(model1[[1]])$adj.r.squared))
# print(paste("McFadden's pseudo-R-squared for multinomial model1:",  pR2(model1[[2]])["r2ML"]))
# 
# model2 <- list(
#   lm(mvpa_total_acc_sqrt ~ modpa_total_sqrt + vigpa_total_sqrt + activity_pattern + source + srvy_yr +
#        age + race + gender + marital + edu + poverty + work + insurance + fitness_access, data = synthpop),
#   multinom(activity_pattern ~ modpa_total_sqrt + vigpa_total_sqrt + mvpa_total_acc_sqrt + source + srvy_yr +
#              age + race + gender + marital + edu + poverty + work + insurance + fitness_access, data = synthpop)
# )
# 
# print(paste("Adjusted R-squared for linear model2:",summary(model2[[1]])$r.squared))
# print(paste("McFadden's pseudo-R-squared for multinomial model2:",  pR2(model2[[2]])["r2ML"]))
# 
# model3 <- list(
#   lm(mvpa_total_acc_sqrt ~ modpa_total_sqrt + vigpa_total_sqrt + activity_pattern + source + srvy_yr +
#        age + race + gender + marital + edu + poverty + work + insurance +
#        self_reported_health + bmi + smoker + alcohol_cat + hypertension + diabetes + heartdiseases + cancers + stroke +
#        fitness_access + health_literacy, data = synthpop),
#   multinom(activity_pattern ~ modpa_total_sqrt + vigpa_total_sqrt + mvpa_total_acc_sqrt + source + srvy_yr +
#              age + race + gender + marital + edu + poverty + work + insurance +
#              self_reported_health + bmi + smoker + alcohol_cat + hypertension + diabetes + heartdiseases + cancers + stroke +
#              fitness_access + health_literacy, data = synthpop))
# print(paste("Adjusted R-squared for linear model3:",summary(model3[[1]])$r.squared))
# print(paste("McFadden's pseudo-R-squared for multinomial model3:",  pR2(model3[[2]])["r2ML"]))

# library(lmtest)
# library(car)
# model<-model3[[1]]
# #heteroscedasticity check
# ncvTest(model)
#
# plot(model$fitted.values, model$residuals,
#      xlab = "Fitted Values",
#      ylab = "Residuals",
#      main = "Residuals vs Fitted Values")
# abline(h = 0, col = "red")
#
# # normality
# qqnorm(residuals(model))
# qqline(residuals(model), col = "red")
#
# vif(model)


#--------Introduce missing data--------#
############################################
exclude_vars <- c("modpa_total", "vigpa_total", "mvpa_total_acc", "cluster", "seqn")
target_vars <- c("mvpa_total_acc_sqrt", "mvpa_total_acc", "activity_pattern")
predictor_vars_mar <- c("srvy_yr", "age", "race", "gender", "marital",
                        "modpa_total", "vigpa_total", "edu", 
                        "poverty", "work", "insurance", "self_reported_health", "bmi", "smoker", "alcohol_cat", 
                        "hypertension", "diabetes", "heartdiseases", "cancers", "stroke", 
                        "fitness_access", "health_literacy")
predictor_vars_mnar <- c("srvy_yr","age", "race", "gender", "marital",
                         "mvpa_total_acc", "activity_pattern")

percents <- c(50, 70, 90)
mech <- c("mar", "mnar")

# Initialize an empty list to store simulated data
simdata_miss_50mar_list <- list()
simdata_miss_50mnar_list <- list()
simdata_miss_70mar_list <- list()
simdata_miss_70mnar_list <- list()
simdata_miss_90mar_list <- list()
simdata_miss_90mnar_list <-list()

for (percent in percents) {
  for (m in mech) {
    # Determine predictor variables based on the mechanism
    if (m == "mar") {
      predictor_vars <- predictor_vars_mar
    } else if (m == "mnar") {
      predictor_vars <- predictor_vars_mnar
    }

    # Loop through each sampled dataset
    for (i in seq_along(sampled_datasets)) { 
      sampled_dataset <- sampled_datasets[[i]]

      # Simulate data using the specified percent and predictor variables
      simdata <- introduce_missingness(sampled_dataset, target_vars, predictor_vars, percent/100)

      # Construct the list name based on the percentage and mechanism
      list_name <- paste0("simdata_miss_", percent, m, "_list")

      # Store the simulated data in the appropriate list
      assign(list_name, append(get(list_name), list(simdata)))
    }
  }
}



rm(synthpop, probabilities, z1, z2)
gc()
############################################
############################################
#--------SRMI imputation --------#
############################################
############################################
# srclib <<- "/Library/srclib/R_copy"  # CHANGE TEST
srclib <<- "/nfs/turbo/isr-bwest1/sw/rhel8/srclib/0.3.1/R" # initialize srclib
source(file.path(srclib, "init.R", fsep=.Platform$file.sep))
# List of source files
source_files <- list("1_srmi_level1_spec1.set", "2_srmi_level2_spec1.set", "3_srmi_level3_spec1.set",
                     "4_srmi_level1_spec2.set", "5_srmi_level2_spec2.set", "6_srmi_level3_spec2.set")
# Function to perform imputation within each dataset's subgroups

num_cores <- detectCores()  
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define the parameters
missing_patterns <- c("50mar", "50mnar", "70mar", "70mnar", "90mar", "90mnar")
data_lists_names <- c("simdata_miss_50mar_list", "simdata_miss_50mnar_list",
                      "simdata_miss_70mar_list", "simdata_miss_70mnar_list",
                      "simdata_miss_90mar_list", "simdata_miss_90mnar_list")

for (i in seq_along(missing_patterns)) {  
  data_list <- get(data_lists_names[i])  
  process_simulation_data_parallel(paste0("miss_", missing_patterns[i]), data_list)
  gc()  # Call garbage collector to free up memory
}

# Combine RDA files in parallel
for (pattern in missing_patterns) {
  combine_rda_files_parallel(file.path(work_dir, paste0("miss_", pattern)))
}

stopCluster(cl)
############################################
############################################



############################################
############################################
#--------PMM imputation--------#
############################################
############################################
# source_files <- list("1_pmm_level1_spec1.R", "2_pmm_level2_spec1.R", "3_pmm_level3_spec1.R",
#                      "4_pmm_level1_spec2.R", "5_pmm_level2_spec2.R", "6_pmm_level3_spec2.R")
# num_cores <- detectCores() 
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# 
# # Define the parameters
# missing_patterns <- c("50mar", "50mnar", "70mar", "70mnar", "90mar", "90mnar")
# 
# data_lists_names <- c("simdata_miss_50mar_list", "simdata_miss_50mnar_list",
#                    "simdata_miss_70mar_list", "simdata_miss_70mnar_list",
#                    "simdata_miss_90mar_list", "simdata_miss_90mnar_list")
# 
# # Loop through each dataset name, using get() to access the dataset directly
# for (i in seq_along(missing_patterns)) {  
#   data_list <- get(data_lists_names[i])  
#   process_simulation_data_parallel_r(paste0("miss_", missing_patterns[i]), data_list)
#   gc()  # Call garbage collector to free up memory
# }
# for (pattern in missing_patterns) {
#   combine_rda_files_parallel_r(file.path(work_dir, paste0("miss_", pattern)))
#   }
# 
# stopCluster(cl)

############################################
############################################



#-------- Combine Estimates --------#
############################################
#--------Bias, SE, RMSE, CR, FMI--------#
#Initialize result sheet
result_sheet_eachrun <- data.frame(
  scenario = character(),
  variable = character(),
  estimate_mi = numeric(),
  mi_vw = numeric(),
  mi_vb = numeric(),
  mi_vt = numeric(),
  bias = numeric(),
  relative_bias = numeric(),
  df = numeric(),
  cr = numeric(),
  lambda = numeric(),
  stringsAsFactors = FALSE
)

variable1 <- "mvpa_total_acc"
benchmark <- 221.7086
variable2 <- "activity_pattern"
benchmarks <- c(0.1124, 0.4469, 0.3254, 0.1153)

# Call the process_folders functions
process_folders_continuous(work_dir, variable1, benchmark)
process_folders_categorical(work_dir, variable2, benchmarks)

write_csv(result_sheet_eachrun, paste0(work_dir, '/result_sheet_eachrun.csv'))

result_sheet_eachrun<-read.csv(paste0(work_dir, '/result_sheet_eachrun.csv'))
result_sheet_eachrun <- result_sheet_eachrun %>%
  mutate(scenario_avg_run = str_replace(scenario, "_run_\\d+", ""))

# variability across samples
# Group by the scenario without the run number and calculate the average of the last five columns
result_sheet <- result_sheet_eachrun %>%
  group_by(scenario_avg_run, variable) %>%
  summarise(
    variance = mean(mi_vt)+ (1 + 1/n_datasets) * var(estimate_mi),
    pooled_estimate = mean(estimate_mi),
    mi_vw =  mean(mi_vw),
    mi_vb =  mean(mi_vb),
    mi_vt =  mean(mi_vt),
    bias = mean(bias),
    relative_bias = mean(relative_bias),
    df = mean(df),
    cr = sum(cr),
    lambda = mean(lambda, na.rm = TRUE)
  ) %>%
  ungroup()

result_sheet$se<-sqrt(result_sheet$variance)
result_sheet$rmse<-sqrt(result_sheet$variance + result_sheet$bias * result_sheet$bias)


result_sheet<- result_sheet %>%
  mutate(
    # Extract miss_mech and miss_rate
    miss_mech = case_when(
      grepl("mar", scenario_avg_run) ~ "MAR",
      grepl("mnar", scenario_avg_run) ~ "MNAR"
    ),
    miss_rate = gsub(".*_(\\d+)(mar|mnar)_levelspec.*", "\\1", scenario_avg_run),

    # Extract levelspec value
    levelspec = gsub(".*levelspec_(\\d+).*", "\\1", scenario_avg_run),

    # Create shared_var and model_spec based on levelspec
    shared_var = case_when(
      levelspec == "1" ~ "Least",
      levelspec == "2" ~ "Moderate",
      levelspec == "3" ~ "Most",
      levelspec == "4" ~ "Least",
      levelspec == "5" ~ "Moderate",
      levelspec == "6" ~ "Most"
    ),
    model_spec = case_when(
      levelspec %in% c("1", "2", "3") ~ "w/_Transf", 
      levelspec %in% c("4", "5", "6") ~ "w/o_Transf"
    ),
    levelspec = NULL,
  )

names(result_sheet)[names(result_sheet) == "scenario_avg_run"] <- "scenario"
result_sheet$scenario <- rep(1:(nrow(result_sheet)/5), each = 5)
result_sheet <- result_sheet[, c(colnames(result_sheet)[1], tail(colnames(result_sheet), 4), colnames(result_sheet)[-c(1, (ncol(result_sheet)-3):ncol(result_sheet))])]
write_csv(result_sheet, paste0(work_dir, '/result_sheet.csv'))



#--------Accuracy--------#
accuracy_df_eachrun_combined <- data.frame(scenario = character(), accuracy = numeric(), f1_score = numeric(), 
                                  stringsAsFactors = FALSE)

folders <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
subfolders <- unlist(lapply(folders, list.dirs, full.names = TRUE, recursive = FALSE))
# subfolders <-subfolders[c(1:600)] # to make sure overburn CPUs
# Set up parallel backend
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, c("work_dir", "subfolders", "sampled_datasets"))

accuracy_df_eachrun_combined <- foreach(subfolder = subfolders, .combine = rbind, .packages = c("caret", "stringr"), .export = c("compute_accuracy_for_subfolder")) %dopar% {
  compute_accuracy_for_subfolder(subfolder)
}
stopCluster(cl)

write_csv(accuracy_df_eachrun_combined, paste0(work_dir, '/accuracy_df_eachrun.csv'))
accuracy_df_eachrun<-read.csv(paste0(work_dir, '/accuracy_df_eachrun.csv'))
accuracy_df_eachrun <- accuracy_df_eachrun %>%
  mutate(scenario_avg_run = str_replace(scenario, "_run_\\d+", ""))
accuracy_df <- accuracy_df_eachrun %>%
  group_by(scenario_avg_run) %>%
  summarise(
    accuracy =  mean(accuracy),
    #f1_score =  mean(f1_score)
  ) %>%
  ungroup()
write_csv(accuracy_df, paste0(work_dir, '/accuracy_df.csv'))



#--------AUC R-squared [LogModel]--------#
# observed data
# Initialize vectors to store results
auc_values <- numeric(length(sampled_datasets))
R_squared_values <- numeric(length(sampled_datasets))

# Loop through each dataset
for (i in seq_along(sampled_datasets)) {
  combined_subset <- sampled_datasets[[i]]  # Select dataset i
  combined_subset$mvpa_total <- combined_subset$modpa_total + combined_subset$vigpa_total
  # Fit the logistic regression model with probit link
  model <- glm(hypertension ~ gender + age + race + edu + mvpa_total, 
               data = combined_subset, 
               family = binomial(link = 'probit'))
  
  # Calculate AUC and R-squared
  auc_values[i] <- roc(combined_subset$hypertension, predict(model, combined_subset, type = "response"))$auc
  R_squared_values[i] <- as.numeric(pR2(model)[4])  # [4] selects McFadden R-squared
}

# Combine the results into a data frame
results <- data.frame(Dataset = seq_along(sampled_datasets), AUC = auc_values, R_squared = R_squared_values)

mean(results$AUC) # 0.7815
mean(results$R_squared) #0.1756

# Initialize the combined results dataframe
log_eachrun_combined <- data.frame(scenario = character(),
                                           auc_values = numeric(),
                                           R_squared_values = numeric(),
                                           stringsAsFactors = FALSE)

# Get the list of folders and subfolders
folders <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
subfolders <- unlist(lapply(folders, list.dirs, full.names = TRUE, recursive = FALSE))
subfolders <-subfolders[c(2001:3000)]
# Set up parallel backend
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, c("work_dir"))

for (i in seq(1, length(subfolders), by = 500)) {
  subfolder_batch <- subfolders[i:min(i + 499, length(subfolders))]
  log_eachrun <- foreach(subfolder = subfolder_batch, .combine = rbind, 
                         .packages = c("stringr", "pscl", "pROC"), 
                         .export = c("logreg_for_subfolder")) %dopar% {
                           # Use tryCatch to handle errors
                           tryCatch({
                             # Call the logistic regression function
                             logreg_for_subfolder(subfolder)
                           }, error = function(e) {
                             # Print a message indicating the failure
                             message(paste("Task failed for subfolder:", subfolder, "->", e$message))
                             # Return a data frame with NA values for failed tasks
                             return(data.frame(scenario = subfolder, mean_auc = NA, mean_R_squared = NA))
                           })
                         }
  
  # Combine the results into the main dataframe
  log_eachrun_combined <- rbind(log_eachrun_combined, log_eachrun)
}
stopCluster(cl)

# Write the combined results to a CSV file
write_csv(log_eachrun_combined, paste0(work_dir, '/log_eachrun.csv'))


log_eachrun <- read.csv(paste0(work_dir, '/log_eachrun.csv'))

# Modify the scenario to get the average run scenario
log_eachrun <- log_eachrun %>%
  mutate(scenario_avg_run = str_replace(scenario, "_run_\\d+", ""))

# Summarise the data
log_result <- log_eachrun %>%
  group_by(scenario_avg_run) %>%
  summarise(
    mean_auc = mean(mean_auc, na.rm = TRUE), # Use correct variable name
    mean_R_squared = mean(mean_R_squared, na.rm = TRUE), # Use correct variable name
  ) %>%
  ungroup()

# Write the summarized data to a CSV file
write_csv(log_result, paste0(work_dir, '/log_result.csv'))
