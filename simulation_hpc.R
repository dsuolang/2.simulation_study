# Via simulation, this study will examine the impacts of 5 key factors on imputation performance: 
#  Missing data mechanism
#  Missing rate
#  Amount of shared auxiliary variables
#  Model misspecification. 
#  Modeling approach 
#  Total scenarios = 36(2 * 3 * 3 * 2) * + 18(2 * 3 * 3 * 1) = 54


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
library(mice)
library(pscl)
library(doParallel)
library(readr)
library(stringr)
library(fastDummies)
library(vcd)
library(magrittr)
library(caret) 
library(pscl)
library(clustMixType)
library(pROC)
options(scipen = 999)

source("functions.R")


#--------Draw samples from a synthetic population--------#
############################################
synthpop <- read_dta("synthpop_sim1.dta") 
synthpop[c("mvpa_total_acc_sqrt", "modpa_total_sqrt", "vigpa_total_sqrt")] <-
  sqrt(synthpop[c("mvpa_total_acc", "modpa_total", "vigpa_total")])
# We only need sqrt transformed mvpa_total_acc, but we prepare all variables for potential use.

synthpop <- synthpop %>%
  mutate(self_reported_health = as.character(self_reported_health)) %>%
  mutate(self_reported_health = recode(self_reported_health,
                                       `1` = "excellent",
                                       `2` = "very good",
                                       `3` = "good",
                                       `4` = "fair",
                                       `5` = "poor"))

synthpop <- synthpop %>% mutate(
  cluster = as.factor(cluster),activity_pattern = as.factor(activity_pattern)
)
synthpop$srvy_yr <- as.numeric(synthpop$srvy_yr) 

summary(synthpop$mvpa_total_acc) # Variable of interest 1: MVPA duration
table(synthpop$activity_pattern)/nrow(synthpop)# Variable of interest 2: activity pattern

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
synthpop$fitness_access<- as.character(cut(probabilities, breaks = quantile(probabilities, probs = seq(0, 1, length.out = 3)), labels = c("yes", "no"), include.lowest = TRUE))
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
synthpop$health_literacy<- as.character(cut(probabilities, breaks = quantile(probabilities, probs = seq(0, 1, length.out = 4)), labels = c("basic", "intermediate", "advanced"), include.lowest = TRUE))
#table(synthpop$health_literacy)

#assocstats(table(synthpop$fitness_access, synthpop$health_literacy))$cramer

variables <- c("cluster","seqn", "modpa_total_sqrt", "vigpa_total_sqrt", 
               "modpa_total", "vigpa_total", "mvpa_total_acc", "mvpa_total_acc_sqrt",
               "activity_pattern", "source", "srvy_yr",
               "age", "race","gender", "marital", "edu","poverty","work","insurance",
               "self_reported_health","bmi", "smoker", "alcohol_cat",
               "hypertension", "diabetes", "heartdiseases", "cancers","stroke", 
               "fitness_access", "health_literacy")
synthpop<- synthpop %>%
  select(all_of(variables))

rm(z1, z2)
save(synthpop, file="synthpop.rda")

load("synthpop.rda")
summary(synthpop) # Ensure the variable type is correct
set.seed(2024)
n_samples <- 10000
n_datasets <- 250

sampled_datasets <- vector("list", n_datasets)

for (i in 1:n_datasets) {
  sampled_datasets[[i]] <- synthpop[sample(nrow(synthpop), n_samples, replace = FALSE), ]
}
names(sampled_datasets) <- paste0("dataset_", 1:n_datasets)
length(sampled_datasets)

# check normality of complete data distribution
# means<-sapply(sampled_datasets, function(dataset) {
#   sum(dataset$activity_pattern == 1, na.rm=T)/nrow(dataset)
# })
# 
# hist(means)
# qqnorm(means)
# qqline(means, col='red')
# shapiro.test(means)




#--------Evaluate R-squared values with varying predictors--------#
############################################
model_results_list <- lapply(sampled_datasets, process_model)
rsquared_results <- do.call(rbind, model_results_list)

#rsquared_results$sample <- sub(".*_(\\d+).*", "\\1", rsquared_results$X)
rsquared_results$sample <- sub(".*_(\\d+).*", "\\1", rsquared_results$dataset_id)
rsquared_results <- rsquared_results %>%
  select(sample, model, variable, rsquared)
rsquared_mean_results <- rsquared_results %>%
  group_by(model, variable) %>%
  summarize(mean_rsquared = mean(rsquared))
rsquared_mean_results
write.csv(rsquared_mean_results, "rsquared_mean_results.csv")


#--------Introduce missing data--------#
############################################
target_vars <- c("mvpa_total_acc_sqrt", "mvpa_total_acc", "activity_pattern")
predictor_vars_mar <- c("srvy_yr", 
                        "age", "gender", "race", "marital", 
                        "edu", "poverty", "work", "insurance",
                        "self_reported_health","bmi", "smoker", "alcohol_cat", 
                        "hypertension", "diabetes", "heartdiseases", "cancers", "stroke",
                        "fitness_access", "health_literacy", "modpa_total", "vigpa_total")
                        
predictor_vars_mnar <- c("srvy_yr", "age", "gender", "race", "marital", 
                         "mvpa_total_acc", "activity_pattern",  "mvpa_total_acc:activity_pattern")

percents <- c(50, 70, 90)
mech <- c("mar", "mnar")

# Initialize an empty list to store simulated data
simdata_miss_50mar_list <- list()
simdata_miss_50mnar_list <- list()
simdata_miss_70mar_list <- list()
simdata_miss_70mnar_list <- list()
simdata_miss_90mar_list <- list()
simdata_miss_90mnar_list <-list()

set.seed(2024)
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

gc()

############################################
############################################
#--------SRMI imputation --------#
############################################
work_dir <- paste0(getwd(), "/srmi")
setwd(work_dir)
############################################ Call IVEware from R, requires IVEware installation.
srclib <<- "/nfs/turbo/isr-bwest1/sw/rhel8/srclib/0.3.1/R" # initialize srclib 
source(file.path(srclib, "init.R", fsep=.Platform$file.sep))
# List of source files
source_files <- list("1_srmi_level1_spec1.set", "2_srmi_level2_spec1.set", "3_srmi_level3_spec1.set",
                     "4_srmi_level1_spec2.set", "5_srmi_level2_spec2.set", "6_srmi_level3_spec2.set")

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define the parameters
missing_patterns <- c("50mar", "50mnar", "70mar", "70mnar", "90mar", "90mnar")

data_lists_names <- c("simdata_miss_50mar_list", "simdata_miss_50mnar_list",
                      "simdata_miss_70mar_list", "simdata_miss_70mnar_list",
                      "simdata_miss_90mar_list", "simdata_miss_90mnar_list")

# Locate the function "process_simulation_data_parallel" and modify it to point to the specific directory where IVEware is installed.
for (i in seq_along(missing_patterns)) {   
  data_list <- get(data_lists_names[i])
  process_simulation_data_parallel(paste0("miss_", missing_patterns[i]), data_list)
  gc()  
}

# Combine RDA files in parallel
for (pattern in missing_patterns) {
  combine_rda_files_parallel(file.path(work_dir, paste0("miss_", pattern)))
}
print("SRMI Imputation Finished")
stopCluster(cl)
############################################
############################################



############################################
#############################################
############################################
#--------PMM imputation--------#
work_dir <- paste0(getwd(), "/pmm")
setwd(work_dir)
############################################
############################################
source_files <- list("1_pmm_level1_spec1.R", "2_pmm_level2_spec1.R", "3_pmm_level3_spec1.R",
                     "4_pmm_level1_spec2.R", "5_pmm_level2_spec2.R", "6_pmm_level3_spec2.R")
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define the parameters
missing_patterns <- c("50mar", "50mnar", "70mar", "70mnar", "90mar", "90mnar")

data_lists_names <- c("simdata_miss_50mar_list", "simdata_miss_50mnar_list", "simdata_miss_70mar_list", 
                   "simdata_miss_70mnar_list","simdata_miss_90mar_list", "simdata_miss_90mnar_list")

# Loop through each dataset name, using get() to access the dataset directly
for (i in seq_along(missing_patterns)) {
  data_list <- get(data_lists_names[i])
  process_simulation_data_parallel_r(paste0("miss_", missing_patterns[i]), data_list)
  gc()  
}
for (pattern in missing_patterns) {
  print(pattern)
  combine_rda_files_parallel_r(file.path(work_dir, paste0("miss_", pattern)))
  }
print("PMM Imputation Finished")
stopCluster(cl)
############################################
############################################

gc()



############################################
# Perform the below analysis separately for files under for SRMI and PMM directory.
############################################
#-------- Combine Estimates --------# 
############################################
# Best possible coverage rate
sampled_datasets_dummy<-lapply(sampled_datasets, function(data){
  dummy_cols(data, select_columns = 'activity_pattern', remove_first_dummy = FALSE)
})
results_var1 <- calculate_ci_with_custom_variance(data_list = sampled_datasets, variable = "mvpa_total_acc")
results_var2 <- lapply(paste0("activity_pattern_", 1:4), function(var){
  calculate_ci_with_custom_variance(data_list = sampled_datasets_dummy, variable = var)
})
results_var2 <- do.call(rbind, results_var2)

# Combine all results into a single data frame
ci_complete_data<- rbind(results_var1, results_var2)


# Use ifelse to assign the benchmark based on the variable name
ci_complete_data$benchmark <- ifelse(ci_complete_data$variable == "mvpa_total_acc", 
                                     209.2866, 
                                     ifelse(ci_complete_data$variable == "activity_pattern_1", 
                                            0.1014771,
                                            ifelse(ci_complete_data$variable == "activity_pattern_2", 
                                                   0.4302855,
                                                   ifelse(ci_complete_data$variable == "activity_pattern_3", 
                                                          0.3392789,
                                                          ifelse(ci_complete_data$variable == "activity_pattern_4", 
                                                                 0.1289586, NA)))))

ci_complete_data$coverage <- ifelse(
  ci_complete_data$benchmark >= ci_complete_data$lower & ci_complete_data$benchmark <= ci_complete_data$upper, 1, 0)
summary_counts <- ci_complete_data %>%
  group_by(variable) %>%
  summarise(total_within_range = sum(coverage == 1)/n_datasets)



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
  t_value = numeric(),
  lower = numeric(),
  upper = numeric(),
  coverage = numeric(),
  lambda = numeric(),
  stringsAsFactors = FALSE
)

variable1 <- "mvpa_total_acc"
benchmark <- 209.2866 #based on the initial population data
variable2 <- "activity_pattern"
benchmarks <- c(0.1014771, 0.4302855, 0.3392789, 0.1289586) #based on the initial population data

# Call the process_folders functions
process_folders_continuous(work_dir, variable1, benchmark)
process_folders_categorical(work_dir, variable2, benchmarks)
write_csv(result_sheet_eachrun, paste0(work_dir, '/result_sheet_eachrun.csv'))

result_sheet_eachrun<-read.csv(paste0(work_dir, '/result_sheet_eachrun.csv'))
result_sheet_eachrun <- result_sheet_eachrun %>%
  mutate(scenario_avg_run = str_replace(scenario, "_run_\\d+", ""))


# variability across samples
# Group by the scenario without the run number and calculate the average of the last five columns
#n_datasets <- 250
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
    t_value = mean(t_value),
    lower = mean(lower),
    upper = mean(upper),
    coverage = sum(coverage)/n_datasets,
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



#--------Accuracy (for activity patterns)--------#
accuracy_df_eachrun_combined <- data.frame(
  scenario = character(),
  accuracy = numeric(),
  stringsAsFactors = FALSE
)

folders <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
subfolders <- unlist(lapply(folders, list.dirs, full.names = TRUE, recursive = FALSE))

batch_size <- 600  # Process 600 subfolders per batch to make sure memory does not exhaust
total_batches <- ceiling(length(subfolders) / batch_size)

for (batch_idx in seq_len(total_batches)) {
  gc()
  # Get the current batch of subfolders
  batch_start <- (batch_idx - 1) * batch_size + 1
  batch_end <- min(batch_idx * batch_size, length(subfolders))
  batch_subfolders <- subfolders[batch_start:batch_end]
  
  # Print progress
  cat("Processing batch", batch_idx, "of", total_batches, "(", batch_start, "-", batch_end, ")\n")
  
  # Set up parallel backend
  num_cores <- detectCores() -1  # Reserve one core for the system
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Export necessary variables and functions
  clusterExport(cl, c("work_dir", "sampled_datasets", "compute_accuracy_for_subfolder"))
  
  # Process the current batch in parallel
  batch_results <- foreach(
    subfolder = batch_subfolders,
    .combine = rbind,
    .packages = c("caret", "stringr", "haven")
  ) %dopar% {
    tryCatch(
      compute_accuracy_for_subfolder(subfolder),
      error = function(e) {
        return(data.frame(scenario = subfolder, accuracy = e$message))
      }
    )
  }
  
  # Stop the cluster after processing the batch
  stopCluster(cl)
  
  # Append the batch results to the combined results
  accuracy_df_eachrun_combined <- rbind(accuracy_df_eachrun_combined, batch_results)
}


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




#-------Logistic regression--------
# Set the outcome variable (change this as needed)
# Self-report
results<-logreg_sampledata(sampled_datasets, outcome_var)
summary(results)
# mean_pr_auc     mean_pseudo_R2   mean_deviance      mean_AIC        mean_BIC 
# Mean   :0.6302   Mean   :0.2118   Mean   : 9770   Mean   : 9784   Mean   : 9835   Hypertension
# Mean   :0.3091   Mean   :0.1707   Mean   :5691   Mean   :5705   Mean   :5756   Diabetes
# Mean   :0.2923   Mean   :0.1351   Mean   :5974   Mean   :5988   Mean   :6038  Heart diesease

# Define all outcome variables
outcome_vars <- c("diabetes", "hypertension", "heartdiseases")

folders <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
subfolders <- unlist(lapply(folders, list.dirs, full.names = TRUE, recursive = FALSE))
length(subfolders)

batch_size <- 9000  # adjust based on resources
total_batches <- 1  

# Loop over outcome variables
for (outcome_var in outcome_vars) {
  cat("\n=== Processing outcome:", outcome_var, "===\n")
  
  # Storage for combined results
  log_eachrun_combined <- data.frame(
    scenario = character(),
    mean_pr_auc = numeric(),
    mean_pseudo_R2 = numeric(),
    mean_deviance = numeric(),
    mean_AIC = numeric(),
    mean_BIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process in batches
  for (batch_idx in seq_len(total_batches)) {
    batch_start <- (batch_idx - 1) * batch_size + 1
    batch_end <- min(batch_idx * batch_size, length(subfolders))
    batch_subfolders <- subfolders[batch_start:batch_end]
    
    cat("Processing batch", batch_idx, "of", total_batches, "(", batch_start, "-", batch_end, ")\n")
    
    # Parallel setup
    num_cores <- detectCores() 
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    clusterExport(cl, c("work_dir", "logreg_for_subfolder", "outcome_var"))
    
    log_eachrun <- foreach(
      subfolder = batch_subfolders,
      .combine = rbind,
      .packages = c("stringr", "pscl", "PRROC")
    ) %dopar% {
      tryCatch(
        logreg_for_subfolder(subfolder, outcome_var = outcome_var),
        error = function(e) {
          data.frame(
            scenario = subfolder,
            mean_pr_auc = NA,
            mean_pseudo_R2 = NA,
            mean_deviance = NA,
            mean_AIC = NA,
            mean_BIC = NA
          )
        }
      )
    }
    
    stopCluster(cl)
    log_eachrun_combined <- rbind(log_eachrun_combined, log_eachrun)
  }
  
  # Save raw batch results
  write_csv(log_eachrun_combined, file = paste0(work_dir, "/log_eachrun_", outcome_var, ".csv"))
  
  # Clean & aggregate results
  log_eachrun <- read.csv(paste0(work_dir, "/log_eachrun_", outcome_var, ".csv"))
  log_eachrun <- log_eachrun %>%
    mutate(across(c(mean_pr_auc, mean_pseudo_R2, mean_deviance, mean_AIC, mean_BIC), as.numeric)) %>%
    drop_na() %>%
    mutate(scenario_avg_run = str_replace(scenario, "_run_\\d+", ""))
  
  log_result <- log_eachrun %>%
    group_by(scenario_avg_run) %>%
    summarise(
      mean_pr_auc = mean(mean_pr_auc, na.rm = TRUE),
      mean_pseudo_R2 = mean(mean_pseudo_R2, na.rm = TRUE),
      mean_deviance = mean(mean_deviance, na.rm = TRUE),
      mean_AIC = mean(mean_AIC, na.rm = TRUE),
      mean_BIC = mean(mean_BIC, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Save aggregated results
  write_csv(log_result, paste0(work_dir, "/log_result_", outcome_var, ".csv"))
}

cat("Finished modeling all outcome variables.\n")