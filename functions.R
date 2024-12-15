#-------Define functions--------#
calculate_squared_se <- function(variable) {
  squared_se <- (sd(variable) / sqrt(length(variable)))^2
  return(squared_se)
}

calculate_squared_se_prop <- function(variable) {
  p<-mean(variable)
  squared_se <- p * (1 - p)/ length(variable)
  return(squared_se)
}

# Define the get_r_squared function
create_binary_variables <- function(df) {
  binary_df <- data.frame(matrix(NA, nrow = nrow(df), ncol = 0))
  for (col in names(df)) {
    df[[col]] <- factor(df[[col]])
    levels(df[[col]]) <- gsub("[ /</>/-]", "_", levels(df[[col]]))
    binary_vars <-
      model.matrix( ~ . - 1, data = df[, col, drop = FALSE])
    colnames(binary_vars) <- paste0(col, "_", levels(df[[col]]))
    binary_df <- cbind(binary_df, binary_vars)
  }
  df <- cbind(df, binary_df)
  return(df)
}

# Function to generate missing data
sim_miss <- function(data, vars, miss_perc) {
  set.seed(2024)
  data <- as.data.frame(data)
  
  logistic_model <- glm(source ~ ., 
                     data = data[, c(vars)], 
                     family = binomial)
  predicted <- 1 - predict(logistic_model, type = "response") # probability of being "not missing" category (NHANES)
    
  data[[paste0("miss_indicator")]] <- as.numeric(ifelse(predicted > quantile(predicted, miss_perc/100), TRUE, FALSE))
  return(data)
}

introduce_missingness <- function(df, cols_to_make_missing, predictors, missing_rate) {
  # Include the source variable as the outcome.
  formula <- as.formula(paste("source", "~", paste(predictors, collapse = " + ")))
  
  # Generate a logistic regression model based on the predictors and source variable.
  logit_model <- glm(formula, data = df, family = binomial(link = "logit"))
  
  # Calculate the probability of missingness for each row. (prob of being in NHANES (not missing))
  df$prob_missing <- 1 - predict(logit_model, type = "response")
  
  # Adjust the predicted probabilities to meet the desired overall missing rate.
  df$adjusted_prob_missing <- df$prob_missing * (missing_rate / mean(df$prob_missing))
  df$adjusted_prob_missing <- pmin(pmax(df$adjusted_prob_missing, 0), 1)
  
  # Ensure reproducibility
  set.seed(2024)
  # Generate a uniform random variable for each row.
  u_random <- runif(nrow(df))
  
  # Create a missing indicator where True means the value should be set to missing (NA).
  df$miss_indicator <- u_random < df$adjusted_prob_missing # T = missing, F = not missing
  
  # Apply the missing mask to the specified columns.
  for (col in cols_to_make_missing) {
    df[[col]] <- ifelse(df$miss_indicator, NA, df[[col]])
  }
  
  # Clean up the temporary columns except miss_indicator and return the modified DataFrame.
  df <- df %>% select(-prob_missing, -adjusted_prob_missing)
  return(df)
}

impute_within_subgroup <- function(dataset, source_file) {
  # Create directory for the dataset
  dataset_dir <- getwd()
  #dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  # Divide dataset into subsets based on "cluster"
  clusters <- split(dataset, dataset$cluster)
  # Perform imputation on each subset
  for (key in names(clusters)) {
    cluster_data <- clusters[[key]]
    print(nrow(cluster_data))
    cluster_dir <- file.path(dataset_dir, key)
    dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Copy the source file to the cluster directory
    file.copy(paste0(work_dir, "/",basename(source_file)),cluster_dir)
    
    # Set working directory to the cluster directory
    setwd(cluster_dir)
    # Save the cluster data
    write_csv(cluster_data, "observed.csv") 
    observed <- read_csv("observed.csv")
    save(observed, file="observed.rda")
    
    # Perform the imputation
    impute(name = sub("\\.set$", "", basename(source_file)))
  }
}


process_simulation_data <- function(miss_type, simdata_list) {
  miss_dir <- paste0(work_dir, "/", miss_type)
  dir.create(miss_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(miss_dir)
  for (i in seq_along(simdata_list)) {
    for (j in seq_along(source_files)+3) {
      impute_within_subgroup(simdata_list[[i]], source_files[[j]])
    }
  }
}


process_simulation_data_parallel <- function(miss_type, simdata_list) {
  clusterExport(cl, varlist = c("impute_within_subgroup", "work_dir", "source_files"))
  
  clusterEvalQ(cl, {
    #srclib <<- "/Library/srclib/R"  # CHANGE TEST
    srclib <<- "/nfs/turbo/isr-bwest1/sw/rhel8/srclib/0.3.1/R" # initialize srclib
    source(file.path(srclib, "init.R", fsep = .Platform$file.sep))
  })
  
  # Define working directory
  miss_dir <- paste0(work_dir, "/", miss_type)
  dir.create(miss_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Parallelize the loops
  tryCatch({
    foreach(i = seq_along(simdata_list), .packages = c('foreach', 'readr')) %:%
      foreach(j = seq_along(source_files), .packages = c('foreach','readr')) %dopar% {
        # Create a local directory for each task
        task_dir <- paste0(miss_dir, "/run_", i, "_levelspec_", j)
        dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Change working directory temporarily within each worker
        setwd(task_dir)
        
        # Perform the imputation
        impute_within_subgroup(simdata_list[[i]], source_files[[j]])
        
        # Restore the original working directory
        setwd(miss_dir)
      }
  }, error = function(e) {
    message("Error: ", e$message)
  })
}

combine_rda_files <- function(work_dir) {
  # List all main folders in the directory
  main_folders <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
  
  for (main_folder in main_folders) {
    # Initialize an empty list to store data frames
    all_data <- list()
    
    # Iterate through each subfolder (1 to 5)
    for (i in 1:5) {
      subfolder_path <- file.path(main_folder, as.character(i))
      
      # List .rda files in the subfolder that start with 'srmi'
      rda_files <- list.files(subfolder_path, pattern = "^srmi.*\\.rda$", full.names = TRUE)
      
      # Iterate through each .rda file
      for (rda_file in rda_files) {
        load(rda_file)
        
        # Check if 'imp' exists in the environment after loading the .rda file
        if (exists("imp")) {
          # Append the 'imp' data frame to the list
          all_data[[length(all_data) + 1]] <- imp
        }
      }
    }
    
    # Combine all data frames if there are any
    if (length(all_data) > 0) {
      combined_data <- do.call(rbind, all_data)
      
      # Write the combined data to an .rda file in the main folder
      save(combined_data, file = file.path(main_folder, "combined_data.rda"))
    }
  }
}


combine_rda_files_parallel <- function(work_dir) {
  # List all main folders in the directory
  main_folders <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
  
  # Process each main folder in parallel
  results <- foreach(main_folder = main_folders, .packages = c("tools")) %dopar% {
    # Initialize an empty list to store data frames
    all_data <- list()
    
    # Iterate through each subfolder (1 to 5)
    for (i in 1:5) {
      subfolder_path <- file.path(main_folder, as.character(i))
      
      # List .rda files in the subfolder that start with 'srmi'
      rda_files <- list.files(subfolder_path, pattern = "^srmi.*\\.rda$", full.names = TRUE)
      
      # Iterate through each .rda file
      for (rda_file in rda_files) {
        # Load the .rda file in a local environment
        local_env <- new.env()
        load(rda_file, envir = local_env)
        
        # Check if 'imp' exists in the local environment
        if (exists("imp", envir = local_env)) {
          # Append the 'imp' data frame to the list
          all_data[[length(all_data) + 1]] <- local_env$imp
        }
      }
    }
    
    # Combine all data frames if there are any
    if (length(all_data) > 0) {
      combined_data <- do.call(rbind, all_data)
      
      # Write the combined data to an .rda file in the main folder
      save(combined_data, file = file.path(main_folder, "combined_data.rda"))
      
      return(TRUE)  # Return success indicator for this iteration
    } else {
      return(NULL)  # Return NULL if no data was combined
    }
  }
  
  return(results)  # Return the result list (TRUE for successful folders, NULL for empty ones)
}


rubin_combining_rule_continuous <- function(variable, benchmark, scenario) {
  # Get the number of imputations
  load('combined_data.rda')
  if (grepl("levelspec_(1|2|3)$", scenario)) {
    variable_sqrt <- paste0(variable, "_sqrt")
    combined_data[[variable]] <- combined_data[[variable_sqrt]]^2
  }
  m <- length(unique(combined_data$`MULT_`))

  estimates_mlevel<-combined_data %>% #each imputed data
    group_by(MULT_) %>%
    summarize(value = mean(get(variable))) 
  estimate_mi<- mean(estimates_mlevel$value) # pooled estimate from MI

  # Calculate the within-imputation variance (mean variance within each group)
  within_variance<-combined_data %>% #each imputed data
    group_by(MULT_) %>%
    summarize(value = calculate_squared_se(get(variable))) 
  # Within-imputation variance (W)
  mi_vw <- mean(within_variance$value) # pooled estimate from MI
  # Between-imputation variance (B)
  mi_vb <- var(estimates_mlevel$value)
  # Total variance (T)
  mi_vt <- mi_vw + (1 + 1/m) * mi_vb
  
  # Fraction of Missing Information (lambda)
  lambda <- (1 + 1/m) * mi_vb / mi_vt
  
  # Calculate bias
  bias <- estimate_mi - benchmark
  relative_bias <- bias/benchmark
  
  df <- (m-1)/lambda^2
  t_value <- qt(0.975, df)
  cr <- as.numeric((benchmark >= estimate_mi- t_value * sqrt(mi_vt)) & (benchmark <= estimate_mi + t_value * sqrt(mi_vt)))
  
  # Append results to result_sheet
  result_sheet_eachrun <<- rbind(result_sheet_eachrun, data.frame(
    scenario = scenario,
    variable = variable,
    estimate_mi = estimate_mi,
    mi_vw = mi_vw,
    mi_vb = mi_vb,
    mi_vt = mi_vt,
    bias = bias,
    relative_bias = relative_bias,
    df = df,
    cr = cr,
    lambda = lambda
  ))
}

rubin_combining_rule_categorical <- function(variable, benchmarks, scenario) {
  # Load the combined data
  load('combined_data.rda')
  
  # Create dummy variables
  df_dummies <- dummy_cols(combined_data, select_columns = variable)
  
  # Define activity pattern columns
  activity_cols <- c("activity_pattern_1", "activity_pattern_2", "activity_pattern_3", "activity_pattern_4")
  
  # Ensure all activity columns are present
  for (col in activity_cols) {
    if (!col %in% colnames(df_dummies)) {
      df_dummies[[col]] <- 0
    }
  }
  m <- length(unique(df_dummies$`MULT_`))
  
  # Loop over the first four columns of interest
  for (i in 1:4) {
    # Calculate the estimate (mean) for each imputation
    estimates_mlevel<-df_dummies  %>% #each imputed data
      group_by(MULT_) %>%
      summarize(value = mean(.data[[activity_cols[[i]]]]))
    estimate_mi<- mean(estimates_mlevel$value) # pooled estimate from MI
    
    # Calculate the within-imputation variance (mean variance within each group)
    within_variance<-df_dummies %>% #each imputed data
      group_by(MULT_) %>%
      summarize(value = calculate_squared_se_prop(.data[[activity_cols[[i]]]]))
    # Within-imputation variance (W)
    mi_vw <- mean(within_variance$value) # pooled estimate from MI
    # Between-imputation variance (B)
    mi_vb <- var(estimates_mlevel$value)
    # Total variance (T)
    mi_vt <- mi_vw + (1 + 1/m) * mi_vb
    
    # Fraction of Missing Information (lambda)
    lambda <- (1 + 1/m) * mi_vb / mi_vt
    
    # Calculate bias
    bias <- estimate_mi - benchmarks[[i]]
    relative_bias <- bias/benchmarks[[i]]
    
    df <- (m-1)/lambda^2
    t_value <- qt(0.95, df)
    cr <- as.numeric((benchmarks[[i]] >= estimate_mi- t_value * sqrt(mi_vt)) & (benchmarks[[i]] <= estimate_mi + t_value * sqrt(mi_vt)))
    
    # Append results to result_sheet
    result_sheet_eachrun <<- rbind(result_sheet_eachrun, data.frame(
      scenario = scenario,
      variable = activity_cols[[i]],
      estimate_mi = estimate_mi,
      mi_vw = mi_vw,
      mi_vb = mi_vb,
      mi_vt = mi_vt,
      bias = bias,
      relative_bias = relative_bias,
      df = df,
      cr = cr,
      lambda = lambda
    ))
  }
}

# Main function to loop through directories and apply the Rubin combining rule
process_folders_continuous <- function(work_dir, variable, benchmark) {
  # List all folders in the main directory
  folders <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
  
  for (folder in folders) {
    # List all subfolders within each folder
    subfolders <- list.dirs(folder, full.names = TRUE, recursive = FALSE)
    
    for (subfolder in subfolders) {
      # Set the subfolder as the working directory
      setwd(subfolder)
      
      # Extract scenario name from folder path
      scenario <- str_replace_all(subfolder, paste0(work_dir, "/"), "")
      scenario <- gsub("/", "_", scenario)
      # Call Rubin combining rule function
      rubin_combining_rule_continuous(variable, benchmark, scenario)
    }
  }
  
  # Return the result sheet at the end of processing
  return('Finished')
}

process_folders_categorical <- function(work_dir, variable, benchmarks) {
  # List all folders in the main directory
  folders <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
  
  for (folder in folders) {
    # List all subfolders within each folder
    subfolders <- list.dirs(folder, full.names = TRUE, recursive = FALSE)
    
    for (subfolder in subfolders) {
      # Set the subfolder as the working directory
      setwd(subfolder)
      
      # Extract scenario name from folder path
      scenario <- str_replace_all(subfolder, paste0(work_dir, "/"), "")
      scenario <- gsub("/", "_", scenario)
      # Call Rubin combining rule function
      rubin_combining_rule_categorical(variable, benchmarks, scenario)
    }
  }
  
  # Return the result sheet at the end of processing
  return('Finished')
}








impute_global_r <- function(dataset, source_file) {
    # Copy the source file to the directory
  file.copy(paste0(work_dir, "/", basename(source_file)), getwd(), overwrite = TRUE)
    # prepare the observed data
  write_csv(dataset, "observed.csv") 
  observed <- read_csv("observed.csv")
  save(observed, file="observed.rda")
    # Perform the imputation
  source(source_file)
}

impute_within_subgroup_r <- function(dataset, source_file) {
  # Create directory for the dataset
  dataset_dir <- getwd()
  #dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  # Divide dataset into subsets based on "cluster"
  clusters <- split(dataset, dataset$cluster)
  # Perform imputation on each subset
  for (key in names(clusters)) {
    cluster_data <- clusters[[key]]
    print(nrow(cluster_data))
    cluster_dir <- file.path(dataset_dir, key)
    dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Copy the source file to the cluster directory
    file.copy(paste0(work_dir, "/",basename(source_file)), cluster_dir)
    
    # Set working directory to the cluster directory
    setwd(cluster_dir)
    # Save the cluster data
    write_csv(cluster_data, "observed.csv") 
    observed <- read_csv("observed.csv")
    save(observed, file="observed.rda")
    
    # Perform the imputation
    source(source_file)
  }
}

combine_rda_files_parallel_r <- function(work_dir) {
  # List all main folders in the directory
  main_folders <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
  
  # Process each main folder in parallel
  results <- foreach(main_folder = main_folders, .packages = c("tools")) %dopar% {
    # Initialize an empty list to store data frames
    all_data <- list()
    
    # Iterate through each subfolder (1 to 5)
    for (i in 1:5) {
      subfolder_path <- file.path(main_folder, as.character(i))
      
      # List .rda files in the subfolder that start with 'pmm'
      rda_files <- list.files(subfolder_path, pattern = "^pmm.*\\.rda$", full.names = TRUE)
      
      # Iterate through each .rda file
      for (rda_file in rda_files) {
        # Load the .rda file in a local environment
        local_env <- new.env()
        load(rda_file, envir = local_env)
        
        # Check if 'imp' exists in the local environment
        if (exists("imp", envir = local_env)) {
          # Append the 'imp' data frame to the list
          all_data[[length(all_data) + 1]] <- local_env$imp
        }
      }
    }
    
    # Combine all data frames if there are any
    if (length(all_data) > 0) {
      combined_data <- do.call(rbind, all_data)
      
      # Write the combined data to an .rda file in the main folder
      save(combined_data, file = file.path(main_folder, "combined_data.rda"))
      
      return(TRUE)  # Return success indicator for this iteration
    } else {
      return(NULL)  # Return NULL if no data was combined
    }
  }
  
  return(results)  # Return the result list (TRUE for successful folders, NULL for empty ones)
}


process_simulation_data_parallel_r <- function(miss_type, simdata_list) {
  clusterExport(cl, varlist = c("impute_within_subgroup_r", "work_dir", "source_files"))
  
  # Define working directory
  miss_dir <- paste0(work_dir, "/", miss_type)
  dir.create(miss_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Parallelize the loops
  tryCatch({
    foreach(i = seq_along(simdata_list), .packages = c('foreach', 'readr', 'mice', "dplyr")) %:% 
      foreach(j = 4:6, .packages = c('foreach','readr', 'mice', "dplyr")) %dopar% {
        # Create a local directory for each task
        task_dir <- paste0(miss_dir, "/run_", i, "_levelspec_", j)
        dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Change working directory temporarily within each worker
        setwd(task_dir)
        
        # Perform the imputation
        impute_within_subgroup_r(simdata_list[[i]], source_files[[j]])
        
        # Restore the original working directory
        setwd(miss_dir)
      }
  }, error = function(e) {
    message("Error: ", e$message)
  })
}



compute_accuracy_for_subfolder <- function(subfolder) {
  # Extract scenario name from folder path
  scenario <- str_replace_all(subfolder, paste0(work_dir, "/"), "")
  scenario <- gsub("/", "_", scenario)
  run_num <- as.numeric(str_extract(scenario, "(?<=run_)\\d+"))
  
  # Load the combined data (use full path instead of changing working directory)
  load(file.path(subfolder, 'combined_data.rda')) # Load the combined data
  
  observed <- sampled_datasets[[run_num]]  # Access the appropriate observed dataset
  observed <- observed[order(observed$seqn), ]
  accuracy_results <- numeric(0)  # Store accuracy results for the current subfolder
  f1_score_results <- numeric(0)
  for (i in 1:5) {
    combined_subset <- subset(combined_data, MULT_ == i)
    combined_subset <- combined_subset[order(combined_subset$seqn), ]
    # Check if row counts match
    if (nrow(observed) == nrow(combined_subset)) {
      observed$activity_pattern <- as.factor(observed$activity_pattern)
      combined_subset$activity_pattern <- as.factor(combined_subset$activity_pattern)
      
      # Add missing levels from observed to combined_subset
      combined_subset$activity_pattern <- factor(combined_subset$activity_pattern,
                                                 levels = levels(observed$activity_pattern))
      
      # Replace missing predictions with 0 counts (optional, depending on how confusionMatrix handles it)
      combined_subset$activity_pattern[is.na(combined_subset$activity_pattern)] <- levels(observed$activity_pattern)[which(is.na(table(combined_subset$activity_pattern)))]
      
      # Now run the confusionMatrix
      confusion_matrix <- confusionMatrix(
        observed$activity_pattern,
        combined_subset$activity_pattern
      )
      accuracy <- confusion_matrix$overall['Accuracy']
      # Store the accuracy result
      accuracy_results <- c(accuracy_results, accuracy)
      precision <- confusion_matrix$byClass[, "Pos Pred Value"]  # Precision per class
      recall <- confusion_matrix$byClass[, "Sensitivity"]        # Recall per class
      f1_score_scores <- 2 * (precision * recall) / (precision + recall)
      support <- table(observed$activity_pattern)
      weighted_f1_score <- sum(f1_score_scores * (support / sum(support)))
      f1_score_results <- c(f1_score_results, weighted_f1_score)
    }
  }
  
  # Calculate and return the mean accuracy for the current subfolder
  mean_accuracy <- mean(accuracy_results)
  mean_f1_score <- mean(f1_score_results)
  return(data.frame(
    scenario = scenario,
    accuracy = mean_accuracy,
    f1_score = mean_f1_score
  ))
}



logreg_for_subfolder <- function(subfolder) {
  # Extract scenario name from folder path
  scenario <- str_replace_all(subfolder, paste0(work_dir, "/"), "")
  scenario <- gsub("/", "_", scenario)
  run_num <- as.numeric(str_extract(scenario, "(?<=run_)\\d+"))
  
  # Load the combined data
  load(file.path(subfolder, 'combined_data.rda')) 
  
  # Initialize vectors to store results
  auc_values <- numeric(5)
  R_squared_values <- numeric(5)
  
  # Loop through subsets
  for (i in 1:5) {
    combined_subset <- subset(combined_data, MULT_ == i)
    combined_subset <- combined_subset[order(combined_subset$seqn), ]
    combined_subset$mvpa_total <- combined_subset$modpa_total + combined_subset$vigpa_total
    
    # Model 2
    model <- glm(hypertension ~ as.factor(gender) + age + as.factor(race) + as.factor(edu) + mvpa_total_acc + as.factor(activity_pattern), 
                  data = combined_subset, 
                  family = binomial(link = 'probit'))
    
    # Calculate AUC and R-squared for Model 2
    auc_values[i] <- roc(combined_subset$hypertension, predict(model, combined_subset, type = "response"))$auc
    R_squared_values[i] <- as.numeric(pR2(model)[4])
  }
  
  # Calculate means for AUC and R-squared values
  mean_auc <- mean(auc_values, na.rm = TRUE)
  mean_R_squared <- mean(R_squared_values, na.rm = TRUE)
  
  # Return the results as a data frame
  return(data.frame(
    scenario = scenario,
    mean_auc = mean_auc,
    mean_R_squared = mean_R_squared
  ))
}



perform_kproto_clustering <- function(df, k) {
  # Step 1: Select and convert categorical variables to dummy variables
  categorical_vars <- df[, c("race", "edu", "marital", "poverty", "smoker", "alcohol_cat", "self_reported_health",
                             "fitness_access", "health_literacy")]
  dummy_vars1 <- as.data.frame(model.matrix(~ . - 1, data = categorical_vars)) %>%
    mutate(across(everything(), as.factor))
  # Step 2: Select and convert additional categorical variables to factors
  dummy_vars2 <- df[, c("hypertension", "heartdiseases", "stroke", "cancers", "diabetes", 
                        "gender", "work", "insurance", "srvy_yr", "source")] %>%
    mutate(across(everything(), as.factor))
  # Step 3: Scale continuous variables
  continuous_vars <- scale(df[, c("age", "bmi", "modpa_total", "vigpa_total")])
  # Step 4: Combine all variables into one data frame
  combined_data <- cbind(dummy_vars1, dummy_vars2, continuous_vars)
  # Step 5: Perform k-prototypes clustering
  kmode <- kproto(combined_data, k)
  # Step 6: Add the cluster assignments back to the original data frame
  df$cluster <- kmode$cluster
  
  return(df)
}
