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
  predicted <- predict(logistic_model, type = "response") # probability of being "not missing" category (NHANES)
  
  data[[paste0("miss_indicator")]] <- as.numeric(ifelse(predicted > quantile(predicted, miss_perc/100), TRUE, FALSE))
  return(data)
}

introduce_missingness <- function(df, cols_to_make_missing, predictors, missing_rate) {
  # Include the source variable as the outcome.
  formula <- as.formula(paste("source", "~", paste(predictors, collapse = " + ")))
  
  # Generate a logistic regression model based on the predictors and source variable.
  logit_model <- glm(formula, data = df, family = binomial(link = "logit"))
  
  # Calculate the probability of missingness for each row. (prob of being in NHANES (not missing))
  df$prob_missing <- predict(logit_model, type = "response")
  
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
  lower = estimate_mi - t_value * sqrt(mi_vt)
  upper = estimate_mi + t_value * sqrt(mi_vt)
  coverage <- as.numeric((benchmark >= lower) & (benchmark <= upper))
  
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
    t_value = t_value,
    lower = lower,
    upper = upper,
    coverage = coverage,
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
    t_value <- qt(0.975, df)
    lower = estimate_mi - t_value * sqrt(mi_vt)
    upper = estimate_mi + t_value * sqrt(mi_vt)
    coverage <- as.numeric((benchmarks[[i]] >= lower) & (benchmarks[[i]] <= upper))
    
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
      t_value = t_value,
      lower = lower,
      upper = upper,
      coverage = coverage,
      lambda = lambda
    ))
  }
}


rubin_combining_rule_categorical_transformation <- function(variable, benchmarks, scenario) {
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
    p_l <- df_dummies %>%
      group_by(MULT_) %>%
      summarise(p_l = mean(.data[[activity_cols[[i]]]]))
    
    # Log-transformed proportion (t_l)
    p_l$t_l <- log(p_l$p_l / (1 - p_l$p_l))
    
    # Calculate the variance for each imputation (v_l)
    p_l$v_l <- 1 / (length(df_dummies[[activity_cols[[i]]]])/m * p_l$p_l * (1 - p_l$p_l))
    
    # Apply Rubin's rule: Pooled estimate (estimate_mi)
    estimate_mi_logit <- mean(p_l$t_l)  # Pooled estimate for proportions
    
    # Within-imputation variance (mi_vw) (average variance within imputations)
    mi_vw_logit <- mean(p_l$v_l)
    
    # Between-imputation variance (mi_vb) (variance of the estimates)
    mi_vb_logit <- var(p_l$t_l)
    
    # Total variance (T)
    mi_vt_logit <- mi_vw_logit + (1 + 1/m) * mi_vb_logit
    
    # Fraction of Missing Information (lambda)
    lambda_logit <- (1 + 1/m) * mi_vb_logit / mi_vt_logit
    
    df_logit <- (m-1)/lambda_logit^2
    t_value_logit <- qt(0.975, df_logit)
    lower_logit <- estimate_mi_logit - t_value_logit * sqrt(mi_vt_logit)
    upper_logit <- estimate_mi_logit + t_value_logit * sqrt(mi_vt_logit)
    lower <- 1 / (1 + exp(-lower_logit))
    upper <- 1 / (1 + exp(-upper_logit))
    
    coverage <- as.numeric((benchmarks[[i]] >= lower) & (benchmarks[[i]] <= upper))
    
    # inverse logit 
    estimate_mi <- 1 / (1 + exp(-estimate_mi_logit))
    bias <- estimate_mi - benchmarks[[i]]
    relative_bias <- bias/benchmarks[[i]]
    
    # Append results to result_sheet
    result_sheet_eachrun <<- rbind(result_sheet_eachrun, data.frame(
      scenario = scenario,
      variable = activity_cols[[i]],
      estimate_mi = estimate_mi,
      mi_vw = 1 / (1 + exp(-mi_vw_logit)),
      mi_vb = 1 / (1 + exp(-mi_vb_logit)),
      mi_vt = 1 / (1 + exp(-mi_vt_logit)),
      bias = bias,
      relative_bias = relative_bias,
      df = df_logit,
      t_value = t_value_logit,
      lower = lower,
      upper = upper,
      coverage = coverage,
      lambda = 1 / (1 + exp(-lambda_logit))
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
      rubin_combining_rule_categorical_transformation(variable, benchmarks, scenario)
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
  scenario <- str_replace_all(subfolder, paste0(work_dir, "/"), "")   
  scenario <- gsub("/", "_", scenario)   
  run_num <- as.numeric(str_extract(scenario, "(?<=run_)\\d+"))      
  
  load(file.path(subfolder, 'combined_data.rda')) # Load the combined data      
  
  observed <- sampled_datasets[[run_num]]  # Access the appropriate observed dataset   
  observed <- observed[order(observed$seqn), ]   
  accuracy_results <- numeric(0)  # Store accuracy results for the current subfolder   
  
  for (i in 1:5) {     
    combined_subset <- subset(combined_data, MULT_ == i)     
    combined_subset <- combined_subset[order(combined_subset$seqn), ]     
    
    # Debugging checks
    print(paste("Processing MULT_ =", i))
    print(paste("Observed rows:", nrow(observed), "Combined rows:", nrow(combined_subset)))
    
    if (nrow(observed) == nrow(combined_subset)) {       
      # Convert haven_labelled to factor
      library(haven)
      observed$activity_pattern <- as_factor(observed$activity_pattern)
      combined_subset$activity_pattern <- as_factor(combined_subset$activity_pattern)
      
      # Align factor levels
      combined_subset$activity_pattern <- factor(
        combined_subset$activity_pattern, 
        levels = levels(observed$activity_pattern)
      )
      
      # Check NA counts
      print(paste("NA in combined_subset activity_pattern:", sum(is.na(combined_subset$activity_pattern))))
      
      # Run confusionMatrix
      confusion_matrix <- confusionMatrix(
        observed$activity_pattern, 
        combined_subset$activity_pattern
      )
      print(confusion_matrix)
      
      accuracy <- confusion_matrix$overall['Accuracy']       
      accuracy_results <- c(accuracy_results, accuracy)     
    } else {
      print("Row counts do not match! Skipping MULT_ =", i)
    }
  }      
  
  mean_accuracy <- mean(accuracy_results, na.rm = TRUE)   
  return(data.frame(scenario = scenario, accuracy = mean_accuracy)) 
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


calculate_ci_with_custom_variance <- function(data_list, variable, alpha = 0.05) {
  ci_results <- lapply(seq_along(data_list), function(i) {
    dataset <- data_list[[i]]
    
    if (variable %in% colnames(dataset)) {
      # Select variance calculation method based on the variable name
      if (variable == "mvpa_total_acc") {
        variance <- calculate_squared_se(dataset[[variable]])
      } else {
        variance <- calculate_squared_se_prop(dataset[[variable]])
      }
      
      # Calculate mean and standard error
      mean_val <- mean(dataset[[variable]])
      se <- sqrt(variance)  
      
      # Calculate sample size
      n <- sum(!is.na(dataset[[variable]]))
      
      # Confidence interval
      ci <- c(
        run = i,
        lower = mean_val - qt(1 - alpha / 2, df = n - 1) * se,
        upper = mean_val + qt(1 - alpha / 2, df = n - 1) * se
      )
    } else {
      # If variable doesn't exist in the dataset, return NA
      ci <- c(run = i, lower = NA, upper = NA)
    }
    
    return(ci)
  })
  
  # Combine results into a data frame
  ci_df <- do.call(rbind, ci_results)
  ci_df <- as.data.frame(ci_df)
  ci_df$variable = variable
  colnames(ci_df) <- c("run", "lower", "upper", "variable")
  
  return(ci_df)
}

logreg_for_subfolder <- function(subfolder, outcome_var) {
  # Extract scenario name from folder path
  scenario <- str_replace_all(subfolder, paste0(work_dir, "/"), "")
  scenario <- gsub("/", "_", scenario)
  run_num <- as.numeric(str_extract(scenario, "(?<=run_)\\d+"))
  scenario_end <- str_sub(scenario, -1)

  # Load the combined data
  load(file.path(subfolder, 'combined_data.rda')) 
  
  if (scenario_end %in% c("1", "2", "3")) {
    combined_data$mvpa_total_acc <- combined_data$mvpa_total_acc_sqrt^2
  }
  # Initialize vectors to store results
  prauc_values <- numeric(5)
  pseudo_R2_values <- numeric(5)
  deviance_values <- numeric(5)
  aic_values <- numeric(5)
  bic_values <- numeric(5)
  
  for (i in 1:5) {
    combined_subset <- subset(combined_data, MULT_ == i)
    combined_subset <- combined_subset[order(combined_subset$seqn), ]
    combined_subset$mvpa_total <- combined_subset$modpa_total + combined_subset$vigpa_total
    
    # Create model formulas dynamically
    formula0 <- as.formula(paste(outcome_var, "~ age + as.factor(race) + bmi + mvpa_total"))
    formula1 <- as.formula(paste(outcome_var, "~ age + as.factor(race) + bmi + mvpa_total_acc + as.factor(activity_pattern)"))
    
    model0 <- glm(formula0, data = combined_subset, family = binomial(link = 'probit'))
    model <- glm(formula1, data = combined_subset, family = binomial(link = 'probit'))
    
    # Prediction probabilities for PR curve
    preds <- predict(model, combined_subset, type = "response")
    labels <- combined_subset[[outcome_var]]
    
    # PRAUC
    pr <- pr.curve(scores.class0 = preds[labels == 1],
                   scores.class1 = preds[labels == 0],
                   curve = FALSE)
    prauc_values[i] <- pr$auc.integral
    
    # Pseudo R2 (McFadden)
    pseudo_R2_values[i] <- as.numeric(pR2(model)["McFadden"])
    
    # Deviance, AIC, BIC
    deviance_values[i] <- deviance(model)
    aic_values[i] <- AIC(model)
    bic_values[i] <- BIC(model)
  }
  
  # Calculate mean metrics
  return(data.frame(
    scenario = scenario,
    mean_pr_auc = mean(prauc_values, na.rm = TRUE),
    mean_pseudo_R2 = mean(pseudo_R2_values, na.rm = TRUE),
    mean_deviance = mean(deviance_values, na.rm = TRUE),
    mean_AIC = mean(aic_values, na.rm = TRUE),
    mean_BIC = mean(bic_values, na.rm = TRUE)
  ))
}

logreg_sampledata <- function(sampled_datasets, outcome_var) {
  # Initialize result storage
  results_list <- vector("list", length(sampled_datasets))
  
  for (i in seq_along(sampled_datasets)) {
    tryCatch({
      data_i <- sampled_datasets[[i]]
      
      # Dynamically calculate outcome and predictors
      data_i$mvpa_total <- data_i$modpa_total + data_i$vigpa_total
      outcome <- data_i[[outcome_var]]
      
      # Fit logistic regression model (Probit link)
      model <- glm(
        formula = as.formula(paste(outcome_var, "~ age + as.factor(race) + bmi + mvpa_total")),
        data = data_i,
        family = binomial(link = 'probit')
      )
      
      # Predictions for metrics
      preds <- predict(model, data_i, type = "response")
      
      # Calculate PR AUC (Precision-Recall AUC) using PRROC package
      pr <- PRROC::pr.curve(scores.class0 = preds[outcome == 1], 
                            scores.class1 = preds[outcome == 0], 
                            curve = FALSE)
      pr_auc <- pr$auc.integral
      
      # Calculate pseudo R²
      pseudo_R2 <- as.numeric(pscl::pR2(model)[4])
      
      # Model evaluation metrics
      deviance_val <- model$deviance
      aic_val <- AIC(model)
      bic_val <- BIC(model)
      
      # Store results
      results_list[[i]] <- data.frame(
        mean_pr_auc = pr_auc,
        mean_pseudo_R2 = pseudo_R2,
        mean_deviance = deviance_val,
        mean_AIC = aic_val,
        mean_BIC = bic_val
      )
    }, error = function(e) {
      # In case of error, store NAs
      results_list[[i]] <- data.frame(
        mean_pr_auc = NA,
        mean_pseudo_R2 = NA,
        mean_deviance = NA,
        mean_AIC = NA,
        mean_BIC = NA
      )
    })
  }
  
  # Combine results into a final dataframe
  results <- do.call(rbind, results_list)
  
  return(results)
}

process_model <- function(dataset) {
  # Define the models
  model1 <- list(
    lm(mvpa_total_acc_sqrt ~ modpa_total + vigpa_total + modpa_indicator + vigpa_indicator + activity_pattern + 
         source + srvy_yr +
         age + race + gender + marital, data = dataset),
    multinom(activity_pattern ~ modpa_total + vigpa_total + modpa_indicator + vigpa_indicator  + mvpa_total_acc_sqrt + 
               source + srvy_yr +
               age + race + gender + marital, data = dataset)
  )
  
  model2 <- list(
    lm(mvpa_total_acc_sqrt ~ modpa_total + vigpa_total + modpa_indicator + vigpa_indicator + activity_pattern + 
         source + srvy_yr +
         age + race + gender + marital + edu + poverty + work + insurance + fitness_access, data = dataset),
    multinom(activity_pattern ~ modpa_total + vigpa_total + modpa_indicator + vigpa_indicator  + mvpa_total_acc_sqrt + 
               source + srvy_yr +
               age + race + gender + marital + edu + poverty + work + insurance + fitness_access, data = dataset)
  )
  
  model3 <- list(
    lm(mvpa_total_acc_sqrt ~ modpa_total + vigpa_total + modpa_indicator + vigpa_indicator + activity_pattern + 
         source + srvy_yr +
         age + race + gender + marital + edu + poverty + work + insurance +
         self_reported_health + bmi + smoker + alcohol_cat + hypertension + diabetes + heartdiseases + cancers + stroke +
         fitness_access + health_literacy, data = dataset),
    multinom(activity_pattern ~ modpa_total + vigpa_total + modpa_indicator + vigpa_indicator + mvpa_total_acc_sqrt + 
               source + srvy_yr +
               age + race + gender + marital + edu + poverty + work + insurance +
               self_reported_health + bmi + smoker + alcohol_cat + hypertension + diabetes + heartdiseases + cancers + stroke +
               fitness_access + health_literacy, data = dataset)
  )
  
  # Collect the results for each model
  model_results <- data.frame(
    model = rep(1:3, each = 2),  # 1, 2, and 3 for each model
    variable = rep(c("mvpa_total_acc_sqrt", "activity_pattern"), 3),
    rsquared = c(
      summary(model1[[1]])$adj.r.squared, pR2(model1[[2]])["r2ML"],
      summary(model2[[1]])$r.squared, pR2(model2[[2]])["r2ML"],
      summary(model3[[1]])$r.squared, pR2(model3[[2]])["r2ML"]
    ),
    dataset_id = rep(deparse(substitute(dataset)), 6)  # Add dataset name to identify it
  )
  
  return(model_results)
}
