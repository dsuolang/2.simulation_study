library(haven)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
#library(devtools)
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
library(clustMixType)
library(pROC)
options(scipen = 999)

path <- "/Users/dsuolang/Desktop/Study2/Data/synthpop.dta"
synthpop <- haven::read_dta(path)

# Continuous variables
cont_vars <- c("age", "bmi", "mvpa_total_acc",
               "modpa_total", "vigpa_total")

# Categorical variables
cat_vars <- c("race","gender","marital","edu","poverty","work","insurance",
              "self_reported_health","smoker","alcohol_cat",
              "hypertension","diabetes","heartdiseases","cancers","stroke",
              "fitness_access","health_literacy","activity_pattern")


# --- Continuous summaries ---
cont_summary <- synthpop %>%
  select(all_of(cont_vars)) %>%
  summarise(across(
    everything(),
    ~ c(mean = mean(.x, na.rm = TRUE),
        se   = sd(.x, na.rm = TRUE)/sqrt(sum(!is.na(.x))))
  )) %>%
  pivot_longer(cols = everything(),
               names_to = c("variable", ".value"),
               names_sep = "_")

# --- Categorical summaries ---
n_total <- nrow(synthpop)

cat_summary <- synthpop %>%
  select(all_of(cat_vars)) %>%
  mutate(across(everything(), as.character)) %>%   # <-- ensure same type
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "level") %>%
  filter(!is.na(level)) %>%
  group_by(variable, level) %>%
  summarise(prop = n()/n_total, .groups = "drop") %>%
  mutate(se = sqrt(prop * (1 - prop) / n_total))

# --- Final results ---
list(
  continuous = cont_summary,
  categorical = cat_summary
)