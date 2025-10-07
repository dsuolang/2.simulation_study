library(ggplot2)
library(gridExtra)
library(grid) 
library(cowplot)
library(dplyr)
library(forcats)
library(scales)
library(stringr)
library(tidyr)
library(ggh4x) 
library(cowplot)

work_dir<-"/Users/dsuolang/Desktop/Study2/hpc_results"
setwd(work_dir)

srmi<- read.csv('result_sheet_srmi.csv')
srmi$method = 'REG'
pmm <- read.csv('result_sheet_pmm.csv')
pmm<-pmm[pmm$model_spec=="w/o_Transf",]
pmm$method = 'PMM'
srmi_pmm<-rbind(srmi, pmm)


datasets <- list(
  srmi_pmm[srmi_pmm$variable == "mvpa_total_acc", ],
  srmi_pmm[srmi_pmm$variable == "activity_pattern_1", ],
  srmi_pmm[srmi_pmm$variable == "activity_pattern_2", ],
  srmi_pmm[srmi_pmm$variable == "activity_pattern_3", ],
  srmi_pmm[srmi_pmm$variable == "activity_pattern_4", ]
)

titles <- c("MVPA Duration (Wearable)",
            "Activity Pattern 1", "Activity Pattern 2", 
            "Activity Pattern 3", "Activity Pattern 4")


# Part 1. Rel_Bias SE RMSE for MVPA
# Function to create plots for a given dataset
create_plots <- function(data, title_text) {
  palette <- scale_color_manual(values = c("lightblue", "blue", "darkblue", "peachpuff", "orange", "orange3"))
  # Plot 1: Relative Bias
  p1 <- ggplot(data, aes(x = miss_rate, y = relative_bias, color = interaction(shared_var, method), 
                         linetype = model_spec, shape = interaction(shared_var, method))) + 
    geom_line() + 
    geom_point(size = 3) + 
    labs(x = "", y = "Rel.Bias") +
    scale_x_continuous(breaks = unique(data$miss_rate)) + 
    scale_y_continuous(
      limits = c(-0.1, 0.15),
      breaks = c(-0.1, -0.05, 0, 0.05, 0.1, 0.15),
      labels = number_format(accuracy = 0.01)
    )+
    #scale_y_continuous(labels = number_format(accuracy = 0.01)) +  # Format y-axis to two decimals
    scale_shape_manual(values = c(
      "Least.REG" = 16, "Moderate.REG" = 16, "Most.REG" = 16,  
      "Least.PMM" = 17, "Moderate.PMM" = 17, "Most.PMM" = 17))+
    scale_linetype_manual(values = c("w/o_Transf" = "solid", "w/_Transf" = "dashed")) +
    palette + 
    theme_minimal() + 
    theme(legend.position = "none", 
          axis.title.x = element_text(size = 12),  # Increase x-axis label size
          axis.title.y = element_text(size = 12),  # Increase y-axis label size
          axis.text.x = element_text(size = 12),   # Increase x-axis tick label size
          axis.text.y = element_text(size = 12),   # Increase y-axis tick label size
          strip.text = element_text(size = 12)) + 
    facet_wrap(~ miss_mech) +
    theme(panel.spacing.x = unit(1, "cm")) 
  # Plot 2: Standard Error
  p2 <- ggplot(data, aes(x = miss_rate, y = se, color = interaction(shared_var, method), 
                         linetype = model_spec, shape = interaction(shared_var, method))) + 
    geom_line() + 
    geom_point(size = 3) + 
    labs(x = "Missing Rate", y = "SE") +
    scale_x_continuous(breaks = unique(data$miss_rate)) + 
    #scale_y_continuous(breaks = c(0.65, 0.70, 0.75)) + 
    scale_shape_manual(values = c(
      "Least.REG" = 16, "Moderate.REG" = 16, "Most.REG" = 16,  
      "Least.PMM" = 17, "Moderate.PMM" = 17, "Most.PMM" = 17)) +
    scale_linetype_manual(values = c("w/o_Transf" = "solid", "w/_Transf" = "dashed")) +
    palette + 
    theme_minimal() + 
    theme(legend.position = "none", 
          axis.title.x = element_text(size = 12),  # Increase x-axis label size
          axis.title.y = element_text(size = 12),  # Increase y-axis label size
          axis.text.x = element_text(size = 12),   # Increase x-axis tick label size
          axis.text.y = element_text(size = 12),   # Increase y-axis tick label size
          strip.text = element_text(size = 12)) + 
    facet_wrap(~ miss_mech) + coord_cartesian(ylim = c(0,20)) +
    theme(panel.spacing.x = unit(1, "cm")) #+ coord_cartesian(ylim = c(0.65, 0.75))  coord_cartesian(ylim = c(0,30))  # Set y-axis limit
  # Plot 3: RMSE
  p3 <- ggplot(data, aes(x = miss_rate, y = rmse, color = interaction(shared_var, method), 
                         linetype = model_spec, shape = interaction(shared_var, method))) + 
    geom_line() + 
    geom_point(size = 3) + 
    labs(x = "", y = "RMSE") +
    scale_x_continuous(breaks = unique(data$miss_rate)) + 
    #scale_y_continuous(breaks = c(0.65, 0.70, 0.75)) + 
    scale_shape_manual(values = c(
      "Least.REG" = 16, "Moderate.REG" = 16, "Most.REG" = 16,  
      "Least.PMM" = 17, "Moderate.PMM" = 17, "Most.PMM" = 17)) +
    scale_linetype_manual(values = c("w/o_Transf" = "solid", "w/_Transf" = "dashed")) +
    palette + 
    theme_minimal() + 
    theme(legend.position = "none", 
          axis.title.x = element_text(size = 12),  # Increase x-axis label size
          axis.title.y = element_text(size = 12),  # Increase y-axis label size
          axis.text.x = element_text(size = 12),   # Increase x-axis tick label size
          axis.text.y = element_text(size = 12),   # Increase y-axis tick label size
          strip.text = element_text(size = 12)) + 
    facet_wrap(~ miss_mech) + coord_cartesian(ylim = c(0,30)) +
    theme(panel.spacing.x = unit(1, "cm")) #+ coord_cartesian(ylim = c(0.65, 0.75))  # Set y-axis limit
  
  # Legend
  legend_plot <- get_legend(p1 + theme(legend.position = "right", legend.title = element_blank(),
                                       legend.text = element_text(size = 12)))
  # Combine Plots
  final_plot <- plot_grid(
    plot_grid(p1, p2, p3, ncol = 3),
    #legend_plot,
    ncol = 1, rel_widths = c(3, 0.4))
  
  # Add Title
  title <- ggdraw() + draw_label(title_text, fontface = 'bold', size = 12)
  plot_grid(title, final_plot, ncol = 1, rel_heights = c(0.1, 1))
}

# Generate and display plots
#plot1 <- lapply(1:length(datasets), function(i) create_plots(datasets[[i]], titles[i]))
plot1 <- create_plots(datasets[[1]], titles[[1]])
plot1 
grid.newpage()
grid.arrange(
  plot1,
  nullGrob(),                    # empty space
  heights = unit.c(unit(1, "null"), unit(1, "in")) # first fills, second = 1 inch
)
legend_plot <- readRDS("legend_plot.rds")
final<-grid.arrange(
  plot1,              # empty space on top
  legend_plot,             # your legend
  heights = c(9, 1)        # push legend to bottom (adjust ratio)
)  
# plots 1050 550
ggsave(
  filename = "mvpa_bias_se_rmse.png",  # file extension determines format
  plot = final,                 # which plot to save
  dpi = 300,               # resolution in dots per inch
  width = 11, height = 5,    # physical size in inches
  units = "in"              # units for width/height
)



# Part 2. Bias,SE for activity pattern
# Function to create plots for a given dataset
create_plots2 <- function(data, title_text) {
  palette <- scale_color_manual(values = c("lightblue", "blue", "darkblue", "peachpuff", "orange", "orange3"))
  # Plot 1: Relative Bias
  p1 <- ggplot(data, aes(x = miss_rate, y = bias, color = interaction(shared_var, method), 
                         linetype = model_spec, shape = interaction(shared_var, method))) + 
    geom_line() + 
    geom_point(size = 3) + 
    labs(x = "Missing Rate", y = "Bias") +
    scale_x_continuous(breaks = unique(data$miss_rate)) + 
    scale_y_continuous(
      limits = c(-0.3, 0.3), 
      breaks = seq(-0.3, 0.3, by = 0.15), 
      labels = number_format(accuracy = 0.01)
    ) +
    #scale_y_continuous(labels = number_format(accuracy = 0.01)) +  # Format y-axis to two decimals
    scale_shape_manual(values = c(
      "Least.REG" = 16, "Moderate.REG" = 16, "Most.REG" = 16,  
      "Least.PMM" = 17, "Moderate.PMM" = 17, "Most.PMM" = 17))+
    scale_linetype_manual(values = c("w/o_Transf" = "solid", "w/_Transf" = "dashed")) +
    palette + 
    theme_minimal() + 
    theme(legend.position = "none", 
          axis.title.x = element_text(size = 12),  # Increase x-axis label size
          axis.title.y = element_text(size = 12),  # Increase y-axis label size
          axis.text.x = element_text(size = 12),   # Increase x-axis tick label size
          axis.text.y = element_text(size = 12),   # Increase y-axis tick label size
          strip.text = element_text(size = 12)) + 
    facet_wrap(~ miss_mech)  +
    theme(panel.spacing.x = unit(1, "cm")) 
  # Plot 2: Standard Error
  p2 <- ggplot(data, aes(x = miss_rate, y = se, color = interaction(shared_var, method), 
                         linetype = model_spec, shape = interaction(shared_var, method))) + 
    geom_line() + 
    geom_point(size = 3) + 
    labs(x = "Missing Rate", y = "SE") +
    scale_x_continuous(breaks = unique(data$miss_rate)) +
    #scale_y_continuous(limits = c(0.6, 0.7), breaks = seq(0.6, 0.75, by = 0.05), 
      #labels = number_format(accuracy = 0.01)) +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +  # Format y-axis to two decimals
    scale_shape_manual(values = c(
      "Least.REG" = 16, "Moderate.REG" = 16, "Most.REG" = 16,  
      "Least.PMM" = 17, "Moderate.PMM" = 17, "Most.PMM" = 17)) +
    scale_linetype_manual(values = c("w/o_Transf" = "solid", "w/_Transf" = "dashed")) +
    palette + 
    theme_minimal() + 
    theme(legend.position = "none", 
          axis.title.x = element_text(size = 12),  # Increase x-axis label size
          axis.title.y = element_text(size = 12),  # Increase y-axis label size
          axis.text.x = element_text(size = 12),   # Increase x-axis tick label size
          axis.text.y = element_text(size = 12),   # Increase y-axis tick label size
          strip.text = element_text(size = 12)) + 
    facet_wrap(~ miss_mech) +  theme(panel.spacing.x = unit(1, "cm")) + coord_cartesian(ylim = c(0.65, 0.75))  # Set y-axis limit
  
  # Combine Plots
  final_plot <- plot_grid(
    plot_grid(p1, p2, ncol = 2),
    #legend_plot,
    ncol = 2, rel_widths = c(3, 0.4))
  
  legend_plot <- get_legend(
    p1 +
      guides(
        shape    = guide_legend(order = 1, nrow = 1, byrow = TRUE),
        color    = guide_legend(order = 1, nrow = 1, byrow = TRUE),
        linetype = guide_legend(order = 1, nrow = 1, byrow = TRUE)
      ) +
      theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.box = "horizontal",          # force everything in one box/row
        legend.text = element_text(size = 12)
      )
  )
  saveRDS(legend_plot, file = "legend_plot.rds", version = 3)  # ensures a modern format
  class(legend_plot) 
  grid.arrange(
    nullGrob(),              # empty space on top
    legend_plot,             # your legend
    heights = c(5, 1)        # push legend to bottom (adjust ratio)
  )  
  # Add Title
  title <- ggdraw() + draw_label(title_text, fontface = 'bold', size = 12)
  plot_grid(title, final_plot, ncol = 1, rel_heights = c(0.1, 1))
}

# List of datasets and titles
datasets <- datasets
titles <- titles

# Generate and display plots
plots <- lapply(2:length(datasets), function(i) create_plots2(datasets[[i]], titles[i]))

# Step 3: Arrange plots and add the legend
plot2<-do.call(grid.arrange, c(plots[1:4], nrow = 2, ncol = 2))
grid.newpage()
grid.arrange(
  plot2,
  nullGrob(),                    # empty space
  heights = unit.c(unit(1, "null"), unit(1, "in")) # first fills, second = 1 inch
)
legend_plot <- readRDS("legend_plot.rds")
final<-grid.arrange(
  plot2,              # empty space on top
  legend_plot,             # your legend
  heights = c(9, 1)        # push legend to bottom (adjust ratio)
)  
# 1300 * 750
ggsave(
  filename = "actpattern_bias_se.png",  # file extension determines format
  plot = final,                 # which plot to save
  dpi = 300,               # resolution in dots per inch
  width = 11, height = 10,    # physical size in inches
  units = "in"              # units for width/height
)


# Part 3 : FMI 
plot_fmi <- function(data, title_text) {
  palette <- scale_color_manual(values = c("orange3", "orange3", "orange3", "lightblue3", "lightblue3", "lightblue3"))
  # Create a new x_label column without breaking the text into separate lines
  data$x_label <- paste(data$shared_var, data$miss_rate, sep = " ")  # Concatenate shared_var and miss_rate
  desired_order <- rev(c(
    "Least 50", "Moderate 50", "Most 50",
    "Least 70", "Moderate 70", "Most 70",
    "Least 90", "Moderate 90", "Most 90"
  ))
  data$group_label <- interaction(data$miss_mech, data$method, data$model_spec, sep = ".")
  
  desired_legend_order <- c(
    "MAR.PMM.w/o_Transf",
    "MAR.REG.w/o_Transf",
    "MAR.REG.w/_Transf",
    "MNAR.PMM.w/o_Transf",
    "MNAR.REG.w/o_Transf",
    "MNAR.REG.w/_Transf"
  )
  
  # Set factor level order
  data$group_label <- factor(data$group_label, levels = desired_legend_order)
  # Step 3: Convert to factor with the desired order
  data$x_label <- factor(data$x_label, levels = desired_order)
  ggplot(data, aes(
    x = lambda,  # X-axis is now lambda (FMI)
    y = x_label,
    shape = group_label, 
    color = group_label
  )) +
    scale_x_continuous(labels = number_format(accuracy = 0.1), 
                       limits = c(min(data$lambda) - 0.1, max(data$lambda) + 0.1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +  # Shift X-Axis slightly
    scale_y_discrete(expand = c(0.05, 0.05)) +  # Shift Y-Axis slightly
    geom_point(size = 4) +  
    labs(
      x =  expression(lambda),
      y = "Shared Variables * %Missing",  # Label for y-axis
      title = title_text
    ) +  
    palette +  
    theme_minimal() + 
    theme(
      axis.title.x = element_text(size = 12),  # Increase x-axis label size
      axis.title.y = element_text(size = 12),  # Increase y-axis label size
      axis.text.x = element_text(size = 12),   # Increase x-axis tick label size
      axis.text.y = element_text(size = 12),   # Increase y-axis tick label size
      strip.text = element_text(size = 12),
      #legend.position = "none",  
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.45, size = 12, face = "bold"),  
      axis.line = element_line(color = "black", linewidth = 0.1)  
    ) +
    scale_shape_manual(values = c(16, 17, 1, 16, 17, 1))
}


# List of datasets and titles
datasets <- datasets
# Generate `p4` plots for each dataset
fmi_plot_mvpa <- plot_fmi(datasets[[1]],"MVPA Duration (Wearable)")
fmi_plot_mvpa

combiend_dataset2_5<-as.data.frame(cbind(
  (datasets[[2]] %>% arrange(scenario))$lambda,
  (datasets[[3]] %>% arrange(scenario))$lambda,
  (datasets[[4]] %>% arrange(scenario))$lambda,
  (datasets[[5]] %>% arrange(scenario))$lambda))
lambda<-rowMeans(combiend_dataset2_5)
datasets[[2]]$lambda<-rowMeans(combiend_dataset2_5)

fmi_plot_activity_pattern <- plot_fmi(datasets[[2]], "Activity Pattern")
fmi_plot_activity_pattern

ggsave(
  filename = "actpattern_bias_se.png",  # file extension determines format
  plot = final,                 # which plot to save
  dpi = 300,               # resolution in dots per inch
  width = 11, height = 10,    # physical size in inches
  units = "in"              # units for width/height
)


# Part 4: Log Models Analysis #800*500

# Load Data
srmi <- read.csv('log_result_hypertension_srmi.csv')
srmi$method <- 'REG'

pmm <- read.csv('log_result_hypertension_pmm.csv')
pmm$method <- 'PMM'

# Combine Data
logresult <- rbind(srmi, pmm)
logresult$mean_pseudo_R2<-round(logresult$mean_pseudo_R2,2)
# Data Transformation
logresult <- logresult %>%
  mutate(
    # Extract Missing Mechanism
    miss_mech = case_when(
      grepl("mar", scenario_avg_run, ignore.case = TRUE) ~ "MAR",
      grepl("mnar", scenario_avg_run, ignore.case = TRUE) ~ "MNAR"
    ),
    # Extract Missing Rate
    miss_rate = gsub(".*_(\\d+)(mar|mnar)_levelspec.*", "\\1", scenario_avg_run),
    # Extract Level Specification
    levelspec = gsub(".*levelspec_(\\d+).*", "\\1", scenario_avg_run),
    # Assign Shared Variables
    shared_var = case_when(
      levelspec == "1" ~ "Least",
      levelspec == "2" ~ "Moderate",
      levelspec == "3" ~ "Most",
      levelspec == "4" ~ "Least",
      levelspec == "5" ~ "Moderate",
      levelspec == "6" ~ "Most"
    ),
    # Assign Model Specification
    model_spec = case_when(
      levelspec %in% c("1", "2", "3") ~ "w/_Transf", 
      levelspec %in% c("4", "5", "6") ~ "w/o_Transf"
    )
  ) %>%
  select(-levelspec)  # Drop the temporary levelspec column

# Filter out PMM results with transformation
logresult <- logresult[!(logresult$model_spec == "w/_Transf" & logresult$method == "PMM"), ]

# Rename for clarity
names(logresult)[names(logresult) == "scenario_avg_run"] <- "scenario"



# Reshape Data to Long Format
logresult_long <- logresult %>%
  pivot_longer(
    cols = c("mean_deviance", "mean_AIC", "mean_BIC","mean_pseudo_R2"),
    names_to = "metric",
    values_to = "value"
  )

# Clean Up Metric Factor Levels
logresult_long$metric <- factor(
  logresult_long$metric, 
  levels = c("mean_deviance", "mean_AIC", "mean_BIC", "mean_pseudo_R2"),
  labels = c("Deviance", "AIC", "BIC", "Pseudo_R2")
)

# Ensure miss_mech is a Factor with Correct Order
logresult_long$miss_mech <- factor(logresult_long$miss_mech, levels = c("MAR", "MNAR"))

# Define Color Palette
palette <- scale_color_manual(values = c(
  "Least.REG" = "peachpuff", 
  "Moderate.REG" = "orange", 
  "Most.REG" = "orange3", 
  "Least.PMM" = "lightblue", 
  "Moderate.PMM" = "blue", 
  "Most.PMM" = "darkblue"
))

# Ensure 'miss_mech' is a factor with desired order
logresult_long$miss_mech <- factor(logresult_long$miss_mech, levels = c("MAR", "MNAR"))
# Prepare self-report lines data
mean_lines <- data.frame(
  metric = c("AIC", "BIC", "Deviance", "Pseudo_R2"),
  mean_value = c(9784, 9835, 9770, 0.212)
  #mean_value = c(5705, 5756, 5691, 0.171)  
  #mean_value = c(5988, 6038, 5974, 0.135)  
)
# mean_pr_auc     mean_pseudo_R2   mean_deviance      mean_AIC        mean_BIC 
# Mean   :0.6302   Mean   :0.2118   Mean   :9770   Mean   :9784   Mean   :9835   Hypertension
# Mean   :0.3091   Mean   :0.1707   Mean   :5691   Mean   :5705   Mean   :5756   Diabetes
# Mean   :0.2923   Mean   :0.1351   Mean   :5974   Mean   :5988   Mean   :6038   Heart diesease
p <- ggplot(
  logresult_long, 
  aes(
    x = as.factor(miss_rate), 
    y = value, 
    color = interaction(shared_var, method),
    linetype = model_spec, 
    shape = interaction(shared_var, method)
  )
) + 
  geom_line(aes(group = interaction(shared_var, method, model_spec))) +
  geom_point(size = 3) + 
  geom_hline(
    data = mean_lines, 
    aes(yintercept = mean_value), 
    color = "red", 
    linetype = "dashed", 
    size = 1
  ) + 
  labs(x = "% Missing", y = "MNAR                                               MAR") +
  scale_shape_manual(values = c(
    "Least.REG" = 16, "Moderate.REG" = 16, "Most.REG" = 16,  
    "Least.PMM" = 17, "Moderate.PMM" = 17, "Most.PMM" = 17
  )) +
  scale_linetype_manual(values = c(
    "w/o_Transf" = "solid", 
    "w/_Transf" = "dashed"
  )) +
  palette +  
  theme_minimal() +
  facet_wrap(
    vars(miss_mech, metric), 
    nrow = 2, 
    scales = "free_y",
    labeller = labeller(miss_mech = c(MAR = "", MNAR = ""))  # Hide facet titles
  )+ 
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )+ 
  ggtitle("Hypertension") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  # Apply y-axis scale only if the metric is Pseudo_R2
  facetted_pos_scales(
    y = list(
      NULL, NULL, NULL, 
      scale_y_continuous(breaks = seq(0.20, 0.23, 0.01), limits = c(0.20, 0.23)),  # 4th facet
      #scale_y_continuous(breaks = seq(0.15, 0.19, 0.01), limits = c(0.15, 0.19)),  
      #scale_y_continuous(breaks = seq(0.12, 0.16, 0.01), limits = c(0.12, 0.16)),
      NULL, NULL, NULL, 
      scale_y_continuous(breaks = seq(0.20, 0.23, 0.01), limits = c(0.20, 0.23))   # 8th facet
      #scale_y_continuous(breaks = seq(0.15, 0.19, 0.01), limits = c(0.15, 0.19))   
      #scale_y_continuous(breaks = seq(0.12, 0.16, 0.01), limits = c(0.12, 0.16))
    )
  )
p




# Part 3 : coverage
plot_cr <- function(data, title_text, hline) {
  palette <- scale_color_manual(values = c("darkblue", "darkblue", "darkblue","lightblue3", "lightblue3", "lightblue3"))
  
  # Create a new x_label column with a newline between miss_rate and shared_var
  data$x_label <- str_replace_all(interaction(data$miss_rate, data$shared_var), "\\.", "\n")
  
  # Desired legend order
  desired_legend_order <- c(
    "MAR.REG.w/_Transf",
    "MAR.REG.w/o_Transf",
    "MAR.PMM.w/o_Transf",
    "MNAR.REG.w/_Transf",
    "MNAR.REG.w/o_Transf",
    "MNAR.PMM.w/o_Transf"
  )
  
  ggplot(data, aes(
    y = coverage, 
    x = fct_reorder(data$x_label, data$miss_rate),
    shape = interaction(miss_mech, method, model_spec), 
    color = interaction(miss_mech, method, model_spec)
  )) +
    geom_point(size = 4) +  
    geom_hline(yintercept = hline, linetype = "dashed", color = "orange", linewidth = 1) +  
    labs(
      x = "Missing Rate and Model",
      y = "CR", 
      title = title_text
    ) +  
    scale_color_manual(
      values = c("darkblue", "darkblue", "darkblue","lightblue3", "lightblue3", "lightblue3"),
      breaks = desired_legend_order
    ) +
    scale_shape_manual(
      values = c(1, 16, 17, 1, 16, 17),
      breaks = desired_legend_order
    ) +
    theme_minimal() + 
    theme(legend.position = "none", 
          axis.title.x = element_text(size = 12),  
          axis.title.y = element_text(size = 12),  
          axis.text.x = element_text(size = 10),   
          axis.text.y = element_text(size = 12),   
          strip.text = element_text(size = 12)) +  
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),  
      limits = c(0, 1)  
    ) +
    theme(
      legend.position = "none",  
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.45, size = 12, face = "bold"),  
      axis.line = element_line(color = "black", linewidth = 0.1)  
    )
}

# List of datasets and titles
datasets <- datasets
titles <- titles
hlines <- c(0.96, 0.94, 0.94, 0.95, 0.95)

# Generate plots
p5_plots <- lapply(1:length(datasets), function(i) plot_cr(datasets[[i]], titles[i], hlines[i]))

# Extract legend from the first plot
legend_plot <- get_legend(
  plot_cr(datasets[[1]], "", hlines[1]) + 
    theme(legend.position = "right", legend.title = element_blank(),
          legend.text = element_text(size = 12))
)

# Final combined plot
final_plot <- plot_grid(
  plot_grid(plotlist = p5_plots, legend_plot, ncol = 2)
)

final_plot
ggsave(
  filename = "coverage.png",  # file extension determines format
  plot = final_plot,                 # which plot to save
  dpi = 300,               # resolution in dots per inch
  width = 11, height = 10,    # physical size in inches
  units = "in"              # units for width/height
)


work_dir <- "/Users/dsuolang/Desktop/study2/"
setwd(work_dir)
synthpop<-read_dta("synthpop.dta")
mean_and_se <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  mu <- mean(x)
  
  se <- if (length(unique(x)) == 2 && all(unique(x) %in% c(0, 1))) {
    sqrt(mu * (1 - mu) / n)         # proportion case
  } else {
    sd(x) / sqrt(n)                 # continuous case
  }
  
  tibble(mean = mu, se = se)
}
synthpop$mvpa_total_acc<-synthpop$mvpa_total_acc_sqrt^2
summary(synthpop$mvpa_total_acc)
# Define variables
continuous_vars <- c("age", "bmi",  "mvpa_total_acc", "modpa_total", "vigpa_total")
categorical_vars <- c("modpa_indicator", "vigpa_indicator", "activity_pattern", "source", "srvy_yr",
                      "race", "gender", "marital", "edu", "poverty", "work", "insurance",
                      "self_reported_health", "smoker", "alcohol_cat", "hypertension", "diabetes", 
                      "heartdiseases", "cancers", "stroke", "fitness_access", "health_literacy"
)

# Function to summarize continuous variable
summarize_continuous <- function(data, var) {
  mean_val <- mean(data[[var]], na.rm = TRUE)
  se_val <- sd(data[[var]], na.rm = TRUE) / sqrt(sum(!is.na(data[[var]])))
  tibble(
    Variable = var,
    `Mean (SE)` = sprintf("%.2f (%.2f)", mean_val, se_val)
  )
}

# Function to summarize categorical variable
summarize_categorical <- function(data, var) {
  tab <- table(data[[var]])
  prop <- prop.table(tab)
  se <- sqrt(prop * (1 - prop) / sum(tab))
  
  tibble(
    Variable = paste0(var, ":", names(prop)),
    `Mean (SE)` = sprintf("%.2f (%.2f)", prop, se)
  )
}

# Apply summary functions
nhanes<-synthpop[synthpop$source==0,]
nhis<-synthpop[synthpop$source==1,]
cont_summary <- map_dfr(continuous_vars, ~summarize_continuous(nhis, .x))
cat_summary <- map_dfr(categorical_vars, ~summarize_categorical(nhis, .x))

# Combine summaries
final_summary <- bind_rows(cont_summary, cat_summary)

# View the result
print(final_summary)





library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# Assuming your data frame is called df

df_clean <- srmi_pmm %>%
  mutate(
    miss_rate = str_extract(scenario_avg_run, "(?<=miss_)\\d+"),  # extract digits after 'miss_'
    miss_mech = ifelse(str_detect(scenario_avg_run, "mar"), "MAR", "MNAR"),
    shared_var = case_when(
      str_detect(scenario_avg_run, "levelspec_1") ~ "least",
      str_detect(scenario_avg_run, "levelspec_2") ~ "moderate",
      str_detect(scenario_avg_run, "levelspec_3") ~ "most",
      TRUE ~ NA_character_
    ),
    model_Spec = ifelse(str_detect(scenario_avg_run, "REG"), "REG", "Other")
  ) %>%
  select(miss_rate, miss_mech, shared_var, model_Spec, accuracy, method)

# Convert miss_rate to numeric for plotting
df_clean$miss_rate <- as.numeric(df_clean$miss_rate)*25
df_clean<-df_clean[!is.na(df_clean$shared_var),]
# Now plot with ggplot2, grouped bars by shared_var, faceted by miss_mech and model_Spec
ggplot(df_clean, aes(x = factor(miss_rate), y = accuracy, fill = shared_var)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(miss_mech ~ model_Spec, scales = "fixed") +
  labs(
    title = "Accuracy Scores by Missing Data Rate, Mechanism, and Model Specification",
    x = "Missing Rate (%)",
    y = "Accuracy Score",
    fill = "Shared Variables"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 10, face = "bold")
  )

