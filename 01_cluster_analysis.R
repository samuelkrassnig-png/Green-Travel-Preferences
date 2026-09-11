## ============================================================
## CLUSTER ANALYSIS — Cross-City Study
## Cities: Vienna | Madrid | Tokyo
##
## Description:
##   Hierarchical cluster analysis (Ward.D2) on aggregated
##   psychological constructs of the extended TPB framework.
##
##   Constructs used for clustering:
##     ATT = Attitude
##     SN  = Subjective Norm
##     PBC = Perceived Behavioral Control
##     INT = Intention
##     HAB = Habit
##     EAP = Early Adopter Profile
##
##   Notes:
##     - PBC is computed from PBC1, PBC2, PBC3. This differs from
##       the SEM specification (where PBC1 is excluded). The
##       decision is intentional and documented in the thesis.
##     - k is chosen based on the hierarchical Ward solution
##       (cutree-based silhouette) rather than k-means.
##     - Nominal demographics (employment, residence) are
##       reported as proportions, not means.
##
## Usage:
##   1. Set the working directory to the repository root.
##   2. Run: source("scripts/01_cluster_analysis.R")
##
## Outputs (per city):
##   outputs/cluster/<city>/*.csv
##   outputs/cluster/<city>/plots/*.png
##   outputs/cluster/<city>/THESIS_SUMMARY_Cluster.txt
## ============================================================

rm(list = ls())
cat("\n=== CLUSTER ANALYSIS (VIENNA | MADRID | TOKYO) ===\n")

## ============================================================
## 1. PACKAGES (only load — do not install)
## ============================================================

required_packages <- c("dplyr", "ggplot2", "cluster", "tidyr",
                       "factoextra", "RColorBrewer")

missing_packages <- required_packages[!required_packages %in% installed.packages()]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       "\nPlease install them before running the analysis.")
}
invisible(lapply(required_packages, library, character.only = TRUE))

## ============================================================
## 2. FILE PATHS
## ============================================================

cities <- list(
  Vienna = list(csv = "data/processed/vienna.csv",
                out = "outputs/cluster/vienna"),
  Madrid = list(csv = "data/processed/madrid.csv",
                out = "outputs/cluster/madrid"),
  Tokyo  = list(csv = "data/processed/tokyo.csv",
                out = "outputs/cluster/tokyo")
)

## Required item variables (must exist — otherwise stop)
required_items <- c(
  "ATT1","ATT2","ATT3",
  "SN1","SN2","SN3",
  "PBC1","PBC2","PBC3",
  "INT1","INT2","INT3","INT4",
  "HAB1","HAB2","HAB3_1","HAB3_2",
  paste0("EAP", 1:11)
)

## Readable labels for plots
construct_labels <- c(
  "ATT" = "Attitude",
  "SN"  = "Subjective Norm",
  "PBC" = "Perceived Behavioral Control",
  "INT" = "Intention",
  "HAB" = "Habit",
  "EAP" = "Early Adopter Profile"
)

## ============================================================
## 3. HELPER FUNCTIONS
## ============================================================

clean_text <- function(x) {
  if (!is.character(x) && !is.factor(x)) return(x)
  trimws(gsub("[^[:print:]]", "", as.character(x)))
}

safe_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

convert_gender <- function(x) {
  if (is.numeric(x)) return(x)
  x <- clean_text(x); out <- rep(NA, length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i]) || x[i] == "" || x[i] == "NA") next
    v <- tolower(x[i])
    if (v %in% c("female","f","weiblich","frau"))           out[i] <- 1
    else if (v %in% c("male","m","männlich","mann"))         out[i] <- 2
    else if (v %in% c("non-binary","nonbinary","diverse","nb")) out[i] <- 3
  }
  out
}

convert_age <- function(x) {
  if (is.numeric(x)) return(x)
  x <- clean_text(x); out <- rep(NA, length(x))
  map <- list("below 18" = 16, "18-25" = 21.5, "26-30" = 28,
              "31-40" = 35.5, "41-50" = 45.5,
              "51-60" = 55.5, "above 60" = 65)
  for (i in seq_along(x)) {
    if (is.na(x[i]) || x[i] == "" || x[i] == "NA") next
    v <- tolower(x[i])
    if (grepl("no answer|keine angabe|prefer not", v)) next
    for (p in names(map)) if (grepl(p, v)) { out[i] <- map[[p]]; break }
  }
  out
}

convert_education <- function(x) {
  if (is.numeric(x)) return(x)
  x <- clean_text(x); out <- rep(NA, length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i]) || x[i] == "" || x[i] == "NA") next
    v <- tolower(x[i])
    if (grepl("primary|lower secondary|grundschule|hauptschule", v)) out[i] <- 1
    else if (grepl("upper secondary|gymnasium|mittelschule|abitur|matura", v)) out[i] <- 2
    else if (grepl("bachelor", v)) out[i] <- 3
    else if (grepl("master|doctoral|phd|doktor|magister|diplom", v)) out[i] <- 4
  }
  out
}

convert_employment <- function(x) {
  if (is.numeric(x)) return(x)
  x <- clean_text(x); out <- rep(NA, length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i]) || x[i] == "" || x[i] == "NA") next
    v <- tolower(x[i])
    if (grepl("student", v)) out[i] <- 1
    else if (grepl("full-time", v)) out[i] <- 2
    else if (grepl("part-time", v)) out[i] <- 3
    else if (grepl("unemployed", v)) out[i] <- 4
    else if (grepl("retired", v)) out[i] <- 5
  }
  out
}

convert_residence <- function(x) {
  if (is.numeric(x)) return(x)
  out <- safe_numeric(clean_text(x))
  out[out < 1 | out > 4] <- NA
  out
}

safe_rowMeans <- function(df, vars) {
  existing <- vars[vars %in% names(df)]
  if (length(existing) == 0) return(rep(NA, nrow(df)))
  r <- rowMeans(df[, existing, drop = FALSE], na.rm = TRUE)
  r[is.nan(r)] <- NA
  r
}

## ============================================================
## 4. MAIN LOOP OVER CITIES
## ============================================================

for (city_name in names(cities)) {
  
  cfg <- cities[[city_name]]
  cat("\n", strrep("=", 65), "\n")
  cat("  CITY:", city_name, "\n")
  cat(strrep("=", 65), "\n")
  
  if (!file.exists(cfg$csv)) stop("File not found: ", cfg$csv)
  if (!dir.exists(cfg$out)) dir.create(cfg$out, recursive = TRUE)
  plot_dir <- file.path(cfg$out, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  ## ---- 4.1 Load data ----
  first_line <- readLines(cfg$csv, n = 1, warn = FALSE)
  data <- if (grepl(";", first_line)) {
    read.csv2(cfg$csv, stringsAsFactors = FALSE,
              na.strings = c("", "NA", "NULL", "N/A"))
  } else {
    read.csv(cfg$csv, stringsAsFactors = FALSE,
             na.strings = c("", "NA", "NULL", "N/A"))
  }
  
  cat("✓ Loaded:", nrow(data), "rows,", ncol(data), "columns\n")
  
  ## ---- 4.2 Required items check ----
  missing_items <- setdiff(required_items, names(data))
  if (length(missing_items) > 0) {
    stop("Missing required item variables for ", city_name, ": ",
         paste(missing_items, collapse = ", "))
  }
  
  ## ---- 4.3 Demographics ----
  if ("Gender" %in% names(data)) {
    data$Gender_numeric  <- convert_gender(data$Gender)
    data$Gender_original <- as.character(data$Gender)
  }
  if ("Age" %in% names(data)) {
    data$Age_numeric  <- convert_age(data$Age)
    data$Age_original <- as.character(data$Age)
  }
  for (col in c("Education.Level","Education Level","Education","Bildung")) {
    if (col %in% names(data)) {
      data$Education_numeric  <- convert_education(data[[col]])
      data$Education_original <- as.character(data[[col]])
      break
    }
  }
  for (col in c("Employment.Status","Employment Status","Employment","Job")) {
    if (col %in% names(data)) {
      data$Employment_numeric  <- convert_employment(data[[col]])
      data$Employment_original <- as.character(data[[col]])
      break
    }
  }
  if ("CC1" %in% names(data)) {
    data$CC1_numeric  <- convert_residence(data$CC1)
    data$CC1_original <- data$CC1
  }
  
  ## ---- 4.4 Items numeric ----
  for (item in required_items) data[[item]] <- safe_numeric(data[[item]])
  
  ## ---- 4.5 Construct aggregation ----
  ## Note: PBC uses PBC1, PBC2, PBC3 (differs from SEM where PBC1 excluded).
  data$ATT <- safe_rowMeans(data, c("ATT1","ATT2","ATT3"))
  data$SN  <- safe_rowMeans(data, c("SN1","SN2","SN3"))
  data$PBC <- safe_rowMeans(data, c("PBC1","PBC2","PBC3"))
  data$INT <- safe_rowMeans(data, c("INT1","INT2","INT3","INT4"))
  data$HAB <- safe_rowMeans(data, c("HAB1","HAB2","HAB3_1","HAB3_2"))
  data$EAP <- safe_rowMeans(data, paste0("EAP", 1:11))
  
  constructs <- c("ATT","SN","PBC","INT","HAB","EAP")
  
  ## ---- 4.6 Complete cases ----
  cluster_vars <- data[, constructs, drop = FALSE]
  cc_idx       <- complete.cases(cluster_vars)
  cluster_data <- cluster_vars[cc_idx, , drop = FALSE]
  data_complete<- data[cc_idx, , drop = FALSE]
  
  cat("✓ Cases for clustering:", nrow(cluster_data),
      "(", round(100 * nrow(cluster_data)/nrow(data), 1), "%)\n")
  
  if (nrow(cluster_data) < 10) stop("Too few complete cases for ", city_name)
  
  ## ---- 4.7 Standardization ----
  CLUSTER.DAT <- as.data.frame(scale(cluster_data))
  d.eucl      <- dist(CLUSTER.DAT, method = "euclidean")
  
  ## ---- 4.8 Hierarchical Ward solution ----
  hc_ward <- hclust(d.eucl, method = "ward.D2")
  
  ## ---- 4.9 Cluster number: based on the Ward solution ----
  max_k <- min(10, nrow(CLUSTER.DAT) - 1)
  
  ## Silhouette computed on Ward cutree solutions
  ward_sil <- sapply(2:max_k, function(kk) {
    cl <- cutree(hc_ward, k = kk)
    mean(silhouette(cl, d.eucl)[, 3])
  })
  
  sil_plot <- ggplot(data.frame(k = 2:max_k, silhouette = ward_sil),
                     aes(k, silhouette)) +
    geom_line(linewidth = 1.2, color = "darkgreen") +
    geom_point(size = 3, color = "darkgreen") +
    labs(title = "Average Silhouette Width — Ward Solution",
         subtitle = paste(city_name, "— silhouette computed on Ward cutree"),
         x = "Number of Clusters (k)", y = "Average Silhouette Width") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  print(sil_plot)
  ggsave(file.path(plot_dir, "02_silhouette_method.png"), sil_plot,
         width = 8, height = 5, dpi = 300)
  
  ## Elbow computed on Ward heights (within-cluster SS of the tree)
  ## Note: elbow and silhouette are both computed on the Ward tree,
  ##       NOT on k-means.
  height_vec <- rev(hc_ward$height)
  wss_ward   <- cumsum(height_vec)[1:max_k]
  
  elbow_plot <- ggplot(data.frame(k = 1:max_k, wss = wss_ward),
                       aes(k, wss)) +
    geom_line(linewidth = 1.2, color = "steelblue") +
    geom_point(size = 3, color = "steelblue") +
    labs(title = "Elbow Plot — Ward Solution",
         subtitle = paste(city_name, "— cumulative Ward heights"),
         x = "Number of Clusters (k)", y = "Cumulative Height") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  print(elbow_plot)
  ggsave(file.path(plot_dir, "01_elbow_method.png"), elbow_plot,
         width = 8, height = 5, dpi = 300)
  
  ## ---- 4.10 Choose k ----
  ## Set k explicitly based on the silhouette plot above.
  ## Example: k <- which.max(ward_sil) + 1
  k <- 3
  
  cat("✓ Chosen k =", k, "\n")
  
  cluster_ward <- as.factor(cutree(hc_ward, k = k))
  data_complete$cluster_ward <- cluster_ward
  
  ## ---- 4.11 Dendrogram ----
  dendro_plot <- fviz_dend(
    hc_ward, k = k, cex = 0.3, rect = TRUE, rect_fill = TRUE,
    main = paste("Dendrogram — Ward's Method (k =", k, ") —", city_name),
    xlab = "Cases", ylab = "Distance",
    palette = "Set2", ggtheme = theme_minimal(base_size = 10)
  ) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(plot_dir, "03_dendrogram.png"), dendro_plot,
         width = 14, height = 8, dpi = 300)
  
  ## ---- 4.12 Cluster profiles (z-scores) ----
  cluster_profiles_z <- data.frame(cluster_ward = data_complete$cluster_ward,
                                   CLUSTER.DAT) %>%
    group_by(cluster_ward) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
  
  profile_long <- cluster_profiles_z %>%
    pivot_longer(cols = -cluster_ward,
                 names_to = "Construct", values_to = "Z_score")
  profile_long$Construct_label <- construct_labels[profile_long$Construct]
  
  profile_plot <- ggplot(profile_long,
                         aes(x = Construct_label, y = Z_score,
                             color = factor(cluster_ward),
                             group = cluster_ward)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    labs(title    = paste("Cluster Profiles —", city_name),
         subtitle = paste("Standardized means (z-scores), k =", k),
         x = "Psychological Construct", y = "Mean Z-score", color = "Cluster") +
    scale_color_brewer(palette = "Set1") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x   = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "bottom") +
    geom_hline(yintercept = 0, linetype = "dashed",
               alpha = 0.5, color = "gray50") +
    ylim(-1.5, 1.5)
  ggsave(file.path(plot_dir, "04_cluster_profiles_line.png"), profile_plot,
         width = 10, height = 6, dpi = 300)
  
  ## ---- 4.13 Boxplots ----
  ward_long <- cbind(CLUSTER.DAT, Cluster = data_complete$cluster_ward) %>%
    pivot_longer(cols = -Cluster, names_to = "Construct", values_to = "Z_score")
  ward_long$Construct_label <- construct_labels[ward_long$Construct]
  
  boxplot_constructs <- ggplot(ward_long,
                               aes(x = Construct_label, y = Z_score,
                                   fill = Construct_label)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.5) +
    facet_wrap(~ Cluster, ncol = 2,
               labeller = labeller(Cluster = function(x) paste("Cluster", x))) +
    labs(title = paste("Construct Distribution Across Clusters —", city_name),
         subtitle = "Standardized values (z-scores)", x = "", y = "Z-score") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x   = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          strip.text      = element_text(face = "bold", size = 11)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
  ggsave(file.path(plot_dir, "05_boxplots_constructs.png"), boxplot_constructs,
         width = 12, height = 8, dpi = 300)
  
  ## ---- 4.14 Silhouette of chosen k ----
  sil <- silhouette(as.numeric(as.character(data_complete$cluster_ward)), d.eucl)
  avg_sil_width <- mean(sil[, 3])
  
  png(file.path(plot_dir, "06_silhouette_plot.png"),
      width = 10, height = 6, units = "in", res = 300)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
  plot(sil,
       main = paste("Silhouette Plot (k =", k, ") —", city_name),
       col = 2:(k + 1), border = NA, cex.names = 0.6)
  text(x = 0, y = par("usr")[4] - 0.5,
       labels = paste("Average silhouette width =", round(avg_sil_width, 3)),
       adj = 0, font = 2, cex = 1.1)
  dev.off()
  
  ## ---- 4.15 Demographic profiling ----
  ## Age (numeric — boxplot)
  if ("Age_numeric" %in% names(data_complete) &&
      sum(!is.na(data_complete$Age_numeric)) > 0) {
    age_boxplot <- ggplot(data_complete %>% filter(!is.na(Age_numeric)),
                          aes(x = factor(cluster_ward),
                              y = Age_numeric, fill = factor(cluster_ward))) +
      geom_boxplot(alpha = 0.7, outlier.size = 1) +
      stat_summary(fun = mean, geom = "point",
                   shape = 18, size = 3, color = "darkred") +
      labs(title = paste("Age Distribution —", city_name),
           x = "Cluster", y = "Age (years)", fill = "Cluster") +
      scale_fill_brewer(palette = "Set1") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "none")
    ggsave(file.path(plot_dir, "07_age_by_cluster.png"), age_boxplot,
           width = 6, height = 5, dpi = 300)
  }
  
  ## Gender (nominal — stacked bar)
  if ("Gender_original" %in% names(data_complete)) {
    gender_counts <- data_complete %>%
      filter(!is.na(Gender_original), Gender_original != "") %>%
      mutate(Gender_clean = case_when(
        grepl("female|weiblich",   tolower(Gender_original)) ~ "Female",
        grepl("male|männlich",     tolower(Gender_original)) ~ "Male",
        grepl("non-binary|diverse",tolower(Gender_original)) ~ "Non-binary",
        TRUE ~ "Other")) %>%
      group_by(cluster_ward, Gender_clean) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(cluster_ward) %>%
      mutate(percentage = count / sum(count) * 100)
    
    if (nrow(gender_counts) > 0) {
      gender_barplot <- ggplot(gender_counts,
                               aes(x = factor(cluster_ward), y = percentage,
                                   fill = Gender_clean)) +
        geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
        geom_text(aes(label = paste0(round(percentage, 1), "%")),
                  position = position_stack(vjust = 0.5), size = 3.5) +
        labs(title = paste("Gender Distribution —", city_name),
             x = "Cluster", y = "Percentage (%)", fill = "Gender") +
        scale_fill_brewer(palette = "Set2") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              legend.position = "bottom")
      ggsave(file.path(plot_dir, "08_gender_by_cluster.png"), gender_barplot,
             width = 6, height = 5, dpi = 300)
    }
  }
  
  ## Residence (nominal — stacked bar)
  if ("CC1_original" %in% names(data_complete)) {
    res_labels <- c("1"="City Centre","2"="Inner Urban",
                    "3"="Suburban","4"="Rural")
    res_counts <- data_complete %>%
      filter(!is.na(CC1_original), CC1_original %in% c(1,2,3,4)) %>%
      mutate(Residence = factor(CC1_original, levels = 1:4, labels = res_labels)) %>%
      group_by(cluster_ward, Residence) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(cluster_ward) %>%
      mutate(percentage = count / sum(count) * 100)
    
    if (nrow(res_counts) > 0) {
      res_barplot <- ggplot(res_counts,
                            aes(x = factor(cluster_ward), y = percentage,
                                fill = Residence)) +
        geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
        geom_text(aes(label = paste0(round(percentage, 1), "%")),
                  position = position_stack(vjust = 0.5), size = 3) +
        labs(title = paste("Residence Type —", city_name),
             x = "Cluster", y = "Percentage (%)", fill = "Residence Type") +
        scale_fill_brewer(palette = "Set3") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              legend.position = "bottom")
      ggsave(file.path(plot_dir, "09_residence_by_cluster.png"), res_barplot,
             width = 6, height = 5, dpi = 300)
    }
  }
  
  ## Education (nominal — stacked bar)
  if ("Education_original" %in% names(data_complete)) {
    edu_counts <- data_complete %>%
      filter(!is.na(Education_original), Education_original != "") %>%
      mutate(Education_short = case_when(
        grepl("primary|lower secondary", tolower(Education_original)) ~ "Primary/Lower Sec.",
        grepl("upper secondary",         tolower(Education_original)) ~ "Upper Secondary",
        grepl("bachelor",                tolower(Education_original)) ~ "Bachelor's",
        grepl("master|doctoral",         tolower(Education_original)) ~ "Master's/Doctoral",
        TRUE ~ "Other")) %>%
      group_by(cluster_ward, Education_short) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(cluster_ward) %>%
      mutate(percentage = count / sum(count) * 100)
    
    if (nrow(edu_counts) > 0) {
      edu_barplot <- ggplot(edu_counts,
                            aes(x = factor(cluster_ward), y = percentage,
                                fill = Education_short)) +
        geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
        labs(title = paste("Education Level —", city_name),
             x = "Cluster", y = "Percentage (%)", fill = "Education") +
        scale_fill_brewer(palette = "Blues") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              legend.position = "bottom",
              legend.text = element_text(size = 8))
      ggsave(file.path(plot_dir, "10_education_by_cluster.png"), edu_barplot,
             width = 7, height = 5, dpi = 300)
    }
  }
  
  ## ---- 4.16 Profile tables ----
  profile_table <- cluster_profiles_z %>%
    mutate(cluster_ward = as.numeric(as.character(cluster_ward))) %>%
    arrange(cluster_ward)
  
  raw_means_table <- data_complete %>%
    group_by(cluster_ward) %>%
    summarise(n   = n(),
              ATT = mean(ATT, na.rm = TRUE),
              SN  = mean(SN,  na.rm = TRUE),
              PBC = mean(PBC, na.rm = TRUE),
              INT = mean(INT, na.rm = TRUE),
              HAB = mean(HAB, na.rm = TRUE),
              EAP = mean(EAP, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(cluster_ward = as.numeric(as.character(cluster_ward))) %>%
    arrange(cluster_ward) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
  
  ## Nominal demographics: proportions (not means)
  ## Build per-cluster proportions for gender, residence, employment.
  prop_table <- function(df, var, labels) {
    if (!var %in% names(df)) return(NULL)
    df %>%
      filter(!is.na(.data[[var]])) %>%
      group_by(cluster_ward, category = .data[[var]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(cluster_ward) %>%
      mutate(proportion = round(100 * count / sum(count), 1)) %>%
      mutate(variable = var) %>%
      rename(level = category)
  }
  
  props_gender <- prop_table(data_complete, "Gender_numeric", NULL)
  props_resid  <- prop_table(data_complete, "CC1_numeric",    NULL)
  props_employ <- prop_table(data_complete, "Employment_numeric", NULL)
  
  nominal_props <- bind_rows(props_gender, props_resid, props_employ)
  
  ## Numeric demographic summary (age, education years)
  numeric_demog <- data_complete %>%
    group_by(cluster_ward) %>%
    summarise(
      n          = n(),
      Age_mean   = mean(Age_numeric, na.rm = TRUE),
      Age_sd     = sd(Age_numeric,   na.rm = TRUE),
      Education_mean = mean(Education_numeric, na.rm = TRUE),
      Education_sd   = sd(Education_numeric,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(cluster_ward = as.numeric(as.character(cluster_ward))) %>%
    arrange(cluster_ward)
  
  ## ---- 4.17 Export ----
  write.csv(data_complete,
            file.path(cfg$out, "clustered_data.csv"), row.names = FALSE)
  write.csv(profile_table,
            file.path(cfg$out, "cluster_profiles_standardized.csv"),
            row.names = FALSE)
  write.csv(raw_means_table,
            file.path(cfg$out, "cluster_profiles_raw_means.csv"),
            row.names = FALSE)
  write.csv(numeric_demog,
            file.path(cfg$out, "demographic_summary_numeric.csv"),
            row.names = FALSE)
  write.csv(nominal_props,
            file.path(cfg$out, "demographic_summary_nominal.csv"),
            row.names = FALSE)
  
  ## ---- 4.18 Thesis summary ----
  sink(file.path(cfg$out, "THESIS_SUMMARY_Cluster.txt"))
  cat(strrep("=", 75), "\n")
  cat("  CLUSTER ANALYSIS —", toupper(city_name), "\n")
  cat(strrep("=", 75), "\n\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
  cat("Number of clusters (k):", k, "\n")
  cat("Cases analyzed:", nrow(data_complete), "\n")
  cat("Average silhouette width:", round(avg_sil_width, 3), "\n\n")
  
  cat(strrep("-", 75), "\n")
  cat("  METHOD NOTES\n")
  cat(strrep("-", 75), "\n\n")
  cat("  - Hierarchical clustering (Ward.D2) on standardized construct scores.\n")
  cat("  - Cluster number chosen based on the silhouette of the Ward solution.\n")
  cat("  - PBC is computed from PBC1, PBC2, PBC3 (differs from SEM where PBC1 excluded).\n")
  cat("  - Nominal demographics (gender, residence, employment) reported as proportions.\n")
  cat("  - EAP is interpreted as 'Early Adopter Profile'.\n\n")
  
  cat(strrep("-", 75), "\n")
  cat("  CLUSTER SIZES\n")
  cat(strrep("-", 75), "\n\n")
  print(table(data_complete$cluster_ward))
  cat("\n")
  
  cat(strrep("-", 75), "\n")
  cat("  CONSTRUCT PROFILES (RAW MEANS)\n")
  cat(strrep("-", 75), "\n\n")
  print(raw_means_table); cat("\n")
  
  cat(strrep("-", 75), "\n")
  cat("  CONSTRUCT PROFILES (STANDARDIZED)\n")
  cat(strrep("-", 75), "\n\n")
  print(profile_table %>% mutate(across(where(is.numeric), ~ round(.x, 2))))
  cat("\n")
  
  cat(strrep("-", 75), "\n")
  cat("  NUMERIC DEMOGRAPHICS\n")
  cat(strrep("-", 75), "\n\n")
  print(numeric_demog); cat("\n")
  
  cat(strrep("-", 75), "\n")
  cat("  NOMINAL DEMOGRAPHICS (PROPORTIONS %)\n")
  cat(strrep("-", 75), "\n\n")
  print(nominal_props); cat("\n")
  
  sink()
  
  cat("\n✓", city_name, "done — outputs in", cfg$out, "\n")
}

cat("\n=== CLUSTER ANALYSIS COMPLETE ===\n")