## ============================================================
## STRUCTURAL EQUATION MODELING (SEM)
## MODEL A — Standard Theory of Planned Behavior (TPB)
## Sustainable Transport Choices — Cross-City Study
##
## Description:
##   Full TPB model with four latent constructs:
##     ATT = Attitude toward sustainable transport
##     SN  = Subjective Norm
##     PBC = Perceived Behavioral Control
##     INT = Behavioral Intention
##
##   PBC is measured with two indicators (PBC2, PBC3).
##   PBC1 is excluded from the measurement model.
##
##   Model A is estimated separately for Vienna and Madrid using
##   an identical specification, so that results are directly
##   comparable across cities.
##
## Hypotheses tested:
##   H1:  ATT → INT  (Direct)        ATT has a positive effect on INT.
##   H2:  SN  → INT  (Direct)        SN has a positive effect on INT.
##   H3:  PBC → INT  (Direct)        PBC has a positive effect on INT.
##   CH1: ATT ↔ SN   (Correlation)   ATT and SN are positively associated.
##   CH2: SN  ↔ PBC  (Correlation)   SN and PBC are positively associated.
##   CH3: ATT ↔ PBC  (Correlation)   ATT and PBC are positively associated.
##
## Admissibility rule:
##   admissible <- converged && post_check
##
##   If admissible == TRUE:
##     - structural paths (H1–H3), correlations (CH1–CH3),
##       R² for INT, and the path diagram are exported.
##
##   If admissible == FALSE:
##     - only fit indices and measurement diagnostics
##       (factor loadings, reliability, CR, AVE, HTMT) are exported.
##     - structural paths, hypotheses, correlations, R² and the
##       path diagram are NOT interpreted and NOT written to disk.
##
## Usage:
##   1. Set the working directory to the repository root.
##   2. Run: source("scripts/02_sem_model_a.R")
##
## Outputs:
##   outputs/model_a/vienna/
##   outputs/model_a/madrid/
## ============================================================

cat("\n=== SEM MODEL A — VIENNA | MADRID ===\n")

## ============================================================
## 1. PACKAGES (only load — do not install)
## ============================================================

packages <- c("lavaan", "semPlot", "tidyverse", "psych",
              "ggplot2", "corrplot")

missing_packages <- packages[!packages %in% installed.packages()]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       "\nPlease install them before running the analysis.")
}
invisible(lapply(packages, library, character.only = TRUE))

## ============================================================
## 2. CITY CONFIGURATION
## ============================================================

cities <- list(
  Vienna = list(csv = "data/processed/vienna.csv",
                out = "outputs/model_a/vienna"),
  Madrid = list(csv = "data/processed/madrid.csv",
                out = "outputs/model_a/madrid")
)

## ============================================================
## 3. VARIABLE DEFINITIONS
## ============================================================

sem_items <- list(
  ATT = c("ATT1", "ATT2", "ATT3"),
  SN  = c("SN1",  "SN2",  "SN3"),
  PBC = c("PBC2", "PBC3"),
  INT = c("INT1", "INT2", "INT3", "INT4")
)

items_all <- unlist(sem_items)

## Named parameters guarantee unambiguous hypothesis assignment,
## independent of lavaan's internal output ordering.
model_A_spec <- '
  # Measurement model
  ATT =~ ATT1 + ATT2 + ATT3
  SN  =~ SN1  + SN2  + SN3
  PBC =~ PBC2 + PBC3
  INT =~ INT1 + INT2 + INT3 + INT4

  # Structural model (named paths)
  INT ~ h1*ATT + h2*SN + h3*PBC

  # Correlations between exogenous constructs (named)
  ATT ~~ ch1*SN
  SN  ~~ ch2*PBC
  ATT ~~ ch3*PBC
'

hyp_defs <- list(
  H1 = list(lhs = "INT", rhs = "ATT", label = "ATT → INT",
            statement = "ATT has a positive effect on INT."),
  H2 = list(lhs = "INT", rhs = "SN",  label = "SN → INT",
            statement = "SN has a positive effect on INT."),
  H3 = list(lhs = "INT", rhs = "PBC", label = "PBC → INT",
            statement = "PBC has a positive effect on INT.")
)

cor_defs <- list(
  CH1 = list(a = "ATT", b = "SN",  label = "ATT ↔ SN",
             statement = "ATT and SN are positively associated."),
  CH2 = list(a = "SN",  b = "PBC", label = "SN ↔ PBC",
             statement = "SN and PBC are positively associated."),
  CH3 = list(a = "ATT", b = "PBC", label = "ATT ↔ PBC",
             statement = "ATT and PBC are positively associated.")
)

## ============================================================
## 4. HELPER FUNCTIONS
## ============================================================

safe_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

load_data <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  first_line <- readLines(path, n = 1, warn = FALSE)
  if (grepl(";", first_line)) {
    read.csv2(path, stringsAsFactors = FALSE,
              na.strings = c("", "NA", "NULL", "N/A"),
              fileEncoding = "UTF-8")
  } else {
    read.csv(path, stringsAsFactors = FALSE,
             na.strings = c("", "NA", "NULL", "N/A"),
             fileEncoding = "UTF-8")
  }
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  ifelse(p < 0.001, "< .001", sprintf("%.4f", p))
}

calculate_cr_ave <- function(loadings) {
  sl  <- sum(loadings)
  cr  <- sl^2 / (sl^2 + sum(1 - loadings^2))
  ave <- mean(loadings^2)
  list(CR = round(cr, 3), AVE = round(ave, 3))
}

## HTMT using absolute correlations and the geometric mean of the
## monotrait block. Matches common HTMT definitions.
calculate_htmt <- function(df, items1, items2) {
  cm  <- abs(cor(df[, c(items1, items2)], use = "pairwise.complete.obs"))
  n1  <- length(items1)
  het <- cm[1:n1, (n1 + 1):ncol(cm)]
  mo1 <- cm[1:n1, 1:n1]
  mo2 <- cm[(n1 + 1):ncol(cm), (n1 + 1):ncol(cm)]
  m_het <- mean(het, na.rm = TRUE)
  m_mo1 <- ifelse(nrow(mo1) > 1, mean(mo1[lower.tri(mo1)], na.rm = TRUE), 1)
  m_mo2 <- ifelse(nrow(mo2) > 1, mean(mo2[lower.tri(mo2)], na.rm = TRUE), 1)
  round(m_het / sqrt(m_mo1 * m_mo2), 3)
}

## ============================================================
## 5. MAIN LOOP OVER CITIES
## ============================================================

for (city_name in names(cities)) {
  
  cfg <- cities[[city_name]]
  cat("\n", strrep("=", 65), "\n")
  cat("  CITY:", city_name, "\n")
  cat(strrep("=", 65), "\n")
  
  if (!file.exists(cfg$csv)) {
    stop("Required input file missing for ", city_name, ": ", cfg$csv)
  }
  if (!dir.exists(cfg$out)) dir.create(cfg$out, recursive = TRUE)
  plot_dir <- file.path(cfg$out, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  ## ---- 5.1 Load data ----
  raw <- load_data(cfg$csv)
  cat("✓ Loaded:", nrow(raw), "rows,", ncol(raw), "columns\n")
  
  ## ---- 5.2 Required variables check ----
  missing_items <- setdiff(items_all, names(raw))
  if (length(missing_items) > 0) {
    stop("Missing required items for ", city_name, ": ",
         paste(missing_items, collapse = ", "))
  }
  
  for (col in items_all) raw[[col]] <- safe_numeric(raw[[col]])
  
  sem_data <- raw[, items_all, drop = FALSE]
  sem_data <- sem_data[complete.cases(sem_data), ]
  cat("✓ Complete cases:", nrow(sem_data), "\n")
  
  if (nrow(sem_data) < 30) {
    stop("Too few complete cases for ", city_name,
         " (N = ", nrow(sem_data), ").")
  }
  
  ## ============================================================
  ## 6. DESCRIPTIVES AND CONSTRUCT CORRELATIONS (always exported)
  ## ============================================================
  
  desc_stats <- data.frame(
    Item = names(sem_data),
    Mean = round(colMeans(sem_data), 2),
    SD   = round(apply(sem_data, 2, sd), 2),
    Min  = apply(sem_data, 2, min),
    Max  = apply(sem_data, 2, max)
  )
  write.csv(desc_stats,
            file.path(cfg$out, "01_descriptives_items.csv"),
            row.names = FALSE)
  
  construct_scores <- data.frame(
    ATT = rowMeans(sem_data[, sem_items$ATT, drop = FALSE]),
    SN  = rowMeans(sem_data[, sem_items$SN,  drop = FALSE]),
    PBC = rowMeans(sem_data[, sem_items$PBC, drop = FALSE]),
    INT = rowMeans(sem_data[, sem_items$INT, drop = FALSE])
  )
  construct_cor <- round(cor(construct_scores), 3)
  write.csv(construct_cor,
            file.path(cfg$out, "02_construct_correlations.csv"))
  
  png(file.path(plot_dir, "01_correlation_matrix.png"),
      width = 7, height = 6, units = "in", res = 300)
  corrplot(construct_cor, method = "color", type = "upper",
           addCoef.col = "black", tl.col = "black",
           title = paste("Construct Correlations —", city_name),
           mar = c(0, 0, 2, 0))
  dev.off()
  
  ## ============================================================
  ## 7. ESTIMATE MODEL A
  ## ============================================================
  
  cat("\n--- ESTIMATING MODEL A (", city_name, ") ---\n")
  
  fit <- sem(model_A_spec, data = sem_data,
             estimator = "MLR", missing = "listwise",
             control = list(iter.max = 1000))
  
  converged  <- lavInspect(fit, "converged")
  post_check <- tryCatch(lavInspect(fit, "post.check"),
                         error = function(e) FALSE)
  cat("  Converged: ", converged, "\n")
  cat("  Post-check:", post_check, "\n")
  
  admissible <- isTRUE(converged) && isTRUE(post_check)
  cat("  Admissible solution:", admissible, "\n")
  
  ## ============================================================
  ## 8. MODEL FIT (always exported)
  ## ============================================================
  
  fit_measures <- fitMeasures(fit, c("chisq", "df", "pvalue",
                                     "cfi", "tli", "rmsea",
                                     "rmsea.ci.lower", "rmsea.ci.upper",
                                     "srmr", "aic", "bic"))
  fit_table <- data.frame(
    Measure = names(fit_measures),
    Value   = round(as.numeric(fit_measures), 3)
  )
  write.csv(fit_table,
            file.path(cfg$out, "03_model_fit_indices.csv"),
            row.names = FALSE)
  
  ## ============================================================
  ## 9. MEASUREMENT DIAGNOSTICS (always exported)
  ## ============================================================
  
  pe_std <- standardizedSolution(fit)
  
  ## ---- 9.1 Factor loadings ----
  loadings <- pe_std[pe_std$op == "=~",
                     c("lhs", "rhs", "est.std", "se", "pvalue")]
  loadings$p_formatted <- format_p(loadings$pvalue)
  loadings$sig <- ifelse(loadings$pvalue < 0.001, "***",
                         ifelse(loadings$pvalue < 0.01, "**",
                                ifelse(loadings$pvalue < 0.05, "*", "")))
  write.csv(loadings,
            file.path(cfg$out, "04_factor_loadings.csv"),
            row.names = FALSE)
  
  ## ---- 9.2 Reliability, CR, AVE ----
  validity_list <- list()
  for (cons in names(sem_items)) {
    items_cons <- sem_items[[cons]]
    items_cons <- items_cons[items_cons %in% names(sem_data)]
    if (length(items_cons) < 2) next
    
    std_loads <- pe_std$est.std[pe_std$op == "=~" & pe_std$lhs == cons]
    std_loads <- std_loads[!is.na(std_loads)]
    
    cr_ave <- calculate_cr_ave(std_loads)
    alpha  <- tryCatch(
      psych::alpha(sem_data[, items_cons],
                   check.keys = FALSE)$total$raw_alpha,
      error = function(e) NA_real_
    )
    validity_list[[cons]] <- data.frame(
      Construct = cons,
      n_items   = length(items_cons),
      Alpha     = round(alpha, 3),
      CR        = cr_ave$CR,
      AVE       = cr_ave$AVE,
      stringsAsFactors = FALSE
    )
  }
  validity_table <- do.call(rbind, validity_list)
  write.csv(validity_table,
            file.path(cfg$out, "05_reliability_cr_ave.csv"),
            row.names = FALSE)
  
  ## ---- 9.3 HTMT ----
  cnames <- names(sem_items)
  htmt_results <- data.frame(
    Pair = character(), HTMT = numeric(),
    Status = character(), stringsAsFactors = FALSE)
  for (i in 1:(length(cnames) - 1)) {
    for (j in (i + 1):length(cnames)) {
      items1 <- sem_items[[cnames[i]]]
      items2 <- sem_items[[cnames[j]]]
      if (length(items1) < 1 || length(items2) < 1) next
      hval <- calculate_htmt(sem_data, items1, items2)
      htmt_results <- rbind(htmt_results, data.frame(
        Pair   = paste(cnames[i], "↔", cnames[j]),
        HTMT   = hval,
        Status = ifelse(hval < 0.85, "✓ Good (<0.85)",
                        ifelse(hval < 0.90, "⚠ Acceptable (<0.90)",
                               "✗ Poor (≥0.90)")),
        stringsAsFactors = FALSE))
    }
  }
  write.csv(htmt_results,
            file.path(cfg$out, "06_htmt.csv"),
            row.names = FALSE)
  
  ## ============================================================
  ## 10. STRUCTURAL RESULTS — ONLY IF ADMISSIBLE
  ## ============================================================
  
  if (admissible) {
    
    cat("\n--- STRUCTURAL RESULTS (admissible) ---\n")
    
    ## ---- 10.1 Direct effects (H1–H3) ----
    hypothesis_results <- data.frame(
      Hypothesis = character(), Type = character(),
      Relationship = character(), Statement = character(),
      Estimate = numeric(), SE = numeric(),
      p_value = numeric(), p_formatted = character(),
      Supported = character(), stringsAsFactors = FALSE)
    
    for (h in names(hyp_defs)) {
      def <- hyp_defs[[h]]
      matched <- pe_std[pe_std$op == "~" &
                          pe_std$lhs == def$lhs &
                          pe_std$rhs == def$rhs, ]
      if (nrow(matched) != 1) {
        stop("Could not uniquely match path for ", h,
             " (", city_name, ").")
      }
      est <- matched$est.std
      se  <- matched$se
      pv  <- matched$pvalue
      supp <- ifelse(!is.na(est) & !is.na(pv) & est > 0 & pv < 0.05,
                     "✓ Yes", "✗ No")
      hypothesis_results <- rbind(hypothesis_results, data.frame(
        Hypothesis = h, Type = "Direct",
        Relationship = def$label, Statement = def$statement,
        Estimate = round(est, 3), SE = round(se, 3),
        p_value = pv, p_formatted = format_p(pv),
        Supported = supp, stringsAsFactors = FALSE))
    }
    
    ## ---- 10.2 Correlations (CH1–CH3) ----
    cor_results <- data.frame(
      Hypothesis = character(), Type = character(),
      Relationship = character(), Statement = character(),
      Estimate = numeric(), SE = numeric(),
      p_value = numeric(), p_formatted = character(),
      Supported = character(), stringsAsFactors = FALSE)
    
    for (h in names(cor_defs)) {
      def <- cor_defs[[h]]
      matched <- pe_std[pe_std$op == "~~" &
                          ((pe_std$lhs == def$a & pe_std$rhs == def$b) |
                             (pe_std$lhs == def$b & pe_std$rhs == def$a)), ]
      if (nrow(matched) != 1) {
        stop("Could not uniquely match correlation for ", h,
             " (", city_name, ").")
      }
      est <- matched$est.std
      se  <- matched$se
      pv  <- matched$pvalue
      supp <- ifelse(!is.na(est) & !is.na(pv) & est > 0 & pv < 0.05,
                     "✓ Yes", "✗ No")
      cor_results <- rbind(cor_results, data.frame(
        Hypothesis = h, Type = "Correlation",
        Relationship = def$label, Statement = def$statement,
        Estimate = round(est, 3), SE = round(se, 3),
        p_value = pv, p_formatted = format_p(pv),
        Supported = supp, stringsAsFactors = FALSE))
    }
    
    ## ---- 10.3 R² ----
    rsq <- inspect(fit, "r2")
    rsq_df <- data.frame(
      Construct        = names(rsq),
      R_squared        = round(as.numeric(rsq), 3),
      Percent_Variance = round(100 * as.numeric(rsq), 1)
    )
    
    ## ---- 10.4 Export hypothesis tables ----
    write.csv(hypothesis_results,
              file.path(cfg$out, "07a_hypotheses_direct.csv"),
              row.names = FALSE)
    write.csv(cor_results,
              file.path(cfg$out, "07b_hypotheses_correlation.csv"),
              row.names = FALSE)
    write.csv(rbind(hypothesis_results, cor_results),
              file.path(cfg$out, "07_hypotheses_all.csv"),
              row.names = FALSE)
    write.csv(rsq_df,
              file.path(cfg$out, "08_rsquared.csv"),
              row.names = FALSE)
    
    ## ---- 10.5 Path diagram ----
    tryCatch({
      png(file.path(plot_dir, "02_path_diagram.png"),
          width = 12, height = 8, units = "in", res = 300)
      semPaths(fit, what = "std", whatLabels = "std",
               layout = "tree2", rotation = 2,
               edge.label.cex = 0.9, style = "lisrel",
               nCharNodes = 8, sizeMan = 5, sizeLat = 8,
               edge.color = "black", fade = FALSE,
               title = TRUE,
               main = paste("Model A —", city_name))
      dev.off()
    }, error = function(e) {
      cat("  ⚠ Path diagram failed:", conditionMessage(e), "\n")
    })
    
    ## ---- 10.6 Structural plots ----
    ggsave(file.path(plot_dir, "03_direct_effects.png"),
           ggplot(hypothesis_results,
                  aes(x = reorder(Relationship, Estimate),
                      y = Estimate,
                      fill = ifelse(p_value < 0.05, "p < 0.05", "n.s."))) +
             geom_bar(stat = "identity", alpha = 0.85) +
             geom_hline(yintercept = 0) +
             labs(title = paste("Direct Effects —", city_name),
                  subtitle = "H1, H2, H3",
                  x = "", y = "Standardized β", fill = "") +
             coord_flip() +
             theme_minimal(base_size = 12) +
             theme(plot.title = element_text(face = "bold", hjust = 0.5),
                   plot.subtitle = element_text(hjust = 0.5),
                   legend.position = "bottom") +
             scale_fill_manual(values = c("p < 0.05" = "steelblue",
                                          "n.s."     = "coral")),
           width = 7, height = 5, dpi = 300)
    
    ggsave(file.path(plot_dir, "04_correlations.png"),
           ggplot(cor_results,
                  aes(x = reorder(Relationship, Estimate),
                      y = Estimate,
                      fill = ifelse(p_value < 0.05, "p < 0.05", "n.s."))) +
             geom_bar(stat = "identity", alpha = 0.85) +
             geom_hline(yintercept = 0) +
             labs(title = paste("Construct Correlations —", city_name),
                  subtitle = "CH1, CH2, CH3",
                  x = "", y = "Standardized r", fill = "") +
             coord_flip() +
             theme_minimal(base_size = 12) +
             theme(plot.title = element_text(face = "bold", hjust = 0.5),
                   plot.subtitle = element_text(hjust = 0.5),
                   legend.position = "bottom") +
             scale_fill_manual(values = c("p < 0.05" = "steelblue",
                                          "n.s."     = "coral")),
           width = 7, height = 5, dpi = 300)
    
    ggsave(file.path(plot_dir, "05_rsquared.png"),
           ggplot(rsq_df,
                  aes(x = reorder(Construct, R_squared),
                      y = R_squared, fill = Construct)) +
             geom_bar(stat = "identity", alpha = 0.85) +
             geom_text(aes(label = paste0(Percent_Variance, "%")),
                       hjust = -0.2, size = 4) +
             labs(title = paste("Variance Explained (R²) —", city_name),
                  x = "", y = "R²") +
             coord_flip() +
             ylim(0, max(rsq_df$R_squared) + 0.1) +
             theme_minimal(base_size = 11) +
             theme(plot.title = element_text(face = "bold", hjust = 0.5),
                   legend.position = "none") +
             scale_fill_brewer(palette = "Set2"),
           width = 7, height = 5, dpi = 300)
    
    ## ---- 10.7 Thesis summary (structural section included) ----
    sink(file.path(cfg$out, "THESIS_SUMMARY_ModelA.txt"))
    cat(strrep("=", 75), "\n")
    cat("  MODEL A —", toupper(city_name), "\n")
    cat(strrep("=", 75), "\n\n")
    cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
    cat("N (complete cases):", nrow(sem_data), "\n")
    cat("Admissible solution:", admissible, "\n\n")
    
    cat(strrep("-", 75), "\n"); cat("  MODEL FIT\n")
    cat(strrep("-", 75), "\n\n")
    print(fit_table); cat("\n")
    
    cat(strrep("-", 75), "\n"); cat("  DIRECT EFFECTS (H1–H3)\n")
    cat(strrep("-", 75), "\n\n")
    for (i in seq_len(nrow(hypothesis_results))) {
      cat(sprintf("  %s: %s\n", hypothesis_results$Hypothesis[i],
                  hypothesis_results$Relationship[i]))
      cat(sprintf("    %s\n", hypothesis_results$Statement[i]))
      cat(sprintf("    β = %6.3f (SE = %5.3f), p = %s\n",
                  hypothesis_results$Estimate[i],
                  hypothesis_results$SE[i],
                  hypothesis_results$p_formatted[i]))
      cat(sprintf("    Supported: %s\n\n",
                  hypothesis_results$Supported[i]))
    }
    
    cat(strrep("-", 75), "\n"); cat("  CORRELATIONS (CH1–CH3)\n")
    cat(strrep("-", 75), "\n\n")
    for (i in seq_len(nrow(cor_results))) {
      cat(sprintf("  %s: %s\n", cor_results$Hypothesis[i],
                  cor_results$Relationship[i]))
      cat(sprintf("    %s\n", cor_results$Statement[i]))
      cat(sprintf("    r = %6.3f (SE = %5.3f), p = %s\n",
                  cor_results$Estimate[i],
                  cor_results$SE[i],
                  cor_results$p_formatted[i]))
      cat(sprintf("    Supported: %s\n\n", cor_results$Supported[i]))
    }
    
    cat(strrep("-", 75), "\n"); cat("  VARIANCE EXPLAINED (R²)\n")
    cat(strrep("-", 75), "\n\n")
    for (i in seq_len(nrow(rsq_df))) {
      cat(sprintf("  %-15s: %5.1f%%\n",
                  rsq_df$Construct[i],
                  rsq_df$Percent_Variance[i]))
    }
    cat("\n")
    
    cat(strrep("-", 75), "\n")
    cat("  MEASUREMENT DIAGNOSTICS\n")
    cat(strrep("-", 75), "\n\n")
    cat("  Reliability / CR / AVE:\n"); print(validity_table); cat("\n")
    cat("  HTMT:\n"); print(htmt_results); cat("\n")
    
    sink()
    
    cat("✓ Full export completed for", city_name, "\n")
    
  } else {
    
    ## ============================================================
    ## 11. INADMISSIBLE — MEASUREMENT DIAGNOSTICS ONLY
    ## ============================================================
    
    cat("\n--- STRUCTURAL RESULTS WITHHELD (inadmissible) ---\n")
    cat("  Structural paths, hypotheses, correlations, R²,\n")
    cat("  and the path diagram are NOT interpreted and NOT written.\n")
    if (!converged)  cat("  Reason: model did not converge.\n")
    if (!post_check) cat("  Reason: post-estimation check failed.\n")
    
    ## Only fit, loadings, reliability and HTMT are already on disk.
    ## No 07*/08* files are created, no structural plots.
    
    sink(file.path(cfg$out, "THESIS_SUMMARY_ModelA.txt"))
    cat(strrep("=", 75), "\n")
    cat("  MODEL A —", toupper(city_name), "\n")
    cat(strrep("=", 75), "\n\n")
    cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
    cat("N (complete cases):", nrow(sem_data), "\n")
    cat("Admissible solution:", admissible, "\n\n")
    
    cat(strrep("-", 75), "\n")
    cat("  MODEL STATUS\n")
    cat(strrep("-", 75), "\n\n")
    cat("  Converged:             ", converged, "\n")
    cat("  Post-estimation check: ", post_check, "\n")
    cat("  Admissible solution:   ", admissible, "\n\n")
    cat("  REASON:\n")
    if (!converged)  cat("    - Model did not converge.\n")
    if (!post_check) cat("    - Post-estimation check failed.\n")
    cat("\n  CONSEQUENCES:\n")
    cat("    - Structural paths (H1–H3) are NOT interpreted.\n")
    cat("    - Construct correlations (CH1–CH3) are NOT interpreted.\n")
    cat("    - R² for INT is NOT reported.\n")
    cat("    - Path diagram and structural plots are NOT produced.\n\n")
    
    cat(strrep("-", 75), "\n")
    cat("  MODEL FIT (diagnostic)\n")
    cat(strrep("-", 75), "\n\n")
    print(fit_table); cat("\n")
    
    cat(strrep("-", 75), "\n")
    cat("  MEASUREMENT DIAGNOSTICS\n")
    cat(strrep("-", 75), "\n\n")
    cat("  Factor loadings:\n"); print(loadings); cat("\n")
    cat("  Reliability / CR / AVE:\n"); print(validity_table); cat("\n")
    cat("  HTMT:\n"); print(htmt_results); cat("\n")
    
    sink()
    
    cat("✓ Diagnostic-only export completed for", city_name, "\n")
  }
}

## ============================================================
## 12. DONE
## ============================================================

cat("\n", strrep("=", 60), "\n")
cat("=== MODEL A ANALYSIS COMPLETE ===\n")
cat(strrep("=", 60), "\n\n")
cat("Outputs:\n")
cat("  outputs/model_a/vienna/\n")
cat("  outputs/model_a/madrid/\n")
