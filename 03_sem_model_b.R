## ============================================================
## STRUCTURAL EQUATION MODELING (SEM)
## MODEL B — Reduced TPB: ATT + PBC combined into one factor
## Sustainable Transport Choices — Tokyo
##
## Description:
##   Tokyo-specific alternative specification of the TPB.
##
##   Model A (standard TPB) did not produce a suitable solution
##   for the Tokyo sample: ATT and PBC were not empirically
##   distinguishable (HTMT ≥ 0.90). Model B therefore combines
##   ATT and PBC into a single latent factor
##   ("Behavioral Readiness" / BR).
##
##   Model B is reported as an alternative specification for
##   Tokyo. It is not framed as a sensitivity analysis, and no
##   formal model comparison (LRT) against Model A is performed.
##
## Relationships tested (Model B):
##   BR → INT   (Direct)        BR has a positive effect on INT.
##   SN → INT   (Direct)        SN has a positive effect on INT.
##   BR ↔ SN    (Correlation)   BR and SN are positively associated.
##
##   Model B uses NAMED parameters (h1, h2, ch1) so that
##   relationship assignment never depends on the row order of
##   the lavaan output.
##
## Usage:
##   1. Set the working directory to the repository root.
##   2. Run: source("scripts/03_sem_model_b.R")
##
## Outputs:
##   outputs/model_b/tokyo/
## ============================================================

cat("\n=== SEM MODEL B — TOKYO ALTERNATIVE SPECIFICATION ===\n")
cat("=== ATT + PBC combined into Behavioral Readiness ===\n\n")

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
## 2. FILE PATHS
## ============================================================

input_csv  <- "data/processed/tokyo.csv"
output_dir <- "outputs/model_b/tokyo"

if (!file.exists(input_csv)) {
  stop("Required input file missing: ", input_csv)
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
plot_dir <- file.path(output_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

cat("✓ Input CSV:   ", input_csv, "\n")
cat("✓ Output dir:  ", output_dir, "\n")

## ============================================================
## 3. MODEL SPECIFICATION
##
##   Model B (Reduced TPB — this script):
##     BR  =~ ATT1 + ATT2 + ATT3 + PBC2 + PBC3
##     SN  =~ SN1  + SN2  + SN3
##     INT =~ INT1 + INT2 + INT3 + INT4
##     INT ~ BR + SN
##     BR ~~ SN
##
##   PBC1 is excluded from the measurement model by design.
## ============================================================

model_B_spec <- '
  # Measurement model
  BR  =~ ATT1 + ATT2 + ATT3 + PBC2 + PBC3
  SN  =~ SN1  + SN2  + SN3
  INT =~ INT1 + INT2 + INT3 + INT4

  # Structural model — named paths for unambiguous extraction
  INT ~ h1*BR + h2*SN

  # Correlation between exogenous constructs — named for CH1
  BR ~~ ch1*SN
'

items_all <- c("ATT1","ATT2","ATT3","SN1","SN2","SN3",
               "PBC2","PBC3","INT1","INT2","INT3","INT4")

## ============================================================
## 4. RELATIONSHIP DEFINITIONS (MODEL B)
## ============================================================

relationships_B <- list(
  BR_INT = list(type      = "direct",
                lhs       = "INT",
                rhs       = "BR",
                label     = "BR → INT",
                statement = "BR has a positive effect on INT."),
  SN_INT = list(type      = "direct",
                lhs       = "INT",
                rhs       = "SN",
                label     = "SN → INT",
                statement = "SN has a positive effect on INT."),
  BR_SN  = list(type      = "correlation",
                lhs       = "BR",
                rhs       = "SN",
                label     = "BR ↔ SN",
                statement = "BR and SN are positively associated.")
)

## ============================================================
## 5. HELPER FUNCTIONS
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
## monotrait block.
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

## ------------------------------------------------------------
## Extract Model B relationships using named parameters
## ------------------------------------------------------------
extract_relationships_B <- function(fit) {
  
  pe_std <- standardizedSolution(fit)
  
  results <- data.frame(
    Relationship = character(),
    Type         = character(),
    Statement    = character(),
    Estimate     = numeric(),
    SE           = numeric(),
    p_value      = numeric(),
    p_formatted  = character(),
    Supported    = character(),
    stringsAsFactors = FALSE
  )
  
  ## --- Direct effects ---
  for (nm in names(relationships_B)) {
    r <- relationships_B[[nm]]
    if (r$type != "direct") next
    matched <- pe_std[pe_std$op == "~" &
                        pe_std$lhs == r$lhs &
                        pe_std$rhs == r$rhs, ]
    if (nrow(matched) == 0) {
      stop("Could not find path: ", r$label)
    }
    if (nrow(matched) > 1) {
      stop("Multiple matches for path: ", r$label)
    }
    est  <- matched$est.std
    se   <- matched$se
    pval <- matched$pvalue
    supp <- ifelse(!is.na(est) & !is.na(pval) & est > 0 & pval < 0.05,
                   "✓ Yes", "✗ No")
    results <- rbind(results, data.frame(
      Relationship = r$label,
      Type         = "Direct",
      Statement    = r$statement,
      Estimate     = round(est, 3),
      SE           = round(se, 3),
      p_value      = pval,
      p_formatted  = format_p(pval),
      Supported    = supp,
      stringsAsFactors = FALSE
    ))
  }
  
  ## --- Correlations ---
  for (nm in names(relationships_B)) {
    r <- relationships_B[[nm]]
    if (r$type != "correlation") next
    matched <- pe_std[pe_std$op == "~~" &
                        ((pe_std$lhs == r$lhs & pe_std$rhs == r$rhs) |
                           (pe_std$lhs == r$rhs & pe_std$rhs == r$lhs)), ]
    if (nrow(matched) == 0) {
      stop("Could not find correlation: ", r$label)
    }
    if (nrow(matched) > 1) {
      stop("Multiple matches for correlation: ", r$label)
    }
    est  <- matched$est.std
    se   <- matched$se
    pval <- matched$pvalue
    supp <- ifelse(!is.na(est) & !is.na(pval) & est > 0 & pval < 0.05,
                   "✓ Yes", "✗ No")
    results <- rbind(results, data.frame(
      Relationship = r$label,
      Type         = "Correlation",
      Statement    = r$statement,
      Estimate     = round(est, 3),
      SE           = round(se, 3),
      p_value      = pval,
      p_formatted  = format_p(pval),
      Supported    = supp,
      stringsAsFactors = FALSE
    ))
  }
  
  results
}

## ============================================================
## 6. DATA IMPORT AND PREPARATION
## ============================================================

cat("\n--- DATA IMPORT ---\n")
raw <- load_data(input_csv)
cat("✓ Loaded:", nrow(raw), "rows,", ncol(raw), "columns\n")

missing_items <- setdiff(items_all, names(raw))
if (length(missing_items) > 0) {
  stop("Missing required items: ", paste(missing_items, collapse = ", "))
}

for (col in items_all) raw[[col]] <- safe_numeric(raw[[col]])

sem_data <- raw[, items_all, drop = FALSE]
sem_data <- sem_data[complete.cases(sem_data), ]
cat("✓ Complete cases:", nrow(sem_data), "\n")

if (nrow(sem_data) < 30) {
  stop("Too few complete cases (N = ", nrow(sem_data), ").")
}

## ============================================================
## 7. ESTIMATE MODEL B
## ============================================================

cat("\n--- Estimating Model B ---\n")

fit_B <- sem(model_B_spec, data = sem_data,
             estimator = "MLR", missing = "listwise",
             control = list(iter.max = 1000))

converged  <- lavInspect(fit_B, "converged")
post_check <- tryCatch(lavInspect(fit_B, "post.check"),
                       error = function(e) FALSE)
cat("  Converged: ", converged, "\n")
cat("  Post-check:", post_check, "\n")

admissible <- isTRUE(converged) && isTRUE(post_check)
cat("  Admissible solution:", admissible, "\n")

## ============================================================
## 8. MODEL FIT
## ============================================================

cat("\n--- Model Fit ---\n")
fit_measures <- fitMeasures(fit_B, c("chisq", "df", "pvalue",
                                     "cfi", "tli", "rmsea",
                                     "rmsea.ci.lower", "rmsea.ci.upper",
                                     "srmr", "aic", "bic"))
fit_table <- data.frame(
  Measure = names(fit_measures),
  Value   = round(as.numeric(fit_measures), 3)
)
print(fit_table)
write.csv(fit_table,
          file.path(output_dir, "01_fit_indices.csv"),
          row.names = FALSE)

## ============================================================
## 9. STRUCTURAL RESULTS — ONLY IF ADMISSIBLE
## ============================================================

pe_std <- standardizedSolution(fit_B)

if (admissible) {
  
  cat("\n--- Structural Results (admissible) ---\n")
  
  ## ---- 9.1 Reported relationships ----
  rel_results <- extract_relationships_B(fit_B)
  print(rel_results[, c("Relationship", "Type", "Estimate",
                        "SE", "p_formatted", "Supported")])
  write.csv(rel_results,
            file.path(output_dir, "02_relationships.csv"),
            row.names = FALSE)
  
  ## ---- 9.2 R² ----
  rsq <- inspect(fit_B, "r2")
  rsq_df <- data.frame(
    Construct        = names(rsq),
    R_squared        = round(as.numeric(rsq), 3),
    Percent_Variance = round(100 * as.numeric(rsq), 1)
  )
  print(rsq_df)
  write.csv(rsq_df,
            file.path(output_dir, "03_rsquared.csv"),
            row.names = FALSE)
  
} else {
  
  cat("\n--- Structural results withheld (inadmissible) ---\n")
  if (!converged)  cat("  Reason: model did not converge.\n")
  if (!post_check) cat("  Reason: post-estimation check failed.\n")
  
  rel_results <- NULL
  rsq_df      <- NULL
}

## ============================================================
## 10. FACTOR LOADINGS
## ============================================================

cat("\n--- Factor Loadings ---\n")
loadings <- pe_std[pe_std$op == "=~",
                   c("lhs", "rhs", "est.std", "se", "pvalue")]
loadings$p_formatted <- format_p(loadings$pvalue)
loadings$sig <- ifelse(loadings$pvalue < 0.001, "***",
                       ifelse(loadings$pvalue < 0.01, "**",
                              ifelse(loadings$pvalue < 0.05, "*", "")))
print(loadings)
write.csv(loadings,
          file.path(output_dir, "04_factor_loadings.csv"),
          row.names = FALSE)

## ============================================================
## 11. RELIABILITY, CR, AVE
## ============================================================

cat("\n--- Reliability & Convergent Validity ---\n")

constructs_B <- list(
  BR  = c("ATT1","ATT2","ATT3","PBC2","PBC3"),
  SN  = c("SN1","SN2","SN3"),
  INT = c("INT1","INT2","INT3","INT4")
)

validity_B <- data.frame(
  Construct = character(),
  n_items   = integer(),
  Alpha     = numeric(),
  CR        = numeric(),
  AVE       = numeric(),
  stringsAsFactors = FALSE
)

for (cons in names(constructs_B)) {
  items_cons <- constructs_B[[cons]]
  items_cons <- items_cons[items_cons %in% names(sem_data)]
  if (length(items_cons) < 2) next
  
  std_loads <- pe_std$est.std[pe_std$op == "=~" & pe_std$lhs == cons]
  std_loads <- std_loads[!is.na(std_loads)]
  
  cr_ave <- calculate_cr_ave(std_loads)
  alpha_val <- tryCatch(
    psych::alpha(sem_data[, items_cons],
                 check.keys = FALSE)$total$raw_alpha,
    error = function(e) NA_real_
  )
  validity_B <- rbind(validity_B, data.frame(
    Construct = cons,
    n_items   = length(std_loads),
    Alpha     = round(alpha_val, 3),
    CR        = cr_ave$CR,
    AVE       = cr_ave$AVE,
    stringsAsFactors = FALSE
  ))
}
print(validity_B)
write.csv(validity_B,
          file.path(output_dir, "05_reliability_cr_ave.csv"),
          row.names = FALSE)

## ============================================================
## 12. HTMT
## ============================================================

cat("\n--- HTMT (Discriminant Validity) ---\n")

htmt_B <- data.frame(Pair = character(), HTMT = numeric(),
                     Status = character(), stringsAsFactors = FALSE)
cnames <- names(constructs_B)
for (i in 1:(length(cnames) - 1)) {
  for (j in (i + 1):length(cnames)) {
    hval <- calculate_htmt(sem_data,
                           constructs_B[[cnames[i]]],
                           constructs_B[[cnames[j]]])
    htmt_B <- rbind(htmt_B, data.frame(
      Pair   = paste(cnames[i], "↔", cnames[j]),
      HTMT   = hval,
      Status = ifelse(hval < 0.85, "✓ Good (<0.85)",
                      ifelse(hval < 0.90, "⚠ Acceptable (<0.90)",
                             "✗ Poor (≥0.90)")),
      stringsAsFactors = FALSE
    ))
  }
}
print(htmt_B)
write.csv(htmt_B,
          file.path(output_dir, "06_htmt.csv"),
          row.names = FALSE)

## ============================================================
## 13. VISUALIZATIONS
## ============================================================

cat("\n--- Generating Visualizations ---\n")

## 13.1 Path diagram
tryCatch({
  png(file.path(plot_dir, "01_path_diagram.png"),
      width = 12, height = 8, units = "in", res = 300)
  semPaths(fit_B,
           what       = "std",
           whatLabels = "std",
           layout     = "tree2",
           rotation   = 2,
           edge.label.cex = 0.9,
           style      = "lisrel",
           nCharNodes = 8,
           sizeMan    = 5,
           sizeLat    = 8,
           edge.color = "black",
           fade       = FALSE,
           title      = TRUE,
           main       = "Model B — Tokyo Alternative Specification")
  dev.off()
  cat("  ✓ Path diagram saved\n")
}, error = function(e) {
  cat("  ⚠ Path diagram failed:", conditionMessage(e), "\n")
})

## 13.2 Relationship barplot
if (!is.null(rel_results) && nrow(rel_results) > 0) {
  rel_plot_data <- rel_results
  rel_plot_data$Significant <- ifelse(rel_plot_data$p_value < 0.05,
                                      "p < 0.05", "n.s.")
  
  rel_plot <- ggplot(rel_plot_data,
                     aes(x = reorder(Relationship, Estimate),
                         y = Estimate, fill = Significant)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "solid") +
    geom_text(aes(label = p_formatted), hjust = -0.1, size = 3.5) +
    labs(title    = "Model B — Tokyo Alternative Specification",
         subtitle = "BR → INT, SN → INT, BR ↔ SN",
         x = "", y = "Standardized Estimate") +
    coord_flip() +
    ylim(min(0, min(rel_plot_data$Estimate) - 0.05),
         max(rel_plot_data$Estimate) + 0.1) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom") +
    scale_fill_manual(values = c("p < 0.05" = "steelblue",
                                 "n.s."     = "coral"))
  
  ggsave(file.path(plot_dir, "02_relationships.png"),
         rel_plot, width = 8, height = 5, dpi = 300)
  cat("  ✓ Relationship barplot saved\n")
}

## 13.3 R² plot
if (!is.null(rsq_df) && nrow(rsq_df) > 0) {
  rsq_plot <- ggplot(rsq_df,
                     aes(x = reorder(Construct, R_squared),
                         y = R_squared, fill = Construct)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    geom_hline(yintercept = 0.3, linetype = "dashed",
               color = "red", alpha = 0.5) +
    geom_text(aes(label = paste0(Percent_Variance, "%")),
              hjust = -0.2, size = 4) +
    labs(title = "Variance Explained (R²) — Model B (Tokyo)",
         x = "", y = "R-squared") +
    coord_flip() +
    ylim(0, max(rsq_df$R_squared) + 0.1) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          legend.position = "none") +
    scale_fill_brewer(palette = "Set2")
  
  ggsave(file.path(plot_dir, "03_rsquared.png"),
         rsq_plot, width = 7, height = 5, dpi = 300)
  cat("  ✓ R² plot saved\n")
}

## 13.4 HTMT barplot
if (!is.null(htmt_B) && nrow(htmt_B) > 0) {
  htmt_plot <- ggplot(htmt_B,
                      aes(x = reorder(Pair, HTMT), y = HTMT,
                          fill = HTMT < 0.85)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    geom_hline(yintercept = 0.85, linetype = "dashed",
               color = "darkgreen") +
    geom_hline(yintercept = 0.90, linetype = "dashed",
               color = "orange") +
    labs(title    = "HTMT Ratios — Model B (Tokyo)",
         subtitle = "Green: 0.85; Orange: 0.90",
         x = "", y = "HTMT") +
    coord_flip() +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none") +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "coral"))
  
  ggsave(file.path(plot_dir, "04_htmt.png"),
         htmt_plot, width = 8, height = 5, dpi = 300)
  cat("  ✓ HTMT barplot saved\n")
}

## ============================================================
## 14. THESIS SUMMARY
## ============================================================

cat("\n--- Generating Thesis Summary ---\n")

sink(file.path(output_dir, "THESIS_SUMMARY_ModelB.txt"))

cat(strrep("=", 75), "\n")
cat("  MODEL B — TOKYO ALTERNATIVE SPECIFICATION\n")
cat("  Reduced TPB: Behavioral Readiness (ATT + PBC) + SN → INT\n")
cat(strrep("=", 75), "\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("N (complete cases):", nrow(sem_data), "\n")
cat("Admissible solution:", admissible, "\n\n")

cat(strrep("-", 75), "\n")
cat("  SPECIFICATION NOTE\n")
cat(strrep("-", 75), "\n\n")
cat("  Model B is the Tokyo-specific alternative specification.\n")
cat("  Model A (standard TPB) did not produce a suitable solution\n")
cat("  for the Tokyo sample, because ATT and PBC were not empirically\n")
cat("  distinguishable (HTMT ≥ 0.90). In Model B, ATT and PBC are\n")
cat("  combined into a single latent factor\n")
cat("  (Behavioral Readiness, BR).\n\n")
cat("  Model B is reported as an alternative specification and is\n")
cat("  not framed as a sensitivity analysis. No formal model\n")
cat("  comparison against Model A is performed.\n\n")

cat(strrep("-", 75), "\n")
cat("  RELATIONSHIPS TESTED\n")
cat(strrep("-", 75), "\n\n")
cat("  BR → INT   (Direct)        BR has a positive effect on INT.\n")
cat("  SN → INT   (Direct)        SN has a positive effect on INT.\n")
cat("  BR ↔ SN    (Correlation)   BR and SN are positively associated.\n\n")
cat("  A relationship is considered supported if the estimate is\n")
cat("  positive and p < .05.\n\n")

cat(strrep("-", 75), "\n")
cat("  MODEL FIT\n")
cat(strrep("-", 75), "\n\n")
print(fit_table); cat("\n")

if (admissible && !is.null(rel_results)) {
  cat(strrep("-", 75), "\n")
  cat("  REPORTED RELATIONSHIPS\n")
  cat(strrep("-", 75), "\n\n")
  for (i in seq_len(nrow(rel_results))) {
    r <- rel_results[i, ]
    cat(sprintf("  %s\n", r$Relationship))
    cat(sprintf("    %s\n", r$Statement))
    cat(sprintf("    Estimate = %6.3f (SE = %5.3f), p = %s\n",
                r$Estimate, r$SE, r$p_formatted))
    cat(sprintf("    Supported: %s\n\n", r$Supported))
  }
} else {
  cat(strrep("-", 75), "\n")
  cat("  REPORTED RELATIONSHIPS\n")
  cat(strrep("-", 75), "\n\n")
  cat("  Not interpreted: the model did not produce an admissible\n")
  cat("  solution.\n\n")
}

if (admissible && !is.null(rsq_df)) {
  cat(strrep("-", 75), "\n")
  cat("  VARIANCE EXPLAINED (R²)\n")
  cat(strrep("-", 75), "\n\n")
  for (i in seq_len(nrow(rsq_df))) {
    cat(sprintf("  %-15s: %5.1f%%\n",
                rsq_df$Construct[i], rsq_df$Percent_Variance[i]))
  }
  cat("\n")
}

cat(strrep("-", 75), "\n")
cat("  FACTOR LOADINGS (STANDARDIZED)\n")
cat(strrep("-", 75), "\n\n")
load_print <- loadings[, c("lhs", "rhs", "est.std", "sig")]
names(load_print) <- c("Construct", "Item", "Std.Loading", "Sig")
print(load_print); cat("\n")

cat(strrep("-", 75), "\n")
cat("  RELIABILITY & CONVERGENT VALIDITY\n")
cat(strrep("-", 75), "\n\n")
cat("  Criteria: Alpha ≥ 0.70, CR ≥ 0.70, AVE ≥ 0.50\n\n")
print(validity_B); cat("\n")

cat(strrep("-", 75), "\n")
cat("  DISCRIMINANT VALIDITY (HTMT)\n")
cat(strrep("-", 75), "\n\n")
cat("  Criteria: HTMT < 0.85 (good), < 0.90 (acceptable)\n\n")
for (i in seq_len(nrow(htmt_B))) {
  cat(sprintf("  %-20s: HTMT = %6.3f  %s\n",
              htmt_B$Pair[i], htmt_B$HTMT[i], htmt_B$Status[i]))
}
cat("\n")

cat(strrep("=", 75), "\n")
cat("  END OF MODEL B SUMMARY\n")
cat(strrep("=", 75), "\n")

sink()

## ============================================================
## 15. DONE
## ============================================================

cat("\n", strrep("=", 60), "\n")
cat("=== MODEL B ANALYSIS COMPLETE (TOKYO) ===\n")
cat(strrep("=", 60), "\n\n")
cat("Outputs saved in:", output_dir, "\n")
cat("Plots saved in:  ", plot_dir, "\n\n")
cat("Files:\n")
cat("  01_fit_indices.csv\n")
cat("  02_relationships.csv\n")
cat("  03_rsquared.csv\n")
cat("  04_factor_loadings.csv\n")
cat("  05_reliability_cr_ave.csv\n")
cat("  06_htmt.csv\n")
cat("  THESIS_SUMMARY_ModelB.txt\n")
