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
## Hypotheses tested:
##   H1:  ATT → INT  (Direct)        ATT has a positive effect on INT.
##   H2:  SN  → INT  (Direct)        SN has a positive effect on INT.
##   H3:  PBC → INT  (Direct)        PBC has a positive effect on INT.
##   CH1: ATT ↔ SN   (Correlation)   ATT and SN are positively associated.
##   CH2: SN  ↔ PBC  (Correlation)   SN and PBC are positively associated.
##   CH3: ATT ↔ PBC  (Correlation)   ATT and PBC are positively associated.
##
## Usage:
##   1. Place the input CSV in the working directory (default: "data.csv").
##   2. Adjust `input_csv` and `output_dir` below if needed.
##   3. Run the entire script.
##
## Outputs:
##   - Model fit indices
##   - Structural path coefficients (H1–H3)
##   - Construct correlations (CH1–CH3)
##   - Factor loadings
##   - Reliability (Cronbach's alpha), CR, AVE
##   - HTMT (discriminant validity)
##   - R² for intention
##   - CSV tables and PNG plots in `output_dir`
##   - A plain-text thesis summary
## ============================================================

rm(list = ls())
cat("\n=== SEM MODEL A — STANDARD TPB ===\n")

## ============================================================
## 1. PACKAGES
## ============================================================

packages <- c("lavaan", "semPlot", "tidyverse", "psych",
              "ggplot2", "corrplot")

installed <- packages %in% installed.packages()
if (any(!installed)) install.packages(packages[!installed])
invisible(lapply(packages, library, character.only = TRUE))

## ============================================================
## 2. FILE PATHS (EDIT HERE)
## ============================================================

input_csv  <- "data.csv"        # <-- change to your CSV file
output_dir <- "output_modelA"   # <-- folder for all outputs

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
plot_dir <- file.path(output_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

cat("✓ Input CSV:   ", input_csv, "\n")
cat("✓ Output dir:  ", output_dir, "\n")

## ============================================================
## 3. VARIABLE DEFINITIONS
## ============================================================

## Items per latent construct (PBC1 excluded by design)
sem_items <- list(
  ATT = c("ATT1", "ATT2", "ATT3"),
  SN  = c("SN1", "SN2", "SN3"),
  PBC = c("PBC2", "PBC3"),
  INT = c("INT1", "INT2", "INT3", "INT4")
)

## Named parameter labels in the model syntax guarantee that
## hypothesis assignment never depends on row order.
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

## ============================================================
## 4. HELPER FUNCTIONS
## ============================================================

safe_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

load_data <- function(path) {
  if (!file.exists(path)) stop(paste("File not found:", path))
  first_line <- readLines(path, n = 1, warn = FALSE)
  if (grepl(";", first_line)) {
    read.csv2(path, stringsAsFactors = FALSE,
              na.strings = c("", "NA", "NULL", "N/A"))
  } else {
    read.csv(path, stringsAsFactors = FALSE,
             na.strings = c("", "NA", "NULL", "N/A"))
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

calculate_htmt <- function(df, items1, items2) {
  cm  <- cor(df[, c(items1, items2)], use = "pairwise.complete.obs")
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
## 5. DATA IMPORT AND PREPARATION
## ============================================================

cat("\n--- DATA IMPORT ---\n")
raw <- load_data(input_csv)
cat("✓ Loaded:", nrow(raw), "rows,", ncol(raw), "columns\n")

items_all <- unlist(sem_items)
for (col in items_all[items_all %in% names(raw)]) {
  raw[[col]] <- safe_numeric(raw[[col]])
}

missing_items <- setdiff(items_all, names(raw))
if (length(missing_items) > 0) {
  cat("⚠ Missing items:", paste(missing_items, collapse = ", "), "\n")
}

sem_data <- raw[, items_all[items_all %in% names(raw)], drop = FALSE]
sem_data <- sem_data[complete.cases(sem_data), ]
cat("✓ Complete cases:", nrow(sem_data), "\n")

## ============================================================
## 6. DESCRIPTIVE STATISTICS
## ============================================================

cat("\n--- DESCRIPTIVE STATISTICS ---\n")

desc_stats <- data.frame(
  Item = names(sem_data),
  Mean = round(colMeans(sem_data), 2),
  SD   = round(apply(sem_data, 2, sd), 2),
  Min  = apply(sem_data, 2, min),
  Max  = apply(sem_data, 2, max)
)
print(desc_stats)
write.csv(desc_stats,
          file.path(output_dir, "01_descriptives_items.csv"),
          row.names = FALSE)

## ============================================================
## 7. CORRELATION MATRIX (CONSTRUCT LEVEL)
## ============================================================

cat("\n--- CONSTRUCT CORRELATIONS ---\n")

construct_scores <- data.frame(
  ATT = rowMeans(sem_data[, sem_items$ATT, drop = FALSE]),
  SN  = rowMeans(sem_data[, sem_items$SN,  drop = FALSE]),
  PBC = rowMeans(sem_data[, sem_items$PBC, drop = FALSE]),
  INT = rowMeans(sem_data[, sem_items$INT, drop = FALSE])
)
construct_cor <- round(cor(construct_scores), 3)
print(construct_cor)
write.csv(construct_cor,
          file.path(output_dir, "02_construct_correlations.csv"))

png(file.path(plot_dir, "01_correlation_matrix.png"),
    width = 7, height = 6, units = "in", res = 300)
corrplot(construct_cor, method = "color", type = "upper",
         addCoef.col = "black", tl.col = "black",
         title = "Construct Correlations — Model A",
         mar = c(0, 0, 2, 0))
dev.off()

## ============================================================
## 8. ESTIMATE MODEL A
## ============================================================

cat("\n--- ESTIMATING MODEL A ---\n")

fit <- sem(model_A_spec, data = sem_data,
           estimator = "MLR", missing = "listwise",
           control = list(iter.max = 1000))

converged  <- lavInspect(fit, "converged")
post_check <- tryCatch(lavInspect(fit, "post.check"),
                       error = function(e) FALSE)
cat("  Converged:", converged, "\n")
cat("  Post-check passed:", post_check, "\n")

## ============================================================
## 9. MODEL FIT
## ============================================================

cat("\n--- MODEL FIT ---\n")
fit_measures <- fitMeasures(fit, c("chisq","df","pvalue",
                                   "cfi","tli","rmsea",
                                   "rmsea.ci.lower","rmsea.ci.upper",
                                   "srmr","aic","bic"))
fit_table <- data.frame(
  Measure = names(fit_measures),
  Value   = round(as.numeric(fit_measures), 3)
)
print(fit_table)
write.csv(fit_table,
          file.path(output_dir, "03_model_fit_indices.csv"),
          row.names = FALSE)

## ============================================================
## 10. STRUCTURAL PATHS (H1, H2, H3)
## ============================================================

cat("\n--- STRUCTURAL PATHS (H1–H3) ---\n")

pe_std <- standardizedSolution(fit)

hyp_defs <- list(
  H1 = list(lhs = "INT", rhs = "ATT", label = "ATT → INT",
            statement = "ATT has a positive effect on INT."),
  H2 = list(lhs = "INT", rhs = "SN",  label = "SN → INT",
            statement = "SN has a positive effect on INT."),
  H3 = list(lhs = "INT", rhs = "PBC", label = "PBC → INT",
            statement = "PBC has a positive effect on INT.")
)

hypothesis_results <- data.frame(
  Hypothesis = character(), Type = character(),
  Relationship = character(), Statement = character(),
  Estimate = numeric(), SE = numeric(),
  p_value = numeric(), p_formatted = character(),
  Supported = character(), stringsAsFactors = FALSE
)

for (h in names(hyp_defs)) {
  def <- hyp_defs[[h]]
  matched <- pe_std[pe_std$op == "~" &
                      pe_std$lhs == def$lhs &
                      pe_std$rhs == def$rhs, ]
  if (nrow(matched) != 1) {
    stop(paste0("Could not uniquely match path for ", h, "."))
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
    Supported = supp, stringsAsFactors = FALSE
  ))
}

print(hypothesis_results[, c("Hypothesis","Relationship",
                             "Estimate","SE","p_formatted","Supported")])
write.csv(hypothesis_results,
          file.path(output_dir, "06a_hypotheses_direct.csv"),
          row.names = FALSE)

## ============================================================
## 11. CORRELATIONS (CH1, CH2, CH3)
## ============================================================

cat("\n--- CORRELATION HYPOTHESES (CH1–CH3) ---\n")

cor_defs <- list(
  CH1 = list(a = "ATT", b = "SN",  label = "ATT ↔ SN",
             statement = "ATT and SN are positively associated."),
  CH2 = list(a = "SN",  b = "PBC", label = "SN ↔ PBC",
             statement = "SN and PBC are positively associated."),
  CH3 = list(a = "ATT", b = "PBC", label = "ATT ↔ PBC",
             statement = "ATT and PBC are positively associated.")
)

cor_results <- data.frame(
  Hypothesis = character(), Type = character(),
  Relationship = character(), Statement = character(),
  Estimate = numeric(), SE = numeric(),
  p_value = numeric(), p_formatted = character(),
  Supported = character(), stringsAsFactors = FALSE
)

for (h in names(cor_defs)) {
  def <- cor_defs[[h]]
  matched <- pe_std[pe_std$op == "~~" &
                      ((pe_std$lhs == def$a & pe_std$rhs == def$b) |
                         (pe_std$lhs == def$b & pe_std$rhs == def$a)), ]
  if (nrow(matched) != 1) {
    stop(paste0("Could not uniquely match correlation for ", h, "."))
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
    Supported = supp, stringsAsFactors = FALSE
  ))
}

print(cor_results[, c("Hypothesis","Relationship",
                      "Estimate","SE","p_formatted","Supported")])
write.csv(cor_results,
          file.path(output_dir, "06b_hypotheses_correlation.csv"),
          row.names = FALSE)

## Combined hypothesis table
all_hypotheses <- rbind(hypothesis_results, cor_results)
write.csv(all_hypotheses,
          file.path(output_dir, "06_hypotheses_all.csv"),
          row.names = FALSE)

## ============================================================
## 12. R²
## ============================================================

cat("\n--- R² (VARIANCE EXPLAINED) ---\n")
rsq <- inspect(fit, "r2")
rsq_df <- data.frame(
  Construct        = names(rsq),
  R_squared        = round(as.numeric(rsq), 3),
  Percent_Variance = round(100 * as.numeric(rsq), 1)
)
print(rsq_df)
write.csv(rsq_df,
          file.path(output_dir, "07_rsquared.csv"),
          row.names = FALSE)

## ============================================================
## 13. FACTOR LOADINGS
## ============================================================

cat("\n--- FACTOR LOADINGS ---\n")
loadings <- pe_std[pe_std$op == "=~",
                   c("lhs","rhs","est.std","se","pvalue")]
loadings$p_formatted <- format_p(loadings$pvalue)
loadings$sig <- ifelse(loadings$pvalue < 0.001, "***",
                       ifelse(loadings$pvalue < 0.01, "**",
                              ifelse(loadings$pvalue < 0.05, "*", "")))
print(loadings)
write.csv(loadings,
          file.path(output_dir, "08_factor_loadings.csv"),
          row.names = FALSE)

## ============================================================
## 14. RELIABILITY, CR, AVE
## ============================================================

cat("\n--- RELIABILITY AND CONVERGENT VALIDITY ---\n")

validity_list <- list()
for (cons in names(sem_items)) {
  items_cons <- sem_items[[cons]]
  items_cons <- items_cons[items_cons %in% names(sem_data)]
  if (length(items_cons) < 2) next
  
  std_loads <- pe_std$est.std[pe_std$op == "=~" & pe_std$lhs == cons]
  std_loads <- std_loads[!is.na(std_loads)]
  
  cr_ave <- calculate_cr_ave(std_loads)
  alpha  <- tryCatch(
    psych::alpha(sem_data[, items_cons], check.keys = TRUE)$total$raw_alpha,
    error = function(e) NA
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
print(validity_table)
write.csv(validity_table,
          file.path(output_dir, "09_reliability_cr_ave.csv"),
          row.names = FALSE)

## ============================================================
## 15. HTMT
## ============================================================

cat("\n--- DISCRIMINANT VALIDITY (HTMT) ---\n")

cnames <- names(sem_items)
htmt_results <- data.frame(
  Pair = character(), HTMT = numeric(),
  Status = character(), stringsAsFactors = FALSE
)
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
      stringsAsFactors = FALSE
    ))
  }
}
print(htmt_results)
write.csv(htmt_results,
          file.path(output_dir, "10_htmt.csv"),
          row.names = FALSE)

## ============================================================
## 16. VISUALIZATIONS
## ============================================================

cat("\n--- VISUALIZATIONS ---\n")

## 16.1 Path diagram
tryCatch({
  png(file.path(plot_dir, "02_path_diagram.png"),
      width = 12, height = 8, units = "in", res = 300)
  semPaths(fit, what = "std", whatLabels = "std",
           layout = "tree2", rotation = 2,
           edge.label.cex = 0.9, style = "lisrel",
           nCharNodes = 8, sizeMan = 5, sizeLat = 8,
           edge.color = "black", fade = FALSE,
           title = TRUE,
           main = "Model A — Standard TPB")
  dev.off()
  cat("  ✓ Path diagram saved\n")
}, error = function(e) cat("  ⚠ Path diagram failed\n"))

## 16.2 Direct effect barplot (H1–H3)
coef_plot_data <- hypothesis_results
coef_plot_data$Significant <- ifelse(coef_plot_data$p_value < 0.05,
                                     "p < 0.05", "n.s.")

coef_barplot <- ggplot(coef_plot_data,
                       aes(x = reorder(Relationship, Estimate),
                           y = Estimate, fill = Significant)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  geom_hline(yintercept = 0) +
  labs(title = "Standardized Path Coefficients — Model A",
       subtitle = "Direct Effects: H1, H2, H3",
       x = "", y = "Standardized β") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("p < 0.05" = "steelblue",
                               "n.s."     = "coral"))
ggsave(file.path(plot_dir, "03_direct_effects.png"),
       coef_barplot, width = 7, height = 5, dpi = 300)

## 16.3 Correlation barplot (CH1–CH3)
cor_plot_data <- cor_results
cor_plot_data$Significant <- ifelse(cor_plot_data$p_value < 0.05,
                                    "p < 0.05", "n.s.")

cor_barplot <- ggplot(cor_plot_data,
                      aes(x = reorder(Relationship, Estimate),
                          y = Estimate, fill = Significant)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  geom_hline(yintercept = 0) +
  labs(title = "Construct Correlations — Model A",
       subtitle = "Correlation Hypotheses: CH1, CH2, CH3",
       x = "", y = "Standardized r") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("p < 0.05" = "steelblue",
                               "n.s."     = "coral"))
ggsave(file.path(plot_dir, "04_correlations.png"),
       cor_barplot, width = 7, height = 5, dpi = 300)

## 16.4 R² plot
rsq_plot <- ggplot(rsq_df,
                   aes(x = reorder(Construct, R_squared),
                       y = R_squared, fill = Construct)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  geom_text(aes(label = paste0(Percent_Variance, "%")),
            hjust = -0.2, size = 4) +
  labs(title = "Variance Explained (R²) — Model A",
       x = "", y = "R²") +
  coord_flip() +
  ylim(0, max(rsq_df$R_squared) + 0.1) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set2")
ggsave(file.path(plot_dir, "05_rsquared.png"),
       rsq_plot, width = 7, height = 5, dpi = 300)

## 16.5 HTMT barplot
htmt_plot <- ggplot(htmt_results,
                    aes(x = reorder(Pair, HTMT),
                        y = HTMT, fill = HTMT < 0.85)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  geom_hline(yintercept = 0.85, linetype = "dashed",
             color = "darkgreen") +
  geom_hline(yintercept = 0.90, linetype = "dashed",
             color = "orange") +
  labs(title = "HTMT Ratios — Model A",
       subtitle = "Green: 0.85; Orange: 0.90",
       x = "", y = "HTMT") +
  coord_flip() +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "coral"))
ggsave(file.path(plot_dir, "06_htmt.png"),
       htmt_plot, width = 8, height = 5, dpi = 300)

cat("  ✓ All plots saved\n")

## ============================================================
## 17. THESIS SUMMARY
## ============================================================

cat("\n--- THESIS SUMMARY ---\n")

sink(file.path(output_dir, "THESIS_SUMMARY_ModelA.txt"))

cat(strrep("=", 75), "\n")
cat("  MODEL A — STANDARD THEORY OF PLANNED BEHAVIOR\n")
cat(strrep("=", 75), "\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("N (complete cases):", nrow(sem_data), "\n\n")

cat(strrep("-", 75), "\n")
cat("  HYPOTHESES\n")
cat(strrep("-", 75), "\n\n")
cat("  H1:  ATT → INT  (Direct)       ATT has a positive effect on INT.\n")
cat("  H2:  SN  → INT  (Direct)       SN has a positive effect on INT.\n")
cat("  H3:  PBC → INT  (Direct)       PBC has a positive effect on INT.\n")
cat("  CH1: ATT ↔ SN   (Correlation)  ATT and SN are positively associated.\n")
cat("  CH2: SN  ↔ PBC  (Correlation)  SN and PBC are positively associated.\n")
cat("  CH3: ATT ↔ PBC  (Correlation)  ATT and PBC are positively associated.\n\n")
cat("  A hypothesis is supported if the estimate is positive and p < .05.\n\n")

cat(strrep("-", 75), "\n")
cat("  MODEL FIT\n")
cat(strrep("-", 75), "\n\n")
for (i in 1:nrow(fit_table)) {
  cat(sprintf("  %-15s: %8.3f\n", fit_table$Measure[i], fit_table$Value[i]))
}
cat("\n")

cat(strrep("-", 75), "\n")
cat("  DIRECT EFFECTS (H1–H3)\n")
cat(strrep("-", 75), "\n\n")
for (i in 1:nrow(hypothesis_results)) {
  cat(sprintf("  %s: %s\n", hypothesis_results$Hypothesis[i],
              hypothesis_results$Relationship[i]))
  cat(sprintf("    %s\n", hypothesis_results$Statement[i]))
  cat(sprintf("    β = %6.3f (SE = %5.3f), p = %s\n",
              hypothesis_results$Estimate[i],
              hypothesis_results$SE[i],
              hypothesis_results$p_formatted[i]))
  cat(sprintf("    Supported: %s\n\n", hypothesis_results$Supported[i]))
}

cat(strrep("-", 75), "\n")
cat("  CORRELATIONS (CH1–CH3)\n")
cat(strrep("-", 75), "\n\n")
for (i in 1:nrow(cor_results)) {
  cat(sprintf("  %s: %s\n", cor_results$Hypothesis[i],
              cor_results$Relationship[i]))
  cat(sprintf("    %s\n", cor_results$Statement[i]))
  cat(sprintf("    r = %6.3f (SE = %5.3f), p = %s\n",
              cor_results$Estimate[i],
              cor_results$SE[i],
              cor_results$p_formatted[i]))
  cat(sprintf("    Supported: %s\n\n", cor_results$Supported[i]))
}

cat(strrep("-", 75), "\n")
cat("  VARIANCE EXPLAINED (R²)\n")
cat(strrep("-", 75), "\n\n")
for (i in 1:nrow(rsq_df)) {
  cat(sprintf("  %-15s: %5.1f%%\n",
              rsq_df$Construct[i], rsq_df$Percent_Variance[i]))
}
cat("\n")

cat(strrep("-", 75), "\n")
cat("  RELIABILITY AND CONVERGENT VALIDITY\n")
cat(strrep("-", 75), "\n\n")
cat("  Criteria: Alpha ≥ 0.70, CR ≥ 0.70, AVE ≥ 0.50\n\n")
print(validity_table)
cat("\n")

cat(strrep("-", 75), "\n")
cat("  DISCRIMINANT VALIDITY (HTMT)\n")
cat(strrep("-", 75), "\n\n")
cat("  Criteria: HTMT < 0.85 (good), < 0.90 (acceptable)\n\n")
for (i in 1:nrow(htmt_results)) {
  cat(sprintf("  %-20s: HTMT = %6.3f  %s\n",
              htmt_results$Pair[i],
              htmt_results$HTMT[i],
              htmt_results$Status[i]))
}
cat("\n")

cat(strrep("-", 75), "\n")
cat("  FACTOR LOADINGS (STANDARDIZED)\n")
cat(strrep("-", 75), "\n\n")
print(loadings[, c("lhs","rhs","est.std","sig")])
cat("\n")

cat(strrep("=", 75), "\n")
cat("  END OF MODEL A SUMMARY\n")
cat(strrep("=", 75), "\n")

sink()

## ============================================================
## 18. DONE
## ============================================================

cat("\n", strrep("=", 60), "\n")
cat("=== MODEL A ANALYSIS COMPLETE ===\n")
cat(strrep("=", 60), "\n\n")
cat("Outputs saved in:", output_dir, "\n")
cat("Plots saved in:  ", plot_dir, "\n\n")
cat("Files:\n")
cat("  01_descriptives_items.csv\n")
cat("  02_construct_correlations.csv\n")
cat("  03_model_fit_indices.csv\n")
cat("  06a_hypotheses_direct.csv\n")
cat("  06b_hypotheses_correlation.csv\n")
cat("  06_hypotheses_all.csv\n")
cat("  07_rsquared.csv\n")
cat("  08_factor_loadings.csv\n")
cat("  09_reliability_cr_ave.csv\n")
cat("  10_htmt.csv\n")
cat("  THESIS_SUMMARY_ModelA.txt\n")