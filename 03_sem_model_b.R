## ============================================================
## STRUCTURAL EQUATION MODELING (SEM)
## MODEL B — Reduced TPB: ATT + PBC combined into one factor
## Sustainable Transport Choices — Cross-City Study
##
## Description:
##   Reduced TPB model in which ATT and PBC are combined into a
##   single latent factor ("Behavioral Readiness" / BR). This is
##   a sensitivity analysis for cases where the standard TPB
##   model (Model A) shows insufficient discriminant validity
##   between ATT and PBC (e.g., HTMT ≥ 0.90).
##
##   Model B is estimated alongside Model A. A likelihood-ratio
##   test (lavTestLRT) determines whether the parsimony of
##   Model B is statistically justified versus Model A.
##
##   Model B uses NAMED parameters (h1, h2, ch1) so that
##   hypothesis assignment never depends on the row order of
##   the lavaan output.
##
## Hypotheses tested (Model B):
##   H1:  BR → INT  (Direct)        BR has a positive effect on INT.
##   H2:  SN → INT  (Direct)        SN has a positive effect on INT.
##   CH1: BR ↔ SN   (Correlation)   BR and SN are positively associated.
##
## Usage:
##   1. Place the input CSV in the working directory (default: "data.csv").
##   2. Adjust `input_csv` and `output_dir` below if needed.
##   3. Run the entire script.
##
## Outputs:
##   - Model fit comparison (Model A vs Model B)
##   - Likelihood-ratio test (LRT)
##   - Structural paths and hypothesis tests (H1, H2, CH1)
##   - R², factor loadings, reliability (alpha), CR, AVE, HTMT
##   - CSV tables and PNG plots in `output_dir`
##   - A plain-text thesis summary
## ============================================================

rm(list = ls())
cat("\n=== SEM MODEL B — REDUCED TPB ===\n")
cat("=== ATT + PBC combined into Behavioral Readiness ===\n\n")

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
output_dir <- "output_modelB"   # <-- folder for all outputs

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
plot_dir <- file.path(output_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

cat("✓ Input CSV:   ", input_csv, "\n")
cat("✓ Output dir:  ", output_dir, "\n")

## ============================================================
## 3. MODEL SPECIFICATIONS
##
##   Model A (Standard TPB — reference):
##     ATT =~ ATT1 + ATT2 + ATT3
##     SN  =~ SN1  + SN2  + SN3
##     PBC =~ PBC2 + PBC3
##     INT =~ INT1 + INT2 + INT3 + INT4
##     INT ~ ATT + SN + PBC
##
##   Model B (Reduced TPB — this script):
##     BR  =~ ATT1 + ATT2 + ATT3 + PBC2 + PBC3
##     SN  =~ SN1  + SN2  + SN3
##     INT =~ INT1 + INT2 + INT3 + INT4
##     INT ~ BR + SN
##
##   PBC1 is excluded from the measurement model by design.
## ============================================================

model_A_spec <- '
  ATT =~ ATT1 + ATT2 + ATT3
  SN  =~ SN1  + SN2  + SN3
  PBC =~ PBC2 + PBC3
  INT =~ INT1 + INT2 + INT3 + INT4

  INT ~ ATT + SN + PBC

  ATT ~~ SN
  ATT ~~ PBC
  SN  ~~ PBC
'

model_B_spec <- '
  # Measurement model
  BR  =~ ATT1 + ATT2 + ATT3 + PBC2 + PBC3
  SN  =~ SN1  + SN2  + SN3
  INT =~ INT1 + INT2 + INT3 + INT4

  # Structural model — named paths for unambiguous hypothesis assignment
  INT ~ h1*BR + h2*SN

  # Correlation between exogenous constructs — named for CH1
  BR ~~ ch1*SN
'

items_all <- c("ATT1","ATT2","ATT3","SN1","SN2","SN3",
               "PBC2","PBC3","INT1","INT2","INT3","INT4")

## ============================================================
## 4. HYPOTHESIS DEFINITIONS (MODEL B)
## ============================================================

hypotheses_B <- list(
  H1  = list(type      = "direct",
             lhs       = "INT",
             rhs       = "BR",
             label     = "BR → INT",
             statement = "BR has a positive effect on INT."),
  H2  = list(type      = "direct",
             lhs       = "INT",
             rhs       = "SN",
             label     = "SN → INT",
             statement = "SN has a positive effect on INT."),
  CH1 = list(type      = "correlation",
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

check_converged <- function(fit) {
  conv <- lavInspect(fit, "converged")
  if (!conv) return(FALSE)
  resid_diag <- diag(lavInspect(fit, "theta"))
  if (any(resid_diag < 0)) return(FALSE)
  pe <- parameterEstimates(fit, standardized = TRUE)
  struct <- pe[pe$op == "~", ]
  if (any(!is.na(struct$std.all) & abs(struct$std.all) > 5)) return(FALSE)
  TRUE
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

format_p <- function(p) {
  if (is.na(p)) return("NA")
  ifelse(p < 0.001, "< .001", sprintf("%.4f", p))
}

## ------------------------------------------------------------
## Extract Model B hypotheses using named parameters
## ------------------------------------------------------------
extract_hypotheses_B <- function(fit) {
  
  pe_std <- standardizedSolution(fit)
  
  results <- data.frame(
    Hypothesis   = character(),
    Type         = character(),
    Relationship = character(),
    Statement    = character(),
    Estimate     = numeric(),
    SE           = numeric(),
    p_value      = numeric(),
    p_formatted  = character(),
    Supported    = character(),
    stringsAsFactors = FALSE
  )
  
  ## --- Direct effects (H1, H2) ---
  direct_defs <- hypotheses_B[sapply(hypotheses_B, function(h) h$type == "direct")]
  for (h_name in names(direct_defs)) {
    h <- direct_defs[[h_name]]
    matched <- pe_std[pe_std$op == "~" &
                        pe_std$lhs == h$lhs &
                        pe_std$rhs == h$rhs, ]
    
    if (nrow(matched) == 0) {
      stop(paste0("Could not find path for ", h_name,
                  " (", h$rhs, " → ", h$lhs, ")."))
    }
    if (nrow(matched) > 1) {
      stop(paste0("Multiple matches for path ", h_name,
                  " (", h$rhs, " → ", h$lhs, ")."))
    }
    
    est  <- matched$est.std
    se   <- matched$se
    pval <- matched$pvalue
    supp <- ifelse(!is.na(est) & !is.na(pval) & est > 0 & pval < 0.05,
                   "✓ Yes", "✗ No")
    
    results <- rbind(results, data.frame(
      Hypothesis   = h_name,
      Type         = "Direct",
      Relationship = h$label,
      Statement    = h$statement,
      Estimate     = round(est, 3),
      SE           = round(se, 3),
      p_value      = pval,
      p_formatted  = format_p(pval),
      Supported    = supp,
      stringsAsFactors = FALSE
    ))
  }
  
  ## --- Correlation (CH1) ---
  corr_defs <- hypotheses_B[sapply(hypotheses_B, function(h) h$type == "correlation")]
  for (h_name in names(corr_defs)) {
    h <- corr_defs[[h_name]]
    matched <- pe_std[pe_std$op == "~~" &
                        ((pe_std$lhs == h$lhs & pe_std$rhs == h$rhs) |
                           (pe_std$lhs == h$rhs & pe_std$rhs == h$lhs)), ]
    
    if (nrow(matched) == 0) {
      stop(paste0("Could not find correlation for ", h_name,
                  " (", h$lhs, " ↔ ", h$rhs, ")."))
    }
    if (nrow(matched) > 1) {
      stop(paste0("Multiple matches for correlation ", h_name,
                  " (", h$lhs, " ↔ ", h$rhs, ")."))
    }
    
    est  <- matched$est.std
    se   <- matched$se
    pval <- matched$pvalue
    supp <- ifelse(!is.na(est) & !is.na(pval) & est > 0 & pval < 0.05,
                   "✓ Yes", "✗ No")
    
    results <- rbind(results, data.frame(
      Hypothesis   = h_name,
      Type         = "Correlation",
      Relationship = h$label,
      Statement    = h$statement,
      Estimate     = round(est, 3),
      SE           = round(se, 3),
      p_value      = pval,
      p_formatted  = format_p(pval),
      Supported    = supp,
      stringsAsFactors = FALSE
    ))
  }
  
  results <- results[order(results$Hypothesis), ]
  rownames(results) <- NULL
  results
}

## ============================================================
## 6. DATA IMPORT AND PREPARATION
## ============================================================

cat("\n--- DATA IMPORT ---\n")
raw <- load_data(input_csv)
cat("✓ Loaded:", nrow(raw), "rows,", ncol(raw), "columns\n")

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
## 7. ESTIMATE MODEL A (REFERENCE)
## ============================================================

cat("\n--- Model A: Standard TPB (reference) ---\n")
fit_A <- tryCatch(
  sem(model_A_spec, data = sem_data, estimator = "MLR",
      missing = "listwise", control = list(iter.max = 1000)),
  error = function(e) { cat("  ✗ Model A error:", conditionMessage(e), "\n"); NULL }
)
conv_A <- if (!is.null(fit_A)) check_converged(fit_A) else FALSE
cat("  Converged:", ifelse(conv_A, "✓ Yes", "✗ No"), "\n")

## ============================================================
## 8. ESTIMATE MODEL B
## ============================================================

cat("\n--- Model B: Reduced TPB (BR = ATT + PBC) ---\n")
fit_B <- tryCatch(
  sem(model_B_spec, data = sem_data, estimator = "MLR",
      missing = "listwise", control = list(iter.max = 1000)),
  error = function(e) { cat("  ✗ Model B error:", conditionMessage(e), "\n"); NULL }
)
conv_B <- if (!is.null(fit_B)) check_converged(fit_B) else FALSE
cat("  Converged:", ifelse(conv_B, "✓ Yes", "✗ No"), "\n")

## ============================================================
## 9. MODEL FIT COMPARISON
## ============================================================

cat("\n--- Model Fit Comparison ---\n")

get_fit <- function(fit, label, conv) {
  if (is.null(fit) || !conv) {
    return(data.frame(Model = label, N = nrow(sem_data),
                      CFI = NA, TLI = NA, RMSEA = NA, SRMR = NA,
                      AIC = NA, BIC = NA, chisq = NA, df = NA,
                      Converged = "✗ No", stringsAsFactors = FALSE))
  }
  fm <- fitMeasures(fit, c("chisq","df","cfi","tli","rmsea","srmr","aic","bic"))
  data.frame(Model = label, N = nrow(sem_data),
             CFI   = round(fm["cfi"],   3),
             TLI   = round(fm["tli"],   3),
             RMSEA = round(fm["rmsea"], 3),
             SRMR  = round(fm["srmr"],  3),
             AIC   = round(fm["aic"],   1),
             BIC   = round(fm["bic"],   1),
             chisq = round(fm["chisq"], 3),
             df    = round(fm["df"],    0),
             Converged = "✓ Yes",
             stringsAsFactors = FALSE)
}

fit_comparison <- rbind(
  get_fit(fit_A, "A (Standard TPB)", conv_A),
  get_fit(fit_B, "B (Reduced TPB)",  conv_B)
)

print(fit_comparison[, c("Model","CFI","TLI","RMSEA","SRMR","AIC","BIC","Converged")])
write.csv(fit_comparison,
          file.path(output_dir, "01_fit_comparison.csv"),
          row.names = FALSE)

## ============================================================
## 10. LRT MODEL COMPARISON (MODEL B vs MODEL A)
## ============================================================

cat("\n--- Model Comparison (LRT) ---\n")

lrt_result <- NULL
if (!is.null(fit_A) && conv_A && !is.null(fit_B) && conv_B) {
  lrt <- tryCatch(lavTestLRT(fit_B, fit_A), error = function(e) NULL)
  if (!is.null(lrt)) {
    chi2_diff <- round(lrt$`Chisq diff`[2], 3)
    df_diff   <- lrt$`Df diff`[2]
    p_lrt     <- lrt$`Pr(>Chisq)`[2]
    preferred <- ifelse(is.na(p_lrt), "Inconclusive",
                        ifelse(p_lrt < 0.05,
                               "Model A — significant fit loss with B",
                               "Model B — no significant fit loss"))
    cat(sprintf("  Δχ² = %.3f, Δdf = %d, p = %s\n",
                chi2_diff, df_diff, format_p(p_lrt)))
    cat("  Preferred model:", preferred, "\n")
    
    lrt_result <- data.frame(
      Chi2_diff   = chi2_diff,
      df_diff     = df_diff,
      p_value     = p_lrt,
      p_formatted = format_p(p_lrt),
      Preferred   = preferred,
      stringsAsFactors = FALSE
    )
    write.csv(lrt_result,
              file.path(output_dir, "02_lrt_results.csv"),
              row.names = FALSE)
  } else {
    cat("  ⚠ LRT could not be computed\n")
  }
} else {
  cat("  ⚠ LRT skipped — one or both models did not converge\n")
}

## ============================================================
## 11. HYPOTHESES MODEL B (H1, H2, CH1)
## ============================================================

cat("\n--- Model B: Hypothesis Testing (H1, H2, CH1) ---\n")

hyp_B <- NULL
if (!is.null(fit_B) && conv_B) {
  hyp_B <- tryCatch(extract_hypotheses_B(fit_B),
                    error = function(e) {
                      cat("  ✗ Hypothesis extraction failed:", conditionMessage(e), "\n")
                      NULL
                    })
  if (!is.null(hyp_B)) {
    print(hyp_B[, c("Hypothesis","Type","Relationship",
                    "Estimate","SE","p_formatted","Supported")])
    write.csv(hyp_B,
              file.path(output_dir, "03_hypotheses.csv"),
              row.names = FALSE)
  }
} else {
  cat("  ⚠ Hypotheses not evaluated — Model B did not converge\n")
}

## ============================================================
## 12. R² MODEL B
## ============================================================

cat("\n--- Model B: R² ---\n")
rsq_df <- NULL
if (!is.null(fit_B) && conv_B) {
  rsq <- tryCatch(inspect(fit_B, "r2"), error = function(e) NULL)
  if (!is.null(rsq)) {
    rsq_df <- data.frame(
      Construct        = names(rsq),
      R_squared        = round(as.numeric(rsq), 3),
      Percent_Variance = round(100 * as.numeric(rsq), 1)
    )
    print(rsq_df)
    write.csv(rsq_df,
              file.path(output_dir, "04_rsquared.csv"),
              row.names = FALSE)
  }
} else {
  cat("  ⚠ R² not available\n")
}

## ============================================================
## 13. FACTOR LOADINGS MODEL B
## ============================================================

cat("\n--- Model B: Factor Loadings ---\n")
load_B <- NULL
if (!is.null(fit_B) && conv_B) {
  pe_B <- parameterEstimates(fit_B, standardized = TRUE)
  load_B <- pe_B[pe_B$op == "=~",
                 c("lhs","rhs","est","se","z","pvalue","std.all")]
  load_B$sig <- ifelse(is.na(load_B$pvalue), "NA",
                       ifelse(load_B$pvalue < 0.001, "***",
                              ifelse(load_B$pvalue < 0.01,  "**",
                                     ifelse(load_B$pvalue < 0.05,  "*", ""))))
  load_B[, c("est","se","z","std.all")] <-
    round(load_B[, c("est","se","z","std.all")], 3)
  print(load_B)
  write.csv(load_B,
            file.path(output_dir, "05_factor_loadings.csv"),
            row.names = FALSE)
}

## ============================================================
## 14. RELIABILITY, CR, AVE MODEL B
## ============================================================

cat("\n--- Model B: Reliability & Convergent Validity ---\n")

validity_B <- data.frame(
  Construct = character(),
  n_items   = integer(),
  Alpha     = numeric(),
  CR        = numeric(),
  AVE       = numeric(),
  stringsAsFactors = FALSE
)

if (!is.null(fit_B) && conv_B) {
  constructs_B <- list(
    BR  = c("ATT1","ATT2","ATT3","PBC2","PBC3"),
    SN  = c("SN1","SN2","SN3"),
    INT = c("INT1","INT2","INT3","INT4")
  )
  pe_B_raw   <- parameterEstimates(fit_B, standardized = TRUE)
  load_B_raw <- pe_B_raw[pe_B_raw$op == "=~", ]
  
  for (cons in names(constructs_B)) {
    items_cons <- constructs_B[[cons]]
    items_cons <- items_cons[items_cons %in% names(sem_data)]
    std_loads  <- load_B_raw[load_B_raw$lhs == cons, "std.all"]
    std_loads  <- std_loads[!is.na(std_loads)]
    
    if (length(std_loads) >= 2) {
      cr_ave <- calculate_cr_ave(std_loads)
      alpha_val <- tryCatch(
        psych::alpha(sem_data[, items_cons], check.keys = TRUE)$total$raw_alpha,
        error = function(e) NA
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
  }
  print(validity_B)
  write.csv(validity_B,
            file.path(output_dir, "06_reliability_cr_ave.csv"),
            row.names = FALSE)
}

## ============================================================
## 15. HTMT MODEL B
## ============================================================

cat("\n--- Model B: HTMT (Discriminant Validity) ---\n")

htmt_B <- NULL
if (!is.null(fit_B) && conv_B) {
  constructs_htmt <- list(
    BR  = c("ATT1","ATT2","ATT3","PBC2","PBC3"),
    SN  = c("SN1","SN2","SN3"),
    INT = c("INT1","INT2","INT3","INT4")
  )
  htmt_B <- data.frame(Pair = character(), HTMT = numeric(),
                       Status = character(), stringsAsFactors = FALSE)
  cnames <- names(constructs_htmt)
  for (i in 1:(length(cnames) - 1)) {
    for (j in (i + 1):length(cnames)) {
      hval <- calculate_htmt(sem_data,
                             constructs_htmt[[cnames[i]]],
                             constructs_htmt[[cnames[j]]])
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
  
  htmt_crit <- htmt_B[htmt_B$HTMT >= 0.90, ]
  if (nrow(htmt_crit) > 0) {
    cat("  ⚠ Remaining HTMT issues in Model B:\n")
    print(htmt_crit)
  } else {
    cat("  ✓ All HTMT values acceptable in Model B\n")
  }
  
  write.csv(htmt_B,
            file.path(output_dir, "07_htmt.csv"),
            row.names = FALSE)
}

## ============================================================
## 16. VISUALIZATIONS
## ============================================================

cat("\n--- Generating Visualizations ---\n")

## 16.1 Path diagram Model B
if (!is.null(fit_B) && conv_B) {
  tryCatch({
    png(file.path(plot_dir, "01_path_diagram_modelB.png"),
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
             main       = "Model B — Reduced TPB (BR = ATT + PBC)")
    dev.off()
    cat("  ✓ Path diagram saved\n")
  }, error = function(e) cat("  ⚠ Path diagram failed\n"))
}

## 16.2 Fit comparison plot (Model A vs Model B)
if (all(!is.na(fit_comparison$CFI))) {
  fit_long <- fit_comparison %>%
    select(Model, CFI, TLI, RMSEA, SRMR) %>%
    pivot_longer(cols = c(CFI, TLI, RMSEA, SRMR),
                 names_to = "Index", values_to = "Value")
  
  fit_plot <- ggplot(fit_long, aes(x = Index, y = Value, fill = Model)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
    geom_hline(yintercept = 0.95, linetype = "dashed",
               color = "darkgreen", alpha = 0.6) +
    geom_hline(yintercept = 0.08, linetype = "dashed",
               color = "orange", alpha = 0.6) +
    labs(title    = "Model Fit Comparison — Model A vs Model B",
         subtitle = "Dashed lines: CFI/TLI ≥ 0.95 (good); RMSEA/SRMR ≤ 0.08 (acceptable)",
         x = "Fit Index", y = "Value") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom") +
    scale_fill_manual(values = c("A (Standard TPB)" = "steelblue",
                                 "B (Reduced TPB)"  = "coral"))
  
  ggsave(file.path(plot_dir, "02_fit_comparison.png"),
         fit_plot, width = 8, height = 5, dpi = 300)
  cat("  ✓ Fit comparison plot saved\n")
}

## 16.3 Hypothesis barplot Model B
if (!is.null(hyp_B) && nrow(hyp_B) > 0) {
  hyp_plot_data <- hyp_B
  hyp_plot_data$Significant <- ifelse(hyp_plot_data$p_value < 0.05,
                                      "p < 0.05", "n.s.")
  
  hyp_plot <- ggplot(hyp_plot_data,
                     aes(x = reorder(Relationship, Estimate),
                         y = Estimate, fill = Significant)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "solid") +
    geom_text(aes(label = p_formatted), hjust = -0.1, size = 3.5) +
    labs(title    = "Model B Hypotheses",
         subtitle = "H1: BR → INT, H2: SN → INT, CH1: BR ↔ SN",
         x = "", y = "Standardized Estimate") +
    coord_flip() +
    ylim(min(0, min(hyp_plot_data$Estimate) - 0.05),
         max(hyp_plot_data$Estimate) + 0.1) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom") +
    scale_fill_manual(values = c("p < 0.05" = "steelblue",
                                 "n.s."     = "coral"))
  
  ggsave(file.path(plot_dir, "03_hypotheses.png"),
         hyp_plot, width = 8, height = 5, dpi = 300)
  cat("  ✓ Hypothesis barplot saved\n")
}

## 16.4 R² plot Model B
if (!is.null(rsq_df) && nrow(rsq_df) > 0) {
  rsq_plot <- ggplot(rsq_df,
                     aes(x = reorder(Construct, R_squared),
                         y = R_squared, fill = Construct)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    geom_hline(yintercept = 0.3, linetype = "dashed",
               color = "red", alpha = 0.5) +
    geom_text(aes(label = paste0(Percent_Variance, "%")),
              hjust = -0.2, size = 4) +
    labs(title = "Variance Explained (R²) — Model B",
         x = "", y = "R-squared") +
    coord_flip() +
    ylim(0, max(rsq_df$R_squared) + 0.1) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          legend.position = "none") +
    scale_fill_brewer(palette = "Set2")
  
  ggsave(file.path(plot_dir, "04_rsquared.png"),
         rsq_plot, width = 7, height = 5, dpi = 300)
  cat("  ✓ R² plot saved\n")
}

## 16.5 HTMT barplot
if (!is.null(htmt_B) && nrow(htmt_B) > 0) {
  htmt_plot <- ggplot(htmt_B,
                      aes(x = reorder(Pair, HTMT), y = HTMT,
                          fill = HTMT < 0.85)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    geom_hline(yintercept = 0.85, linetype = "dashed",
               color = "darkgreen") +
    geom_hline(yintercept = 0.90, linetype = "dashed",
               color = "orange") +
    labs(title    = "HTMT Ratios — Model B",
         subtitle = "Green: 0.85; Orange: 0.90",
         x = "", y = "HTMT") +
    coord_flip() +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none") +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "coral"))
  
  ggsave(file.path(plot_dir, "05_htmt.png"),
         htmt_plot, width = 8, height = 5, dpi = 300)
  cat("  ✓ HTMT barplot saved\n")
}

## ============================================================
## 17. THESIS SUMMARY
## ============================================================

cat("\n--- Generating Thesis Summary ---\n")

sink(file.path(output_dir, "THESIS_SUMMARY_ModelB.txt"))

cat(strrep("=", 75), "\n")
cat("  MODEL B — REDUCED TPB (BR = ATT + PBC) + SN → INT\n")
cat(strrep("=", 75), "\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("N (complete cases):", nrow(sem_data), "\n\n")

cat(strrep("-", 75), "\n")
cat("  RATIONALE\n")
cat(strrep("-", 75), "\n\n")
cat("  Model B is a sensitivity analysis for cases where the standard\n")
cat("  TPB model (Model A) shows insufficient discriminant validity\n")
cat("  between ATT and PBC (HTMT ≥ 0.90). In Model B, ATT and PBC are\n")
cat("  combined into a single latent factor (Behavioral Readiness, BR).\n")
cat("  A likelihood-ratio test determines whether the parsimony of\n")
cat("  Model B is statistically justified versus Model A.\n\n")

cat(strrep("-", 75), "\n")
cat("  HYPOTHESES TESTED (MODEL B)\n")
cat(strrep("-", 75), "\n\n")
cat("  H1:  BR → INT  (Direct)       BR has a positive effect on INT.\n")
cat("  H2:  SN → INT  (Direct)       SN has a positive effect on INT.\n")
cat("  CH1: BR ↔ SN   (Correlation)  BR and SN are positively associated.\n\n")
cat("  A hypothesis is considered supported if the estimate is positive\n")
cat("  and p < .05.\n\n")

cat(strrep("-", 75), "\n")
cat("  MODEL FIT COMPARISON\n")
cat(strrep("-", 75), "\n\n")
print(fit_comparison[, c("Model","CFI","TLI","RMSEA","SRMR","AIC","BIC","Converged")])
cat("\n")

if (!is.null(lrt_result)) {
  cat(strrep("-", 75), "\n")
  cat("  LRT MODEL COMPARISON (Model B vs Model A)\n")
  cat(strrep("-", 75), "\n\n")
  cat("  H0: Model B fits as well as Model A (no significant fit loss).\n")
  cat("  p < 0.05: Model A preferred (combining ATT+PBC costs fit).\n")
  cat("  p ≥ 0.05: Model B acceptable (combining ATT+PBC is justified).\n\n")
  cat(sprintf("  Δχ² = %.3f, Δdf = %d, p = %s\n",
              lrt_result$Chi2_diff, lrt_result$df_diff,
              lrt_result$p_formatted))
  cat("  Preferred model:", lrt_result$Preferred, "\n\n")
}

if (!is.null(hyp_B) && nrow(hyp_B) > 0) {
  cat(strrep("-", 75), "\n")
  cat("  MODEL B HYPOTHESIS TESTING (H1, H2, CH1)\n")
  cat(strrep("-", 75), "\n\n")
  for (i in seq_len(nrow(hyp_B))) {
    cat(sprintf("  %s: %s\n", hyp_B$Hypothesis[i], hyp_B$Relationship[i]))
    cat(sprintf("    %s\n", hyp_B$Statement[i]))
    cat(sprintf("    Estimate = %6.3f (SE = %5.3f), p = %s\n",
                hyp_B$Estimate[i], hyp_B$SE[i], hyp_B$p_formatted[i]))
    cat(sprintf("    Supported: %s\n\n", hyp_B$Supported[i]))
  }
}

if (!is.null(rsq_df) && nrow(rsq_df) > 0) {
  cat(strrep("-", 75), "\n")
  cat("  VARIANCE EXPLAINED (R²)\n")
  cat(strrep("-", 75), "\n\n")
  for (i in seq_len(nrow(rsq_df))) {
    cat(sprintf("  %-15s: %5.1f%%\n",
                rsq_df$Construct[i], rsq_df$Percent_Variance[i]))
  }
  cat("\n")
}

if (!is.null(load_B) && nrow(load_B) > 0) {
  cat(strrep("-", 75), "\n")
  cat("  FACTOR LOADINGS (STANDARDIZED)\n")
  cat(strrep("-", 75), "\n\n")
  load_print <- load_B[, c("lhs","rhs","std.all","sig")]
  names(load_print) <- c("Construct","Item","Std.Loading","Sig")
  print(load_print)
  cat("\n")
}

if (!is.null(validity_B) && nrow(validity_B) > 0) {
  cat(strrep("-", 75), "\n")
  cat("  RELIABILITY & CONVERGENT VALIDITY\n")
  cat(strrep("-", 75), "\n\n")
  cat("  Criteria: Alpha ≥ 0.70, CR ≥ 0.70, AVE ≥ 0.50\n\n")
  print(validity_B)
  cat("\n")
}

if (!is.null(htmt_B) && nrow(htmt_B) > 0) {
  cat(strrep("-", 75), "\n")
  cat("  DISCRIMINANT VALIDITY (HTMT)\n")
  cat(strrep("-", 75), "\n\n")
  cat("  Criteria: HTMT < 0.85 (good), < 0.90 (acceptable)\n\n")
  for (i in seq_len(nrow(htmt_B))) {
    cat(sprintf("  %-20s: HTMT = %6.3f  %s\n",
                htmt_B$Pair[i], htmt_B$HTMT[i], htmt_B$Status[i]))
  }
  cat("\n")
}

cat(strrep("=", 75), "\n")
cat("  END OF MODEL B SUMMARY\n")
cat(strrep("=", 75), "\n")

sink()

## ============================================================
## 18. DONE
## ============================================================

cat("\n", strrep("=", 60), "\n")
cat("=== MODEL B ANALYSIS COMPLETE ===\n")
cat(strrep("=", 60), "\n\n")
cat("Outputs saved in:", output_dir, "\n")
cat("Plots saved in:  ", plot_dir, "\n\n")
cat("Files:\n")
cat("  01_fit_comparison.csv\n")
cat("  02_lrt_results.csv\n")
cat("  03_hypotheses.csv\n")
cat("  04_rsquared.csv\n")
cat("  05_factor_loadings.csv\n")
cat("  06_reliability_cr_ave.csv\n")
cat("  07_htmt.csv\n")
cat("  THESIS_SUMMARY_ModelB.txt\n")