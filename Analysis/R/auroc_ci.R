# AUROC with confidence intervals for Table 3.7
# Answers the "we can do this in R? just tell me what to run it on" comment.
#
# WHAT TO RUN IT ON
# The same data frame that produced the current Table 3.7: one row per predicted
# structure for the native streptogramin dataset (n = 900), with
#   - a logical/0-1 column marking whether that structure met all three
#     attack-geometry criteria   (the current table's positive class)
#   - one column per confidence metric: nucleophile pLDDT, pTM, ipTM, ranking score
#   - the domain each structure came from, for the within-domain rows
#
# If your existing script writes stats/Table_3.6.csv (the AUROC table), the input
# is whatever data frame it computes roc() on. Point `df` at that.

library(pROC)
library(dplyr)

df <- readr::read_csv("/path/to/your/Experimental/Bioinformatics/prisnucleophile_distances_long.csv")

# --- adjust these four names to match your columns -------------------------
outcome_col <- "near_attack"        # logical or 0/1
domain_col  <- "te_domain"
metrics     <- c(nucleophile_pLDDT = "nuc_plddt",
                 pTM               = "ptm",
                 ipTM              = "iptm",
                 ranking_score     = "ranking_score")
# ---------------------------------------------------------------------------

auroc_ci <- function(data, metric_col, label) {
  # returns AUROC with a 95% DeLong interval, or NA if one class is absent
  y <- data[[outcome_col]]
  if (length(unique(y[!is.na(y)])) < 2) {
    return(tibble::tibble(metric = label, n = nrow(data), n_pos = sum(y, na.rm = TRUE),
                          auc = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  r  <- pROC::roc(response = y, predictor = data[[metric_col]],
                  quiet = TRUE, direction = "auto")
  ci <- pROC::ci.auc(r, method = "delong", conf.level = 0.95)
  tibble::tibble(metric = label, n = nrow(data), n_pos = sum(y, na.rm = TRUE),
                 auc = as.numeric(ci[2]), lower = as.numeric(ci[1]), upper = as.numeric(ci[3]))
}

# pooled across all six domains
pooled <- purrr::map2_dfr(metrics, names(metrics),
                          ~ auroc_ci(df, .x, .y)) |>
  mutate(stratum = "Pooled, six domains", .before = 1)

# within each domain separately
within <- df |>
  group_split(.data[[domain_col]]) |>
  purrr::map_dfr(function(sub) {
    purrr::map2_dfr(metrics, names(metrics), ~ auroc_ci(sub, .x, .y)) |>
      mutate(stratum = unique(sub[[domain_col]]), .before = 1)
  })

results <- bind_rows(pooled, within) |>
  mutate(across(c(auc, lower, upper), ~ round(.x, 3)),
         reported = ifelse(is.na(auc), "n/a",
                           sprintf("%.2f (%.2f\u2013%.2f)", auc, lower, upper)))

print(results, n = Inf)
readr::write_csv(results, "stats/Table_3.7_auroc_ci.csv")


# ---------------------------------------------------------------------------
# WHAT THE ANSWER DECIDES
#
# The thesis currently says every metric "fell to chance" within a single
# domain. Three did. Nucleophile pLDDT sits at 0.61 and 0.59, which is the
# number this script puts an interval around.
#
#   If the interval for nucleophile pLDDT INCLUDES 0.50
#       -> "indistinguishable from chance" is defensible; keep the wording but
#          report the interval so the reader can see why.
#
#   If the interval EXCLUDES 0.50
#       -> reword to: pooled 0.91 falls to 0.61 and 0.59 within domain, while
#          pTM, ipTM and ranking score fall to chance. This is the stronger
#          result anyway, and it is the honest one.
#
# Four places need the same treatment: the Abstract, Section 3.6, Section 4.4
# and the Conclusions.
#
# Note on independence: structures sampled from one seed are not independent,
# which Section 4.5 already acknowledges. If you want the conservative version,
# run the same function on one structure per seed (e.g. the top-ranked model)
# and report that interval instead; it will be wider but not open to the
# objection.
# ---------------------------------------------------------------------------
