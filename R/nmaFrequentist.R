# nmaFrequentist.R
# Frequentist random-effects network meta-analysis of smoking cessation
# interventions (Hasselblad 1998; 24 trials, 4 interventions) using netmeta.
#
# Run from the repo root:  Rscript R/nmaFrequentist.R

library(netmeta)

out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Data 
# One row/trial, wide format
# Outcome: successful cessation at 6-12 months
data(smokingcessation)

trt_labels <- c(
  A = "No intervention",
  B = "Self-help",
  C = "Individual counselling",
  D = "Group counselling"
)

smk <- smokingcessation
for (v in c("treat1", "treat2", "treat3")) {
  smk[[v]] <- unname(trt_labels[as.character(smk[[v]])])
}
smk$study <- paste("Study", seq_len(nrow(smk)))

# 2. Arm-level counts to pairwise contrasts (log odds ratios)
# pairwise() handles the two three-arm trials, and netmeta() then accounts for
# the correlation between contrasts from the same trial
pw <- pairwise(
  treat   = list(treat1, treat2, treat3),
  event   = list(event1, event2, event3),
  n       = list(n1, n2, n3),
  studlab = study,
  data    = smk,
  sm      = "OR"
)

# 3. Fit NMA 
nma <- netmeta(
  pw,
  common          = FALSE,
  random          = TRUE,
  reference.group = "No intervention",
  small.values    = "undesirable",   # cessation is good
  method.tau      = "REML"
)

sink(file.path(out_dir, "nma_summary.txt"))
print(summary(nma))
sink()

# 4. Network plot 
png(file.path(out_dir, "fig1_network.png"), width = 2000, height = 1600, res = 300)
netgraph(
  nma,
  plastic          = FALSE,
  thickness        = "number.of.studies",
  number.of.studies = TRUE,
  pos.number.of.studies = 0.35,
  points           = TRUE,
  cex.points       = 3,
  col.points       = "grey20",
  offset           = 0.05,
  col              = "grey60",
  multiarm         = FALSE
)
dev.off()

# 5. Relative effects versus the reference
png(file.path(out_dir, "fig2_forest.png"), width = 2400, height = 1100, res = 300)
forest(
  nma,
  sortvar    = -TE,
  xlim       = c(0.5, 8),
  at         = c(0.5, 1, 2, 4, 8),
  smlab      = "Odds ratio of cessation\nvs. no intervention",
  label.left  = "Favours no intervention",
  label.right = "Favours intervention",
  leftcols   = c("studlab", "k"),
  leftlabs   = c("Intervention", "Direct\ncomparisons")
)
dev.off()

# 6. All pairwise ORs and random effects
league <- netleague(nma, common = FALSE, digits = 2, bracket = "(", separator = " to ")
write.csv(league$random, file.path(out_dir, "league_table.csv"), row.names = FALSE)

# 7. Treatment ranking (P-scores)
ranking <- netrank(nma, small.values = "undesirable")
pscores <- data.frame(
  treatment = names(ranking$ranking.random),
  p_score   = round(unname(ranking$ranking.random), 3)
)
pscores <- pscores[order(-pscores$p_score), ]
write.csv(pscores, file.path(out_dir, "p_scores.csv"), row.names = FALSE)

# 8. Heterogeneity, inconsistencies
# 8a. Global: design-by-treatment interaction, Q decomposition
decomp <- decomp.design(nma)
sink(file.path(out_dir, "inconsistency_global.txt"))
print(decomp)
sink()

# 8b. Direct vs indirect evidence for each comparison
split <- netsplit(nma)
sink(file.path(out_dir, "inconsistency_nodesplit.txt"))
print(split)
sink()

png(file.path(out_dir, "fig3_nodesplit.png"), width = 2400, height = 2400, res = 300)
forest(split, show = "both")
dev.off()

# 9. Small-study effects: comparison-adjusted funnel plot
png(file.path(out_dir, "fig4_funnel.png"), width = 2200, height = 1800, res = 300)
funnel(
  nma,
  order = c("No intervention", "Self-help",
            "Individual counselling", "Group counselling"),
  pch   = 19,
  col   = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"),
  method.bias = "Egger",
  pos.legend  = "bottomright",
  digits.pval = 2
)
dev.off()

# 10. Tidy results table for the README and for Bayesian script
ref  <- "No intervention"
trts <- setdiff(nma$trts, ref)
results <- data.frame(
  treatment = trts,
  or        = exp(nma$TE.random[trts, ref]),
  lower     = exp(nma$lower.random[trts, ref]),
  upper     = exp(nma$upper.random[trts, ref])
)
results <- merge(results, pscores, by = "treatment")
results <- results[order(-results$or), ]
results[, c("or", "lower", "upper")] <- round(results[, c("or", "lower", "upper")], 2)
write.csv(results, file.path(out_dir, "results_frequentist.csv"), row.names = FALSE)

cat("\nRandom-effects NMA, odds ratio of cessation vs no intervention\n")
print(results, row.names = FALSE)
cat(sprintf("\ntau^2 = %.3f, I^2 = %.0f%%\n", nma$tau2, 100 * nma$I2))
# The common-effects Q is inflated by heterogeneity, so report both versions.
inc_re <- decomp$Q.inc.random
cat(sprintf("Inconsistency between designs, common effects: Q = %.2f, df = %d, p = %.3f\n",
            nma$Q.inconsistency, nma$df.Q.inconsistency, nma$pval.Q.inconsistency))
cat(sprintf("Inconsistency between designs, random effects: Q = %.2f, df = %d, p = %.3f\n",
            inc_re$Q, inc_re$df, inc_re$pval))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info_01.txt"))

cat("\nOutputs written to:", normalizePath(out_dir), "\n")