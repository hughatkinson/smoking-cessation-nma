# nmaBayesian.R
# Bayesian random-effects NMA of the same smoking cessation network using
# multinma (Stan), with model comparison and inconsistency checks.
#
# Run after  nmaFrequentist.R:  Rscript R/nmaBayesian.R

library(multinma)
library(ggplot2)

options(mc.cores = parallel::detectCores())
out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Network
# multinma::smoking, same Hasselblad (1998) dataset
net <- set_agd_arm(
  smoking,
  study   = studyn,
  trt     = trtc,
  r       = r,
  n       = n,
  trt_ref = "No intervention"
)
print(net)

# 2. Models, Binomial likelihood, logit link.
fit_nma <- function(...) {
  nma(
    net,
    prior_intercept = normal(scale = 100),
    prior_trt       = normal(scale = 100),
    seed            = 2026,
    ...
  )
}

fit_fe  <- fit_nma(trt_effects = "fixed")
fit_re  <- fit_nma(trt_effects = "random", prior_het = half_normal(scale = 5))
fit_ume <- fit_nma(trt_effects = "random", prior_het = half_normal(scale = 5),
                   consistency = "ume")

# 3. Convergence
sink(file.path(out_dir, "bayes_model_summary.txt"))
print(fit_re, pars = c("d", "tau"))
sink()

ggsave(file.path(out_dir, "fig5_prior_posterior.png"),
       plot_prior_posterior(fit_re, prior = c("trt", "het")),
       width = 7, height = 4, dpi = 300)

# 4. Model fit (fixed vs random effects, consistency vs inconsistency)
dic_fe  <- dic(fit_fe)
dic_re  <- dic(fit_re)
dic_ume <- dic(fit_ume)

model_fit <- data.frame(
  model     = c("Fixed effects", "Random effects", "Random effects, UME (inconsistency)"),
  resdev    = c(dic_fe$resdev, dic_re$resdev, dic_ume$resdev),
  pD        = c(dic_fe$pd,     dic_re$pd,     dic_ume$pd),
  DIC       = c(dic_fe$dic,    dic_re$dic,    dic_ume$dic)
)
model_fit[, -1] <- round(model_fit[, -1], 1)
write.csv(model_fit, file.path(out_dir, "bayes_model_fit.csv"), row.names = FALSE)
print(model_fit, row.names = FALSE)

# Dev-dev plot
ggsave(file.path(out_dir, "fig6_dev_dev.png"),
       plot(dic_re, dic_ume, show_uncertainty = FALSE) +
         labs(x = "Residual deviance, consistency model",
              y = "Residual deviance, UME model"),
       width = 5, height = 5, dpi = 300)

# 5. Relative effects vs no intervention
rel <- as.data.frame(relative_effects(fit_re))
rel$treatment <- sub("^d\\[(.*)\\]$", "\\1", rel$parameter)

bayes <- data.frame(
  treatment   = rel$treatment,
  or_bayes    = round(exp(rel[["50%"]]), 2),
  lower_bayes = round(exp(rel[["2.5%"]]), 2),
  upper_bayes = round(exp(rel[["97.5%"]]), 2)
)

or_ticks <- c(0.5, 1, 2, 4, 8)
ggsave(file.path(out_dir, "fig7_forest_bayes.png"),
       plot(relative_effects(fit_re), ref_line = 0) +
         scale_x_continuous(breaks = log(or_ticks), labels = or_ticks) +
         labs(x = "Odds ratio of cessation vs no intervention (log scale)"),
       width = 7, height = 3, dpi = 300)

# 6. Ranking 
ranks <- as.data.frame(posterior_ranks(fit_re, lower_better = FALSE, sucra = TRUE))
ranks$treatment <- sub("^rank\\[(.*)\\]$", "\\1", ranks$parameter)
bayes <- merge(bayes, ranks[, c("treatment", "sucra")], by = "treatment")
bayes$sucra <- round(bayes$sucra, 3)

ggsave(file.path(out_dir, "fig8_rank_probs.png"),
       plot(posterior_rank_probs(fit_re, lower_better = FALSE)),
       width = 8, height = 3, dpi = 300)

# 7. Heterogeneity 
tau <- as.data.frame(summary(fit_re, pars = "tau"))
cat(sprintf("\nBetween-study SD (tau): %.2f (95%% CrI %.2f to %.2f)\n",
            tau[["50%"]], tau[["2.5%"]], tau[["97.5%"]]))

# 8. With the frequentist results
write.csv(bayes, file.path(out_dir, "results_bayesian.csv"), row.names = FALSE)

if (file.exists(file.path(out_dir, "results_frequentist.csv"))) {
  freq <- read.csv(file.path(out_dir, "results_frequentist.csv"))
  names(freq)[-1] <- c("or_freq", "lower_freq", "upper_freq", "p_score")
  comparison <- merge(freq, bayes, by = "treatment")
  comparison <- comparison[order(-comparison$or_freq), ]
  write.csv(comparison, file.path(out_dir, "results_comparison.csv"), row.names = FALSE)
  cat("\nFrequentist (netmeta) vs Bayesian (multinma), OR vs no intervention\n")
  print(comparison, row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info_02.txt"))

cat("\nOutputs written to:", normalizePath(out_dir), "\n")