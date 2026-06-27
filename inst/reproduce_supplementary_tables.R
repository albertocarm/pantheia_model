# Reproduces Supplementary Table S1 (collinearity diagnostics: GVIF and
# redundancy analysis) and Supplementary Table S2 (pairwise associations among
# the cachexia, ECOG and SIRI predictors) of the PANTHEIA-SIRI validation paper,
# from the de-identified derivation analysis data shipped with this package.
# Run from the package root:  Rscript inst/reproduce_supplementary_tables.R
suppressMessages({ library(survival); library(car); library(Hmisc) })

find_dir <- function(sub) {
  p <- system.file(sub, package = "pantheia_model"); if (nzchar(p)) return(p)
  for (cand in c(file.path("inst", sub), sub)) if (dir.exists(cand)) return(cand); stop("Run from package root.")
}
extd <- find_dir("extdata")
imp <- read.csv(file.path(extd, "derivation_os_imputed.csv"))
d <- imp[imp$imp == 1, ]                         # one multiply-imputed dataset (n = 593)
d <- d[d$os_time > 0, ]
for (cc in c("diam3", "regimen_cat", "ecog_cat_3", "CACS")) d[[cc]] <- factor(d[[cc]])

cat("=== Supplementary Table S1: collinearity diagnostics (derivation) ===\n")
m <- lm(log(os_time) ~ diam3 + logsiri + regimen_cat + ecog_cat_3 + CACS, data = d)
cat("(a) Generalized variance inflation factors:\n"); print(round(car::vif(m), 3))
cat("\n(b) Redundancy analysis (Hmisc::redun, R2 = 0.9):\n")
print(redun(~ diam3 + logsiri + regimen_cat + ecog_cat_3 + CACS, data = d, r2 = 0.9, nk = 4))

cat("\n=== Supplementary Table S2: pairwise associations (cachexia / ECOG / SIRI) ===\n")
cramerV <- function(x, y) { tb <- table(x, y); chi <- suppressWarnings(chisq.test(tb)); k <- min(dim(tb))
  c(V = sqrt(as.numeric(chi$statistic) / (sum(tb) * (k - 1))), p = chi$p.value) }
ce <- cramerV(d$CACS, d$ecog_cat_3)
cat(sprintf("  Cachexia x ECOG (Cramer's V)        = %.3f   (p = %.3f)\n", ce["V"], ce["p"]))
cat(sprintf("  log-SIRI by cachexia (Wilcoxon)     :         p = %.3f\n", wilcox.test(logsiri ~ CACS, data = d)$p.value))
cat(sprintf("  log-SIRI by ECOG (Kruskal-Wallis)   :         p = %.3f\n", kruskal.test(logsiri ~ ecog_cat_3, data = d)$p.value))
cat(sprintf("  log-SIRI by regimen (Kruskal-Wallis):         p = %.3f\n", kruskal.test(logsiri ~ regimen_cat, data = d)$p.value))
cat(sprintf("  log-SIRI by tumour burden (Kruskal) :         p = %.3f\n", kruskal.test(logsiri ~ diam3, data = d)$p.value))
