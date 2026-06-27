# Reproduces Table 2 (OS model coefficients) and Table 3 (discrimination and
# calibration) of the PANTHEIA-SIRI validation paper.
#
# Derivation performance is computed over the 10 multiply-imputed datasets of the
# survival-analysis population (n = 593) and pooled with Rubin's rules -- i.e. the
# missing tumour-burden split is imputed, NOT dropped. External validation uses the
# independent cohort (n = 62). Run from the package root:
#   Rscript inst/reproduce_tables.R
suppressMessages({ library(survival); library(splines) })

find_dir <- function(sub) {
  p <- system.file(sub, package = "pantheia_model")
  if (nzchar(p)) return(p)
  for (cand in c(file.path("inst", sub), sub)) if (dir.exists(cand)) return(cand)
  stop("Could not locate '", sub, "'. Run from the package root.")
}
app  <- find_dir("app"); extd <- find_dir("extdata")
m_os <- readRDS(file.path(app, "FINAL_OS_RUBIN.rds"))
knots_logsiri <- m_os$knots

## ---- TABLE 2: OS model coefficients (Weibull AFT) ----
se <- sqrt(diag(m_os$vcov))[seq_along(m_os$coef)]
cat("=== TABLE 2 (OS model coefficients) ===\n")
print(data.frame(coef = round(m_os$coef, 3), SE = round(se, 3), TR = round(exp(m_os$coef), 2)))
cat("Weibull scale (sigma) =", round(m_os$scale, 3), "\n\n")

## ---- helpers ----
rhs <- ~ diam3 + ns(logsiri, knots = knots_logsiri[2], Boundary.knots = knots_logsiri[c(1, 3)]) +
  regimen_cat + ecog_cat_3 + CACS + regimen_cat:logsiri
prep <- function(x) {
  for (cc in c("diam3", "regimen_cat", "ecog_cat_3", "CACS")) x[[cc]] <- factor(x[[cc]], levels = m_os$xlevels[[cc]])
  x[complete.cases(x[, c("logsiri","diam3","regimen_cat","ecog_cat_3","CACS","os_time","os_event")]), ]
}
lp_of <- function(x) { X <- model.matrix(rhs, x, xlev = m_os$xlevels); cf <- m_os$coef[colnames(X)]; cf[is.na(cf)] <- 0; as.numeric(X %*% cf) }
wsurv <- function(t, lp, s) exp(-(t / exp(lp))^(1 / s))
brier_ipcw <- function(time, event, ps, te) {
  kmc <- survfit(Surv(time, 1 - event) ~ 1); G <- function(ti){i<-max(which(kmc$time<=ti),0);if(i==0)1 else kmc$surv[i]}
  Gt <- G(te); bs <- 0
  for (i in seq_along(time)) { if (time[i]<=te && event[i]==1){g<-G(time[i]);if(g>.01)bs<-bs+ps[i]^2/g}
    else if (time[i]>te){if(Gt>.01)bs<-bs+(1-ps[i])^2/Gt} }
  bs / length(time)
}
km_at <- function(time, event, te){km<-survfit(Surv(time,event)~1);i<-max(which(km$time<=te),0);if(i>0)km$surv[i] else 1}
metrics <- function(d) {
  lp <- lp_of(d); sc <- m_os$scale; cc <- concordance(Surv(d$os_time, d$os_event) ~ lp)
  out <- c(C = cc$concordance, Cvar = cc$var)
  for (te in c(6, 12)) {
    ps <- wsurv(te, lp, sc); b <- brier_ipcw(d$os_time, d$os_event, ps, te)
    s0 <- km_at(d$os_time, d$os_event, te); b0 <- brier_ipcw(d$os_time, d$os_event, rep(s0, nrow(d)), te)
    out[paste0("Brier", te)] <- b; out[paste0("IPA", te)] <- 1 - b / b0
  }
  out
}

cat("=== TABLE 3 (discrimination + calibration) ===\n")
## Derivation: pool over the 10 imputed datasets (Rubin's rules)
imp <- read.csv(file.path(extd, "derivation_os_imputed.csv"))
ms  <- t(sapply(sort(unique(imp$imp)), function(k) metrics(prep(imp[imp$imp == k, ]))))
m_  <- colMeans(ms); within <- mean(ms[, "Cvar"]); between <- var(ms[, "C"]); mImp <- nrow(ms)
Ctot <- within + (1 + 1/mImp) * between; Clo <- m_["C"] - 1.96*sqrt(Ctot); Chi <- m_["C"] + 1.96*sqrt(Ctot)
n1 <- nrow(prep(imp[imp$imp == 1, ])); ev1 <- sum(prep(imp[imp$imp == 1, ])$os_event == 1)
cat(sprintf("Derivation (n=%d, pooled over %d imputations)  events=%d (%.1f%%)\n", n1, mImp, ev1, 100*ev1/n1))
cat(sprintf("   C-index = %.3f (%.3f-%.3f)   6mo Brier/IPA = %.3f/%.3f   12mo = %.3f/%.3f\n",
            m_["C"], Clo, Chi, m_["Brier6"], m_["IPA6"], m_["Brier12"], m_["IPA12"]))
## Validation: independent cohort (single dataset)
v <- prep(read.csv(file.path(extd, "validation_os.csv"))); mv <- metrics(v)
cat(sprintf("Validation (n=%d)  events=%d (%.1f%%)\n", nrow(v), sum(v$os_event==1), 100*mean(v$os_event==1)))
cat(sprintf("   C-index = %.3f (%.3f-%.3f)   6mo Brier/IPA = %.3f/%.3f   12mo = %.3f/%.3f\n",
            mv["C"], mv["C"]-1.96*sqrt(mv["Cvar"]), mv["C"]+1.96*sqrt(mv["Cvar"]), mv["Brier6"], mv["IPA6"], mv["Brier12"], mv["IPA12"]))
