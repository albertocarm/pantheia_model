# Reproduces the figures of the PANTHEIA-SIRI validation paper from the fitted
# model and the de-identified analysis datasets shipped with this package:
#   Figure 1     - internal-validation OS by risk group (KM) + calibration (6/12 mo)
#   Supp Fig S1  - variable clustering of the model predictors (derivation)
#   Supp Fig S3  - probability density of death times (validation cohort)
# Run from the package root:  Rscript inst/reproduce_figures.R
suppressMessages({ library(survival); library(splines); library(ggplot2) })

find_dir <- function(sub) {
  p <- system.file(sub, package = "pantheia_model"); if (nzchar(p)) return(p)
  for (cand in c(file.path("inst", sub), sub)) if (dir.exists(cand)) return(cand); stop("Run from package root.")
}
app <- find_dir("app"); extd <- find_dir("extdata"); outdir <- "figures"; dir.create(outdir, showWarnings = FALSE)
m <- readRDS(file.path(app, "FINAL_OS_RUBIN.rds")); knots_logsiri <- m$knots; scale_param <- m$scale
rhs <- ~ diam3 + ns(logsiri, knots = knots_logsiri[2], Boundary.knots = knots_logsiri[c(1,3)]) +
  regimen_cat + ecog_cat_3 + CACS + regimen_cat:logsiri
lp_of <- function(x){ for(cc in c("diam3","regimen_cat","ecog_cat_3","CACS")) x[[cc]]<-factor(x[[cc]],levels=m$xlevels[[cc]])
  x<-x[complete.cases(x[,c("logsiri","diam3","regimen_cat","ecog_cat_3","CACS","os_time","os_event")]),]
  X<-model.matrix(rhs,x,xlev=m$xlevels); cf<-m$coef[colnames(X)]; cf[is.na(cf)]<-0; x$lp<-as.numeric(X%*%cf); x }

## ===== Figure 1: validation KM (A,B) + calibration (C,D) =====
v <- lp_of(read.csv(file.path(extd, "validation_os.csv")))
df <- data.frame(time = v$os_time, event = v$os_event, lp = v$lp)
wsurv <- function(t, lp, s) exp(-(t/exp(lp))^(1/s))
mk_km <- function(grp, labs, pal, title){
  d2 <- df; d2$g <- grp; d2$g <- factor(d2$g, levels = labs)
  fit <- survfit(Surv(time, event) ~ g, data = d2)
  survminer::ggsurvplot(fit, data = d2, pval = TRUE, pval.size = 3.5, conf.int = TRUE, conf.int.alpha = 0.15,
    risk.table = TRUE, risk.table.height = 0.25, risk.table.fontsize = 3, risk.table.y.text = FALSE,
    xlab = "Time (months)", ylab = "Overall survival", title = title, legend.title = "", legend.labs = labs,
    palette = unname(pal), ggtheme = ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=11)))$plot
}
col2 <- c("Low risk"="#2E86C1","High risk"="#E74C3C")
col3 <- c("Low risk (T3)"="#2E86C1","Intermediate (T2)"="#F39C12","High risk (T1)"="#E74C3C")
pA <- mk_km(ifelse(df$lp <= median(df$lp),"High risk","Low risk"), c("Low risk","High risk"), col2, "A.  OS by risk group (median)")
q3 <- quantile(df$lp, c(1/3,2/3))
pB <- mk_km(as.character(cut(df$lp, c(-Inf,q3[1],q3[2],Inf), labels=c("High risk (T1)","Intermediate (T2)","Low risk (T3)"))),
            c("Low risk (T3)","Intermediate (T2)","High risk (T1)"), col3, "B.  OS by risk group (tertiles)")
calib <- function(t_cal, lab, col){
  df$pred <- wsurv(t_cal, df$lp, scale_param)
  df$grp <- cut(df$pred, quantile(df$pred, seq(0,1,length.out=5), na.rm=TRUE), include.lowest=TRUE, labels=FALSE)
  cal <- do.call(rbind, lapply(sort(unique(df$grp)), function(g){ s<-df[df$grp==g,]; km<-survfit(Surv(time,event)~1,data=s)
    i<-max(which(km$time<=t_cal),0); os<-if(i>0)km$surv[i] else 1; lo<-if(i>0)max(km$lower[i],0) else NA; hi<-if(i>0)min(km$upper[i],1) else NA
    data.frame(n=nrow(s), pred=mean(s$pred), obs=os, lo=lo, hi=hi) }))
  ggplot(cal, aes(pred, obs)) + geom_abline(linetype="dashed", colour="grey50") +
    geom_errorbar(aes(ymin=lo, ymax=hi), width=.03, colour=col, linewidth=.6) + geom_point(size=3.5, colour=col) +
    geom_text(aes(label=paste0("n=",n)), hjust=-.3, vjust=-.5, size=3) + coord_equal(xlim=c(0,1), ylim=c(0,1)) +
    labs(x="Predicted survival", y="Observed survival (KM)", title=paste0(lab,".  Calibration at ",t_cal," months")) +
    theme_classic(base_size=10) + theme(plot.title=element_text(face="bold", size=11))
}
fig1 <- cowplot::plot_grid(cowplot::plot_grid(pA,pB,ncol=2), cowplot::plot_grid(calib(6,"C","#2E86C1"),calib(12,"D","#E74C3C"),ncol=2), nrow=2)
ggsave(file.path(outdir,"Figure_1.pdf"), fig1, width=7.96, height=5.2, dpi=300)

## ===== Supp Fig S3: density of death times (validation) =====
deaths <- sort(df$time[df$event==1]); nd <- length(deaths); med <- median(deaths)
sw <- function(fr){k<-ceiling(fr*nd);best<-Inf;bi<-1;for(i in 1:(nd-k+1)){w<-deaths[i+k-1]-deaths[i];if(w<best){best<-w;bi<-i}};c(lo=deaths[bi],hi=deaths[bi+k-1],width=best)}
w50<-sw(.5); w75<-sw(.75); dens<-density(deaths,from=0,to=max(deaths)+2); dfd<-data.frame(t=dens$x,d=dens$y); dd<-data.frame(time=deaths)
ps1 <- ggplot() +
  annotate("rect",xmin=w75["lo"],xmax=w75["hi"],ymin=0,ymax=Inf,fill="#F5B7B1",alpha=.30) +
  annotate("rect",xmin=w50["lo"],xmax=w50["hi"],ymin=0,ymax=Inf,fill="#E6B0AA",alpha=.50) +
  geom_area(data=dfd,aes(t,d),fill="#C0392B",alpha=.35) + geom_line(data=dfd,aes(t,d),colour="#922B21",linewidth=1.1) +
  geom_rug(data=dd,aes(x=time),sides="b",colour="#7B241C",alpha=.6) + geom_vline(xintercept=med,colour="#1F618D",linewidth=.9) +
  geom_vline(xintercept=c(6,12,18),linetype="dashed",colour="grey45") +
  labs(x="Time to death (months)", y="Probability density",
       title=sprintf("Validation cohort: density of death times (n=%d deaths)", nd)) +
  theme_classic(base_size=11) + theme(plot.title=element_text(face="bold", size=11))
ggsave(file.path(outdir,"Supp_Figure_S3.pdf"), ps1, width=8.5, height=4.6)

## ===== Supp Fig S1: variable clustering (derivation) =====
imp <- read.csv(file.path(extd,"derivation_os_imputed.csv")); d1 <- imp[imp$imp==1,]
for(cc in c("diam3","regimen_cat","ecog_cat_3","CACS")) d1[[cc]] <- factor(d1[[cc]])
vc <- Hmisc::varclus(~ diam3 + logsiri + regimen_cat + ecog_cat_3 + CACS, data = d1)
pdf(file.path(outdir,"Supp_Figure_S1.pdf"), width=7, height=5); plot(vc); title("Variable clustering of the model predictors (derivation)"); dev.off()

cat("Figures written to", outdir, ":\n  Figure_1.pdf, Supp_Figure_S1.pdf, Supp_Figure_S3.pdf\n")
