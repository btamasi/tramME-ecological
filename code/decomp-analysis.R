## ----report-info--------------------------------------------------------------
##' **Flexible regression for correlated and censored ecological data with
##' mixed-effects additive transformation models**
##'
##' Balint Tamasi and Torsten Hothorn
##'
##' This document reproduces the results of the main text of the article


## ----packages, message=FALSE, warning=FALSE, results='hide'-------------------
##' Check packages, install if necessary.
##' Check version of tramME.
pkgs <- c("tramME", "xtable", "multcomp", "survival", "mgcv")
inst <- pkgs %in% installed.packages()
if (any(!inst)) install.packages(pkgs[!inst])
if (packageDescription("tramME")$Version < "1.0.3") {
  install.packages("tramME")
}
pkgs_out <- lapply(pkgs, library, character.only = TRUE, quietly = TRUE)


## ----session-info, eval=do_report, include=do_report--------------------------
##' Additional information about the session.
sessionInfo()
date()


## ----load-rat, message=FALSE--------------------------------------------------
##' The data is available from the Supplementary Information
##' of Englmeier et al. (2022)
##' https://doi.org/10.1007/s10021-022-00764-7
carrion <- read.csv("../data/carrion.csv")
carrion$Time <- with(carrion, Surv(time1, time2, type = "interval2"))
carrion$Time_rc <- with(carrion,
  Surv(ifelse(is.finite(time2), time2, time1), event = is.finite(time2)))
carrion$Insects <- factor(carrion$Insects, labels = c("no", "yes"))
carrion$Habitat <- factor(carrion$Habitat)
carrion$Habitat <- relevel(carrion$Habitat, ref = "forest")
carrion$Landscape <- factor(carrion$Landscape)
carrion$Landscape <- relevel(carrion$Landscape, ref = "seminatural")
carrion$PlotID <- factor(carrion$PlotID)


## ----rat1, cache=TRUE---------------------------------------------------------
##' The proportional hazards mixed-effects additive model
mc1 <- CoxphME(Time ~ Insects + Habitat + Landscape
               + s(Temperature, k = 20) + s(Elevation100, k = 20)
               + (1 | PlotID), data = carrion,
               log_first = TRUE, order = 6)
stopifnot(mc1$opt$convergence == 0)


## ----rat1-summary, eval=do_report, include=do_report--------------------------
summary(mc1)


## ----rat-nonprop, cache=TRUE--------------------------------------------------
##' Mixed effects additive model with time-varying
##' effect of instects
mc_np <- CoxphME(Time | Insects ~ Habitat + Landscape
                 + s(Temperature, k = 20) + s(Elevation100, k = 20)
                 + (1 | PlotID), data = carrion,
                 log_first = TRUE, order = 6)
stopifnot(mc_np$opt$convergence == 0)


## ----rat1-smooth, dependson="rat1"--------------------------------------------
##' Evaluate smooth terms
sm <- smooth_terms(mc1)


## ----rat-ph-vs-tve, results="hide", warning=FALSE, dependson=c("rat1", "rat-nonprop")----
##' Evaluate the effect of insect in the proportional and non-proportional model
cf <- coef(mc_np, with_baseline = TRUE)
vc <- vcov(mc_np, pargroup = "baseline")
cf <- cf[grep("Insectsyes", names(cf), fixed = TRUE)]
idx <- grep("Insectsyes", colnames(vc), fixed = TRUE)
vc <- vc[idx, idx]
ns <- 200
nd <- model.frame(mc_np)[rep(1, ns), ]
nd[[variable.names(mc_np, "response")]] <- seq(1, 100, length.out = ns)
X <- model.matrix(mc_np, data = nd, type = "Y", simplify = TRUE)$Ye
idx <- grep("Insectsyes", colnames(X), fixed = TRUE)
X <- X[, idx]
ci <- confint(multcomp::glht(multcomp::parm(cf, vc), linfct = X),
              calpha = multcomp::univariate_calpha(), level = 0.95)$confint
## use multcomp::adjusted_calpha() for multiplicity adjustments
ci2 <- confint(mc1, parm = "Insectsyes", estimate = TRUE)
ci2 <- ci2[rep(1, ns), ]


## ----fig-rat-effects, fig.width=7.5, fig.height=7.5, results="hide", out.width="\\textwidth"----
##' Plots
op <- par(mfrow = c(2, 2))
## -- Plot 1
par(cex = 0.9, mar = c(4.5, 4.1, 2, 1), las = 1)
plot(sm[1], col = 1, fill = grey(0.5, 0.2), lwd = 2, panel.first = {
  grid()
  abline(h = 0, col = "grey", lty = 1, lwd = 2)
  rug(mc1$data$Temperature, col = grey(0.5, 0.5), quiet = TRUE)
},
xlab = "Average temperature $(^{\\circ}\\text{C})$",
ylab = "$\\widehat{f}(\\texttt{temperature})$")
blx <- grconvertX(0.10, from = "nfc", to = "user")
bly <- grconvertY(0.98, from = "nfc", to = "user")
text(blx, bly, labels = "(a)", xpd = TRUE, cex = 1.2)
## -- Plot 2
par(cex = 0.9, mar = c(4.5, 4.1, 2, 1), las = 1)
plot(sm[2], col = 1, fill = grey(0.5, 0.2), lwd = 2, panel.first = {
  grid()
  abline(h = 0, col = "grey", lty = 1, lwd = 2)
  rug(mc1$data$Elevation100, col = grey(0.5, 0.5), quiet = TRUE)
},
xlab = "Elevation (100 meters)", ylab = "$\\widehat{f}(\\texttt{elevation})$")
blx <- grconvertX(0.10, from = "nfc", to = "user")
bly <- grconvertY(0.98, from = "nfc", to = "user")
text(blx, bly, labels = "(b)", xpd = TRUE, cex = 1.2)
## -- Plot 3
par(cex = 0.9, mar = c(4.5, 4, 2, 1), las = 1)
plot(nd$Time, ci[, 1], type = "l", col = 1, lwd = 2, ylim = c(-1, 3),
     panel.first = {grid(); abline(h = 0, lwd = 2, col = "lightgrey")},
     ylab = "Log-hazard ratio", xlab = "")
title(xlab = "Time (days)", line = 2.3)
polygon(c(nd$Time, rev(nd$Time)), c(ci[, 2], rev(ci[, 3])), border = NA,
        col = grey(0.5, 0.2))
matlines(nd$Time, ci2, col = 1, lwd = c(1, 1, 2), lty = c(3, 3, 2))
legend("bottomright", c("Time-varying effect", "Proportional effect"),
       lty = 1:2, lwd = 2, bty = "n", cex = 1)
mx <- grconvertX(1, "npc", "user")
rr <- model.frame(mc_np)
op <- par(xpd = NA)
for (i in 1:nrow(rr)) {
  x0 <- rr$Time[i][, 1]
  x1 <- min(mx, rr$Time[i][, 2])
  if (x1 == 1) x1 <- mx
  rect(xleft = x0, xright = x1, ytop = grconvertY(0.04, "nfc", "user"),
       ybottom = grconvertY(0.0, "nfc", "user"), col = grey(0.3, 0.015),
       border = NA)
}
blx <- grconvertX(0.10, from = "nfc", to = "user")
bly <- grconvertY(0.98, from = "nfc", to = "user")
text(blx, bly, labels = "(c)", xpd = TRUE, cex = 1.2)
par(op)


sv <- survfit(Time ~ Insects, data = carrion)
par(cex = 0.9, mar = c(4.5, 4, 2, 1), las = 1)
plot(sv, lty = c(1, 2), lwd = 2, xlab = "",
     ylab = "Probability")
title(xlab = "Time (days)", line = 2.3)
grid()
legend("topright", c("Insect = no", "Insect = yes"),
       lty = c(1, 2), lwd = 2, bty = "n")
mx <- grconvertX(1, "npc", "user")
op <- par(xpd = NA)
for (i in 1:nrow(carrion)) {
  x0 <- carrion$Time[i][, 1]
  x1 <- min(mx, carrion$Time[i][, 2])
  if (x1 == 1) x1 <- mx
    rect(xleft = x0, xright = x1, ytop = grconvertY(0.04, "nfc", "user"),
         ybottom = grconvertY(0, "nfc", "user"), col = grey(0.3, 0.015),
         border = NA)
}
blx <- grconvertX(0.10, from = "nfc", to = "user")
bly <- grconvertY(0.98, from = "nfc", to = "user")
text(blx, bly, labels = "(d)", xpd = TRUE, cex = 1.2)
par(op)


## ----rat2, cache=TRUE---------------------------------------------------------
##' Proportional hazards model with linear effects for temperature and time
mc2 <- CoxphME(Time ~ Insects + Habitat + Landscape
               + Temperature + Elevation100
               + (1 | PlotID), data = carrion,
               log_first = TRUE, order = 6)
stopifnot(mc2$opt$convergence == 0)


## ----rat2-summary, eval=do_report, include=do_report--------------------------
summary(mc2)


## ----rats2-rc, cache=TRUE-----------------------------------------------------
##' Ignoring interval-censoring, same as reported by Englmeier et al. (2022)
carrion$time <- carrion$Time_rc[, 1]
carrion$event <- carrion$Time_rc[, 2]
mc2_rc <- gam(time ~ Insects + Habitat + Landscape + Temperature
               + Elevation100 + s(PlotID, bs = "re"), data = carrion,
               family = cox.ph(), weights = event)


## ----rats2-rc-summary, eval=do_report, include=do_report----------------------
summary(mc2_rc)


## ----rats-tbl, results="asis"-------------------------------------------------
b <- coef(mc2)
tau <- formatC(sqrt(varcov(mc2)[[1]][1, 1]), format = "f", digits = 2)
ll <- formatC(logLik(mc2), format = "f", digits = 2)
tbl <- formatC(b, format = "f", digits = 2)
se1 <- formatC(sqrt(diag(vcov(mc2, pargroup = "shift"))), format = "f", digits = 2)
tbl <- paste0("$", tbl, "$", " ($", se1, "$)")
tbl <- c(tbl, tau, ll)
hr <- formatC(exp(b), format = "f", digits = 2)
hr <- c(hr, NA, NA)
ci <- rbind(exp(confint(mc2, pargroup = "shift")))#,
ci <- formatC(ci, format = "f", digits = 2)
ci <- paste0("$", ci[, 1], "$", " -- ", "$", ci[, 2], "$")
ci <- c(ci, NA, NA)
pi <- 1/(1+exp(b))
pi <- formatC(pi, format = "f", digits = 2)
pi <- c(pi, NA, NA)
pici <- 1 / (1 + exp(confint(mc2, pargroup = "shift")))
pici <- formatC(pici, format = "f", digits = 2)
pici <- paste0("$", pici[, 2], "$", " -- ", "$", pici[, 1], "$")
pici <- c(pici, NA, NA)
b2 <- formatC(coef(mc2_rc)[1:8], format = "f", digits = 2)
se2 <- formatC(sqrt(diag(vcov(mc2_rc)[1:8, 1:8])), format = "f", digits = 2)
oo <- capture.output(tau2 <- gam.vcomp(mc2_rc)[1])
tau2 <- formatC(tau2, format = "f", digits = 2)
ll2 <- formatC(logLik(mc2_rc), format = "f", digits = 2)
tbl2 <- formatC(b2, format = "f", digits = 2)
tbl2 <- paste0("$", tbl2, "$", " ($", se2, "$)")
tbl2 <- c(tbl2, tau2, ll2)
tbl <- cbind(tbl, hr, ci, pi, pici, tbl2)
tbl <- tbl[c(1, 1, 2, 2:5, 5:10), ]
tbl[c(1, 3, 7), ] <- c(rep("$0$", 3), rep("$1$", 3), rep(NA, 3),
                       rep("$0.50$", 3), rep(NA, 3), rep("$0$", 3))
rownames(tbl) <- c("Insect = no", "Insect = yes", "Habitat = forest", "Habitat = arable field",
                   "Habitat = meadow", "Habitat = settlement",
                   "Landscape = seminatural", "Landscape = agricultural",
                   "Landscape = urban", "Average temperature $(^{\\circ}\\text{C})$",
                   "Elevation (100 m)",
                   "$\\widehat\\tau$",
                   "Log-likelihood")
add <- list()
add$pos <- list(0)
add$command <- paste(c("& \\multicolumn{5}{c}{Interval-censored} & ",
                       "\\multicolumn{1}{c}{Right-censored} \\\\\n",
                       "\\cmidrule(lr){2-6}\\cmidrule(lr){7-7}",
                       "& $\\widehat\\beta$ (SE) & HR & HR 95\\% CI ",
                       "& PI & PI 95\\% CI & $\\widehat\\beta$ (SE) \\\\\n"),
                     collapse = "")
print(xtable(tbl, align = "lrrrrrr"), floating = FALSE,
      sanitize.text.function = function(x){x},
      hline.after = c(-1, 0, nrow(tbl)-2, nrow(tbl)), add.to.row = add,
      NA.string = "", include.colnames = FALSE, booktabs = TRUE)


