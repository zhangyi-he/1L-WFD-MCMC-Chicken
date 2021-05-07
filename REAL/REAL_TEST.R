#' @title Inferring the timing and strength of natural selection and gene migration in the evolution of chicken from ancient DNA data
#' @author Wenyang Lyu, Xiaoyang Dai, Mark Beaumont, Feng Yu, Zhangyi He

#' TSHR and BCDO2 (the starting time of gene migration is given)

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggsci")
library("ggsci")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

#install.packages("readr")
library("readr")

#install.packages("xtable")
library("xtable")

#' call R functions
source("./REAL_RFUN.R")

################################################################################

#' Analysis of TSHR

# TSHR
TSHR <- read_csv("./Data/groupedTSHR.csv")
TSHR <- as.data.frame(TSHR)
TSHR$smp_siz <- TSHR$D.count + TSHR$A.count
TSHR <- cbind(aggregate(subset(TSHR, select = c(`Mean(yearsAD)`)), by = list(TSHR$Grouped.year), FUN = mean)[, -1],
              aggregate(subset(TSHR, select = c(smp_siz, D.count)), by = list(TSHR$Grouped.year), FUN = sum)[, -1])
TSHR[, 1] <- as.integer(round(TSHR[, 1] - max(TSHR[, 1])))
TSHR[, 2] <- as.integer(round(TSHR[, 2]))
TSHR[, 3] <- as.integer(round(TSHR[, 3]))
rownames(TSHR) <- NULL
colnames(TSHR) <- c("smp_gen", "smp_siz", "smp_cnt")

smp_gen <- TSHR$smp_gen
smp_siz <- TSHR$smp_siz
smp_cnt <- matrix(NA, nrow = 2, ncol = length(smp_gen))
smp_cnt[1, ] <- TSHR$smp_cnt
smp_cnt[2, length(smp_gen)] <- round(0.15 * smp_siz[length(smp_gen)])

############################################################

# Joint estimation of the selection coefficient, selection timing and migration rate

# The case of population size N = 26000
sel_cof <- 5e-03
dom_par <- 1e-00
mig_rat <- 0.15 / 250
pop_siz <- 26000
sel_gen <- -920
mig_gen <- -250
ext_frq <- 0.99
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH_1L(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_TSHR_PMMH1.rda")

load("./REAL_TSHR_PMMH1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
frq_pth_chn <- PMMH$frq_pth_chn[1, , ] + PMMH$frq_pth_chn[3, , ]

pdf(file = "./REAL_TSHR_PMMH1_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.5)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
frq_pth_chn <- frq_pth_chn[, brn_num:dim(frq_pth_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, (1:floor(dim(frq_pth_chn)[2] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn)
colnames(smp) <- c("selection coefficient", "selection timing", "migration rate")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]

# sel_cof_est <- mean(sel_cof_chn)
# sel_gen_est <- mean(sel_gen_chn)
# mig_rat_est <- mean(mig_rat_chn)
frq_pth_est <- rowMeans(frq_pth_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
frq_pth_hpd <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_chn)[1]) {
  frq_pth_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_chn[k, ]), prob = 0.95)
}

pdf(file = "./REAL_TSHR_PMMH1_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Migration rate",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 2], col = "blue", lty = 1, lwd = 2)
dev.off()

####################

# The case of population size N = 180000
sel_cof <- 5e-03
dom_par <- 1e-00
mig_rat <- 0.15 / 250
pop_siz <- 180000
sel_gen <- -920
mig_gen <- -250
ext_frq <- 0.99
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH_1L(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_TSHR_PMMH2.rda")

load("./REAL_TSHR_PMMH2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
frq_pth_chn <- PMMH$frq_pth_chn[1, , ] + PMMH$frq_pth_chn[3, , ]

pdf(file = "./REAL_TSHR_PMMH2_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.5)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
frq_pth_chn <- frq_pth_chn[, brn_num:dim(frq_pth_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, (1:floor(dim(frq_pth_chn)[2] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn)
colnames(smp) <- c("selection coefficient", "selection timing", "migration rate")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]

# sel_cof_est <- mean(sel_cof_chn)
# sel_gen_est <- mean(sel_gen_chn)
# mig_rat_est <- mean(mig_rat_chn)
frq_pth_est <- rowMeans(frq_pth_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
frq_pth_hpd <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_chn)[1]) {
  frq_pth_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_chn[k, ]), prob = 0.95)
}

pdf(file = "./REAL_TSHR_PMMH2_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Migration rate",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 2], col = "blue", lty = 1, lwd = 2)
dev.off()

####################

# The case of population size N = 460000
sel_cof <- 5e-03
dom_par <- 1e-00
mig_rat <- 0.15 / 250
pop_siz <- 460000
sel_gen <- -920
mig_gen <- -250
ext_frq <- 0.99
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH_1L(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_TSHR_PMMH3.rda")

load("./REAL_TSHR_PMMH3.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
frq_pth_chn <- PMMH$frq_pth_chn[1, , ] + PMMH$frq_pth_chn[3, , ]

pdf(file = "./REAL_TSHR_PMMH3_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.5)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
frq_pth_chn <- frq_pth_chn[, brn_num:dim(frq_pth_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, (1:floor(dim(frq_pth_chn)[2] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn)
colnames(smp) <- c("selection coefficient", "selection timing", "migration rate")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]

# sel_cof_est <- mean(sel_cof_chn)
# sel_gen_est <- mean(sel_gen_chn)
# mig_rat_est <- mean(mig_rat_chn)
frq_pth_est <- rowMeans(frq_pth_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
frq_pth_hpd <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_chn)[1]) {
  frq_pth_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_chn[k, ]), prob = 0.95)
}

pdf(file = "./REAL_TSHR_PMMH3_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Migration rate",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 2], col = "blue", lty = 1, lwd = 2)
dev.off()

################################################################################

#' Analysis of BCDO2

#' BCDO2
BCDO2 <- read_csv("./Data/groupedBCDO2.csv")
BCDO2 <- as.data.frame(BCDO2)
BCDO2$smp_siz <- BCDO2$D.count + BCDO2$A.count
BCDO2 <- cbind(aggregate(subset(BCDO2, select = c(`Mean(yearsAD)`)), by = list(BCDO2$Grouped.year), FUN = mean)[, -1],
              aggregate(subset(BCDO2, select = c(smp_siz, D.count)), by = list(BCDO2$Grouped.year), FUN = sum)[, -1])
BCDO2[, 1] <- as.integer(round(BCDO2[, 1] - max(BCDO2[, 1])))
BCDO2[, 2] <- as.integer(round(BCDO2[, 2]))
BCDO2[, 3] <- as.integer(round(BCDO2[, 3]))
rownames(BCDO2) <- NULL
colnames(BCDO2) <- c("smp_gen", "smp_siz", "smp_cnt")

smp_gen <- BCDO2$smp_gen
smp_siz <- BCDO2$smp_siz
smp_cnt <- matrix(NA, nrow = 2, ncol = length(smp_gen))
smp_cnt[1, ] <- BCDO2$smp_cnt
smp_cnt[2, length(smp_gen)] <- round(0.15 * smp_siz[length(smp_gen)])

############################################################

# Joint estimation of the selection coefficient, selection timing and migration rate

# The case of population size N = 26000
sel_cof <- 5e-03
dom_par <- 1e-00
mig_rat <- 0.15 / 250
pop_siz <- 26000
sel_gen <- -920
mig_gen <- -250
ext_frq <- 0.99
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH_1L(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_BCDO2_PMMH1.rda")

load("./REAL_BCDO2_PMMH1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
frq_pth_chn <- PMMH$frq_pth_chn[1, , ] + PMMH$frq_pth_chn[3, , ]

pdf(file = "./REAL_BCDO2_PMMH1_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.5)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
frq_pth_chn <- frq_pth_chn[, brn_num:dim(frq_pth_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, (1:floor(dim(frq_pth_chn)[2] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn)
colnames(smp) <- c("selection coefficient", "selection timing", "migration rate")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]

# sel_cof_est <- mean(sel_cof_chn)
# sel_gen_est <- mean(sel_gen_chn)
# mig_rat_est <- mean(mig_rat_chn)
frq_pth_est <- rowMeans(frq_pth_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
frq_pth_hpd <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_chn)[1]) {
  frq_pth_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_chn[k, ]), prob = 0.95)
}

pdf(file = "./REAL_BCDO2_PMMH1_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Migration rate",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 2], col = "blue", lty = 1, lwd = 2)
dev.off()

####################

# The case of population size N = 180000
sel_cof <- 5e-03
dom_par <- 1e-00
mig_rat <- 0.15 / 250
pop_siz <- 180000
sel_gen <- -920
mig_gen <- -250
ext_frq <- 0.99
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH_1L(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_BCDO2_PMMH2.rda")

load("./REAL_BCDO2_PMMH2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
frq_pth_chn <- PMMH$frq_pth_chn[1, , ] + PMMH$frq_pth_chn[3, , ]

pdf(file = "./REAL_BCDO2_PMMH2_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.5)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
frq_pth_chn <- frq_pth_chn[, brn_num:dim(frq_pth_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, (1:floor(dim(frq_pth_chn)[2] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn)
colnames(smp) <- c("selection coefficient", "selection timing", "migration rate")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]

# sel_cof_est <- mean(sel_cof_chn)
# sel_gen_est <- mean(sel_gen_chn)
# mig_rat_est <- mean(mig_rat_chn)
frq_pth_est <- rowMeans(frq_pth_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
frq_pth_hpd <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_chn)[1]) {
  frq_pth_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_chn[k, ]), prob = 0.95)
}

pdf(file = "./REAL_BCDO2_PMMH2_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Migration rate",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 2], col = "blue", lty = 1, lwd = 2)
dev.off()

####################

# The case of population size N = 460000
sel_cof <- 5e-03
dom_par <- 1e-00
mig_rat <- 0.15 / 250
pop_siz <- 460000
sel_gen <- -920
mig_gen <- -250
ext_frq <- 0.99
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH_1L(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_BCDO2_PMMH3.rda")

load("./REAL_BCDO2_PMMH3.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
frq_pth_chn <- PMMH$frq_pth_chn[1, , ] + PMMH$frq_pth_chn[3, , ]

pdf(file = "./REAL_BCDO2_PMMH3_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.5)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
frq_pth_chn <- frq_pth_chn[, brn_num:dim(frq_pth_chn)[2]]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, (1:floor(dim(frq_pth_chn)[2] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn)
colnames(smp) <- c("selection coefficient", "selection timing", "migration rate")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]

# sel_cof_est <- mean(sel_cof_chn)
# sel_gen_est <- mean(sel_gen_chn)
# mig_rat_est <- mean(mig_rat_chn)
frq_pth_est <- rowMeans(frq_pth_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
frq_pth_hpd <- matrix(NA, nrow = dim(frq_pth_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_chn)[1]) {
  frq_pth_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_chn[k, ]), prob = 0.95)
}

pdf(file = "./REAL_BCDO2_PMMH3_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Migration rate",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn), max(frq_pth_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[, 2], col = "blue", lty = 1, lwd = 2)
dev.off()

################################################################################

#' Joint analysis of TSHR and BCDO2

# TSHR
TSHR <- read_csv("./Data/groupedTSHR.csv")
TSHR <- as.data.frame(TSHR)
TSHR$smp_siz <- TSHR$D.count + TSHR$A.count
TSHR <- cbind(aggregate(subset(TSHR, select = c(`Mean(yearsAD)`)), by = list(TSHR$Grouped.year), FUN = mean)[, -1],
              aggregate(subset(TSHR, select = c(smp_siz, D.count)), by = list(TSHR$Grouped.year), FUN = sum)[, -1])
TSHR[, 1] <- as.integer(round(TSHR[, 1] - max(TSHR[, 1])))
TSHR[, 2] <- as.integer(round(TSHR[, 2]))
TSHR[, 3] <- as.integer(round(TSHR[, 3]))
rownames(TSHR) <- NULL
colnames(TSHR) <- c("smp_gen", "smp_siz", "smp_cnt")

smp_gen_A <- TSHR$smp_gen
smp_siz_A <- TSHR$smp_siz
smp_cnt_A <- matrix(NA, nrow = 2, ncol = length(smp_gen_A))
smp_cnt_A[1, ] <- TSHR$smp_cnt
smp_cnt_A[2, length(smp_gen_A)] <- round(0.15 * smp_siz_A[length(smp_gen_A)])

####################

#' BCDO2
BCDO2 <- read_csv("./Data/groupedBCDO2.csv")
BCDO2 <- as.data.frame(BCDO2)
BCDO2$smp_siz <- BCDO2$D.count + BCDO2$A.count
BCDO2 <- cbind(aggregate(subset(BCDO2, select = c(`Mean(yearsAD)`)), by = list(BCDO2$Grouped.year), FUN = mean)[, -1],
               aggregate(subset(BCDO2, select = c(smp_siz, D.count)), by = list(BCDO2$Grouped.year), FUN = sum)[, -1])
BCDO2[, 1] <- as.integer(round(BCDO2[, 1] - max(BCDO2[, 1])))
BCDO2[, 2] <- as.integer(round(BCDO2[, 2]))
BCDO2[, 3] <- as.integer(round(BCDO2[, 3]))
rownames(BCDO2) <- NULL
colnames(BCDO2) <- c("smp_gen", "smp_siz", "smp_cnt")

smp_gen_B <- BCDO2$smp_gen
smp_siz_B <- BCDO2$smp_siz
smp_cnt_B <- matrix(NA, nrow = 2, ncol = length(smp_gen_B))
smp_cnt_B[1, ] <- BCDO2$smp_cnt
smp_cnt_B[2, length(smp_gen_B)] <- round(0.15 * smp_siz_B[length(smp_gen_B)])

############################################################

# Joint estimation of the selection coefficient, selection timing and migration rate

# The case of population size N = 26000
sel_cof <- c(5e-03, 5e-03)
dom_par <- c(1e-00, 1e-00)
mig_rat <- 0.15 / 250
pop_siz <- 26000
sel_gen <- c(-920, -920)
mig_gen <- -250
ext_frq <- c(0.99, 0.99)
smp_gen_A
smp_siz_A
smp_cnt_A
smp_gen_B
smp_siz_B
smp_cnt_B
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH_2L(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen_A, smp_siz_A, smp_cnt_A, smp_gen_B, smp_siz_B, smp_cnt_B, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen_A, smp_siz_A, smp_cnt_A, smp_gen_B, smp_siz_B, smp_cnt_B, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_PMMH1.rda")

load("./REAL_PMMH1.rda")

sel_cof_A_chn <- PMMH$sel_cof_A_chn
sel_gen_A_chn <- PMMH$sel_gen_A_chn
frq_pth_A_chn <- PMMH$frq_pth_A_chn[1, , ] + PMMH$frq_pth_A_chn[3, , ]
sel_cof_B_chn <- PMMH$sel_cof_B_chn
sel_gen_B_chn <- PMMH$sel_gen_B_chn
frq_pth_B_chn <- PMMH$frq_pth_B_chn[1, , ] + PMMH$frq_pth_B_chn[3, , ]
mig_rat_chn <- PMMH$mig_rat_chn

pdf(file = "./REAL_PMMH1_Traceplot.pdf", width = 16, height = 18)
par(mfrow = c(3, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_A_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient of TSHR")

plot(1:itn_num, sel_gen_A_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing of TSHR")

plot(1:itn_num, sel_cof_B_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient of BCDO2")

plot(1:itn_num, sel_gen_B_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing of BCDO2")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(mig_rat_chn) * 0.5)
sel_cof_A_chn <- sel_cof_A_chn[brn_num:length(sel_cof_A_chn)]
sel_gen_A_chn <- sel_gen_A_chn[brn_num:length(sel_gen_A_chn)]
frq_pth_A_chn <- frq_pth_A_chn[, brn_num:dim(frq_pth_A_chn)[2]]
sel_cof_B_chn <- sel_cof_B_chn[brn_num:length(sel_cof_B_chn)]
sel_gen_B_chn <- sel_gen_B_chn[brn_num:length(sel_gen_B_chn)]
frq_pth_B_chn <- frq_pth_B_chn[, brn_num:dim(frq_pth_B_chn)[2]]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]

thn_num <- 5e+00
sel_cof_A_chn <- sel_cof_A_chn[(1:floor(length(sel_cof_A_chn) / thn_num)) * thn_num]
sel_gen_A_chn <- sel_gen_A_chn[(1:floor(length(sel_gen_A_chn) / thn_num)) * thn_num]
frq_pth_A_chn <- frq_pth_A_chn[, (1:floor(dim(frq_pth_A_chn)[2] / thn_num)) * thn_num]
sel_cof_B_chn <- sel_cof_B_chn[(1:floor(length(sel_cof_B_chn) / thn_num)) * thn_num]
sel_gen_B_chn <- sel_gen_B_chn[(1:floor(length(sel_gen_B_chn) / thn_num)) * thn_num]
frq_pth_B_chn <- frq_pth_B_chn[, (1:floor(dim(frq_pth_B_chn)[2] / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]

smp <- data.frame(sel_cof_A_chn, sel_gen_A_chn, sel_cof_B_chn, sel_gen_B_chn, mig_rat_chn)
colnames(smp) <- c("selection coefficient A", "selection timing A", "selection coefficient B", "selection timing B", "migration rate")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_A_est <- est[1]
sel_gen_A_est <- est[2]
sel_cof_B_est <- est[3]
sel_gen_B_est <- est[4]
mig_rat_est <- est[5]

# sel_cof_A_est <- mean(sel_cof_A_chn)
# sel_gen_A_est <- mean(sel_gen_A_chn)
frq_pth_A_est <- rowMeans(frq_pth_A_chn)
# sel_cof_B_est <- mean(sel_cof_B_chn)
# sel_gen_B_est <- mean(sel_gen_B_chn)
frq_pth_B_est <- rowMeans(frq_pth_B_chn)
# mig_rat_est <- mean(mig_rat_chn)

sel_cof_A_hpd <- HPDinterval(as.mcmc(sel_cof_A_chn), prob = 0.95)
sel_gen_A_hpd <- HPDinterval(as.mcmc(sel_gen_A_chn), prob = 0.95)
frq_pth_A_hpd <- matrix(NA, nrow = dim(frq_pth_A_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_A_chn)[1]) {
  frq_pth_A_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_A_chn[k, ]), prob = 0.95)
}
sel_cof_B_hpd <- HPDinterval(as.mcmc(sel_cof_B_chn), prob = 0.95)
sel_gen_B_hpd <- HPDinterval(as.mcmc(sel_gen_B_chn), prob = 0.95)
frq_pth_B_hpd <- matrix(NA, nrow = dim(frq_pth_B_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_B_chn)[1]) {
  frq_pth_B_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_B_chn[k, ]), prob = 0.95)
}
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)

pdf(file = "./REAL_PMMH1_Posterior.pdf", width = 24, height = 18)
par(mfrow = c(3, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_A_chn, breaks = seq(min(sel_cof_A_chn), max(sel_cof_A_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient of TSHR")
lines(density(sel_cof_A_chn), lwd = 2, col = 'black')
abline(v = sel_cof_A_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_A_chn, breaks = seq(min(sel_gen_A_chn), max(sel_gen_A_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing of TSHR")
lines(density(sel_gen_A_chn), lwd = 2, col = 'black')
abline(v = sel_gen_A_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_A_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_A_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_A_chn), max(frq_pth_A_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory of TSHR")
for (i in 1:dim(frq_pth_A_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_A_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_A_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_A_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_A_hpd[, 2], col = "blue", lty = 1, lwd = 2)

hist(sel_cof_B_chn, breaks = seq(min(sel_cof_B_chn), max(sel_cof_B_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient of BCDO2")
lines(density(sel_cof_B_chn), lwd = 2, col = 'black')
abline(v = sel_cof_B_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_B_chn, breaks = seq(min(sel_gen_B_chn), max(sel_gen_B_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing of BCDO2")
lines(density(sel_gen_B_chn), lwd = 2, col = 'black')
abline(v = sel_gen_B_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_B_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_B_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_B_chn), max(frq_pth_B_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory of BCDO2")
for (i in 1:dim(frq_pth_B_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_B_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_B_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_B_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_B_hpd[, 2], col = "blue", lty = 1, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Migration rate",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

####################

# The case of population size N = 180000
sel_cof <- c(5e-03, 5e-03)
dom_par <- c(1e-00, 1e-00)
mig_rat <- 0.15 / 250
pop_siz <- 180000
sel_gen <- c(-920, -920)
mig_gen <- -250
ext_frq <- c(0.99, 0.99)
smp_gen_A
smp_siz_A
smp_cnt_A
smp_gen_B
smp_siz_B
smp_cnt_B
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH_2L(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen_A, smp_siz_A, smp_cnt_A, smp_gen_B, smp_siz_B, smp_cnt_B, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen_A, smp_siz_A, smp_cnt_A, smp_gen_B, smp_siz_B, smp_cnt_B, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_PMMH2.rda")

load("./REAL_PMMH2.rda")

sel_cof_A_chn <- PMMH$sel_cof_A_chn
sel_gen_A_chn <- PMMH$sel_gen_A_chn
frq_pth_A_chn <- PMMH$frq_pth_A_chn[1, , ] + PMMH$frq_pth_A_chn[3, , ]
sel_cof_B_chn <- PMMH$sel_cof_B_chn
sel_gen_B_chn <- PMMH$sel_gen_B_chn
frq_pth_B_chn <- PMMH$frq_pth_B_chn[1, , ] + PMMH$frq_pth_B_chn[3, , ]
mig_rat_chn <- PMMH$mig_rat_chn

pdf(file = "./REAL_PMMH2_Traceplot.pdf", width = 16, height = 18)
par(mfrow = c(3, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_A_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient of TSHR")

plot(1:itn_num, sel_gen_A_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing of TSHR")

plot(1:itn_num, sel_cof_B_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient of BCDO2")

plot(1:itn_num, sel_gen_B_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing of BCDO2")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(mig_rat_chn) * 0.5)
sel_cof_A_chn <- sel_cof_A_chn[brn_num:length(sel_cof_A_chn)]
sel_gen_A_chn <- sel_gen_A_chn[brn_num:length(sel_gen_A_chn)]
frq_pth_A_chn <- frq_pth_A_chn[, brn_num:dim(frq_pth_A_chn)[2]]
sel_cof_B_chn <- sel_cof_B_chn[brn_num:length(sel_cof_B_chn)]
sel_gen_B_chn <- sel_gen_B_chn[brn_num:length(sel_gen_B_chn)]
frq_pth_B_chn <- frq_pth_B_chn[, brn_num:dim(frq_pth_B_chn)[2]]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]

thn_num <- 5e+00
sel_cof_A_chn <- sel_cof_A_chn[(1:floor(length(sel_cof_A_chn) / thn_num)) * thn_num]
sel_gen_A_chn <- sel_gen_A_chn[(1:floor(length(sel_gen_A_chn) / thn_num)) * thn_num]
frq_pth_A_chn <- frq_pth_A_chn[, (1:floor(dim(frq_pth_A_chn)[2] / thn_num)) * thn_num]
sel_cof_B_chn <- sel_cof_B_chn[(1:floor(length(sel_cof_B_chn) / thn_num)) * thn_num]
sel_gen_B_chn <- sel_gen_B_chn[(1:floor(length(sel_gen_B_chn) / thn_num)) * thn_num]
frq_pth_B_chn <- frq_pth_B_chn[, (1:floor(dim(frq_pth_B_chn)[2] / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]

smp <- data.frame(sel_cof_A_chn, sel_gen_A_chn, sel_cof_B_chn, sel_gen_B_chn, mig_rat_chn)
colnames(smp) <- c("selection coefficient A", "selection timing A", "selection coefficient B", "selection timing B", "migration rate")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_A_est <- est[1]
sel_gen_A_est <- est[2]
sel_cof_B_est <- est[3]
sel_gen_B_est <- est[4]
mig_rat_est <- est[5]

# sel_cof_A_est <- mean(sel_cof_A_chn)
# sel_gen_A_est <- mean(sel_gen_A_chn)
frq_pth_A_est <- rowMeans(frq_pth_A_chn)
# sel_cof_B_est <- mean(sel_cof_B_chn)
# sel_gen_B_est <- mean(sel_gen_B_chn)
frq_pth_B_est <- rowMeans(frq_pth_B_chn)
# mig_rat_est <- mean(mig_rat_chn)

sel_cof_A_hpd <- HPDinterval(as.mcmc(sel_cof_A_chn), prob = 0.95)
sel_gen_A_hpd <- HPDinterval(as.mcmc(sel_gen_A_chn), prob = 0.95)
frq_pth_A_hpd <- matrix(NA, nrow = dim(frq_pth_A_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_A_chn)[1]) {
  frq_pth_A_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_A_chn[k, ]), prob = 0.95)
}
sel_cof_B_hpd <- HPDinterval(as.mcmc(sel_cof_B_chn), prob = 0.95)
sel_gen_B_hpd <- HPDinterval(as.mcmc(sel_gen_B_chn), prob = 0.95)
frq_pth_B_hpd <- matrix(NA, nrow = dim(frq_pth_B_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_B_chn)[1]) {
  frq_pth_B_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_B_chn[k, ]), prob = 0.95)
}
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)

pdf(file = "./REAL_PMMH2_Posterior.pdf", width = 24, height = 18)
par(mfrow = c(3, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_A_chn, breaks = seq(min(sel_cof_A_chn), max(sel_cof_A_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient of TSHR")
lines(density(sel_cof_A_chn), lwd = 2, col = 'black')
abline(v = sel_cof_A_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_A_chn, breaks = seq(min(sel_gen_A_chn), max(sel_gen_A_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing of TSHR")
lines(density(sel_gen_A_chn), lwd = 2, col = 'black')
abline(v = sel_gen_A_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_A_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_A_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_A_chn), max(frq_pth_A_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory of TSHR")
for (i in 1:dim(frq_pth_A_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_A_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_A_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_A_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_A_hpd[, 2], col = "blue", lty = 1, lwd = 2)

hist(sel_cof_B_chn, breaks = seq(min(sel_cof_B_chn), max(sel_cof_B_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient of BCDO2")
lines(density(sel_cof_B_chn), lwd = 2, col = 'black')
abline(v = sel_cof_B_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_B_chn, breaks = seq(min(sel_gen_B_chn), max(sel_gen_B_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing of BCDO2")
lines(density(sel_gen_B_chn), lwd = 2, col = 'black')
abline(v = sel_gen_B_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_B_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_B_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_B_chn), max(frq_pth_B_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory of BCDO2")
for (i in 1:dim(frq_pth_B_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_B_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_B_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_B_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_B_hpd[, 2], col = "blue", lty = 1, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Migration rate",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

####################

# The case of population size N = 460000
sel_cof <- c(5e-03, 5e-03)
dom_par <- c(1e-00, 1e-00)
mig_rat <- 0.15 / 250
pop_siz <- 460000
sel_gen <- c(-920, -920)
mig_gen <- -250
ext_frq <- c(0.99, 0.99)
smp_gen_A
smp_siz_A
smp_cnt_A
smp_gen_B
smp_siz_B
smp_cnt_B
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(PMMH <- cmprunPMMH_2L(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen_A, smp_siz_A, smp_cnt_A, smp_gen_B, smp_siz_B, smp_cnt_B, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen_A, smp_siz_A, smp_cnt_A, smp_gen_B, smp_siz_B, smp_cnt_B, ptn_num, pcl_num, itn_num, PMMH,
     file = "./REAL_PMMH3.rda")

load("./REAL_PMMH3.rda")

sel_cof_A_chn <- PMMH$sel_cof_A_chn
sel_gen_A_chn <- PMMH$sel_gen_A_chn
frq_pth_A_chn <- PMMH$frq_pth_A_chn[1, , ] + PMMH$frq_pth_A_chn[3, , ]
sel_cof_B_chn <- PMMH$sel_cof_B_chn
sel_gen_B_chn <- PMMH$sel_gen_B_chn
frq_pth_B_chn <- PMMH$frq_pth_B_chn[1, , ] + PMMH$frq_pth_B_chn[3, , ]
mig_rat_chn <- PMMH$mig_rat_chn

pdf(file = "./REAL_PMMH3_Traceplot.pdf", width = 16, height = 18)
par(mfrow = c(3, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_A_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient of TSHR")

plot(1:itn_num, sel_gen_A_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing of TSHR")

plot(1:itn_num, sel_cof_B_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient of BCDO2")

plot(1:itn_num, sel_gen_B_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection timing",
     main = "Trace plot of the selection timing of BCDO2")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(mig_rat_chn) * 0.5)
sel_cof_A_chn <- sel_cof_A_chn[brn_num:length(sel_cof_A_chn)]
sel_gen_A_chn <- sel_gen_A_chn[brn_num:length(sel_gen_A_chn)]
frq_pth_A_chn <- frq_pth_A_chn[, brn_num:dim(frq_pth_A_chn)[2]]
sel_cof_B_chn <- sel_cof_B_chn[brn_num:length(sel_cof_B_chn)]
sel_gen_B_chn <- sel_gen_B_chn[brn_num:length(sel_gen_B_chn)]
frq_pth_B_chn <- frq_pth_B_chn[, brn_num:dim(frq_pth_B_chn)[2]]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]

thn_num <- 5e+00
sel_cof_A_chn <- sel_cof_A_chn[(1:floor(length(sel_cof_A_chn) / thn_num)) * thn_num]
sel_gen_A_chn <- sel_gen_A_chn[(1:floor(length(sel_gen_A_chn) / thn_num)) * thn_num]
frq_pth_A_chn <- frq_pth_A_chn[, (1:floor(dim(frq_pth_A_chn)[2] / thn_num)) * thn_num]
sel_cof_B_chn <- sel_cof_B_chn[(1:floor(length(sel_cof_B_chn) / thn_num)) * thn_num]
sel_gen_B_chn <- sel_gen_B_chn[(1:floor(length(sel_gen_B_chn) / thn_num)) * thn_num]
frq_pth_B_chn <- frq_pth_B_chn[, (1:floor(dim(frq_pth_B_chn)[2] / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]

smp <- data.frame(sel_cof_A_chn, sel_gen_A_chn, sel_cof_B_chn, sel_gen_B_chn, mig_rat_chn)
colnames(smp) <- c("selection coefficient A", "selection timing A", "selection coefficient B", "selection timing B", "migration rate")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_A_est <- est[1]
sel_gen_A_est <- est[2]
sel_cof_B_est <- est[3]
sel_gen_B_est <- est[4]
mig_rat_est <- est[5]

# sel_cof_A_est <- mean(sel_cof_A_chn)
# sel_gen_A_est <- mean(sel_gen_A_chn)
frq_pth_A_est <- rowMeans(frq_pth_A_chn)
# sel_cof_B_est <- mean(sel_cof_B_chn)
# sel_gen_B_est <- mean(sel_gen_B_chn)
frq_pth_B_est <- rowMeans(frq_pth_B_chn)
# mig_rat_est <- mean(mig_rat_chn)

sel_cof_A_hpd <- HPDinterval(as.mcmc(sel_cof_A_chn), prob = 0.95)
sel_gen_A_hpd <- HPDinterval(as.mcmc(sel_gen_A_chn), prob = 0.95)
frq_pth_A_hpd <- matrix(NA, nrow = dim(frq_pth_A_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_A_chn)[1]) {
  frq_pth_A_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_A_chn[k, ]), prob = 0.95)
}
sel_cof_B_hpd <- HPDinterval(as.mcmc(sel_cof_B_chn), prob = 0.95)
sel_gen_B_hpd <- HPDinterval(as.mcmc(sel_gen_B_chn), prob = 0.95)
frq_pth_B_hpd <- matrix(NA, nrow = dim(frq_pth_B_chn)[1], ncol = 2)
for (k in 1:dim(frq_pth_B_chn)[1]) {
  frq_pth_B_hpd[k, ] = HPDinterval(as.mcmc(frq_pth_B_chn[k, ]), prob = 0.95)
}
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)

pdf(file = "./REAL_PMMH3_Posterior.pdf", width = 24, height = 18)
par(mfrow = c(3, 3), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_A_chn, breaks = seq(min(sel_cof_A_chn), max(sel_cof_A_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient of TSHR")
lines(density(sel_cof_A_chn), lwd = 2, col = 'black')
abline(v = sel_cof_A_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_A_chn, breaks = seq(min(sel_gen_A_chn), max(sel_gen_A_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing of TSHR")
lines(density(sel_gen_A_chn), lwd = 2, col = 'black')
abline(v = sel_gen_A_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_A_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_A_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_A_chn), max(frq_pth_A_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory of TSHR")
for (i in 1:dim(frq_pth_A_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_A_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_A_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_A_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_A_hpd[, 2], col = "blue", lty = 1, lwd = 2)

hist(sel_cof_B_chn, breaks = seq(min(sel_cof_B_chn), max(sel_cof_B_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient of BCDO2")
lines(density(sel_cof_B_chn), lwd = 2, col = 'black')
abline(v = sel_cof_B_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_B_chn, breaks = seq(min(sel_gen_B_chn), max(sel_gen_B_chn), length.out = 50), freq = FALSE,
     xlab = "Selection timing",
     main = "Posterior for the selection timing of BCDO2")
lines(density(sel_gen_B_chn), lwd = 2, col = 'black')
abline(v = sel_gen_B_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_B_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_B_hpd[2], col = 'blue', lty = 2, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_B_chn), max(frq_pth_B_chn)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated mutant allele frequency trajectory of BCDO2")
for (i in 1:dim(frq_pth_B_chn)[2]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_B_chn[, i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), frq_pth_B_est, col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_B_hpd[, 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_B_hpd[, 2], col = "blue", lty = 1, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Migration rate",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
