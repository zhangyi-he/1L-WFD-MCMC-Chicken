#' @title Inferring natural selection and gene migration in the evolution of chickens from ancient DNA data
#' @author Zhangyi He, Wenyang Lyu, Xiaoyang Dai, Mark Beaumont, Feng Yu

#' version 1.2

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2018/HE2021-WFM-1L-DiffusApprox-PMMH-MolEcolResour")

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

################################################################################

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

# Joint estimation of the selection coefficient, selection timing, migration rate and migration timing

# call R functions
source("./Code/Code v1.0/Code v1.2/RFUN.R")

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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N1.rda")

load("./Output/Output v1.0/REAL_TSHR_PMMH_4N1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N1_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N1_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N1_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N2.rda")

load("./Output/Output v1.0/REAL_TSHR_PMMH_4N2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N2_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N2_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N2_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N3.rda")

load("./Output/Output v1.0/REAL_TSHR_PMMH_4N3.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N3_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N3_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_4N3_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
dev.off()

########################################

# Joint estimation of the selection coefficient, selection timing and migration rate

# call R functions
source("./Code/Code v1.0/Code v1.2/RFUN_REAL.R")

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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N1.rda")

load("./Output/Output v1.0/REAL_TSHR_PMMH_3N1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N1_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N1_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N1_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N2.rda")

load("./Output/Output v1.0/REAL_TSHR_PMMH_3N2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N2_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N2_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N2_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N3.rda")

load("./Output/Output v1.0/REAL_TSHR_PMMH_3N3.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N3_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N3_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_TSHR_PMMH_3N3_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
dev.off()

################################################################################

#' BCD02
BCD02 <- read_csv("./Data/groupedBCD02.csv")
BCD02 <- as.data.frame(BCD02)
BCD02$smp_siz <- BCD02$D.count + BCD02$A.count
BCD02 <- cbind(aggregate(subset(BCD02, select = c(`Mean(yearsAD)`)), by = list(BCD02$Grouped.year), FUN = mean)[, -1],
              aggregate(subset(BCD02, select = c(smp_siz, D.count)), by = list(BCD02$Grouped.year), FUN = sum)[, -1])
BCD02[, 1] <- as.integer(round(BCD02[, 1] - max(BCD02[, 1])))
BCD02[, 2] <- as.integer(round(BCD02[, 2]))
BCD02[, 3] <- as.integer(round(BCD02[, 3]))
rownames(BCD02) <- NULL
colnames(BCD02) <- c("smp_gen", "smp_siz", "smp_cnt")

smp_gen <- BCD02$smp_gen
smp_siz <- BCD02$smp_siz
smp_cnt <- matrix(NA, nrow = 2, ncol = length(smp_gen))
smp_cnt[1, ] <- BCD02$smp_cnt
smp_cnt[2, length(smp_gen)] <- round(0.15 * smp_siz[length(smp_gen)])

############################################################

# Joint estimation of the selection coefficient, selection timing, migration rate and migration timing

# call R functions
source("./Code/Code v1.0/Code v1.2/RFUN.R")

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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N1.rda")

load("./Output/Output v1.0/REAL_BCD02_PMMH_4N1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N1_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N1_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N1_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N2.rda")

load("./Output/Output v1.0/REAL_BCD02_PMMH_4N2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N2_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N2_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N2_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N3.rda")

load("./Output/Output v1.0/REAL_BCD02_PMMH_4N3.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N3_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N3_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_4N3_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
dev.off()

########################################

# Joint estimation of the selection coefficient, selection timing and migration rate

# call R functions
source("./Code/Code v1.0/Code v1.2/RFUN_REAL.R")

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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N1.rda")

load("./Output/Output v1.0/REAL_BCD02_PMMH_3N1.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N1_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N1_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N1_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N2.rda")

load("./Output/Output v1.0/REAL_BCD02_PMMH_3N2.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N2_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N2_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N2_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
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
bck_smp <- TRUE # blockwise sampling

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp, PMMH,
     file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N3.rda")

load("./Output/Output v1.0/REAL_BCD02_PMMH_3N3.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
frq_pth_chn <- frq_pth_chn[, (0:(max(smp_gen) - min(smp_gen))) * ptn_num + 1, ]

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N3_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
dev.off()

# brn_num <- 1e+04
brn_num <- floor(length(sel_cof_chn) * 0.8)
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
frq_pth_chn <- frq_pth_chn[, , brn_num:dim(frq_pth_chn)[3]]

thn_num <- 4e+00
sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]
frq_pth_chn <- frq_pth_chn[, , (1:floor(dim(frq_pth_chn)[3] / thn_num)) * thn_num]

smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
pdf <- kepdf(smp, bwtype = "adaptive")
plot(pdf, main = "Posterior for strength and timing of natural selection and gene migration", col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32),
     props = c(100, 100, 100, 100), method = "perspective", gap = 0.5, phi = 30, theta = 10)
est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
sel_cof_est <- est[1]
sel_gen_est <- est[2]
mig_rat_est <- est[3]
mig_gen_est <- est[4]

# sel_cof_est <- median(sel_cof_chn)
# sel_gen_est <- median(sel_gen_chn)
# mig_rat_est <- median(mig_rat_chn)
# mig_gen_est <- median(mig_gen_chn)
frq_pth_est <- matrix(NA, nrow = 4, ncol = dim(frq_pth_chn)[2])
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_est[i, k] <- mean(frq_pth_chn[i, k, ])
  }
}

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)
frq_pth_hpd <- array(NA, dim = c(4, dim(frq_pth_chn)[2], 2))
for (k in 1:dim(frq_pth_chn)[2]) {
  for (i in 1:4) {
    frq_pth_hpd[i, k, ] = HPDinterval(as.mcmc(frq_pth_chn[i, k, ]), prob = 0.95)
  }
}

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N3_Posterior.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the selection time")
lines(density(sel_gen_chn), lwd = 2, col = 'black')
abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for the migration time")
lines(density(mig_gen_chn), lwd = 2, col = 'black')
abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

pdf(file = "./Output/Output v1.0/REAL_BCD02_PMMH_3N3_Trajectory.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[1, , ]), max(frq_pth_chn[1, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[1, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[1, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[1, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[1, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[2, , ]), max(frq_pth_chn[2, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated island ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[2, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[2, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[2, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[2, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[3, , ]), max(frq_pth_chn[3, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent mutant allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[3, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[3, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[3, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[3, , 2], col = "blue", lty = 1, lwd = 2)

plot(0, type = "n", xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(frq_pth_chn[4, , ]), max(frq_pth_chn[4, , ])),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Estimated continent ancestral allele frequency trajectory")
for (i in 1:dim(frq_pth_chn)[3]) {
  lines(min(smp_gen):max(smp_gen), frq_pth_chn[4, , i], col = "grey", lty = 1, lwd = 2)
}
lines(min(smp_gen):max(smp_gen), pop_ale_frq[4, ], col = "red", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_est[4, ], col = "black", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 1], col = "blue", lty = 1, lwd = 2)
lines(min(smp_gen):max(smp_gen), frq_pth_hpd[4, , 2], col = "blue", lty = 1, lwd = 2)
dev.off()

################################################################################
