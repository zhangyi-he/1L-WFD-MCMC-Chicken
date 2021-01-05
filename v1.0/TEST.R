#' @title Inferring natural selection and gene migration in the evolution of chickens from ancient DNA data
#' @author Zhangyi He, Wenyang Lyu, Xiaoyang Dai, Sile Hu, Mark Beaumont, Feng Yu

#' version 1.0

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2018/HE2021-WFM-1L-DiffusApprox-PMMH-Chicken-MolEcol")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

# call R functions
source("./Code/Code v1.0/Code v1.0/RFUN.R")

################################################################################

#' Simulate the allele frequency trajectories according to the one-locus Wright-Fisher model with selection and migration
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the island population
#' @param ext_frq the mutant allele frequency (of the continent population)
#' @param int_frq the initial allele frequencies (of the island population)
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory

# sel_cof <- 1e-02
# dom_par <- 5e-01
# mig_rat <- 1e-03
# pop_siz <- 5e+03
# ext_frq <- 9e-01
# int_frq <- c(1e-01, 9e-01, 0e-00, 0e-00)
# int_gen <- 0
# lst_gen <- 500
#
# frq_pth <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen)
#
# k <- int_gen:lst_gen
# plot(k, frq_pth[1, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Allele frequency",
#      main = "WFM: the island mutant allele frequency trajectory")
# plot(k, frq_pth[2, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Allele frequency",
#      main = "WFM: the island ancestral allele frequency trajectory")
# plot(k, frq_pth[3, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Allele frequency",
#      main = "WFM: the continent mutant allele frequency trajectory")
# plot(k, frq_pth[4, ], type = "l", lwd = 1.5,
#      xlab = "Generation", ylab = "Allele frequency",
#      main = "WFM: the continent ancestral allele frequency trajectory")

########################################

#' Simulate the allele frequency trajectories under the one-locus Wright-Fisher diffusion with selection migration using the Euler-Maruyama method
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the island population
#' @param ext_frq the mutant allele frequency (of the continent population)
#' @param int_frq the initial allele frequencies (of the island population)
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

# sel_cof <- 1e-02
# dom_par <- 5e-01
# mig_rat <- 1e-03
# pop_siz <- 5e+03
# ext_frq <- 9e-01
# int_frq <- c(1e-01, 9e-01, 0e-00, 0e-00)
# int_gen <- 0
# lst_gen <- 500
# ptn_num <- 5e+00
#
# frq_pth <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE)
#
# t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
# plot(t, frq_pth[1, ], type = "l", lwd = 1.5,
#      xlab = "Time", ylab = "Allele frequency",
#      main = "WFD: the island mutant allele frequency trajectory")
# plot(t, frq_pth[2, ], type = "l", lwd = 1.5,
#      xlab = "Time", ylab = "Allele frequency",
#      main = "WFD: the island ancestral allele frequency trajectory")
# plot(t, frq_pth[3, ], type = "l", lwd = 1.5,
#      xlab = "Time", ylab = "Allele frequency",
#      main = "WFD: the continent mutant allele frequency trajectory")
# plot(t, frq_pth[4, ], type = "l", lwd = 1.5,
#      xlab = "Time", ylab = "Allele frequency",
#      main = "WFD: the continent ancestral allele frequency trajectory")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
# sel_cof <- 1e-02
# dom_par <- 5e-01
# mig_rat <- 1e-03
# pop_siz <- 5e+03
# ext_frq <- 9e-01
# int_frq <- c(1e-01, 9e-01, 0e-00, 0e-00)
# int_gen <- 0
# lst_gen <- 500
# ptn_num <- 5e+00
# sim_num <- 1e+06
#
# smp_WFM <- matrix(NA, nrow = 4, ncol = sim_num)
# smp_WFD <- matrix(NA, nrow = 4, ncol = sim_num)
# for (i in 1:sim_num) {
#   print(i)
#   smp_WFM[, i] <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
#   smp_WFD[, i] <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num, dat_aug = FALSE)[, (lst_gen - int_gen) + 1]
# }
#
# hist(smp_WFM[1, ], breaks = seq(min(smp_WFM[1, ], smp_WFD[1, ]), max(smp_WFM[1, ], smp_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
#      xlim = c(min(smp_WFM[1, ], smp_WFD[1, ]), max(smp_WFM[1, ], smp_WFD[1, ])),
#      xlab = "Allele frequency", main = paste("Histogram of the island mutant allele at generation", lst_gen))
# hist(smp_WFD[1, ], breaks = seq(min(smp_WFM[1, ], smp_WFD[1, ]), max(smp_WFM[1, ], smp_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
# hist(smp_WFM[2, ], breaks = seq(min(smp_WFM[2, ], smp_WFD[2, ]), max(smp_WFM[2, ], smp_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
#      xlim = c(min(smp_WFM[2, ], smp_WFD[2, ]), max(smp_WFM[2, ], smp_WFD[2, ])),
#      xlab = "Allele frequency", main = paste("Histogram of the island ancestral allele at generation", lst_gen))
# hist(smp_WFD[2, ], breaks = seq(min(smp_WFM[2, ], smp_WFD[2, ]), max(smp_WFM[2, ], smp_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
# hist(smp_WFM[3, ], breaks = seq(min(smp_WFM[3, ], smp_WFD[3, ]), max(smp_WFM[3, ], smp_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
#      xlim = c(min(smp_WFM[3, ], smp_WFD[3, ]), max(smp_WFM[3, ], smp_WFD[3, ])),
#      xlab = "Allele frequency", main = paste("Histogram of the continent mutant allele at generation", lst_gen))
# hist(smp_WFD[3, ], breaks = seq(min(smp_WFM[3, ], smp_WFD[3, ]), max(smp_WFM[3, ], smp_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
# hist(smp_WFM[4, ], breaks = seq(min(smp_WFM[4, ], smp_WFD[4, ]), max(smp_WFM[4, ], smp_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
#      xlim = c(min(smp_WFM[4, ], smp_WFD[4, ]), max(smp_WFM[4, ], smp_WFD[4, ])),
#      xlab = "Allele frequency", main = paste("Histogram of the continent ancestral allele at generation", lst_gen))
# hist(smp_WFD[4, ], breaks = seq(min(smp_WFM[4, ], smp_WFD[4, ]), max(smp_WFM[4, ], smp_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

################################################################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param sel_gen the starting time of natural selection
#' @param mig_gen the starting time of gene migration
#' @param ext_frq the mutant allele frequency (of the continent population)
#' @param int_frq the initial allele frequencies (of the island population)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param cna_gen the sampling time points that the count of the continent alleles are not available in the sample
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model
# model <- "WFM"
# sel_cof <- 1e-02
# dom_par <- 5e-01
# mig_rat <- 1e-03
# pop_siz <- 5e+03
# sel_gen <- 120
# mig_gen <- 320
# ext_frq <- 9e-01
# int_frq <- c(1e-01, 9e-01, 0e-00, 0e-00)
# smp_gen <- (0:10) * 50
# smp_siz <- rep(100, 11)
# cna_gen <- NULL
#
# sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, int_frq, smp_gen, smp_siz, cna_gen)
# smp_gen <- sim_HMM_WFM$smp_gen
# smp_siz <- sim_HMM_WFM$smp_siz
# smp_cnt <- sim_HMM_WFM$smp_cnt
# smp_frq <- sim_HMM_WFM$smp_frq
# pop_frq <- sim_HMM_WFM$pop_frq
# # smp_ale_cnt <- sim_HMM_WFM$smp_ale_cnt
# # smp_ale_frq <- sim_HMM_WFM$smp_ale_frq
# # pop_ale_frq <- sim_HMM_WFM$pop_ale_frq
#
# k <- min(smp_gen):max(smp_gen)
# plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[1, ], pop_frq[1, ], na.rm = TRUE), max(smp_frq[1, ], pop_frq[1, ], na.rm = TRUE)),
#      xlab = "Generation", ylab = "Allele frequency",
#      main = "WFM-HMM: the mutant allele")
# points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[2, ], pop_frq[2, ], na.rm = TRUE), max(smp_frq[2, ], pop_frq[2, ], na.rm = TRUE)),
#      xlab = "Generation", ylab = "Allele frequency",
#      main = "WFM-HMM: the continent allele")
# points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
# model <- "WFD"
# sel_cof <- 1e-02
# dom_par <- 5e-01
# mig_rat <- 1e-03
# pop_siz <- 5e+03
# sel_gen <- 120
# mig_gen <- 320
# ext_frq <- 9e-01
# int_frq <- c(1e-01, 9e-01, 0e-00, 0e-00)
# smp_gen <- (0:10) * 50
# smp_siz <- rep(100, 11)
# cna_gen <- NULL
# ptn_num <- 5e+00
#
# sim_HMM_WFD <- cmpsimulateHMM(model, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, int_frq, smp_gen, smp_siz, cna_gen, ptn_num)
# smp_gen <- sim_HMM_WFD$smp_gen
# smp_siz <- sim_HMM_WFD$smp_siz
# smp_cnt <- sim_HMM_WFD$smp_cnt
# smp_frq <- sim_HMM_WFD$smp_frq
# pop_frq <- sim_HMM_WFD$pop_frq
# # smp_ale_cnt <- sim_HMM_WFD$smp_ale_cnt
# # smp_ale_frq <- sim_HMM_WFD$smp_ale_frq
# # pop_ale_frq <- sim_HMM_WFD$pop_ale_frq
#
# k <- min(smp_gen):max(smp_gen)
# plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[1, ], pop_frq[1, ], na.rm = TRUE), max(smp_frq[1, ], pop_frq[1, ], na.rm = TRUE)),
#      xlab = "Generation", ylab = "Allele frequency",
#      main = "WFD-HMM: the mutant allele")
# points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)
# plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
#      xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq[2, ], pop_frq[2, ], na.rm = TRUE), max(smp_frq[2, ], pop_frq[2, ], na.rm = TRUE)),
#      xlab = "Generation", ylab = "Allele frequency",
#      main = "WFD-HMM: the continent allele")
# points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model
test_seed <- 19
set.seed(test_seed)

#' Simulate the dataset under the Wright-Fisher model
model <- "WFM"
sel_cof <- 8e-03
dom_par <- 5e-01
mig_rat <- 5e-03
pop_siz <- 5e+03
sel_gen <- 180
mig_gen <- 360
ext_frq <- 9e-01
int_frq <- c(4e-01, 6e-01, 0e-00, 0e-00)
smp_gen <- (0:10) * 50
smp_siz <- rep(100, 11)
# cna_gen <- sample(smp_gen, 5)
cna_gen <- NULL

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, int_frq, smp_gen, smp_siz, cna_gen)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_cnt
smp_frq <- sim_HMM_WFM$smp_frq
pop_frq <- sim_HMM_WFM$pop_frq
smp_ale_cnt <- sim_HMM_WFM$smp_ale_cnt
smp_ale_frq <- sim_HMM_WFM$smp_ale_frq
pop_ale_frq <- sim_HMM_WFM$pop_ale_frq

save(model, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, int_frq, smp_gen, smp_siz, cna_gen, smp_cnt, smp_frq, pop_frq, smp_ale_cnt, smp_ale_frq, pop_ale_frq,
     file = "./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_SimData.pdf", width = 16, height = 6)
par(mfrow = c(1, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq[1, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_frq[1, ], smp_frq[1, ], na.rm = TRUE), max(pop_frq[1, ], smp_frq[1, ], na.rm = TRUE)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Mutant alleles")
points(smp_gen, smp_frq[1, ], col = 'red', pch = 17, cex = 1)

plot(k, pop_frq[2, ], type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_frq[2, ], smp_frq[2, ], na.rm = TRUE), max(pop_frq[2, ], smp_frq[2, ], na.rm = TRUE)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Continent alleles")
points(smp_gen, smp_frq[2, ], col = 'red', pch = 17, cex = 1)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the two-locus Wright-Fisher diffusion with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param sel_gen the starting time of natural selection
#' @param mig_gen the starting time of gene migration
#' @param ext_frq the mutant allele frequency (of the continent population)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and continent alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
mig_rat
pop_siz
sel_gen
mig_gen
ext_frq
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+05

system.time(BPF <- cmprunBPF(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, BPF,
     file = "./Output/Output v1.0/Test v1.0/TEST_BPF.rda")

load("./Output/Output v1.0/Test v1.0/TEST_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_BPF_Likelihood.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:pcl_num, log(lik), type = 'l',
     xlab = "Number of particles", ylab = "Log likelihood",
     main = "Log likelihood through the bootstrap particle filter")
dev.off()

pop_frq_pre_resmp <- BPF$pop_frq_pre_resmp
pop_frq_pst_resmp <- BPF$pop_frq_pst_resmp

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_BPF_Particle.pdf", width = 32, height = 66)
par(mfrow = c(11, 4), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
for (k in 1:length(smp_gen)) {
  hist_pst_resmp <- hist(pop_frq_pst_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  if (!any(is.nan(c(hist_pst_resmp$density, hist_pre_resmp$density)))) {
    hist(pop_frq_pst_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
         xlim = c(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k], smp_ale_frq[1, k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k], smp_ale_frq[1, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
         xlab = "Allele frequency",
         main = paste("Island mutant allele at generation", smp_gen[k]))
    hist(pop_frq_pre_resmp[1, , k], breaks = seq(min(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), max(pop_frq_pst_resmp[1, , k], pop_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
    abline(v = smp_ale_frq[1, k], col = 'red', lty = 2, lwd = 2)
  } else {
    plot(0, type = "n", axes = FALSE, ann = FALSE)
  }

  hist_pst_resmp <- hist(pop_frq_pst_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  if (!any(is.nan(c(hist_pst_resmp$density, hist_pre_resmp$density)))) {
    hist(pop_frq_pst_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
         xlim = c(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k], smp_ale_frq[2, k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k], smp_ale_frq[2, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
         xlab = "Allele frequency",
         main = paste("Island ancestral allele at generation", smp_gen[k]))
    hist(pop_frq_pre_resmp[2, , k], breaks = seq(min(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), max(pop_frq_pst_resmp[2, , k], pop_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
    abline(v = smp_ale_frq[2, k], col = 'red', lty = 2, lwd = 2)
  } else {
    plot(0, type = "n", axes = FALSE, ann = FALSE)
  }

  hist_pst_resmp <- hist(pop_frq_pst_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  if (!any(is.nan(c(hist_pst_resmp$density, hist_pre_resmp$density)))) {
    hist(pop_frq_pst_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
         xlim = c(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k], smp_ale_frq[3, k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k], smp_ale_frq[3, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
         xlab = "Allele frequency",
         main = paste("Continent mutant allele at generation", smp_gen[k]))
    hist(pop_frq_pre_resmp[3, , k], breaks = seq(min(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), max(pop_frq_pst_resmp[3, , k], pop_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
    abline(v = smp_ale_frq[3, k], col = 'red', lty = 2, lwd = 2)
  } else {
    plot(0, type = "n", axes = FALSE, ann = FALSE)
  }

  hist_pst_resmp <- hist(pop_frq_pst_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  if (!any(is.nan(c(hist_pst_resmp$density, hist_pre_resmp$density)))) {
    hist(pop_frq_pst_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
         xlim = c(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k], smp_ale_frq[4, k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k], smp_ale_frq[4, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
         xlab = "Allele frequency",
         main = paste("Continent ancestral allele at generation", smp_gen[k]))
    hist(pop_frq_pre_resmp[4, , k], breaks = seq(min(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), max(pop_frq_pst_resmp[4, , k], pop_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
    abline(v = smp_ale_frq[4, k], col = 'red', lty = 2, lwd = 2)
  } else {
    plot(0, type = "n", axes = FALSE, ann = FALSE)
  }
}
dev.off()

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param sel_gen the starting time of natural selection
#' @param mig_gen the starting time of gene migration
#' @param ext_frq the mutant allele frequency (of the continent population)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and continent alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
mig_rat
pop_siz
sel_gen
mig_gen
ext_frq
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
gap_num <- 1e+02

system.time(OptNum <- calculateOptimalParticleNum(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num, OptNum,
     file = "./Output/Output v1.0/Test v1.0/TEST_OptNum.rda")

load("./Output/Output v1.0/Test v1.0/TEST_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_OptNum.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(opt_pcl_num, log_lik_sdv, type = 'b', lwd = 2,
     xlab = "Particle number", ylab = "Log-likelihood standard deviation",
     main = "Optimal particle number in the PMMH")
abline(h = 1.7, col = 'red', lty = 2, lwd = 2)
abline(h = 1.0, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param sel_gen the starting time of natural selection
#' @param mig_gen the starting time of gene migration
#' @param ext_frq the mutant allele frequency (of the continent population)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and continent alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

set.seed(test_seed)

sel_cof <- 0e-00
dom_par
mig_rat <- 0e-00
pop_siz
sel_gen <- min(smp_gen)
mig_gen <- ifelse(smp_gen[which(smp_cnt[2, ] > 0)[1]] - min(smp_gen) > 0, smp_gen[which(smp_cnt[2, ] > 0)[1] - 1], smp_gen[which(smp_cnt[2, ] > 0)[1]])
ext_frq
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, PMMH,
     file = "./Output/Output v1.0/Test v1.0/TEST_PMMH.rda")

load("./Output/Output v1.0/Test v1.0/TEST_PMMH.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_PMMH_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")
abline(h = sel_cof, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")
abline(h = sel_gen, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
abline(h = mig_rat, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
abline(h = mig_gen, col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 1e+04
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]

thn_num <- 8e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:round(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:round(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:round(length(mig_gen_chn) / thn_num)) * thn_num]

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

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_PMMH_Posterior.pdf", width = 16, height = 12)
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

########################################

#' Run the particle marginal Metropolis-Hastings within Gibbs (PMMHwGibbs)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param sel_gen the starting time of natural selection
#' @param mig_gen the starting time of gene migration
#' @param ext_frq the mutant allele frequency (of the continent population)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and continent alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH within Gibbs

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

set.seed(test_seed)

sel_cof <- 0e-00
dom_par
mig_rat <- 0e-00
pop_siz
sel_gen <- min(smp_gen)
mig_gen <- ifelse(smp_gen[which(smp_cnt[2, ] > 0)[1]] - min(smp_gen) > 0, smp_gen[which(smp_cnt[2, ] > 0)[1] - 1], smp_gen[which(smp_cnt[2, ] > 0)[1]])
ext_frq
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04

system.time(PMMHwGibbs <- cmprunPMMHwGibbs(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, PMMHwGibbs,
     file = "./Output/Output v1.0/Test v1.0/TEST_PMMHwGibbs.rda")

load("./Output/Output v1.0/Test v1.0/TEST_PMMHwGibbs.rda")

sel_cof_chn <- PMMHwGibbs$sel_cof_chn
sel_gen_chn <- PMMHwGibbs$sel_gen_chn
mig_rat_chn <- PMMHwGibbs$mig_rat_chn
mig_gen_chn <- PMMHwGibbs$mig_gen_chn

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_PMMHwGibbs_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")
abline(h = sel_cof, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")
abline(h = sel_gen, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
abline(h = mig_rat, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
abline(h = mig_gen, col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 1e+04
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]

thn_num <- 8e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:round(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:round(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:round(length(mig_gen_chn) / thn_num)) * thn_num]

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

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_PMMHwGibbs_Posterior.pdf", width = 16, height = 12)
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

########################################

#' Run the adaptive particle marginal Metropolis-Hastings (AdaptPMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param sel_gen the starting time of natural selection
#' @param mig_gen the starting time of gene migration
#' @param ext_frq the mutant allele frequency (of the continent population)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and continent alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

set.seed(test_seed)

sel_cof <- 0e-00
dom_par
mig_rat <- 0e-00
pop_siz
sel_gen <- min(smp_gen)
mig_gen <- ifelse(smp_gen[which(smp_cnt[2, ] > 0)[1]] - min(smp_gen) > 0, smp_gen[which(smp_cnt[2, ] > 0)[1] - 1], smp_gen[which(smp_cnt[2, ] > 0)[1]])
ext_frq
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(AdaptPMMH <- cmprunAdaptPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, AdaptPMMH,
     file = "./Output/Output v1.0/Test v1.0/TEST_AdaptPMMH.rda")

load("./Output/Output v1.0/Test v1.0/TEST_AdaptPMMH.rda")

sel_cof_chn <- AdaptPMMH$sel_cof_chn
sel_gen_chn <- AdaptPMMH$sel_gen_chn
mig_rat_chn <- AdaptPMMH$mig_rat_chn
mig_gen_chn <- AdaptPMMH$mig_gen_chn

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_AdaptPMMH_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")
abline(h = sel_cof, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")
abline(h = sel_gen, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
abline(h = mig_rat, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
abline(h = mig_gen, col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 1e+04
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]

thn_num <- 8e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:round(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:round(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:round(length(mig_gen_chn) / thn_num)) * thn_num]

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

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_AdaptPMMH_Posterior.pdf", width = 16, height = 12)
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

########################################

#' Run the adaptive particle marginal Metropolis-Hastings within Gibbs (AdaptPMMHwGibbs)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param sel_gen the starting time of natural selection
#' @param mig_gen the starting time of gene migration
#' @param ext_frq the mutant allele frequency (of the continent population)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and continent alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH within Gibbs
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

set.seed(test_seed)

sel_cof <- 0e-00
dom_par
mig_rat <- 0e-00
pop_siz
sel_gen <- min(smp_gen)
mig_gen <- ifelse(smp_gen[which(smp_cnt[2, ] > 0)[1]] - min(smp_gen) > 0, smp_gen[which(smp_cnt[2, ] > 0)[1] - 1], smp_gen[which(smp_cnt[2, ] > 0)[1]])
ext_frq
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04
stp_siz <- (1:itn_num)^(-2 / 3)
apt_rto <- 4e-01

system.time(AdaptPMMHwGibbs <- cmprunAdaptPMMHwGibbs(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto))

load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto, AdaptPMMHwGibbs,
     file = "./Output/Output v1.0/Test v1.0/TEST_AdaptPMMHwGibbs.rda")

load("./Output/Output v1.0/Test v1.0/TEST_AdaptPMMHwGibbs.rda")

sel_cof_chn <- AdaptPMMHwGibbs$sel_cof_chn
sel_gen_chn <- AdaptPMMHwGibbs$sel_gen_chn
mig_rat_chn <- AdaptPMMHwGibbs$mig_rat_chn
mig_gen_chn <- AdaptPMMHwGibbs$mig_gen_chn

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_AdaptPMMHwGibbs_Traceplot.pdf", width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot of the selection coefficient")
abline(h = sel_cof, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection time",
     main = "Trace plot of the selection time")
abline(h = sel_gen, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration rate",
     main = "Trace plot of the migration rate")
abline(h = mig_rat, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Migration time",
     main = "Trace plot of the migration time")
abline(h = mig_gen, col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 1e+04
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]

thn_num <- 8e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]
sel_gen_chn <- sel_gen_chn[(1:round(length(sel_gen_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:round(length(mig_rat_chn) / thn_num)) * thn_num]
mig_gen_chn <- mig_gen_chn[(1:round(length(mig_gen_chn) / thn_num)) * thn_num]

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

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)

pdf(file = "./Output/Output v1.0/Test v1.0/TEST_AdaptPMMHwGibbs_Posterior.pdf", width = 16, height = 12)
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

########################################

#' Run the Bayesian procedure for the inference of natural selection and gene migration
#' Parameter settings
#' @param method = "PMMH"/"PMMHwGibss" (return the results obtained from with the PMMH/PMMHwGibbs)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param sel_gen the starting time of natural selection
#' @param mig_gen the starting time of gene migration
#' @param ext_frq the mutant allele frequency (of the continent population)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and continent alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH/PMMH within Gibbs
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param adp_set = TRUE/FALSE (return the result with the adaptive setting or not)
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

# load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")
#
# set.seed(test_seed)
#
# method <- "PMMHwGibbs"
# sel_cof <- 0e-00
# dom_par
# mig_rat <- 0e-00
# pop_siz
# sel_gen <- min(smp_gen)
# mig_gen <- ifelse(smp_gen[which(smp_cnt[2, ] > 0)[1]] - min(smp_gen) > 0, smp_gen[which(smp_cnt[2, ] > 0)[1] - 1], smp_gen[which(smp_cnt[2, ] > 0)[1]])
# ext_frq
# smp_gen
# smp_siz
# smp_cnt
# ptn_num <- 5e+00
# pcl_num <- 1e+03
# itn_num <- 2e+04
# brn_num <- 5e+03
# thn_num <- 5e+00
# adp_set <- TRUE
# stp_siz <- (1:itn_num)^(-2 / 3)
# apt_rto <- 4e-01
#
# system.time(BayesianProcedure <- cmprunBayesianProcedure(method, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, adp_set, stp_siz, apt_rto))
#
# load("./Output/Output v1.0/Test v1.0/TEST_SimData.rda")
#
# save(method, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, adp_set, stp_siz, apt_rto, BayesianProcedure,
#      file = "./Output/Output v1.0/Test v1.0/TEST_BayesianProcedure.rda")
#
# load("./Output/Output v1.0/Test v1.0/TEST_BayesianProcedure.rda")
#
# sel_cof_chn <- BayesianProcedure$sel_cof_chn
# sel_gen_chn <- BayesianProcedure$sel_gen_chn
# mig_rat_chn <- BayesianProcedure$mig_rat_chn
# mig_gen_chn <- BayesianProcedure$mig_gen_chn
#
# sel_cof_est <- BayesianProcedure$sel_cof_est
# sel_gen_est <- BayesianProcedure$sel_gen_est
# mig_rat_est <- BayesianProcedure$mig_rat_est
# mig_gen_est <- BayesianProcedure$mig_gen_est
#
# mig_rat_hpd <- BayesianProcedure$mig_rat_hpd
# mig_gen_hpd <- BayesianProcedure$mig_gen_hpd
# sel_cof_hpd <- BayesianProcedure$sel_cof_hpd
# sel_geh_hpd <- BayesianProcedure$sel_gen_hpd
#
# pdf(file = "./Output/Output v1.0/Test v1.0/TEST_BayesianProcedure_Posterior.pdf", width = 16, height = 12)
# par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
# hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
#      xlab = "Selection coefficient",
#      main = "Posterior for the selection coefficient")
# lines(density(sel_cof_chn), lwd = 2, col = 'black')
# abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
# abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
# abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
# abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
#
# hist(sel_gen_chn, breaks = seq(min(sel_gen_chn), max(sel_gen_chn), length.out = 50), freq = FALSE,
#      xlab = "Selection coefficient",
#      main = "Posterior for the selection time")
# lines(density(sel_gen_chn), lwd = 2, col = 'black')
# abline(v = sel_gen, col = 'red', lty = 2, lwd = 2)
# abline(v = sel_gen_est, col = 'black', lty = 2, lwd = 2)
# abline(v = sel_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
# abline(v = sel_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
#
# hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE,
#      xlab = "Selection coefficient",
#      main = "Posterior for the migration rate")
# lines(density(mig_rat_chn), lwd = 2, col = 'black')
# abline(v = mig_rat, col = 'red', lty = 2, lwd = 2)
# abline(v = mig_rat_est, col = 'black', lty = 2, lwd = 2)
# abline(v = mig_rat_hpd[1], col = 'blue', lty = 2, lwd = 2)
# abline(v = mig_rat_hpd[2], col = 'blue', lty = 2, lwd = 2)
#
# hist(mig_gen_chn, breaks = seq(min(mig_gen_chn), max(mig_gen_chn), length.out = 50), freq = FALSE,
#      xlab = "Selection coefficient",
#      main = "Posterior for the migration time")
# lines(density(mig_gen_chn), lwd = 2, col = 'black')
# abline(v = mig_gen, col = 'red', lty = 2, lwd = 2)
# abline(v = mig_gen_est, col = 'black', lty = 2, lwd = 2)
# abline(v = mig_gen_hpd[1], col = 'blue', lty = 2, lwd = 2)
# abline(v = mig_gen_hpd[2], col = 'blue', lty = 2, lwd = 2)
# dev.off()

################################################################################
