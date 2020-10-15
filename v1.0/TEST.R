#' @title Inferring natural selection and gene migration in the evolution of chickens from ancient DNA data
#' @author Zhangyi He, Wenyang Lyu, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' version 2.0

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2018/HE2020-WFM-1L-DiffusApprox-PMMHwGibbs-Chicken")

source("./Output/Output v2.1/Test v2.1/Code v2-2.0chicken/RFUN.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")  
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection and migration
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param ext_frq the external frequency of the mutant allele of the continent population
#' @param int_frq the initial frequency of the mutant allele of the island population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory

sel_cof <- 5e-03
dom_par <- 5e-01
mig_rat <- 5e-03
pop_siz <- 5e+03
ext_frq <- 9e-01
int_frq <- c(2e-01, 8e-01, 0e-01, 0e-01)
int_gen <- 0
lst_gen <- 500

frq_pth <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A1^i haplotype generated with the Wright-Fisher model")
plot(k, frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A2^i haplotype generated with the Wright-Fisher model")
plot(k, frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A1^c haplotype generated with the Wright-Fisher model")
plot(k, frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Generation", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A2^c haplotype generated with the Wright-Fisher model")

########################################

#' Simulate the mutant allele frequency trajectory under the one-locus Wright-Fisher diffusion with selection migration using the Euler-Maruyama method
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param ext_frq the external frequency of the mutant allele of the continent population
#' @param int_frq the initial frequency of the mutant allele of the island population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method
#' @param data_augmentation = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- 5e-03
dom_par <- 5e-01
mig_rat <- 5e-03
pop_siz <- 5e+03
ext_frq <- 9e-01
int_frq <- c(2e-01, 8e-01, 0e-01, 0e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
plot(t, frq_pth[1, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A1^i haplotype generated with the Wright-Fisher diffusion")
plot(t, frq_pth[2, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A2^i haplotype generated with the Wright-Fisher diffusion")
plot(t, frq_pth[3, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A1^c haplotype generated with the Wright-Fisher diffusion")
plot(t, frq_pth[4, ], type = "l", lwd = 1.5,
     xlab = "Time", ylab = "Haplotype frequency",
     main = "A frequency trajectory of the A2^c haplotype generated with the Wright-Fisher diffusion")
########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- 1e-02
dom_par <- 5e-01
mig_rat <- 5e-03
pop_siz <- 5e+03
ext_frq <- 9e-01
int_frq <- c(2e-01, 8e-01, 0e-01, 0e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00
sim_num <- 1e+04

sim_frq_WFM <- matrix(NA, nrow = 4, ncol = sim_num)
sim_frq_WFD <- matrix(NA, nrow = 4, ncol = sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[, i] <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
  sim_frq_WFD[, i] <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)[, (lst_gen - int_gen) + 1]
}

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_frq_WFM, sim_frq_WFD, 
     file = "./Output/Output v2.1/Test v2.1/TEST_2L_WFM_vs_WFD.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_WFM_vs_WFD.rda")

pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_WFM_vs_WFD.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ])), 
     xlab = "Haplotype frequency", main = "Haplotype A1^i")
hist(sim_frq_WFD[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ])), 
     xlab = "Haplotype frequency", main = "Haplotype A2^i")
hist(sim_frq_WFD[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ])), 
     xlab = "Haplotype frequency", main = "Haplotype A1^c")
hist(sim_frq_WFD[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ])), 
     xlab = "Haplotype frequency", main = "Haplotype A2^c")
hist(sim_frq_WFD[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
title(paste("Histograms of the haplotype frequencies in generation", lst_gen, "under the Wright-Fisher model and the Wright-Fisher diffusion"), outer = TRUE)
dev.off()

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
#' @param ext_frq the external frequency of the mutant allele of the continent population
#' @param int_frq the initial haplotype frequencies of the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model 
model <- "WFM"
sel_cof <- 1e-02
dom_par <- 5e-01
mig_rat <- 5e-03
pop_siz <- 5e+03
ext_frq <- 9e-01
int_frq <- c(2e-01, 8e-01, 0e-01, 0e-01)
smp_gen <- (0:10) * 50
smp_siz <- rep(100, 11)
sel_gen <- 150
mig_gen <- 350

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, int_frq, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_hap_cnt <- sim_HMM_WFM$smp_hap_cnt
pop_hap_frq <- sim_HMM_WFM$pop_hap_frq
smp_ale_cnt <- sim_HMM_WFM$smp_ale_cnt
pop_ale_frq <- sim_HMM_WFM$pop_ale_frq

k <- min(smp_gen):max(smp_gen)
smp_ale_frq <- smp_ale_cnt %*% diag(1 / smp_siz)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[1, ], pop_ale_frq[1, ]), max(smp_ale_frq[1, ], pop_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant alleles generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[2, ], pop_ale_frq[2, ]), max(smp_ale_frq[2, ], pop_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the continent alleles generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion 
model <- "WFD"
sel_cof <- 1e-02
dom_par <- 5e-01
mig_rat <- 5e-03
pop_siz <- 5e+03
ext_frq <- 9e-01
int_frq <- c(2e-01, 8e-01, 0e-01, 0e-01)
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)
ptn_num <- 5e+00

sim_HMM_WFD <- cmpsimulateHMM(model, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, int_frq, smp_gen, smp_siz, ptn_num)
smp_gen <- sim_HMM_WFD$smp_gen
smp_siz <- sim_HMM_WFD$smp_siz
smp_hap_cnt <- sim_HMM_WFD$smp_hap_cnt
pop_hap_frq <- sim_HMM_WFD$pop_hap_frq
smp_ale_cnt <- sim_HMM_WFD$smp_ale_cnt
pop_ale_frq <- sim_HMM_WFD$pop_ale_frq

k <- min(smp_gen):max(smp_gen)
smp_ale_frq <- smp_ale_cnt %*% diag(1 / smp_siz)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[1, ], pop_ale_frq[1, ]), max(smp_ale_frq[1, ], pop_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant alleles generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[2, ], pop_ale_frq[2, ]), max(smp_ale_frq[2, ], pop_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the continent alleles generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model 
test_seed <- 7
set.seed(test_seed)

model <- "WFM"
sel_cof <- 1e-02
dom_par <- 5e-01
mig_rat <- 5e-03
pop_siz <- 5e+03
ext_frq <- 9e-01
int_frq <- c(2e-01, 8e-01, 0e-01, 0e-01)
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)
sel_gen <- 150
mig_gen <- 350

SimData <- cmpsimulateHMM(model, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, int_frq, smp_gen, smp_siz)
smp_gen <- SimData$smp_gen
smp_siz <- SimData$smp_siz
smp_hap_cnt <- SimData$smp_hap_cnt
pop_hap_frq <- SimData$pop_hap_frq
smp_ale_cnt <- SimData$smp_ale_cnt
pop_ale_frq <- SimData$pop_ale_frq

save(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_hap_cnt, pop_hap_frq, smp_ale_cnt, pop_ale_frq, 
     file = "./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

k <- min(smp_gen):max(smp_gen)
smp_ale_frq <- smp_ale_cnt %*% diag(1 / smp_siz)

pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_SimData.pdf", width = 20, height = 10)
par(mfrow = c(1, 2), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_ale_frq[1, ], smp_ale_frq[1, ]), max(pop_ale_frq[1, ], smp_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "Mutant alleles")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)

plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_ale_frq[2, ], smp_ale_frq[2, ]), max(pop_ale_frq[2, ], smp_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "Continent alleles")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)
title("A simulated dataset without missing values generated with the Wright-Fisher model", outer = TRUE)
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
#' @param ext_frq the external frequency of the mutant allele of the continent population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and ancestral alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

load("./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

set.seed(test_seed)

sel_cof <- 1e-02
dom_par <- 5e-01
mig_rat <- 5e-03
pop_siz <- 5e+03
ext_frq <- 9e-01
int_frq <- c(2e-01, 8e-01, 0e-01, 0e-01)
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)
sel_gen
mig_gen 
smp_ale_cnt
ptn_num <- 5e+00
pcl_num <- 5e+03

system.time(BPF <- cmprunBPF(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, BPF, 
     file = "./Output/Output v2.1/Test v2.1/TEST_2L_BPF.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}

pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_BPF_Likelihood.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:pcl_num, log(lik), type = 'l', 
     xlab = "Number of particles", ylab = "Log likelihood", 
     main = "Log likelihood estimated with the bootstrap particle filter")
dev.off()

smp_hap_frq <- smp_hap_cnt %*% diag(1 / smp_siz)
pop_hap_frq_pre_resmp <- BPF$pop_frq_pre_resmp
pop_hap_frq_pst_resmp <- BPF$pop_frq_pst_resmp

pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_BPF_Particle.pdf", width = 20, height = 55)
par(mfrow = c(11, 4), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
for (k in 1:length(smp_gen)) {
  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist(pop_hap_frq_pst_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k], smp_hap_frq[1, k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k], smp_hap_frq[1, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)), 
       xlab = "Haplotype frequency", 
       main = paste("Haplotype A1^i in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[1, k], col = 'red', lty = 2, lwd = 2)
  
  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)  
  hist(pop_hap_frq_pst_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k], smp_hap_frq[2, k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k], smp_hap_frq[2, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)), 
       xlab = "Haplotype frequency", 
       main = paste("Haplotype A2^i in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[2, k], col = 'red', lty = 2, lwd = 2)
  
  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist(pop_hap_frq_pst_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k], smp_hap_frq[3, k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k], smp_hap_frq[3, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)), 
       xlab = "Haplotype frequency", 
       main = paste("Haplotype A1^c in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[3, k], col = 'red', lty = 2, lwd = 2)
  
  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist(pop_hap_frq_pst_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k], smp_hap_frq[4, k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k], smp_hap_frq[4, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)), 
       xlab = "Haplotype frequency", 
       main = paste("Haplotype A2^c in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[4, k], col = 'red', lty = 2, lwd = 2)
}
title("Histograms of the pre- and post-resampling particles", outer = TRUE)
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
#' @param ext_frq the external frequency of the mutant allele of the continent population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

load("./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
rec_rat
pop_siz
smp_gen
smp_siz
smp_ale_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
gap_num <- 1e+02

system.time(OptNum <- calculateOptimalParticleNum(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, gap_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, gap_num, OptNum,
     file = "./Output/Output v2.1/Test v2.1/TEST_2L_OptNum.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_OptNum.pdf", width = 12, height = 9)
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
#' @param ext_frq the external frequency of the mutant allele of the continent population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings

load("./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

set.seed(test_seed)

sel_cof <- 1e-02
dom_par <- 5e-01
mig_rat <- 5e-03
pop_siz <- 5e+03
ext_frq <- 9e-01
int_frq <- c(2e-01, 8e-01, 0e-01, 0e-01)
sel_gen 
mig_gen 
smp_gen
smp_siz
smp_ale_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 1e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, itn_num, PMMH, 
     file = "./Output/Output v2.1/Test v2.1/TEST_2L_PMMH.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_PMMH.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

sel_cof_chn <- PMMH$sel_cof_chn
sel_gen_chn <- PMMH$sel_gen_chn
mig_rat_chn <- PMMH$mig_rat_chn
mig_gen_chn <- PMMH$mig_gen_chn
frq_pth_chn <- PMMH$frq_pth_chn
pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_PMMH_Traceplot.pdf", width = 20, height = 10)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", 
     main = "Trace plot of the selection coefficient")
abline(h = sel_cof, col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, sel_gen_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection generation", 
     main = "Trace plot of the selection generation")
abline(h = sel_gen, col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, mig_rat_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Migration rate", 
     main = "Trace plot of the migration rate")
abline(h = mig_rat, col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, mig_gen_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Migration generation", 
     main = "Trace plot of the migration generation")
abline(h = mig_gen, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, cumsum(sel_cof_chn[1:itn_num])/(1:itn_num), type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", 
     main = "Mean ")
abline(h = sel_cof, col = 'red', lty = 2, lwd = 2)

plot(1:itn_num, cumsum(mig_rat_chn[1:itn_num])/(1:itn_num), type = 'l', 
     xlab = "Iteration", ylab = "Migration rate", 
     main = "Mean ")
abline(h = mig_rat, col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 1e+03
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]

thn_num <- 8e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]
mig_rat_chn <- mig_rat_chn[(1:round(length(mig_rat_chn) / thn_num)) * thn_num]

grd_num <- 1e+03
sel_cof_mig_rat_pdf <- kde2d(sel_cof_chn, mig_rat_chn, n = grd_num)
pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_PMMH_Posterior.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
image(sel_cof_mig_rat_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32), 
      xlab = "Selection coefficient", ylab = "Migration", 
      main = "Posterior for the selection coefficient and the migration")
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(h = mig_rat, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the Bayesian procedure for the inference of natural selection
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mig_rat the migration rate
#' @param pop_siz the number of the diploid individuals in the population
#' @param sel_gen the starting time of natural selection
#' @param mig_gen the starting time of gene migration
#' @param ext_frq the external frequency of the mutant allele of the continent population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param grd_num the number of the grids in the kernel density estimation

load("./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

set.seed(test_seed)

sel_cof <- 1e-02
dom_par <- 5e-01
mig_rat <- 5e-03
pop_siz <- 5e+03
ext_frq <- 9e-01
int_frq <- c(2e-01, 8e-01, 0e-01, 0e-01)
sel_gen 
mig_gen 
smp_gen
smp_siz
smp_ale_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04
brn_num <- 1e+04
thn_num <- 8e+00
grd_num <- 1e+03

system.time(BayesianProcedure <- cmprunBayesianProcedure(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num))

save(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num, BayesianProcedure, 
     file = "./Output/Output v2.1/Test v2.1/TEST_2L_BayesianProcedure.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_BayesianProcedure.rda")

load("./Output/Output v2.1/Test v2.1/TEST_2L_SimData.rda")

sel_cof_chn <- BayesianProcedure$sel_cof_chn
mig_rat_chn <- BayesianProcedure$mig_rat_chn

sel_cof_mig_rat_pdf <- BayesianProcedure$sel_cof_mig_rat_pdf

sel_cof_map <- BayesianProcedure$sel_cof_map
mig_rat_map <- BayesianProcedure$mig_rat_map

sel_cof_mmse <- BayesianProcedure$sel_cof_mmse
mig_rat_mmse <- BayesianProcedure$mig_rat_mmse

mig_rat_hpd <- BayesianProcedure$mig_rat_hpd
sel_cof_hpd <- BayesianProcedure$sel_cof_hpd

pdf(file = "./Output/Output v2.1/Test v2.1/TEST_2L_BayesianProcedure_Posterior.pdf", width = 20, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
image(sel_cof_mig_rat_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32), 
      xlab = "Selection coefficient", ylab = "Migration rate", 
      main = "Joint posterior for the selection coefficient and the migration rate")
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(h = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_map, col = 'black', lty = 4, lwd = 2)
abline(h = mig_rat_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_mmse, col = 'black', lty = 2, lwd = 2)
abline(h = mig_rat_mmse, col = 'black', lty = 2, lwd = 2)

hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE, 
     xlab = "Selection coefficient", 
     main = "Marginal posterior for the selection coefficients")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd, col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd, col = 'blue', lty = 2, lwd = 2)

hist(mig_rat_chn, breaks = seq(min(mig_rat_chn), max(mig_rat_chn), length.out = 50), freq = FALSE, 
     xlab = "Migration rate", 
     main = "Marginal posterior for the migration rate")
lines(density(mig_rat_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = mig_rat_map, col = 'black', lty = 4, lwd = 2)
abline(v = mig_rat_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = mig_rat_hpd, col = 'blue', lty = 2, lwd = 2)
abline(v = mig_rat_hpd, col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################