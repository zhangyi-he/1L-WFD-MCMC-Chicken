#' @title Inferring natural selection and gene migration in the evolution of chickens from ancient DNA data
#' @author Zhangyi He, Wenyang Lyu, Xiaoyang Dai, Sile Hu, Mark Beaumont, Feng Yu

#' version 1.1

#' R functions

# install.packages("pdfCluster")
library("pdfCluster")

#install.packages("MASS")
library("MASS")

#install.packages("coda")
library("coda")

#install.packages("inline")
library("inline")
#install.packages("Rcpp")
library("Rcpp")
#install.packages("RcppArmadillo")
library("RcppArmadillo")

#install.packages("compiler")
library("compiler")
#enableJIT(1)

# call C++ functions
sourceCpp("./Code/Code v1.0/Code v1.1/CFUN.cpp")

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

#' Standard version
simulateWFM <- function(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen) {
  fts_mat <- calculateFitnessMat_arma(sel_cof, dom_par)
  frq_pth <- simulateWFM_arma(fts_mat, sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen)
  frq_pth <- as.matrix(frq_pth)

  return(frq_pth)
}
#' Compiled version
cmpsimulateWFM <- cmpfun(simulateWFM)

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

#' Standard version
simulateWFD <- function(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  frq_pth <- simulateWFD_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num)
  frq_pth <- as.matrix(frq_pth)

  if (dat_aug == FALSE) {
    return(frq_pth[, (0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(frq_pth)
  }
}
#' Compiled version
cmpsimulateWFD <- cmpfun(simulateWFD)

########################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the samples obtained from the underlying population evolving according to the WFM or the WFD)
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

#' Standard version
simulateHMM <- function(model = "WFM", sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, int_frq, smp_gen, smp_siz, cna_gen = NULL, ...) {
  int_gen <- min(sel_gen, mig_gen, smp_gen)
  lst_gen <- max(sel_gen, mig_gen, smp_gen)

  # generate the population allele frequency trajectories
  if (model == "WFM") {
    if (sel_gen < mig_gen) {
      if (sel_gen > int_gen && mig_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen)
        pop_frq_stg_2 <- cmpsimulateWFM(sel_cof, dom_par, 0, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], sel_gen, mig_gen)
        pop_frq_stg_3 <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], mig_gen, lst_gen)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen == int_gen && mig_gen < lst_gen) {
        pop_frq_stg_2 <- cmpsimulateWFM(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, sel_gen, mig_gen)
        pop_frq_stg_3 <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], mig_gen, lst_gen)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen > int_gen && mig_gen == lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen)
        pop_frq_stg_2 <- cmpsimulateWFM(sel_cof, dom_par, 0, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], sel_gen, mig_gen)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)])
      }
      if (sel_gen == int_gen && mig_gen == lst_gen) {
        pop_ale_frq <- cmpsimulateWFM(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, sel_gen, mig_gen)
      }
    } else if (sel_gen > mig_gen) {
      if (mig_gen > int_gen && sel_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen)
        pop_frq_stg_2 <- cmpsimulateWFM(0, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, sel_gen)
        pop_frq_stg_3 <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], sel_gen, lst_gen)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(mig_gen - int_gen + 1)], pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - sel_gen + 1)])
      }
      if (mig_gen == int_gen && sel_gen < lst_gen) {
        pop_frq_stg_2 <- cmpsimulateWFM(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, mig_gen, sel_gen)
        pop_frq_stg_3 <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], sel_gen, lst_gen)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - sel_gen + 1)])
      }
      if (mig_gen > int_gen && sel_gen == lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen)
        pop_frq_stg_2 <- cmpsimulateWFM(0, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, sel_gen)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(mig_gen - int_gen + 1)], pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)])
      }
      if (mig_gen == int_gen && sel_gen == lst_gen) {
        pop_ale_frq <- cmpsimulateWFM(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, mig_gen, sel_gen)
      }
    } else {
      if (sel_gen > int_gen && mig_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen)
        pop_frq_stg_2 <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, lst_gen)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen == int_gen) {
        pop_ale_frq <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, sel_gen, lst_gen)
      }
      if (mig_gen == lst_gen) {
        pop_ale_frq <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen)
      }
    }

    pop_frq <- matrix(NA, nrow = 2, ncol = (lst_gen - int_gen) + 1)
    pop_frq[1, ] <- pop_ale_frq[1, ] + pop_ale_frq[3, ]
    pop_frq[2, ] <- pop_ale_frq[3, ] + pop_ale_frq[4, ]
  }
  if (model == "WFD") {
    if (sel_gen < mig_gen) {
      if (sel_gen > int_gen && mig_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen, ptn_num, dat_aug = FALSE)
        pop_frq_stg_2 <- cmpsimulateWFD(sel_cof, dom_par, 0, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], sel_gen, mig_gen, ptn_num, dat_aug = FALSE)
        pop_frq_stg_3 <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], mig_gen, lst_gen, ptn_num, dat_aug = FALSE)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen == int_gen && mig_gen < lst_gen) {
        pop_frq_stg_2 <- cmpsimulateWFD(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, sel_gen, mig_gen, ptn_num, dat_aug = FALSE)
        pop_frq_stg_3 <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], mig_gen, lst_gen, ptn_num, dat_aug = FALSE)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen > int_gen && mig_gen == lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen, ptn_num, dat_aug = FALSE)
        pop_frq_stg_2 <- cmpsimulateWFD(sel_cof, dom_par, 0, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], sel_gen, mig_gen, ptn_num, dat_aug = FALSE)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)])
      }
      if (sel_gen == int_gen && mig_gen == lst_gen) {
        pop_ale_frq <- cmpsimulateWFD(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, sel_gen, mig_gen, ptn_num, dat_aug = FALSE)
      }
    } else if (sel_gen > mig_gen) {
      if (mig_gen > int_gen && sel_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen, ptn_num, dat_aug = FALSE)
        pop_frq_stg_2 <- cmpsimulateWFD(0, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, sel_gen, ptn_num, dat_aug = FALSE)
        pop_frq_stg_3 <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], sel_gen, lst_gen, ptn_num, dat_aug = FALSE)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(mig_gen - int_gen + 1)], pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - sel_gen + 1)])
      }
      if (mig_gen == int_gen && sel_gen < lst_gen) {
        pop_frq_stg_2 <- cmpsimulateWFD(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, mig_gen, sel_gen, ptn_num, dat_aug = FALSE)
        pop_frq_stg_3 <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], sel_gen, lst_gen, ptn_num, dat_aug = FALSE)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - sel_gen + 1)])
      }
      if (mig_gen > int_gen && sel_gen == lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen, ptn_num, dat_aug = FALSE)
        pop_frq_stg_2 <- cmpsimulateWFD(0, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, sel_gen, ptn_num, dat_aug = FALSE)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(mig_gen - int_gen + 1)], pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)])
      }
      if (mig_gen == int_gen && sel_gen == lst_gen) {
        pop_ale_frq <- cmpsimulateWFD(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, mig_gen, sel_gen, ptn_num, dat_aug = FALSE)
      }
    } else {
      if (sel_gen > int_gen && mig_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen, ptn_num, dat_aug = FALSE)
        pop_frq_stg_2 <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, lst_gen, ptn_num, dat_aug = FALSE)

        pop_ale_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen == int_gen) {
        pop_ale_frq <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, sel_gen, lst_gen, ptn_num, dat_aug = FALSE)
      }
      if (mig_gen == lst_gen) {
        pop_ale_frq <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen, ptn_num, dat_aug = FALSE)
      }
    }

    pop_frq <- matrix(NA, nrow = 2, ncol = (lst_gen - int_gen) + 1)
    pop_frq[1, ] <- pop_ale_frq[1, ] + pop_ale_frq[3, ]
    pop_frq[2, ] <- pop_ale_frq[3, ] + pop_ale_frq[4, ]
  }

  # generate the sample allele counts at all sampling time points
  smp_ale_cnt <- matrix(NA, nrow = 4, ncol = length(smp_gen))
  smp_ale_frq <- matrix(NA, nrow = 4, ncol = length(smp_gen))
  smp_cnt <- matrix(NA, nrow = 2, ncol = length(smp_gen))
  smp_frq <- matrix(NA, nrow = 2, ncol = length(smp_gen))
  for (k in 1:length(smp_gen)) {
    smp_ale_cnt[, k] <- rmultinom(1, size = smp_siz[k], prob = pop_ale_frq[, which(int_gen:lst_gen %in% smp_gen)[k]])
    smp_ale_frq[, k] <- smp_ale_cnt[, k] / smp_siz[k]
    smp_cnt[1, k] <- smp_ale_cnt[1, k] + smp_ale_cnt[3, k]
    smp_cnt[2, k] <- smp_ale_cnt[3, k] + smp_ale_cnt[4, k]
    smp_frq[, k] <- smp_cnt[, k] / smp_siz[k]
  }

  if (is.null(cna_gen)) {
    return(list(smp_gen = smp_gen,
                smp_siz = smp_siz,
                smp_cnt = smp_cnt,
                smp_frq = smp_frq,
                pop_frq = pop_frq,
                smp_ale_cnt = smp_ale_cnt,
                smp_ale_frq = smp_ale_frq,
                pop_ale_frq = pop_ale_frq))
  } else {
    smp_cnt[2, which(smp_gen %in% cna_gen)] <- NA
    smp_frq[2, which(smp_gen %in% cna_gen)] <- NA

    return(list(smp_gen = smp_gen,
                smp_siz = smp_siz,
                smp_cnt = smp_cnt,
                smp_frq = smp_frq,
                pop_frq = pop_frq,
                smp_ale_cnt = smp_ale_cnt,
                smp_ale_frq = smp_ale_frq,
                pop_ale_frq = pop_ale_frq))
  }
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

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

#' Standard version
runBPF <- function(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num) {
  if (any(is.na(smp_cnt[2, ]))) {
    smp_cnt[2, which(is.na(smp_cnt[2, ]))] <- -1
  }

  # run the BPF
  BPF <- runBPF_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num)

  return(list(lik = BPF$lik,
              wght = BPF$wght,
              pop_frq_pre_resmp = BPF$part_pre_resmp,
              pop_frq_pst_resmp = BPF$part_pst_resmp))
}
#' Compiled version
cmprunBPF <- cmpfun(runBPF)

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

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num) {
  if (any(is.na(smp_cnt[2, ]))) {
    smp_cnt[2, which(is.na(smp_cnt[2, ]))] <- -1
  }

  # calculate the optimal particle number
  OptNum <- calculateOptimalParticleNum_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num)

  return(list(opt_pcl_num = as.vector(OptNum$opt_pcl_num),
              log_lik_sdv = as.vector(OptNum$log_lik_sdv)))
}
#' Compiled version
cmpcalculateOptimalParticleNum <- cmpfun(calculateOptimalParticleNum)

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
#' @param bck_smp = TRUE/FALSE (return the result with the blockwise sampling or not)

#' Standard version
runPMMH <- function(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, bck_smp = TRUE) {
  if (any(is.na(smp_cnt[2, ]))) {
    smp_cnt[2, which(is.na(smp_cnt[2, ]))] <- -1
  }

  # run the PMMH
  if (bck_smp == TRUE) {
    PMMH <- runPMMH_BlockSmp_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  } else {
    PMMH <- runPMMH_ComponentSmp_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  }

  return(list(sel_cof_chn = as.vector(PMMH$sel_cof_chn),
              sel_gen_chn = as.vector(PMMH$sel_gen_chn),
              mig_rat_chn = as.vector(PMMH$mig_rat_chn),
              mig_gen_chn = as.vector(PMMH$mig_gen_chn)))
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

########################################

#' Run the Bayesian procedure for the inference of natural selection and gene migration
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
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param bck_smp = TRUE/FALSE (return the result with the blockwise sampling or not)

#' Standard version
runBayesianProcedure <- function(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, bck_smp = TRUE) {
  if (any(is.na(smp_cnt[2, ]))) {
    smp_cnt[2, which(is.na(smp_cnt[2, ]))] <- -1
  }

  # run the PMMH
  if (bck_smp == TRUE) {
    PMMH <- runPMMH_BlockSmp_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  } else {
    PMMH <- runPMMH_ComponentSmp_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  }

  sel_cof_chn <- as.vector(PMMH$sel_cof_chn)
  sel_gen_chn <- as.vector(PMMH$sel_gen_chn)
  mig_rat_chn <- as.vector(PMMH$mig_rat_chn)
  mig_gen_chn <- as.vector(PMMH$mig_gen_chn)

  # burn-in and thinning
  # brn_num <- floor(0.8 * length(sel_cof_chn))
  sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
  sel_cof_chn <- sel_cof_chn[(1:floor(length(sel_cof_chn) / thn_num)) * thn_num]
  # brn_num <- floor(0.8 * length(sel_gen_chn))
  sel_gen_chn <- sel_gen_chn[brn_num:length(sel_gen_chn)]
  sel_gen_chn <- sel_gen_chn[(1:floor(length(sel_gen_chn) / thn_num)) * thn_num]
  # brn_num <- floor(0.8 * length(mig_rat_chn))
  mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
  mig_rat_chn <- mig_rat_chn[(1:floor(length(mig_rat_chn) / thn_num)) * thn_num]
  # brn_num <- floor(0.8 * length(mig_gen_chn))
  mig_gen_chn <- mig_gen_chn[brn_num:length(mig_gen_chn)]
  mig_gen_chn <- mig_gen_chn[(1:floor(length(mig_gen_chn) / thn_num)) * thn_num]

  # MAP estimates for the population genetic parameters
  smp <- data.frame(sel_cof_chn, sel_gen_chn, mig_rat_chn, mig_gen_chn)
  # colnames(smp) <- c("selection coefficient", "selection time", "migration rate", "migration time")
  pdf <- kepdf(smp, bwtype = "adaptive")
  est <- pdf@eval.points[which(pdf@estimate == max(pdf@estimate))[1], ]
  sel_cof_est <- est[1]
  sel_gen_est <- est[2]
  mig_rat_est <- est[3]
  mig_gen_est <- est[4]

  # # MMSE estimates for the population genetic parameters
  # sel_cof_est <- mean(sel_cof_chn)
  # sel_gen_est <- mean(sel_gen_chn)
  # mig_rat_est <- mean(mig_rat_chn)
  # mig_gen_est <- mean(mig_gen_chn)

  # 95% HPD intervals for the population genetic parameters
  sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
  sel_gen_hpd <- HPDinterval(as.mcmc(sel_gen_chn), prob = 0.95)
  mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
  mig_gen_hpd <- HPDinterval(as.mcmc(mig_gen_chn), prob = 0.95)

  return(list(sel_cof_est = sel_cof_est,
              sel_cof_hpd = sel_cof_hpd,
              sel_cof_chn = sel_cof_chn,
              sel_gen_est = sel_gen_est,
              sel_gen_hpd = sel_gen_hpd,
              sel_gen_chn = sel_gen_chn,
              mig_rat_est = mig_rat_est,
              mig_rat_hpd = mig_rat_hpd,
              mig_rat_chn = mig_rat_chn,
              mig_gen_est = mig_gen_est,
              mig_gen_hpd = mig_gen_hpd,
              mig_gen_chn = mig_gen_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
