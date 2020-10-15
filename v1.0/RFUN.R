#' @title Inferring natural selection and gene migration in the evolution of chickens from ancient DNA data
#' @author Zhangyi He, Wenyang Lyu, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' version 2.0

#' R functions

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
sourceCpp("./Output/Output v2.1/Test v2.1/Code v2-2.0chicken/CFUN.cpp")

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

#' Standard version
simulateWFM <- function(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen) {
  fts_mat <- calculateFitnessMat_2L_arma(sel_cof, dom_par)
  frq_pth <- simulateWFM_2L_arma(fts_mat, sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen)
  frq_pth <- as.matrix(frq_pth)
  
  return(frq_pth)
}
#' Compiled version
cmpsimulateWFM <- cmpfun(simulateWFM)

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

#' Standard version
simulateWFD <- function(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE) {
  frq_pth <- simulateWFD_2L_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num)
  frq_pth <- as.matrix(frq_pth)
  
  if (data_augmentation == FALSE) {
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

#' Standard version
simulateHMM <- function(model, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, int_frq, smp_gen, smp_siz, ...) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)
  
  # generate the population haplotype frequency trajectories and the population allele frequency trajectories
  if (model == "WFM") {
    if (sel_gen < mig_gen) {
      if (sel_gen > int_gen && mig_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen)
        pop_frq_stg_2 <- cmpsimulateWFM(sel_cof, dom_par, 0, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], sel_gen, mig_gen)
        pop_frq_stg_3 <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], mig_gen, lst_gen)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen == int_gen && mig_gen < lst_gen) {
        pop_frq_stg_2 <- cmpsimulateWFM(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, sel_gen, mig_gen)
        pop_frq_stg_3 <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], mig_gen, lst_gen)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen > int_gen && mig_gen == lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen)
        pop_frq_stg_2 <- cmpsimulateWFM(sel_cof, dom_par, 0, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], sel_gen, mig_gen)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)])
      }
      if (sel_gen == int_gen && mig_gen == lst_gen) {
        pop_frq <- cmpsimulateWFM(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, sel_gen, mig_gen)
      }
    }
    if (mig_gen < sel_gen) {
      if (mig_gen > int_gen && sel_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen)
        pop_frq_stg_2 <- cmpsimulateWFM(0, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, sel_gen)
        pop_frq_stg_3 <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], sel_gen, lst_gen)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(mig_gen - int_gen + 1)], pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - sel_gen + 1)])
      }
      if (mig_gen == int_gen && sel_gen < lst_gen) {
        pop_frq_stg_2 <- cmpsimulateWFM(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, mig_gen, sel_gen)
        pop_frq_stg_3 <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], sel_gen, lst_gen)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - sel_gen + 1)])
      }
      if (mig_gen > int_gen && sel_gen == lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen)
        pop_frq_stg_2 <- cmpsimulateWFM(0, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, sel_gen)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(mig_gen - int_gen + 1)], pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)])
      }
      if (mig_gen == int_gen && sel_gen == lst_gen) {
        pop_frq <- cmpsimulateWFM(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, mig_gen, sel_gen)
      }
    }
    if (sel_gen == mig_gen) {
      if (sel_gen > int_gen && mig_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen)
        pop_frq_stg_2 <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, lst_gen)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen == int_gen) {
        pop_frq <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, sel_gen, lst_gen)
      }
      if (mig_gen == lst_gen) {
        pop_frq <- cmpsimulateWFM(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen)
      }
    }
  }
  if (model == "WFD") {
    if (sel_gen < mig_gen) {
      if (sel_gen > int_gen && mig_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen, ptn_num, data_augmentation = FALSE)
        pop_frq_stg_2 <- cmpsimulateWFD(sel_cof, dom_par, 0, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], sel_gen, mig_gen, ptn_num, data_augmentation = FALSE)
        pop_frq_stg_3 <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], mig_gen, lst_gen, ptn_num, data_augmentation = FALSE)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen == int_gen && mig_gen < lst_gen) {
        pop_frq_stg_2 <- cmpsimulateWFD(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, sel_gen, mig_gen, ptn_num, data_augmentation = FALSE)
        pop_frq_stg_3 <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], mig_gen, lst_gen, ptn_num, data_augmentation = FALSE)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen > int_gen && mig_gen == lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen, ptn_num, data_augmentation = FALSE)
        pop_frq_stg_2 <- cmpsimulateWFD(sel_cof, dom_par, 0, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], sel_gen, mig_gen, ptn_num, data_augmentation = FALSE)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(mig_gen - sel_gen + 1)])
      }
      if (sel_gen == int_gen && mig_gen == lst_gen) {
        pop_frq <- cmpsimulateWFD(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, sel_gen, mig_gen, ptn_num, data_augmentation = FALSE)
      }
    }
    if (mig_gen < sel_gen) {
      if (mig_gen > int_gen && sel_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen, ptn_num, data_augmentation = FALSE)
        pop_frq_stg_2 <- cmpsimulateWFD(0, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, sel_gen, ptn_num, data_augmentation = FALSE)
        pop_frq_stg_3 <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], sel_gen, lst_gen, ptn_num, data_augmentation = FALSE)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(mig_gen - int_gen + 1)], pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - sel_gen + 1)])
      }
      if (mig_gen == int_gen && sel_gen < lst_gen) {
        pop_frq_stg_2 <- cmpsimulateWFD(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, mig_gen, sel_gen, ptn_num, data_augmentation = FALSE)
        pop_frq_stg_3 <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_2[, ncol(pop_frq_stg_2)], sel_gen, lst_gen, ptn_num, data_augmentation = FALSE)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)], pop_frq_stg_3[, 2:(lst_gen - sel_gen + 1)])
      }
      if (mig_gen > int_gen && sel_gen == lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen, ptn_num, data_augmentation = FALSE)
        pop_frq_stg_2 <- cmpsimulateWFD(0, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, sel_gen, ptn_num, data_augmentation = FALSE)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(mig_gen - int_gen + 1)], pop_frq_stg_2[, 2:(sel_gen - mig_gen + 1)])
      }
      if (mig_gen == int_gen && sel_gen == lst_gen) {
        pop_frq <- cmpsimulateWFD(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, mig_gen, sel_gen, ptn_num, data_augmentation = FALSE)
      }
    }
    if (sel_gen == mig_gen) {
      if (sel_gen > int_gen && mig_gen < lst_gen) {
        pop_frq_stg_1 <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen, ptn_num, data_augmentation = FALSE)
        pop_frq_stg_2 <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, pop_frq_stg_1[, ncol(pop_frq_stg_1)], mig_gen, lst_gen, ptn_num, data_augmentation = FALSE)
        
        pop_frq <- cbind(int_frq, pop_frq_stg_1[, 2:(sel_gen - int_gen + 1)], pop_frq_stg_2[, 2:(lst_gen - mig_gen + 1)])
      }
      if (sel_gen == int_gen) {
        pop_frq <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, sel_gen, lst_gen, ptn_num, data_augmentation = FALSE)
      }
      if (mig_gen == lst_gen) {
        pop_frq <- cmpsimulateWFD(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen, ptn_num, data_augmentation = FALSE)
      }
    }
  }
  #if (model == "WFM") {
  #  pop_hap_frq <- cmpsimulateWFM(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen)
  #}
  #if (model == "WFD") {
  #  pop_hap_frq <- cmpsimulateWFD(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)
  #}
  pop_ale_frq <- matrix(NA, nrow = 2, ncol = (lst_gen - int_gen) + 1)
  pop_ale_frq[1, ] <- pop_frq[1, ] + pop_frq[3, ]
  pop_ale_frq[2, ] <- pop_frq[3, ] + pop_frq[4, ]
  
  # generate the sample haplotype counts and the sample allele counts at all sampling time points
  smp_hap_cnt <- matrix(NA, nrow = 4, ncol = length(smp_gen))
  smp_ale_cnt <- matrix(NA, nrow = 2, ncol = length(smp_gen))
  for (k in 1:length(smp_gen)) {
    smp_hap_cnt[, k] <- rmultinom(1, size = smp_siz[k], prob = pop_frq[, smp_gen[k] - int_gen + 1])
    smp_ale_cnt[1, k] <- smp_hap_cnt[1, k] + smp_hap_cnt[3, k]
    smp_ale_cnt[2, k] <- smp_hap_cnt[3, k] + smp_hap_cnt[4, k]
  }
  
  return(list(smp_gen = smp_gen, 
              smp_siz = smp_siz, 
              smp_hap_cnt = smp_hap_cnt, 
              smp_ale_cnt = smp_ale_cnt, 
              pop_hap_frq = pop_frq, 
              pop_ale_frq = pop_ale_frq))
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
#' @param ext_frq the external frequency of the mutant allele of the continent population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' Standard version
runBPF <- function(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num) {
  
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
#' @param ext_frq the external frequency of the mutant allele of the continent population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num) {

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
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings

#' Standard version
runPMMH <- function(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num) {
  
  # run the PMMH
  PMMH <- runPMMH_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  
  return(list(sel_cof_chn = as.vector(PMMH$sel_cof_chn),
              sel_gen_chn = as.vector(PMMH$sel_gen_chn),
              mig_rat_chn = as.vector(PMMH$mig_rat_chn),
              mig_gen_chn = as.vector(PMMH$mig_gen_chn),
              log_lik_chn = as.vector(PMMH$log_lik_chn),
              frq_pth_chn = as.array(PMMH$frq_pth_chn)))
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

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

#' Standard version
runBayesianProcedure <- function(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num) {

  # run the PMMH
  PMMH <- runPMMH_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  
  # burn-in and thinning
  sel_cof_chn <- as.vector(PMMH$sel_cof_chn)
  sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
  sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]
  mig_rat_chn <- as.vector(PMMH$mig_rat_chn)
  mig_rat_chn <- mig_rat_chn[brn_num:length(mig_rat_chn)]
  mig_rat_chn <- mig_rat_chn[(1:round(length(mig_rat_chn) / thn_num)) * thn_num]
  
  # MAP estimates for the selection coefficients
  if (length(sel_cof_chn) < 1e+05) {
    sel_cof_mig_rat_pdf <- kde2d(sel_cof_chn, mig_rat_chn, n = grd_num)
    sel_cof_grd <- sel_cof_mig_rat_pdf$x
    mig_rat_grd <- sel_cof_mig_rat_pdf$y
  } else {
    sel_cof_mig_rat_pdf <- kde2d(tail(sel_cof_grd, 1e+05), tail(mig_rat_grd, 1e+05), n = grd_num)
    sel_cof_grd <- sel_cof_mig_rat_pdf$x
    mig_rat_grd <- sel_cof_mig_rat_pdf$y
  }
  sel_cof_map <- sel_cof_grd[which(sel_cof_mig_rat_pdf$z == max(sel_cof_mig_rat_pdf$z), arr.ind = TRUE)[1]]
  mig_rat_map <- mig_rat_grd[which(sel_cof_mig_rat_pdf$z == max(sel_cof_mig_rat_pdf$z), arr.ind = TRUE)[2]]
  
  # MMSE estimates for the selection coefficients
  sel_cof_mmse <- mean(sel_cof_chn)
  mig_rat_mmse <- mean(mig_rat_chn)
  
  # 95% HPD intervals for the selection coefficients
  sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
  mig_rat_hpd <- HPDinterval(as.mcmc(mig_rat_chn), prob = 0.95)
  
  return(list(sel_cof_mig_rat_pdf = sel_cof_mig_rat_pdf, 
              sel_cof_map = sel_cof_map, 
              mig_rat_map = mig_rat_map, 
              sel_cof_mmse = sel_cof_mmse, 
              mig_rat_mmse = mig_rat_mmse, 
              sel_cof_hpd = sel_cof_hpd, 
              mig_rat_hpd = mig_rat_hpd, 
              sel_cof_chn = sel_cof_chn, 
              mig_rat_chn = mig_rat_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
