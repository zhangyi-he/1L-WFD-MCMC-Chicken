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

#' call R functions
source("./Code/Code v1.0/Code v1.2/RFUN.R")

####################################################################################################

#' @section Section 1. Introduction

################################################################################

#' @section Section 2. Materials and Methods

#' @section Section 2.1. Wright-Fisher model

#' @section Section 2.1.1. Wright-Fisher model with selection and migration

#' @section Section 2.1.2. Wright-Fisher diffusion with selection and migration

############################################################

#' @section Section 2.2. Bayesian inference of natural selection and gene migration

#' @section Section 2.2.1. Hidden Markov model

#' @section Section 2.2.2. Particle marginal Metropolis-Hastings

################################################################################

#' @section Section 3. Results

#' @section Section 3.1. Robustness and performance

#' Figure 1. Simulation studies for the case that the continent allele frequencies of the sample are not available at 3 sampling time points
df_sel_cof <- NULL
df_sel_gen <- NULL
df_mig_rat <- NULL
df_mig_gen <- NULL

sel_cof <- c(3e-03, 6e-03, 9e-03)
mig_gen <- c(90, 360)
pop_siz <- c(5e+03, 5e+04)
for (s in 1:3) {
  for (N in 1:2) {
    for (t in 1:2) {
      load(paste0("./Output/Output v1.0/SIMU_u3_s", s, "N", N, "t", t, ".rda"))
      sel_cof_data <- cbind(sel_cof_data, rep(sel_cof[s], length.out = nrow(sel_cof_data)), rep(pop_siz[N], length.out = nrow(sel_cof_data)), rep(mig_gen[t], length.out = nrow(sel_cof_data)))
      df_sel_cof <- rbind(df_sel_cof, sel_cof_data)

      sel_gen_data <- cbind(sel_gen_data, rep(sel_cof[s], length.out = nrow(sel_gen_data)), rep(pop_siz[N], length.out = nrow(sel_gen_data)), rep(mig_gen[t], length.out = nrow(sel_gen_data)))
      df_sel_gen <- rbind(df_sel_gen, sel_gen_data)

      mig_rat_data <- cbind(mig_rat_data, rep(sel_cof[s], length.out = nrow(mig_rat_data)), rep(pop_siz[N], length.out = nrow(mig_rat_data)), rep(mig_gen[t], length.out = nrow(mig_rat_data)))
      df_mig_rat <- rbind(df_mig_rat, mig_rat_data)

      mig_gen_data <- cbind(mig_gen_data, rep(sel_cof[s], length.out = nrow(mig_gen_data)), rep(pop_siz[N], length.out = nrow(mig_gen_data)), rep(mig_gen[t], length.out = nrow(mig_gen_data)))
      df_mig_gen <- rbind(df_mig_gen, mig_gen_data)
    }
  }
}
df_sel_cof <- na.omit(df_sel_cof)
df_sel_cof <- df_sel_cof[, c(5, 6, 7, 1, 2, 3, 4)]
df_sel_cof <- as.data.frame(df_sel_cof)
rownames(df_sel_cof) <- NULL
colnames(df_sel_cof) <- c("sel_cof", "pop_siz", "mig_gen", "tru_val", "map_est", "low_hpd", "upp_hpd")
df_sel_cof$sel_cof <- factor(df_sel_cof$sel_cof, levels = c(3e-03, 6e-03, 9e-03), labels = c("3e-03", "6e-03", "9e-03"))
df_sel_cof$mig_gen <- factor(df_sel_cof$mig_gen, levels = c(90, 360), labels = c("90", "360"))
df_sel_cof$pop_siz <- factor(df_sel_cof$pop_siz, levels = c(5e+03, 5e+04), labels = c("5e+03", "5e+04"))

df_sel_gen <- na.omit(df_sel_gen)
df_sel_gen <- df_sel_gen[, c(5, 6, 7, 1, 2, 3, 4)]
df_sel_gen <- as.data.frame(df_sel_gen)
rownames(df_sel_gen) <- NULL
colnames(df_sel_gen) <- c("sel_cof", "pop_siz", "mig_gen", "tru_val", "map_est", "low_hpd", "upp_hpd")
df_sel_gen$sel_cof <- factor(df_sel_gen$sel_cof, levels = c(3e-03, 6e-03, 9e-03), labels = c("3e-03", "6e-03", "9e-03"))
df_sel_gen$mig_gen <- factor(df_sel_gen$mig_gen, levels = c(90, 360), labels = c("90", "360"))
df_sel_gen$pop_siz <- factor(df_sel_gen$pop_siz, levels = c(5e+03, 5e+04), labels = c("5e+03", "5e+04"))

df_mig_rat <- na.omit(df_mig_rat)
df_mig_rat <- df_mig_rat[, c(5, 6, 7, 1, 2, 3, 4)]
df_mig_rat <- as.data.frame(df_mig_rat)
rownames(df_mig_rat) <- NULL
colnames(df_mig_rat) <- c("sel_cof", "pop_siz", "mig_gen", "tru_val", "map_est", "low_hpd", "upp_hpd")
df_mig_rat$sel_cof <- factor(df_mig_rat$sel_cof, levels = c(3e-03, 6e-03, 9e-03), labels = c("3e-03", "6e-03", "9e-03"))
df_mig_rat$mig_gen <- factor(df_mig_rat$mig_gen, levels = c(90, 360), labels = c("90", "360"))
df_mig_rat$pop_siz <- factor(df_mig_rat$pop_siz, levels = c(5e+03, 5e+04), labels = c("5e+03", "5e+04"))

df_mig_gen <- na.omit(df_mig_gen)
df_mig_gen <- df_mig_gen[, c(5, 6, 7, 1, 2, 3, 4)]
df_mig_gen <- as.data.frame(df_mig_gen)
rownames(df_mig_gen) <- NULL
colnames(df_mig_gen) <- c("sel_cof", "pop_siz", "mig_gen", "tru_val", "map_est", "low_hpd", "upp_hpd")
df_mig_gen$sel_cof <- factor(df_mig_gen$sel_cof, levels = c(3e-03, 6e-03, 9e-03), labels = c("3e-03", "6e-03", "9e-03"))
df_mig_gen$mig_gen <- factor(df_mig_gen$mig_gen, levels = c(90, 360), labels = c("90", "360"))
df_mig_gen$pop_siz <- factor(df_mig_gen$pop_siz, levels = c(5e+03, 5e+04), labels = c("5e+03", "5e+04"))

df_sel_cof$bias <- df_sel_cof$map_est - df_sel_cof$tru_val
ggplot(df_sel_cof, aes(x = sel_cof, y = bias, fill = pop_siz)) +
  geom_boxplot(position = position_dodge(0.9)) +
  stat_summary(fun = "median", fun.args = list(mult = 1), position = position_dodge(0.9), geom = "pointrange", color = "black") +
  scale_fill_npg() +
  theme_bw() +
  facet_wrap(~mig_gen, scale = "free")

df_sel_gen$bias <- df_sel_gen$map_est - df_sel_gen$tru_val
ggplot(df_sel_gen, aes(x = sel_cof, y = bias, fill = pop_siz)) +
  geom_boxplot(position = position_dodge(0.9)) +
  stat_summary(fun = "median", fun.args = list(mult = 1), position = position_dodge(0.9), geom = "pointrange", color = "black") +
  scale_fill_npg() +
  theme_bw() +
  facet_wrap(~mig_gen, scale = "free")

df_mig_rat$bias <- df_mig_rat$map_est - df_mig_rat$tru_val
ggplot(df_mig_rat, aes(x = sel_cof, y = bias, fill = pop_siz)) +
  geom_boxplot(position = position_dodge(0.9)) +
  stat_summary(fun = "median", fun.args = list(mult = 1), position = position_dodge(0.9), geom = "pointrange", color = "black") +
  scale_fill_npg() +
  theme_bw() +
  facet_wrap(~mig_gen, scale = "free")

df_mig_gen$bias <- df_mig_gen$map_est - df_mig_gen$tru_val
ggplot(df_mig_gen, aes(x = sel_cof, y = bias, fill = pop_siz)) +
  geom_boxplot(position = position_dodge(0.9)) +
  stat_summary(fun = "median", fun.args = list(mult = 1), position = position_dodge(0.9), geom = "pointrange", color = "black") +
  scale_fill_npg() +
  theme_bw() +
  facet_wrap(~mig_gen, scale = "free")

############################################################

#' @section Section 3.2. Application to ancient chicken samples

#' @section Section 3.2.1. TSHR

#' @section Section 3.2.2. BCOD2

################################################################################

#' @section Section 4. Discussion

################################################################################

#' Supplemental Material

#' File S1. Additional results for the analysis of simulated data

#' Figure S1. Simulation studies for the case that the continent allele frequencies of the sample are not available at 7 sampling time points
df_sel_cof <- NULL
df_sel_gen <- NULL
df_mig_rat <- NULL
df_mig_gen <- NULL

sel_cof <- c(3e-03, 6e-03, 9e-03)
mig_gen <- c(90, 360)
pop_siz <- c(5e+03, 5e+04)
for (s in 1:3) {
  for (N in 1:2) {
    for (t in 1:2) {
      load(paste0("./Output/Output v1.0/SIMU_u7_s", s, "N", N, "t", t, ".rda"))
      sel_cof_data <- cbind(sel_cof_data, rep(sel_cof[s], length.out = nrow(sel_cof_data)), rep(pop_siz[N], length.out = nrow(sel_cof_data)), rep(mig_gen[t], length.out = nrow(sel_cof_data)))
      df_sel_cof <- rbind(df_sel_cof, sel_cof_data)

      sel_gen_data <- cbind(sel_gen_data, rep(sel_cof[s], length.out = nrow(sel_gen_data)), rep(pop_siz[N], length.out = nrow(sel_gen_data)), rep(mig_gen[t], length.out = nrow(sel_gen_data)))
      df_sel_gen <- rbind(df_sel_gen, sel_gen_data)

      mig_rat_data <- cbind(mig_rat_data, rep(sel_cof[s], length.out = nrow(mig_rat_data)), rep(pop_siz[N], length.out = nrow(mig_rat_data)), rep(mig_gen[t], length.out = nrow(mig_rat_data)))
      df_mig_rat <- rbind(df_mig_rat, mig_rat_data)

      mig_gen_data <- cbind(mig_gen_data, rep(sel_cof[s], length.out = nrow(mig_gen_data)), rep(pop_siz[N], length.out = nrow(mig_gen_data)), rep(mig_gen[t], length.out = nrow(mig_gen_data)))
      df_mig_gen <- rbind(df_mig_gen, mig_gen_data)
    }
  }
}

ind <- which(df_sel_cof[, 2] < 0.1 & df_sel_cof[, 2] > -0.1)
df_sel_cof <- df_sel_cof[ind, ]
df_sel_cof <- na.omit(df_sel_cof)
df_sel_cof <- df_sel_cof[, c(5, 6, 7, 1, 2, 3, 4)]
df_sel_cof <- as.data.frame(df_sel_cof)
rownames(df_sel_cof) <- NULL
colnames(df_sel_cof) <- c("sel_cof", "pop_siz", "mig_gen", "tru_val", "map_est", "low_hpd", "upp_hpd")
df_sel_cof$sel_cof <- factor(df_sel_cof$sel_cof, levels = c(3e-03, 6e-03, 9e-03), labels = c("3e-03", "6e-03", "9e-03"))
df_sel_cof$mig_gen <- factor(df_sel_cof$mig_gen, levels = c(90, 360), labels = c("90", "360"))
df_sel_cof$pop_siz <- factor(df_sel_cof$pop_siz, levels = c(5e+03, 5e+04), labels = c("5e+03", "5e+04"))

df_sel_gen <- df_sel_gen[ind, ]
df_sel_gen <- na.omit(df_sel_gen)
df_sel_gen <- df_sel_gen[, c(5, 6, 7, 1, 2, 3, 4)]
df_sel_gen <- as.data.frame(df_sel_gen)
rownames(df_sel_gen) <- NULL
colnames(df_sel_gen) <- c("sel_cof", "pop_siz", "mig_gen", "tru_val", "map_est", "low_hpd", "upp_hpd")
df_sel_gen$sel_cof <- factor(df_sel_gen$sel_cof, levels = c(3e-03, 6e-03, 9e-03), labels = c("3e-03", "6e-03", "9e-03"))
df_sel_gen$mig_gen <- factor(df_sel_gen$mig_gen, levels = c(90, 360), labels = c("90", "360"))
df_sel_gen$pop_siz <- factor(df_sel_gen$pop_siz, levels = c(5e+03, 5e+04), labels = c("5e+03", "5e+04"))

df_mig_rat <- df_mig_rat[ind, ]
df_mig_rat <- na.omit(df_mig_rat)
df_mig_rat <- df_mig_rat[, c(5, 6, 7, 1, 2, 3, 4)]
df_mig_rat <- as.data.frame(df_mig_rat)
rownames(df_mig_rat) <- NULL
colnames(df_mig_rat) <- c("sel_cof", "pop_siz", "mig_gen", "tru_val", "map_est", "low_hpd", "upp_hpd")
df_mig_rat$sel_cof <- factor(df_mig_rat$sel_cof, levels = c(3e-03, 6e-03, 9e-03), labels = c("3e-03", "6e-03", "9e-03"))
df_mig_rat$mig_gen <- factor(df_mig_rat$mig_gen, levels = c(90, 360), labels = c("90", "360"))
df_mig_rat$pop_siz <- factor(df_mig_rat$pop_siz, levels = c(5e+03, 5e+04), labels = c("5e+03", "5e+04"))

df_mig_gen <- df_mig_gen[ind, ]
df_mig_gen <- na.omit(df_mig_gen)
df_mig_gen <- df_mig_gen[, c(5, 6, 7, 1, 2, 3, 4)]
df_mig_gen <- as.data.frame(df_mig_gen)
rownames(df_mig_gen) <- NULL
colnames(df_mig_gen) <- c("sel_cof", "pop_siz", "mig_gen", "tru_val", "map_est", "low_hpd", "upp_hpd")
df_mig_gen$sel_cof <- factor(df_mig_gen$sel_cof, levels = c(3e-03, 6e-03, 9e-03), labels = c("3e-03", "6e-03", "9e-03"))
df_mig_gen$mig_gen <- factor(df_mig_gen$mig_gen, levels = c(90, 360), labels = c("90", "360"))
df_mig_gen$pop_siz <- factor(df_mig_gen$pop_siz, levels = c(5e+03, 5e+04), labels = c("5e+03", "5e+04"))


df_sel_cof$bias <- df_sel_cof$map_est - df_sel_cof$tru_val
ggplot(df_sel_cof, aes(x = sel_cof, y = bias, fill = pop_siz)) +
  geom_boxplot(position = position_dodge(0.9)) +
  stat_summary(fun = "median", fun.args = list(mult = 1), position = position_dodge(0.9), geom = "pointrange", color = "black") +
  scale_fill_npg() +
  theme_bw() +
  facet_wrap(~mig_gen, scale = "free")

df_sel_gen$bias <- df_sel_gen$map_est - df_sel_gen$tru_val
ggplot(df_sel_gen, aes(x = sel_cof, y = bias, fill = pop_siz)) +
  geom_boxplot(position = position_dodge(0.9)) +
  stat_summary(fun = "median", fun.args = list(mult = 1), position = position_dodge(0.9), geom = "pointrange", color = "black") +
  scale_fill_npg() +
  theme_bw() +
  facet_wrap(~mig_gen, scale = "free")

df_mig_rat$bias <- df_mig_rat$map_est - df_mig_rat$tru_val
ggplot(df_mig_rat, aes(x = sel_cof, y = bias, fill = pop_siz)) +
  geom_boxplot(position = position_dodge(0.9)) +
  stat_summary(fun = "median", fun.args = list(mult = 1), position = position_dodge(0.9), geom = "pointrange", color = "black") +
  scale_fill_npg() +
  theme_bw() +
  facet_wrap(~mig_gen, scale = "free")

df_mig_gen$bias <- df_mig_gen$map_est - df_mig_gen$tru_val
ggplot(df_mig_gen, aes(x = sel_cof, y = bias, fill = pop_siz)) +
  geom_boxplot(position = position_dodge(0.9)) +
  stat_summary(fun = "median", fun.args = list(mult = 1), position = position_dodge(0.9), geom = "pointrange", color = "black") +
  scale_fill_npg() +
  theme_bw() +
  facet_wrap(~mig_gen, scale = "free")

############################################################

#' File S2. Additional results for the analysis of real data

####################################################################################################
