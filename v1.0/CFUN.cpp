// Inferring natural selection and gene migration in the evolution of chickens from ancient DNA data
// Zhangyi He, Wenyang Lyu, Xiaoyang Dai, Mark Beaumont and Feng Yu

// version 2.0

// C functions

#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
#include <math.h>

using namespace Rcpp;
using namespace std;
// using namespace arma;

/********** WFM **********/
// Calculate the fitness matrix for the one-locus Wright-Fisher model with selection and migration
// [[Rcpp::export]]
arma::dmat calculateFitnessMat_2L_arma(const double& sel_cof, const double& dom_par) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // calculate the fitness matrix
  arma::dmat fts_mat = arma::ones<arma::dmat>(4, 4);
  // fts_mat(0, 0) = 1;
  fts_mat(1, 0) = (1 - dom_par * sel_cof);
  //fts_mat(2, 0) = 1;
  fts_mat(3, 0) = (1 - dom_par * sel_cof);
  fts_mat(0, 1) = fts_mat(1, 0);
  fts_mat(1, 1) = (1 - sel_cof);
  fts_mat(2, 1) = (1 - dom_par * sel_cof);
  fts_mat(3, 1) = (1 - sel_cof);
  fts_mat(0, 2) = fts_mat(2, 0);
  fts_mat(1, 2) = fts_mat(2, 1);
  //fts_mat(2, 2) = 1;
  fts_mat(3, 2) = (1 - dom_par * sel_cof);
  fts_mat(0, 3) = fts_mat(3, 0);
  fts_mat(1, 3) = fts_mat(3, 1);
  fts_mat(2, 3) = fts_mat(3, 2);
  fts_mat(3, 3) = (1 - sel_cof);
  
  // return the fitness matrix for the Wright-Fisher model
  return fts_mat;
}

// Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection and migration
// [[Rcpp::export]]
arma::dmat simulateWFM_2L_arma(const arma::dmat& fts_mat, const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const double& ext_frq, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the haplotype frequency trajectories
  arma::dmat frq_pth(4, arma::uword(lst_gen - int_gen) + 1);
  
  // initialise the haplotype frequencies in generation 0
  frq_pth.col(0) = int_frq;
  
  // declare eta
  arma::dcolvec eta = {0.0, 0.0, ext_frq, 1.0 - ext_frq};
  
  // simulate the haplotype frequency trajectories
  arma::dcolvec hap_frq = int_frq;
  for(arma::uword k = 1; k < arma::uword(lst_gen - int_gen) + 1; k++) {
    // calculate the sampling probabilities
    arma::dcolvec prob = hap_frq;
    prob = hap_frq % (fts_mat * hap_frq) / arma::as_scalar(hap_frq.t() * fts_mat * hap_frq);
    prob = (1 - mig_rat) * prob + eta * mig_rat;
    
    // proceed the Wright-Fisher sampling
    IntegerVector hap_cnt(4);
    R::rmultinom(2 * pop_siz, prob.begin(), 4, hap_cnt.begin());
    hap_frq = as<arma::dcolvec>(hap_cnt) / 2 / pop_siz;
    frq_pth.col(k) = hap_frq;
  }
  
  // return the haplotype frequency trajectories under the Wright-Fisher model
  return frq_pth;
}
/*************************/


/********** WFD **********/
// Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher diffusion with selection migration using the Euler-Maruyama method
// [[Rcpp::export]]
arma::dmat simulateWFD_2L_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const double& ext_frq, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // rescale the selection coefficients and the recombination rate
  double scl_sel_cof = 2 * pop_siz * sel_cof;
  double scl_mig_rat = 4 * pop_siz * mig_rat;
  
  // calculate delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  // generate delta W
  arma::dmat dW = pow(dt, 0.5) * arma::randn<arma::dmat>(6, arma::uword(lst_gen - int_gen) * ptn_num);
  
  // declare the haplotype frequency trajectories
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);
  
  // initialise the haplotype frequencies in generation 0
  frq_pth.col(0) = int_frq;
  
  // simulate the haplotype frequency trajectories
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient vector
    arma::dcolvec mu = arma::zeros<arma::dcolvec>(4);
    mu(0) =  scl_sel_cof * frq_pth(0, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par)) - 0.5 * scl_mig_rat * (frq_pth(0, t - 1));
    mu(1) = -scl_sel_cof * frq_pth(1, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par)) - 0.5 * scl_mig_rat * (frq_pth(1, t - 1));
    mu(2) = scl_sel_cof * frq_pth(2, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par))  - 0.5 * scl_mig_rat * (frq_pth(2, t - 1) - ext_frq);
    mu(3) = -scl_sel_cof * frq_pth(3, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par)) - 0.5 * scl_mig_rat * (frq_pth(3, t - 1) + ext_frq - 1);
    
    // calculate the diffusion coefficient matrix
    arma::dmat sigma = arma::zeros<arma::dmat>(4, 6);
    sigma(0, 0) = pow(frq_pth(0, t - 1) * frq_pth(1, t - 1), 0.5);
    sigma(0, 1) = pow(frq_pth(0, t - 1) * frq_pth(2, t - 1), 0.5);
    sigma(0, 2) = pow(frq_pth(0, t - 1) * frq_pth(3, t - 1), 0.5);
    // sigma(0, 3) = 0;
    // sigma(0, 4) = 0;
    // sigma(0, 5) = 0;
    sigma(1, 0) = -pow(frq_pth(1, t - 1) * frq_pth(0, t - 1), 0.5);
    // sigma(1, 1) = 0;
    // sigma(1, 2) = 0;
    sigma(1, 3) = pow(frq_pth(1, t - 1) * frq_pth(2, t - 1), 0.5);
    sigma(1, 4) = pow(frq_pth(1, t - 1) * frq_pth(3, t - 1), 0.5);
    // sigma(1, 5) = 0;
    // sigma(2, 0) = 0;
    sigma(2, 1) = -pow(frq_pth(2, t - 1) * frq_pth(0, t - 1), 0.5);
    // sigma(2, 2) = 0;
    sigma(2, 3) = -pow(frq_pth(2, t - 1) * frq_pth(1, t - 1), 0.5);
    sigma(2, 4) = 0;
    sigma(2, 5) = pow(frq_pth(2, t - 1) * frq_pth(3, t - 1), 0.5);
    // sigma(3, 0) = 0;
    // sigma(3, 1) = 0;
    sigma(3, 2) = -pow(frq_pth(3, t - 1) * frq_pth(0, t - 1), 0.5);
    // sigma(3, 3) = 0;
    sigma(3, 4) = -pow(frq_pth(3, t - 1) * frq_pth(1, t - 1), 0.5);
    sigma(3, 5) = -pow(frq_pth(3, t - 1) * frq_pth(2, t - 1), 0.5);
    
    // proceed the Euler-Maruyama scheme
    frq_pth.col(t) = frq_pth.col(t - 1) + mu * dt + sigma * dW.col(t - 1);
    
    // remove the noise from the numerical techniques
    for(arma::uword i = 0; i < 4; i++) {
      if(frq_pth(i, t) < 0) {
        frq_pth(i, t) = 0;
      }
      if(frq_pth(i, t) > 1) {
        frq_pth(i, t) = 1;
      }
    }
    frq_pth.col(t) = frq_pth.col(t) / sum(frq_pth.col(t));
  }
  
  // return the haplotype frequency trajectories under the Wright-Fisher diffusion
  return frq_pth;
}
/*************************/


/********** BPF **********/
// Calculate the possible haplotype counts in the sample
// [[Rcpp::export]]
arma::imat calculateHaploCnt_arma(const int& smp_siz, const arma::icolvec& smp_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  if (smp_cnt.n_elem == 2) {
    // smp_cnt = (c_k, d_k)
    arma::imat smp_hap_cnt = arma::zeros<arma::imat>(4, 1);
    
    for (int i = 0; i <= min(smp_cnt(0), smp_cnt(1)); i++) {
      int j = smp_cnt(0) - i;
      int k = smp_cnt(1) - i;
      if (i + j + k <= smp_siz) {
        smp_hap_cnt(0, 0) = j;
        smp_hap_cnt(1, 0) = smp_siz - i - j - k;
        smp_hap_cnt(2, 0) = i;
        smp_hap_cnt(3, 0) = k;
        smp_hap_cnt.insert_cols(0, 1);
      }
    }
    smp_hap_cnt.shed_cols(0, 0);
    
    return smp_hap_cnt;
  } else {
    arma::imat smp_hap_cnt = arma::zeros<arma::imat>(4, 1);
    smp_hap_cnt.col(0) = smp_cnt;
    
    return smp_hap_cnt;
  }
}

// Calculate the multinomial probabilities
// [[Rcpp::export]]
double calculateMultinomProb_arma(const arma::icolvec& smp_cnt, const int& smp_siz, const arma::dcolvec& pop_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  if (arma::any(pop_frq == 0)) {
    if (arma::any(smp_cnt.elem(arma::find(pop_frq == 0)) != 0)) {
      return 0;
    }
    
    arma::icolvec smp_cnt_nonzero = smp_cnt.elem(arma::find(pop_frq != 0));
    arma::dcolvec pop_frq_nonzero = pop_frq.elem(arma::find(pop_frq != 0));
    
    return exp(lgamma(smp_siz + 1) + sum(arma::conv_to<arma::dcolvec>::from(smp_cnt_nonzero) % log(pop_frq_nonzero) - lgamma(arma::conv_to<arma::dcolvec>::from(smp_cnt_nonzero) + 1)));
  } else {
    return exp(lgamma(smp_siz + 1) + sum(arma::conv_to<arma::dcolvec>::from(smp_cnt) % log(pop_frq) - lgamma(arma::conv_to<arma::dcolvec>::from(smp_cnt) + 1)));
  }
}

// Initialise the particles in the particle filter (uniform generation from the flat Dirichlet distribution)
// [[Rcpp::export]]
arma::dmat initialiseParticle(const arma::uword& pcl_num){
  // ensure RNG gets set/reset
  RNGScope scope;
  
  NumericMatrix part(pcl_num, 4);
  for(int j = 0; j < 4; j++){
    part(_, j) = rgamma(pcl_num, 1.0, 1.0);
  }
  for(int i = 0; i < pcl_num; i++){
    part(i, _) = part(i, _) / sum(part(i, _));
  }
  
  return as<arma::dmat>(transpose(part));
}

// Initialise the particles in the particle filter (uniform generation from the flat Dirichlet distribution)
// [[Rcpp::export]]
arma::dmat initialiseParticleNoMigration(const arma::uword& pcl_num){
  // ensure RNG gets set/reset
  RNGScope scope;
  
  NumericMatrix part(4, pcl_num);
  part(0, _) = runif(pcl_num, 0.0, 1.0);
  part(1, _) = 1 - part(0, _);
  
  return as<arma::dmat>(part);
}
/********** BPF **********/
// Run the bootstrap particle filter for two consecutive time points (from one time point to the sampling time point)
// [[Rcpp::export]]
arma::dmat runBPF_to_Smp_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const double& ext_frq, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num, const bool& non_sel, const bool& non_mig) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::dmat path = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);
  
  if (non_sel == true && non_mig == true) {
      path = simulateWFD_2L_arma(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);

  } else if (non_sel == true && non_mig == false) {
      path = simulateWFD_2L_arma(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);
    
  } else if (non_sel == false && non_mig == true) {
      path = simulateWFD_2L_arma(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);

  } else {
      path = simulateWFD_2L_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);
  }
  return (path);
}

// Run the bootstrap particle filter for two consecutive time points (from one time point to the non-sampling time point)
// [[Rcpp::export]]
arma::dmat runBPF_to_NonSmp_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const double& ext_frq, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num, const bool& non_sel, const bool& non_mig) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::dmat path = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);
  
  if (non_sel == true && non_mig == true) {
    path = simulateWFD_2L_arma(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);

  } else if (non_sel == true && non_mig == false) {
      path = simulateWFD_2L_arma(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);
    
  } else if (non_sel == false && non_mig == true) {
      path = simulateWFD_2L_arma(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);
  } else {
      path = simulateWFD_2L_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);
  }
  return (path);
}


// Run the bootstrap particle filter
// [[Rcpp::export]]
List runBPF_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const int& sel_gen, const int& mig_gen, const double& ext_frq, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::irowvec all_gen(smp_gen.n_elem + 2);
  all_gen.subvec(1, smp_gen.n_elem) = smp_gen;
  all_gen.head(1) = sel_gen;
  all_gen.tail(1) = mig_gen;
  all_gen = arma::unique(all_gen);
  all_gen = arma::sort(all_gen);
  
  double lik = 1;
  
  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, all_gen.n_elem);
  arma::dcube part_pre = arma::zeros<arma::dcube>(4, pcl_num, all_gen.n_elem);
  arma::dcube part_pst = arma::zeros<arma::dcube>(4, pcl_num, all_gen.n_elem);
  
  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_tmp = arma::zeros<arma::dmat>(4, pcl_num);
  
  bool non_sel = true;
  bool non_mig = true;
  
  arma::uword smp_ind = 0;
  
  // initialise the particles
  cout << "generation: " << all_gen(0) << endl;
  if (mig_gen != 0 ){
    part_tmp = initialiseParticleNoMigration(pcl_num);
  } else{
    part_tmp = initialiseParticle(pcl_num);
  }
  arma::imat smp_hap_cnt = calculateHaploCnt_arma(smp_siz(0), smp_cnt.col(0));
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
      wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(0), part_tmp.col(i));
    }
  }
  
  if (arma::sum(wght_tmp) > 0) {
    arma::dcolvec prob = arma::normalise(wght_tmp, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
    
    lik = lik * arma::mean(wght_tmp);
    wght.col(0) = wght_tmp;
    part_pre.slice(0) = part_tmp;
    part_pst.slice(0) = part_tmp.cols(indx);
  } else {
    lik = 0;
    wght.shed_cols(0, all_gen.n_elem - 1);
    part_pre.shed_slices(0, all_gen.n_elem - 1);
    part_pst.shed_slices(0, all_gen.n_elem - 1);
    
    return List::create(Named("lik", lik), 
                        Named("wght", wght), 
                        Named("part_pre_resmp", part_pre), 
                        Named("part_pst_resmp", part_pst));
  }
  
  // run the bootstrap particle filter
  for (arma::uword k = 1; k < all_gen.n_elem; k++) {
    if (sel_gen > all_gen(k - 1)) {
      non_sel = true;
    } else {
      non_sel = false;
    }
    
    if (mig_gen > all_gen(k - 1)) {
      non_mig = true;
    } else {
      non_mig = false;
    }
    
    cout << "generation: " << all_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    part_tmp = part_pst.slice(k - 1);
    if (arma::any(smp_gen == all_gen(k))) {
      smp_ind = smp_ind + 1;
      
      arma::imat smp_hap_cnt = calculateHaploCnt_arma(smp_siz(smp_ind), smp_cnt.col(smp_ind));
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::dmat path = runBPF_to_Smp_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, part_tmp.col(i), all_gen(k - 1), all_gen(k), ptn_num, non_sel, non_mig);
        part_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
        for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
          wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(smp_ind), part_tmp.col(i));
        }
      }
      if (arma::sum(wght_tmp) > 0) {
        arma::dcolvec prob = arma::normalise(wght_tmp, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
        
        lik = lik * arma::mean(wght_tmp);
        wght.col(k) = wght_tmp;
        part_pre.slice(k) = part_tmp;
        part_pst.slice(k) = part_tmp.cols(indx);
      } else {
        lik = 0;
        wght.shed_cols(k, all_gen.n_elem - 1);
        part_pre.shed_slices(k, all_gen.n_elem - 1);
        part_pst.shed_slices(k, all_gen.n_elem - 1);
        break;
      }
      
    } else {
      smp_ind = smp_ind + 0;
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::dmat path = runBPF_to_NonSmp_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, part_tmp.col(i), all_gen(k - 1), all_gen(k), ptn_num, non_sel, non_mig);
        part_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
        wght_tmp(i) = 1;
      }
        
        lik = lik * arma::mean(wght_tmp);
        wght.col(k) = wght_tmp;
        part_pre.slice(k) = part_tmp;
        part_pst.slice(k) = part_tmp;
    }
    
  }
  
  return List::create(Named("lik", lik),
                      Named("wght", wght), 
                      Named("part_pre_resmp", part_pre), 
                      Named("part_pst_resmp", part_pst));
}
/*************************/

/*** PMMH within Gibbs ***/
// Run the bootstrap particle filter for two consecutive sampling time points
// [[Rcpp::export]]
arma::dmat runBPF_Smp_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const int& sel_gen, const int& mig_gen, const double& ext_frq, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  bool chg_sel;
  if (sel_gen <= int_gen || sel_gen > lst_gen) {
    chg_sel = false;
  } else {
    chg_sel = true;
  }
  bool chg_mig;
  if (mig_gen <= int_gen || mig_gen > lst_gen) {
    chg_mig = false;
  } else {
    chg_mig = true;
  }
  
  bool non_sel;
  if (sel_gen <= int_gen) {
    non_sel = false;
  } else {
    non_sel = true;
  }
  bool non_mig;
  if (mig_gen <= int_gen) {
    non_mig = false;
  } else {
    non_mig = true;
  }
  
  arma::dcolvec part = arma::ones<arma::dcolvec>(4);
  
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);
  
  if (chg_sel == false && chg_mig == false) {
    if (non_sel == true && non_mig == true) {
      arma::dmat path = simulateWFD_2L_arma(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);
      part = (path.tail_cols(1));
      
      frq_pth = path;
    } else if (non_sel == true && non_mig == false) {
      arma::dmat path = simulateWFD_2L_arma(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);
      part = (path.tail_cols(1));
      
      frq_pth = path;
    } else if (non_sel == false && non_mig == true) {
      arma::dmat path = simulateWFD_2L_arma(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);
      part = (path.tail_cols(1));
      
      frq_pth = path;
    } else {
      arma::dmat path = simulateWFD_2L_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, lst_gen, ptn_num);
      part = (path.tail_cols(1));
      
      frq_pth = path;
    }
  } else if (chg_sel == true && chg_mig == false) {
    if (non_mig == true) {
      arma::dmat path_int = simulateWFD_2L_arma(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen, ptn_num);
      part = (path_int.tail_cols(1));
      
      arma::dmat path_lst = simulateWFD_2L_arma(sel_cof, dom_par, 0, pop_siz, ext_frq, part, sel_gen, lst_gen, ptn_num);
      part = (path_lst.tail_cols(1));
      
      frq_pth.cols(0, (sel_gen - int_gen) * ptn_num) = path_int;
      frq_pth.cols((sel_gen - int_gen) * ptn_num, (lst_gen - int_gen) * ptn_num) = path_lst;
    } else {
      arma::dmat path_int = simulateWFD_2L_arma(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq, int_gen, sel_gen, ptn_num);
      part = (path_int.tail_cols(1));
      
      arma::dmat path_lst = simulateWFD_2L_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, part, sel_gen, lst_gen, ptn_num);
      part = (path_lst.tail_cols(1));
      
      frq_pth.cols(0, (sel_gen - int_gen) * ptn_num) = path_int;
      frq_pth.cols((sel_gen - int_gen) * ptn_num, (lst_gen - int_gen) * ptn_num) = path_lst;
    }
  } else if (chg_sel == false && chg_mig == true) {
    if (non_sel == true) {
      arma::dmat path_int = simulateWFD_2L_arma(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen, ptn_num);
      part = (path_int.tail_cols(1));
      
      arma::dmat path_lst = simulateWFD_2L_arma(0, dom_par, mig_rat, pop_siz, ext_frq, part, mig_gen, lst_gen, ptn_num);
      part = (path_lst.tail_cols(1));
 
      frq_pth.cols(0, (mig_gen - int_gen) * ptn_num) = path_int;
      frq_pth.cols((mig_gen - int_gen) * ptn_num, (lst_gen - int_gen) * ptn_num) = path_lst;
    } else {
      arma::dmat path_int = simulateWFD_2L_arma(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen, ptn_num);
      part = (path_int.tail_cols(1));
      
      arma::dmat path_lst = simulateWFD_2L_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, part, mig_gen, lst_gen, ptn_num);
      part = (path_lst.tail_cols(1));
      
      frq_pth.cols(0, (mig_gen - int_gen) * ptn_num) = path_int;
      frq_pth.cols((mig_gen - int_gen) * ptn_num, (lst_gen - int_gen) * ptn_num) = path_lst;
    }
  } else {
    if (sel_gen > mig_gen) {
      arma::dmat path_int = simulateWFD_2L_arma(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, mig_gen, ptn_num);
      part = (path_int.tail_cols(1));
      
      arma::dmat path_mid = simulateWFD_2L_arma(0, dom_par, mig_rat, pop_siz, ext_frq, part, mig_gen, sel_gen, ptn_num);
      part = (path_mid.tail_cols(1));
      
      arma::dmat path_lst = simulateWFD_2L_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, part, sel_gen, lst_gen, ptn_num);
      part = (path_lst.tail_cols(1));
      
      frq_pth.cols(0, (mig_gen - int_gen) * ptn_num) = path_int;
      frq_pth.cols((mig_gen - int_gen) * ptn_num, (sel_gen - int_gen) * ptn_num) = path_mid;
      frq_pth.cols((sel_gen - int_gen) * ptn_num, (lst_gen - int_gen) * ptn_num) = path_lst;
    } else if (sel_gen < mig_gen) {
      arma::dmat path_int = simulateWFD_2L_arma(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen, ptn_num);
      part = (path_int.tail_cols(1));
      
      arma::dmat path_mid = simulateWFD_2L_arma(sel_cof, dom_par, 0, pop_siz, ext_frq, part, sel_gen, mig_gen, ptn_num);
      part = (path_mid.tail_cols(1));
      
      arma::dmat path_lst = simulateWFD_2L_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, part, mig_gen, lst_gen, ptn_num);
      part = (path_lst.tail_cols(1));
      
      frq_pth.cols(0, (sel_gen - int_gen) * ptn_num) = path_int;
      frq_pth.cols((sel_gen - int_gen) * ptn_num, (mig_gen - int_gen) * ptn_num) = path_mid;
      frq_pth.cols((mig_gen - int_gen) * ptn_num, (lst_gen - int_gen) * ptn_num) = path_lst;
    } else {
      arma::dmat path_int = simulateWFD_2L_arma(0, dom_par, 0, pop_siz, ext_frq, int_frq, int_gen, sel_gen, ptn_num);
      part = (path_int.tail_cols(1));
      
      arma::dmat path_lst = simulateWFD_2L_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, part, mig_gen, lst_gen, ptn_num);
      part = (path_lst.tail_cols(1));

      
      frq_pth.cols(0, (sel_gen - int_gen) * ptn_num) = path_int;
      frq_pth.cols((mig_gen - int_gen) * ptn_num, (lst_gen - int_gen) * ptn_num) = path_lst;
    }
  }
  return (frq_pth);
}


/********** PMMH **********/
// Calculate the log-likelihood using the bootstrap particle filter
// [[Rcpp::export]]
void calculateLogLikelihood_arma(double& log_lik, arma::dcube& frq_pth, const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const int& sel_gen, const int& mig_gen, const double& ext_frq, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::field<arma::imat>& ptl_hap_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  log_lik = 0;
  
  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_pre = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat part_pst = arma::zeros<arma::dmat>(4, pcl_num);
  frq_pth = arma::zeros<arma::dcube>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1, pcl_num);
  
  // initialise the particles
  if (mig_gen != 0 ){
    part_pre = initialiseParticleNoMigration(pcl_num);
  } else{
    part_pre = initialiseParticle(pcl_num);
  }
  arma::imat smp_hap_cnt = ptl_hap_cnt(0);
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
      wght(i) = wght(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(0), part_pre.col(i));
    }
  }
  
  if (arma::mean(wght) > 0) {
    log_lik = log_lik + log(arma::mean(wght));
    arma::dcolvec prob = arma::normalise(wght, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
    part_pst = part_pre.cols(indx);
  } else {
    log_lik = -(arma::datum::inf);
    return;
  }
  
  // run the bootstrap particle filter
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat smp_hap_cnt = ptl_hap_cnt(k);
    arma::dcube frq_pth_smp_k = arma::zeros<arma::dcube>(4, arma::uword(smp_gen(k) - smp_gen(k - 1)) * ptn_num + 1, pcl_num);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = runBPF_Smp_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, part_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      part_pre.col(i) = arma::vectorise(path.tail_cols(1), 0);
      frq_pth_smp_k.slice(i) = path;
      for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(k), part_pre.col(i));
      }
    }
    
    if (arma::mean(wght) > 0) {
      log_lik = log_lik + log(arma::mean(wght));
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
      frq_pth.subcube(0, (smp_gen(k - 1) - smp_gen(0)) * ptn_num, 0, 3, (smp_gen(k) - smp_gen(0)) * ptn_num, pcl_num - 1) = frq_pth_smp_k;
      frq_pth = frq_pth.slices(indx);
    } else {
      log_lik = -(arma::datum::inf);
      return;
    }
  }
  
}

// Calculate the optimal particle number in the particle marginal Metropolis-Hastings
// [[Rcpp::export]]
List calculateOptimalParticleNum_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const int& sel_gen, const int& mig_gen, const double& ext_frq, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::field<arma::imat> ptl_hap_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_hap_cnt(k) = calculateHaploCnt_arma(smp_siz(k), smp_cnt.col(k));
  }
  
  arma::drowvec log_lik(300);
  for (arma::uword i = 0; i < 300; i++) {
    double log_lik_i = 0;
    arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1, pcl_num);
    calculateLogLikelihood_arma(log_lik_i, frq_pth, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, pcl_num);
    log_lik(i) = log_lik_i;
  }
  
  arma::drowvec log_lik_sdv(1);
  log_lik_sdv(0) = arma::stddev(log_lik);
  log_lik_sdv.print();
  arma::urowvec opt_pcl_num(1);
  opt_pcl_num(0) = pcl_num;
  
  if (log_lik_sdv(0) > 1.7) {
    while (log_lik_sdv(0) > 1.0) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) + gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        double log_lik_i = 0;
        arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1, pcl_num);
        calculateLogLikelihood_arma(log_lik_i, frq_pth, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, opt_pcl_num(0));
        log_lik(i) = log_lik_i;
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
    opt_pcl_num = arma::reverse(opt_pcl_num);
    log_lik_sdv = arma::reverse(log_lik_sdv);
  } else if (log_lik_sdv(0) < 1.0) {
    while (log_lik_sdv(0) < 1.7 && opt_pcl_num(0) > gap_num) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) - gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        double log_lik_i = 0;
        arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1, pcl_num);
        calculateLogLikelihood_arma(log_lik_i, frq_pth, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, opt_pcl_num(0));
        log_lik(i) = log_lik_i;
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
  } else {
    while (log_lik_sdv(0) > 1.0) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) + gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        double log_lik_i = 0;
        arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1, pcl_num);
        calculateLogLikelihood_arma(log_lik_i, frq_pth, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, opt_pcl_num(0));
        log_lik(i) = log_lik_i;
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
    opt_pcl_num = arma::reverse(opt_pcl_num);
    log_lik_sdv = arma::reverse(log_lik_sdv);
    
    while (log_lik_sdv(0) < 1.7 && opt_pcl_num(0) > gap_num) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) - gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        log_lik(i) = 0.5;
        double log_lik_i = 0;
        arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1, pcl_num);
        calculateLogLikelihood_arma(log_lik_i, frq_pth, sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, opt_pcl_num(0));
        log_lik(i) = log_lik_i;
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
  }
  
  return List::create(Named("opt_pcl_num", opt_pcl_num), 
                      Named("log_lik_sdv", log_lik_sdv));
}

// Run the particle marginal Metropolis-Hastings
//[[Rcpp::export]]
List runPMMH_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const int& sel_gen, const int& mig_gen, const double& ext_frq, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::field<arma::imat> ptl_hap_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_hap_cnt(k) = calculateHaploCnt_arma(smp_siz(k), smp_cnt.col(k));
  }
  
  arma::drowvec sel_cof_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::irowvec sel_gen_chn = arma::zeros<arma::irowvec>(itn_num);
  arma::drowvec mig_rat_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::irowvec mig_gen_chn = arma::zeros<arma::irowvec>(itn_num);
  arma::dcube frq_pth_chn = arma::zeros<arma::dcube>(4, arma::uword(max(smp_gen) - min(smp_gen)) * ptn_num + 1, itn_num);
  
  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);
  double sel_cof_sd = 5e-03;
  //double sel_gen_sd = 2e+01;
  double mig_rat_sd = 1e-03;
  //double mig_gen_sd = 2e+01;
  
  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn(0) = sel_cof;
  sel_gen_chn(0) = sel_gen;
  mig_rat_chn(0) = mig_rat;
  mig_gen_chn(0) = mig_gen;
  double apt_rate = 0;
  
  arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1, pcl_num);
  double log_lik = 0;
  calculateLogLikelihood_arma(log_lik, frq_pth, sel_cof_chn(0), dom_par, mig_rat_chn(0), pop_siz, sel_gen_chn(0), mig_gen_chn(0), ext_frq, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, pcl_num);
  log_lik_chn(0) = log_lik;
  frq_pth_chn.slice(0) = frq_pth.slice(0);
  double apt_rto = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    
    // draw the candidates of the selection coefficients from the random walk proposal
    apt_rto = 1;
    sel_cof_chn(i) = sel_cof_chn(i - 1) + sel_cof_sd * arma::randn();
    //sel_cof_chn(i) = sel_cof;
    mig_rat_chn(i) = mig_rat_chn(i - 1) + mig_rat_sd * arma::randn();
    //mig_rat_chn(i) = mig_rat;
    if (sel_cof_chn(i) > 1 || mig_rat_chn(i) > 1 ) {
      apt_rto = 0;
    }

    //if (sel_cof_chn(i) < 0 || mig_rat_chn(i) < 0 ) {
    //  apt_rto = 0;
    //}    
    
    // draw the candidate of the selection time from a random walk proposal
    //sel_gen_chn(i) = sel_gen_chn(i - 1) + int(round(sel_gen_sd * arma::randn()));
    sel_gen_chn(i) = sel_gen;
    if (sel_gen_chn(i) < smp_gen.min() || sel_gen_chn(i) > smp_gen.max()) {
      apt_rto = 0;
    }
    //mig_gen_chn(i) = mig_gen_chn(i - 1) + int(round(mig_gen_sd * arma::randn()));
    mig_gen_chn(i) = mig_gen;
    if (mig_gen_chn(i) < smp_gen.min() || mig_gen_chn(i) > smp_gen.max()) {
      apt_rto = 0;
    }
    
    if (apt_rto == 0)  {
      sel_cof_chn(i) = sel_cof_chn(i - 1);
      sel_gen_chn(i) = sel_gen_chn(i - 1);
      mig_rat_chn(i) = mig_rat_chn(i - 1);
      mig_gen_chn(i) = mig_gen_chn(i - 1);
      frq_pth_chn(i) = frq_pth_chn(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      // calculate the proposal
      //double log_psl_old_new = log(arma::normpdf(sel_cof_A_chn(i - 1), sel_cof_A_chn(i), sel_cof_sd)) + log(arma::normpdf(sel_cof_B_chn(i - 1), sel_cof_B_chn(i), sel_cof_sd));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_A_chn(i), sel_cof_A_chn(i - 1), sel_cof_sd)) + log(arma::normpdf(sel_cof_B_chn(i), sel_cof_B_chn(i - 1), sel_cof_sd));
      
      // calculate the likelihood
      arma::dcube frq_pth = arma::zeros<arma::dcube>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1, pcl_num);
      calculateLogLikelihood_arma(log_lik, frq_pth, sel_cof_chn(i), dom_par, mig_rat_chn(i), pop_siz, sel_gen_chn(i), mig_gen_chn(i), ext_frq, smp_gen, smp_siz, ptl_hap_cnt, ptn_num, pcl_num);
      log_lik_chn(i) = log_lik;
      frq_pth_chn.slice(i) = frq_pth.slice(0);
      // calculate the acceptance ratio
      apt_rto = exp(log_lik_chn(i) - log_lik_chn(i - 1));
      //apt_rto = exp((log_pri_chn(i) + log_lik_chn(i) + log_psl_old_new) - (log_pri_chn(i - 1) + log_lik_chn(i - 1) + log_psl_new_old));
      
      apt_rate = apt_rate + 1;
      if (arma::randu() > apt_rto) {
        sel_cof_chn(i) = sel_cof_chn(i - 1);
        sel_gen_chn(i) = sel_gen_chn(i - 1);
        mig_rat_chn(i) = mig_rat_chn(i - 1);
        mig_gen_chn(i) = mig_gen_chn(i - 1);
        frq_pth_chn(i) = frq_pth_chn(i - 1);
        log_lik_chn(i) = log_lik_chn(i - 1);
        apt_rate = apt_rate - 1;
      }
    }
  cout << "apt_rate: " << apt_rate << endl;
  }
  
  return List::create(Named("sel_cof_chn", sel_cof_chn),
                      Named("sel_gen_chn", sel_gen_chn),
                      Named("mig_rat_chn", mig_rat_chn),
                      Named("mig_gen_chn", mig_gen_chn),
                      Named("log_lik_chn", log_lik_chn),
                      Named("frq_pth_chn", frq_pth_chn));
}
/*************************/
