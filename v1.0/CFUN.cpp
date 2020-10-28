// Inferring natural selection and gene migration in the evolution of chickens from ancient DNA data
// Wenyang Lyu, Xiaoyang Dai, Mark Beaumont, Feng Yu and Zhangyi He

// version 1.0

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
// Calculate the fitness matrix for the Wright-Fisher model with selection and migration
// [[Rcpp::export]]
arma::dmat calculateFitnessMat_arma(const double& sel_cof, const double& dom_par) {
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

// Simulate the allele frequency trajectories according to the Wright-Fisher model with selection and migration
// [[Rcpp::export]]
arma::dmat simulateWFM_arma(const arma::dmat& fts_mat, const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const double& ext_frq, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // declare the allele frequency trajectories
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) + 1);

  // initialise the island allele frequencies in generation 0
  frq_pth.col(0) = int_frq;

  // initialise the continent allele frequencies
  arma::dcolvec con_frq = {0.0, 0.0, ext_frq, 1.0 - ext_frq};

  // simulate the allele frequency trajectories
  arma::dcolvec isl_frq = int_frq;
  for (arma::uword k = 1; k <= arma::uword(lst_gen - int_gen); k++) {
    // calculate the sampling probabilities
    isl_frq = isl_frq % (fts_mat * isl_frq) / as_scalar(isl_frq.t() * fts_mat * isl_frq);
    arma::dcolvec prob = (1 - mig_rat) * isl_frq + mig_rat * con_frq;

    // proceed the Wright-Fisher sampling
    IntegerVector isl_cnt(4);
    R::rmultinom(2 * pop_siz, prob.begin(), 4, isl_cnt.begin());
    isl_frq = as<arma::dcolvec>(isl_cnt) / 2 / pop_siz;
    frq_pth.col(k) = isl_frq;
  }

  // return the allele frequency trajectories under the Wright-Fisher model
  return frq_pth;
}
/*************************/


/********** WFD **********/
// Simulate the allele frequency trajectories according to the Wright-Fisher diffusion with selection and migration using the Euler-Maruyama method
// [[Rcpp::export]]
arma::dmat simulateWFD_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const double& ext_frq, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  // rescale the selection coefficient and the migration rate
  double scl_sel_cof = 2 * pop_siz * sel_cof;
  double scl_mig_rat = 4 * pop_siz * mig_rat;

  // calculate delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  // generate delta W
  arma::dmat dW = pow(dt, 0.5) * arma::randn<arma::dmat>(6, arma::uword(lst_gen - int_gen) * ptn_num);

  // declare the allele frequency trajectories
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);

  // initialise the allele frequencies in generation 0
  frq_pth.col(0) = int_frq;

  // simulate the allele frequency trajectories
  for (arma::uword t = 1; t <= arma::uword(lst_gen - int_gen) * ptn_num; t++) {
    // calculate the drift coefficient vector
    arma::dcolvec mu = arma::zeros<arma::dcolvec>(4);
    mu(0) =  scl_sel_cof * frq_pth(0, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par)) - 0.5 * scl_mig_rat * (frq_pth(0, t - 1));
    mu(1) = -scl_sel_cof * frq_pth(1, t - 1) * (frq_pth(0, t - 1) + frq_pth(2, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par)) - 0.5 * scl_mig_rat * (frq_pth(1, t - 1));
    mu(2) =  scl_sel_cof * frq_pth(2, t - 1) * (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * ((frq_pth(0, t - 1) + frq_pth(2, t - 1)) * dom_par + (frq_pth(1, t - 1) + frq_pth(3, t - 1)) * (1 - dom_par)) - 0.5 * scl_mig_rat * (frq_pth(2, t - 1) - ext_frq);
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
    for (arma::uword i = 0; i < 4; i++) {
      if (frq_pth(i, t) < 0) {
        frq_pth(i, t) = 0;
      }
      if (frq_pth(i, t) > 1) {
        frq_pth(i, t) = 1;
      }
    }
    frq_pth.col(t) = frq_pth.col(t) / sum(frq_pth.col(t));
  }

  // return the allele frequency trajectories under the Wright-Fisher diffusion
  return frq_pth;
}
/*************************/


/********** BPF **********/
// Calculate the possible allele counts in the sample
// [[Rcpp::export]]
arma::imat calculateAlleleCnt_arma(const int& smp_siz, const arma::icolvec& smp_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::imat ptl_cnt = arma::zeros<arma::imat>(4, 1);
  if (smp_cnt.n_elem == 2) {
    for (int i = 0; i <= min(smp_cnt(0), smp_cnt(1)); i++) {
      int j = smp_cnt(0) - i;
      int k = smp_cnt(1) - i;
      if (i + j + k <= smp_siz) {
        ptl_cnt(0, 0) = j;
        ptl_cnt(1, 0) = smp_siz - i - j - k;
        ptl_cnt(2, 0) = i;
        ptl_cnt(3, 0) = k;
        ptl_cnt.insert_cols(0, 1);
      }
    }
    ptl_cnt.shed_cols(0, 0);
  } else {
    ptl_cnt.col(0) = smp_cnt;
  }

  return ptl_cnt;
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

// Initialise the particles in the bootstrap particle filter at the initial time point
// [[Rcpp::export]]
arma::dmat initialiseParticle_arma(const int& mig_gen, const arma::irowvec& smp_gen, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  NumericMatrix part(4, pcl_num);
  if (mig_gen > smp_gen.min()) {
    // gene migration begins before the initial sampling time point
    part(0, _) = runif(pcl_num, 0.0, 1.0);
    part(1, _) = 1 - part(0, _);
  } else {
    // gene migration begins after the initial sampling time point
    for (int i = 0; i < 4; i++) {
      part(i, _) = rgamma(pcl_num, 1.0, 1.0);
    }
    for (int j = 0; j < pcl_num; j++) {
      part(_, j) = part(_, j) / sum(part(_, j));
    }
  }

  return as<arma::dmat>(part);
}

// Generate the particles in the bootstrap particle filter between two consecutive time points
// [[Rcpp::export]]
arma::dmat generateParticle_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const bool& non_sel, const bool& non_mig, const double& ext_frq, const arma::dmat& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::dcube path = arma::zeros<arma::dcube>(4, arma::uword(lst_gen - int_gen) * ptn_num + 1, pcl_num);
  if (non_sel == true && non_mig == true) {
    // no natural selection and gene migration
    for(arma::uword i = 0; i < pcl_num; i++) {
      path.slice(i) = simulateWFD_arma(0, dom_par, 0, pop_siz, ext_frq, int_frq.col(i), int_gen, lst_gen, ptn_num);
    }
  } else if (non_sel == true && non_mig == false) {
    // no natural selection but gene migration
    for(arma::uword i = 0; i < pcl_num; i++) {
      path.slice(i) = simulateWFD_arma(0, dom_par, mig_rat, pop_siz, ext_frq, int_frq.col(i), int_gen, lst_gen, ptn_num);
    }
  } else if (non_sel == false && non_mig == true) {
    // natural selection but no gene migration
    for(arma::uword i = 0; i < pcl_num; i++) {
      path.slice(i) = simulateWFD_arma(sel_cof, dom_par, 0, pop_siz, ext_frq, int_frq.col(i), int_gen, lst_gen, ptn_num);
    }
  } else {
    // natural selection and gene migration
    for(arma::uword i = 0; i < pcl_num; i++) {
      path.slice(i) = simulateWFD_arma(sel_cof, dom_par, mig_rat, pop_siz, ext_frq, int_frq.col(i), int_gen, lst_gen, ptn_num);
    }
  }

  arma::dmat part = path.col(arma::uword(lst_gen - int_gen) * ptn_num);

  return part;
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

  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dcube part_pre = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  arma::dcube part_pst = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);

  // initialise the particles
  cout << "generation: " << all_gen(0) << endl;
  arma::uword smp_ind = 0;
  arma::imat ptl_cnt = calculateAlleleCnt_arma(smp_siz(smp_ind), smp_cnt.col(smp_ind));
  arma::dmat part_tmp = initialiseParticle_arma(mig_gen, smp_gen, pcl_num);
  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < ptl_cnt.n_cols; j++) {
      wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(ptl_cnt.col(j), smp_siz(smp_ind), part_tmp.col(i));
    }
  }

  if (arma::sum(wght_tmp) > 0) {
    arma::dcolvec prob = arma::normalise(wght_tmp, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);

    lik = lik * arma::mean(wght_tmp);
    wght.col(smp_ind) = wght_tmp;
    part_pre.slice(smp_ind) = part_tmp;
    part_pst.slice(smp_ind) = part_tmp.cols(indx);

    part_tmp = part_pst.slice(smp_ind);
  } else {
    lik = 0;
    wght.shed_cols(smp_ind, smp_gen.n_elem - 1);
    part_pre.shed_slices(smp_ind, smp_gen.n_elem - 1);
    part_pst.shed_slices(smp_ind, smp_gen.n_elem - 1);

    return List::create(Named("lik", lik),
                        Named("wght", wght),
                        Named("part_pre_resmp", part_pre),
                        Named("part_pst_resmp", part_pst));
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < all_gen.n_elem; k++) {
    bool non_sel = (sel_gen > all_gen(k - 1)) ? true : false;
    bool non_mig = (mig_gen > all_gen(k - 1)) ? true : false;

    cout << "generation: " << all_gen(k) << endl;
    if (arma::any(smp_gen == all_gen(k))) {
      smp_ind = smp_ind + 1;
      arma::imat ptl_cnt = calculateAlleleCnt_arma(smp_siz(smp_ind), smp_cnt.col(smp_ind));
      part_tmp = generateParticle_arma(sel_cof, dom_par, mig_rat, pop_siz, non_sel, non_mig, ext_frq, part_tmp, all_gen(k - 1), all_gen(k), ptn_num, pcl_num);
      wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
      for (arma::uword i = 0; i < pcl_num; i++) {
        for (arma::uword j = 0; j < ptl_cnt.n_cols; j++) {
          wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(ptl_cnt.col(j), smp_siz(smp_ind), part_tmp.col(i));
        }
      }
      if (arma::sum(wght_tmp) > 0) {
        arma::dcolvec prob = arma::normalise(wght_tmp, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);

        lik = lik * arma::mean(wght_tmp);
        wght.col(smp_ind) = wght_tmp;
        part_pre.slice(smp_ind) = part_tmp;
        part_pst.slice(smp_ind) = part_tmp.cols(indx);

        part_tmp = part_pst.slice(smp_ind);
      } else {
        lik = 0;
        wght.shed_cols(smp_ind, smp_gen.n_elem - 1);
        part_pre.shed_slices(smp_ind, smp_gen.n_elem - 1);
        part_pst.shed_slices(smp_ind, smp_gen.n_elem - 1);

        break;
      }
    } else {
      part_tmp = generateParticle_arma(sel_cof, dom_par, mig_rat, pop_siz, non_sel, non_mig, ext_frq, part_tmp, all_gen(k - 1), all_gen(k), ptn_num, pcl_num);
    }
  }

  return List::create(Named("lik", lik),
                      Named("wght", wght),
                      Named("part_pre_resmp", part_pre),
                      Named("part_pst_resmp", part_pst));
}
/*************************/


/****** PMMHwGibbs *******/
// Calculate the log-likelihood using the bootstrap particle filter
// [[Rcpp::export]]
double calculateLogLikelihood_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const int& sel_gen, const int& mig_gen, const double& ext_frq, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::field<arma::imat>& ptl_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::irowvec all_gen(smp_gen.n_elem + 2);
  all_gen.subvec(1, smp_gen.n_elem) = smp_gen;
  all_gen.head(1) = sel_gen;
  all_gen.tail(1) = mig_gen;
  all_gen = arma::unique(all_gen);
  all_gen = arma::sort(all_gen);

  double log_lik = 0;

  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_pre = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat part_pst = arma::zeros<arma::dmat>(4, pcl_num);

  // initialise the particles
  arma::uword smp_ind = 0;
  arma::imat ale_cnt = ptl_cnt(smp_ind);
  part_pre = initialiseParticle_arma(mig_gen, smp_gen, pcl_num);
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < ale_cnt.n_cols; j++) {
      wght(i) = wght(i) + calculateMultinomProb_arma(ale_cnt.col(j), smp_siz(0), part_pre.col(i));
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

    return log_lik;
  }

  // run the bootstrap particle filter
  for (arma::uword k = 1; k < all_gen.n_elem; k++) {
    bool non_sel = (sel_gen > all_gen(k - 1)) ? true : false;
    bool non_mig = (mig_gen > all_gen(k - 1)) ? true : false;

    if (arma::any(smp_gen == all_gen(k))) {
      smp_ind = smp_ind + 1;
      arma::imat ale_cnt = ptl_cnt(smp_ind);
      part_pre = generateParticle_arma(sel_cof, dom_par, mig_rat, pop_siz, non_sel, non_mig, ext_frq, part_pst, all_gen(k - 1), all_gen(k), ptn_num, pcl_num);
      wght = arma::zeros<arma::dcolvec>(pcl_num);
      for (arma::uword i = 0; i < pcl_num; i++) {
        for (arma::uword j = 0; j < ale_cnt.n_cols; j++) {
          wght(i) = wght(i) + calculateMultinomProb_arma(ale_cnt.col(j), smp_siz(smp_ind), part_pre.col(i));
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

        break;
      }
    } else {
      part_pre = generateParticle_arma(sel_cof, dom_par, mig_rat, pop_siz, non_sel, non_mig, ext_frq, part_pst, all_gen(k - 1), all_gen(k), ptn_num, pcl_num);
      part_pst = part_pre;
    }
  }

  return log_lik;
}

// Calculate the optimal particle number in the particle marginal Metropolis-Hastings
// [[Rcpp::export]]
List calculateOptimalParticleNum_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const int& sel_gen, const int& mig_gen, const double& ext_frq, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateAlleleCnt_arma(smp_siz(k), smp_cnt.col(k));
  }

  arma::drowvec log_lik = arma::zeros<arma::drowvec>(300);
  arma::dmat frq_pth = arma::zeros<arma::dmat>(4, arma::uword(arma::max(smp_gen) - arma::min(smp_gen)) * ptn_num + 1);
  for (arma::uword i = 0; i < 300; i++) {
    log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, mig_rat, pop_siz, sel_gen, mig_gen, ext_frq, smp_gen, smp_siz, ptl_cnt, ptn_num, opt_pcl_num(0));
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

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateAlleleCnt_arma(smp_siz(k), smp_cnt.col(k));
  }

  arma::drowvec sel_cof_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::irowvec sel_gen_chn = arma::zeros<arma::irowvec>(itn_num);
  arma::drowvec mig_rat_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::irowvec mig_gen_chn = arma::zeros<arma::irowvec>(itn_num);

  //arma::drowvec log_pri = arma::zeros<arma::drowvec>(2);
  arma::drowvec log_lik = arma::zeros<arma::drowvec>(2);

  double sel_cof_sd = 5e-03;
  double sel_gen_sd = 5e+01;
  double mig_rat_sd = 5e-03;
  double mig_gen_sd = 5e+01;

  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the population genetic parameters to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn(0) = sel_cof;
  sel_gen_chn(0) = sel_gen;
  mig_rat_chn(0) = mig_rat;
  mig_gen_chn(0) = mig_gen;

  log_lik(0) = calculateLogLikelihood_arma(sel_cof_chn(0), dom_par, mig_rat_chn(0), pop_siz, sel_gen_chn(0), mig_gen_chn(0), ext_frq, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    double apt_rto = 0;
    // draw the candidate of the selection coefficient from the random walk proposal
    sel_cof_chn(i) = sel_cof_chn(i - 1) + sel_cof_sd * arma::randn();
    if (sel_cof_chn(i) > 1) {
      sel_cof_chn(i) = sel_cof_chn(i - 1);
      sel_gen_chn(i) = sel_gen_chn(i - 1);
      mig_rat_chn(i) = mig_rat_chn(i - 1);
      mig_gen_chn(i) = mig_gen_chn(i - 1);
      log_lik(1) = log_lik(0);
    } else {
      apt_rto = 1;
    }
    // draw the candidate of the selection time from the random walk proposal
    sel_gen_chn(i) = sel_gen_chn(i - 1) + int(round(sel_gen_sd * arma::randn()));
    if (sel_gen_chn(i) < smp_gen.min() || sel_gen_chn(i) > smp_gen.max()) {
      sel_cof_chn(i) = sel_cof_chn(i - 1);
      sel_gen_chn(i) = sel_gen_chn(i - 1);
      mig_rat_chn(i) = mig_rat_chn(i - 1);
      mig_gen_chn(i) = mig_gen_chn(i - 1);
      log_lik(1) = log_lik(0);
    } else {
      apt_rto = 1;
    }
    // draw the candidate of the migration rate from the random walk proposal
    mig_rat_chn(i) = mig_rat_chn(i - 1) + mig_rat_sd * arma::randn();
    if (mig_rat_chn(i) > 1) {
      sel_cof_chn(i) = sel_cof_chn(i - 1);
      sel_gen_chn(i) = sel_gen_chn(i - 1);
      mig_rat_chn(i) = mig_rat_chn(i - 1);
      mig_gen_chn(i) = mig_gen_chn(i - 1);
      log_lik(1) = log_lik(0);
    } else {
      apt_rto = 1;
    }
    // draw the candidate of the migration time from the random walk proposal
    mig_gen_chn(i) = mig_gen_chn(i - 1) + int(round(mig_gen_sd * arma::randn()));
    if (mig_gen_chn(i) < smp_gen.min() || mig_gen_chn(i) > smp_gen(arma::find(smp_gen > 0)).min()) {
      sel_cof_chn(i) = sel_cof_chn(i - 1);
      sel_gen_chn(i) = sel_gen_chn(i - 1);
      mig_rat_chn(i) = mig_rat_chn(i - 1);
      mig_gen_chn(i) = mig_gen_chn(i - 1);
      log_lik(1) = log_lik(0);
    } else {
      apt_rto = 1;
    }

    if (apt_rto == 1)  {
      // calculate the proposal
      //arma::drowvec log_psl = arma::zeros<arma::drowvec>(2);

      // calculate the likelihood
      log_lik(1) = calculateLogLikelihood_arma(sel_cof_chn(i), dom_par, mig_rat_chn(i), pop_siz, sel_gen_chn(i), mig_gen_chn(i), ext_frq, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

      // calculate the acceptance ratio
      apt_rto = exp(log_lik(1) - log_lik(0));
      //apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      if (arma::randu() > apt_rto) {
        sel_cof_chn(i) = sel_cof_chn(i - 1);
        sel_gen_chn(i) = sel_gen_chn(i - 1);
        mig_rat_chn(i) = mig_rat_chn(i - 1);
        mig_gen_chn(i) = mig_gen_chn(i - 1);
        log_lik(1) = log_lik(0);
      } else {
        log_lik(0) = log_lik(1);
      }
    }
  }

  return List::create(Named("sel_cof_chn", sel_cof_chn),
                      Named("sel_gen_chn", sel_gen_chn),
                      Named("mig_rat_chn", mig_rat_chn),
                      Named("mig_gen_chn", mig_gen_chn));
}

// Run the particle marginal Metropolis-Hastings within Gibbs
//[[Rcpp::export]]
List runPMMHwGibbs_arma(const double& sel_cof, const double& dom_par, const double& mig_rat, const int& pop_siz, const int& sel_gen, const int& mig_gen, const double& ext_frq, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;

  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    ptl_cnt(k) = calculateAlleleCnt_arma(smp_siz(k), smp_cnt.col(k));
  }

  arma::drowvec sel_cof_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::irowvec sel_gen_chn = arma::zeros<arma::irowvec>(itn_num);
  arma::drowvec mig_rat_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::irowvec mig_gen_chn = arma::zeros<arma::irowvec>(itn_num);

  //arma::drowvec log_pri = arma::zeros<arma::drowvec>(2);
  arma::drowvec log_lik = arma::zeros<arma::drowvec>(2);

  double sel_cof_sd = 5e-03;
  double sel_gen_sd = 5e+01;
  double mig_rat_sd = 5e-03;
  double mig_gen_sd = 5e+01;

  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the population genetic parameters to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn(0) = sel_cof;
  sel_gen_chn(0) = sel_gen;
  mig_rat_chn(0) = mig_rat;
  mig_gen_chn(0) = mig_gen;

  log_lik(0) = calculateLogLikelihood_arma(sel_cof_chn(0), dom_par, mig_rat_chn(0), pop_siz, sel_gen_chn(0), mig_gen_chn(0), ext_frq, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

  double apt_rto = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;

    // update natural selection related parameters in the Gibbs step
    apt_rto = 0;
    // draw the candidate of the selection coefficient from the random walk proposal
    sel_cof_chn(i) = sel_cof_chn(i - 1) + sel_cof_sd * arma::randn();
    if (sel_cof_chn(i) > 1) {
      sel_cof_chn(i) = sel_cof_chn(i - 1);
      sel_gen_chn(i) = sel_gen_chn(i - 1);
      log_lik(1) = log_lik(0);
    } else {
      apt_rto = 1;
    }
    // draw the candidate of the selection time from the random walk proposal
    sel_gen_chn(i) = sel_gen_chn(i - 1) + int(round(sel_gen_sd * arma::randn()));
    if (sel_gen_chn(i) < smp_gen.min() || sel_gen_chn(i) > smp_gen.max()) {
      sel_cof_chn(i) = sel_cof_chn(i - 1);
      sel_gen_chn(i) = sel_gen_chn(i - 1);
      log_lik(1) = log_lik(0);
    } else {
      apt_rto = 1;
    }

    if (apt_rto == 1)  {
      // calculate the proposal
      //arma::drowvec log_psl = arma::zeros<arma::drowvec>(2);

      // calculate the likelihood
      log_lik(1) = calculateLogLikelihood_arma(sel_cof_chn(i), dom_par, mig_rat_chn(i - 1), pop_siz, sel_gen_chn(i), mig_gen_chn(i - 1), ext_frq, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

      // calculate the acceptance ratio
      apt_rto = exp(log_lik(1) - log_lik(0));
      //apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      if (arma::randu() > apt_rto) {
        sel_cof_chn(i) = sel_cof_chn(i - 1);
        sel_gen_chn(i) = sel_gen_chn(i - 1);
        log_lik(1) = log_lik(0);
      }
    }

    // update gene migration related parameters in the Gibbs step
    apt_rto = 0;
    // draw the candidate of the migration rate from the random walk proposal
    mig_rat_chn(i) = mig_rat_chn(i - 1) + mig_rat_sd * arma::randn();
    if (mig_rat_chn(i) > 1) {
      mig_rat_chn(i) = mig_rat_chn(i - 1);
      mig_gen_chn(i) = mig_gen_chn(i - 1);
      log_lik(0) = log_lik(1);
    } else {
      apt_rto = 1;
    }
    // draw the candidate of the migration time from the random walk proposal
    mig_gen_chn(i) = mig_gen_chn(i - 1) + int(round(mig_gen_sd * arma::randn()));
    if (mig_gen_chn(i) < smp_gen.min() || mig_gen_chn(i) > smp_gen(arma::find(smp_gen > 0)).min()) {
      mig_rat_chn(i) = mig_rat_chn(i - 1);
      mig_gen_chn(i) = mig_gen_chn(i - 1);
      log_lik(0) = log_lik(1);
    } else {
      apt_rto = 1;
    }

    if (apt_rto == 1)  {
      // calculate the proposal
      //arma::drowvec log_psl = arma::zeros<arma::drowvec>(2);

      // calculate the likelihood
      log_lik(0) = calculateLogLikelihood_arma(sel_cof_chn(i), dom_par, mig_rat_chn(i), pop_siz, sel_gen_chn(i), mig_gen_chn(i), ext_frq, smp_gen, smp_siz, ptl_cnt, ptn_num, pcl_num);

      // calculate the acceptance ratio
      apt_rto = exp(log_lik(0) - log_lik(1));
      //apt_rto = exp((log_pri(1) + log_lik(1) + log_psl(1)) - (log_pri(0) + log_lik(0) + log_psl(0)));

      if (arma::randu() > apt_rto) {
        mig_rat_chn(i) = mig_rat_chn(i - 1);
        mig_gen_chn(i) = mig_gen_chn(i - 1);
        log_lik(0) = log_lik(1);
      }
    }
  }

  return List::create(Named("sel_cof_chn", sel_cof_chn),
                      Named("sel_gen_chn", sel_gen_chn),
                      Named("mig_rat_chn", mig_rat_chn),
                      Named("mig_gen_chn", mig_gen_chn));
}
/*************************/
