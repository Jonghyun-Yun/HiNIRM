#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <omp.h>

#define nITEM    72
#define nSCHOOL  62
#define nDIM      2
#define nCOV      2
#define nMAX     81

static int **dataset;
static int totalsize, MM;
static int iter, niter, nburn, repeat, thin, print;
static int *ncount, **jump_index;
static double *oldsigma, *oldtau;
static double *olddelta, **oldmu;
static double *oldgamma, *oldvarphi;
static double *jump_Z, jump_beta, jump_theta, jump_mu, jump_W;
static double prior_a = 0.001, prior_b = 0.001;
static double pr_mean_beta = 0.0, pr_var_beta = 10.0;
static double pr_mean_theta = 0.0, pr_var_theta = 10.0;
static double pr_mean_Z = 0.0, pr_var_Z;
static double pr_mean_mu = 0.0, pr_var_mu = 10.0;
static double pr_mean_delta = 0.0, pr_var_delta = 10.0;
static double pr_mean_gamma = 0.0, pr_var_gamma = 10.0;

typedef struct school_type{
	int cbsize;
	int equal;
	int *count_samp, *count_item;
	int **dataset, ***Y, ***U;
	double *oldbeta, *newbeta, *oldtheta, *newtheta;
	double **old_Zsamp, **new_Zsamp, **old_Zmean, **new_Zmean, **old_Zitem, **new_Zitem;
	double **sample_beta,  *sum_beta,  *acc_beta,  *var_beta;
	double **old_item_mat, **new_item_mat;
	double oldsigma, *mean_Z;
	double post_a, post_b;
	double *sample_sigma, sum_sigma, var_sigma;
	double **sample_theta, *sum_theta, *acc_theta, *var_theta;
	double ***sample_Zsamp, **sum_Zsamp, **var_Zsamp, *acc_Zsamp;
	double ***sample_Zitem, **sum_Zitem, **var_Zitem, *acc_Zitem;
	double **sample_item_mat, **sum_item_mat, **var_item_mat;
}YEWON;

static YEWON *SCHOOL;
/*
Define Variables for Dataset
  nITEM: the number of item;
  nSAMPLE: the number of respondent;
  nDIM: the number of dimension;
  nCOV: the number of covariates; generally assumed as 0.
  **dataset: item response dataset (nSAMPLE X nITEM)
Define Variables for MCMC
  iter: current iterations at each run;
  niter: total number of iterations;
  nburn: burn-in point;
  repeat: total number of MCMC runs;
  thin: save one sample per thin iterations;
  print: print outcome values per print iterations;
  jump_Z: jumping rule for updating Z;
  jump_beta: jumping rule for updating beta;
  jump_sigma: jumping rule for updating log(sigma);
Define Variables for Prior
  pr_mean_beta: prior mean of beta;
  pr_var_beta: prior variance of beta;
  pr_mean_beta: prior mean of Z;
Define Variables for Parameters
  **old_Z: current values of latent positions for respondents (nSAMPLE * nDIM);
  **new_Z: proposed values of latent positions for respondents (nSAMPLE * nDIM);
  oldsigma: current log-transformed values of variance for latent positions;
  newsigma: proposed log-transformed values of variance for latent positions;
Define Variables in typedef network
  totalsize: memory size of typedef network;
  **network: construct network based on item response dataset (nSAMPLE * nSAMPLE);
  oldbeta: current values of difficulty parameters;
  newbeta: proposed values of difficulty paramters;
*/
