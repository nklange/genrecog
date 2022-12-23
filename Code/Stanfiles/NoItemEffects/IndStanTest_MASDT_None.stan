functions {
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
  /* cumulative-probit log-PDF for a single response
   * Args:
   *   y: response category
   *   mu: latent mean parameter
   *   disc: discrimination parameter
   *   thres: ordinal thresholds
   * Returns:
   *   a scalar to be added to the log posterior
   */
   vector cumulative_probit_pmf(real mu, real disc, vector thres, real mix, real nonattmu) {
     int nthres = num_elements(thres);
     vector[(nthres+1)] p;
     
     if(mix){
       
       for(y in 1:(nthres+1)){
         if (y == 1) {
           p[y] = disc * Phi(thres[1] - mu) + (1 - disc) * Phi(thres[1] - nonattmu);
         } else if (y == nthres + 1) {
           p[y] = disc * (1 - Phi(thres[nthres] - mu)) + (1 - disc) * (1 - Phi(thres[nthres] - nonattmu));
         } else {
           p[y] = disc * (Phi(thres[y] - mu) -
           Phi(thres[y - 1] - mu)) + 
           (1 - disc) * (Phi(thres[y] - nonattmu) -
           Phi(thres[y - 1] - nonattmu));
         }
       }
       
     } else {
       
       for(y in 1:(nthres+1)){
         if (y == 1) {
           p[y] = Phi(thres[1] - mu);
         } else if (y == nthres + 1) {
           p[y] = 1 - Phi(thres[nthres] - mu);
         } else {
           p[y] = Phi(thres[y] - mu) -
           Phi(thres[y - 1] - mu);
         }
       }
     }
     
     
     return p;
   }
}
data {
  int<lower=1> N;  // total number of observations (Conditions x ID)
  int<lower=2> C; // total number of response/rating categories
  int<lower=2> nthres;  // number of thresholds (C - 1)
  int Y[N, C];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_disc;  // number of population-level effects
  matrix[N, K_disc] X_disc;  // population-level design matrix
  int<lower=1> K_nonattmu;  // number of population-level effects
  matrix[N, K_nonattmu] X_nonattmu;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels / participants
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  int<lower=1> G; // number of group-level effects
  int<lower=1> G_disc; // number of group-level effects
  int<lower=1> G_nonattmu; // number of group-level effects
  matrix[N, G] Z;
  matrix[N, G_disc] Z_disc;
  matrix[N, G_nonattmu] Z_nonattmu;
  //int<lower=1> NC_1;  // number of group-level correlations
  // data for group-level effects of ID 2
  int prior_only;  // should the likelihood be ignored?
}
parameters {
  vector[K] b;  // population-level effects
    // temporary thresholds for centered predictors
  vector[K_disc] b_disc;  // population-level effects
  vector[K_nonattmu] b_nonattmu;
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  ordered[nthres] Intercept[N_1]; // only per-person not per-condition intercept
  ordered[nthres] mu_cr;
  real<lower=0> sigma_cr[nthres];
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  matrix[N_1, G] r_1_mu;
  matrix[N_1, G_disc] r_1_disc;
  matrix[N_1, G_nonattmu] r_1_nonattmu;
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_mu = r_1[,1:G];
  r_1_disc = r_1[,(G + 1):(G + G_disc)];
  r_1_nonattmu = r_1[,(G + G_disc + 1):M_1];
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = X * b;
    // initialize linear predictor term
    // contrast matrix is multiplied over b_disc to make sure that recollection
    // can vary across conditions, i.e., calculate rec cond 1, delta rec cond 2 etc
    // it is coded so that when new: disc = 0
    // problem: as Phi(disc = 0) = 0.5, recollection != 0 for new items
    vector[N] disc = X_disc * b_disc;
    vector[N] nonattmu = X_nonattmu * b_nonattmu;
    matrix[N,C] theta;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += dot_product(r_1_mu[J_1[n],], Z[n,]);
      disc[n] += dot_product(r_1_disc[J_1[n],],Z_disc[n,]);
      nonattmu[n] += dot_product(r_1_nonattmu[J_1[n],], Z_nonattmu[n,]);
    }
    //for (n in 1:N) {
      // apply the inverse link function
      // make sure than disc[n] = recollection= 0 for new items
      // by multiplying with isold main effect
      // now: old items = 1 * Phi(disc[n]); new items = 0 * Phi(disc[n])

    //}
    for (n in 1:N){
        disc[n] = X_disc[n,1] * Phi(disc[n]);
        nonattmu[n] = X_disc[n,1] * Phi(nonattmu[n]) ;
        theta[n,] = to_row_vector(cumulative_probit_pmf(mu[n], disc[n], Intercept[J_1[n]], X_disc[n,1],nonattmu[n]*mu[n]));
    }
    for (n in 1:N) {
      target += multinomial_lpmf( to_array_1d(Y[n,]) | to_vector(theta[n,]));
    }
  }
  // priors including constants

  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - M_1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(to_vector(z_1));
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
  target += student_t_lpdf(b[1] | 4, 1, 1.5);
  //target += student_t_lpdf(b[2:K] | 4, 0, 1/sqrt(2));
  target += normal_lpdf(b_disc[1] | 0, 1);
  target += normal_lpdf(b_nonattmu[1] |-0.4,0.7);
  //target += student_t_lpdf(b_disc[2:K] | 4, 0, 0.4);
  target += student_t_lpdf(sigma_cr | 3, 0, 2)  - 1 * student_t_lccdf(0 | 3, 0, 2);
  for (i in 1:nthres) {
  target += student_t_lpdf(mu_cr[i] | 4,
                            -0.75 + (3.0 / (nthres - 1)) * (i - 1),
                           0.5);
  Intercept[,i] ~ normal(mu_cr[i], sigma_cr[i]);
};
}
