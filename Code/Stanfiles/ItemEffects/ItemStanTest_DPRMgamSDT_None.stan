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
     real cumulative_probit_pmf(int y, real mu, real disc, vector thres, vector s_m) {
     int nthres = num_elements(thres);
     int C = nthres + 1;
     
     real p;
     
     if (y == 1) {
       p = (1.0 - disc) * Phi(thres[1] - mu);
     } else if (y == C) {
       p = disc * s_m[1] + (1.0 - disc) * (1 - Phi(thres[nthres] - mu));
     } else if (y == C - 1){
       p = disc * s_m[2]  + (1.0 - disc) *(Phi(thres[y] - mu) - Phi(thres[y - 1] - mu));
     } else {
      p = (1.0 - disc) * (Phi(thres[y] - mu) - Phi(thres[y - 1] - mu));
     }
    
     return log(p);
   }
}

data {
  int<lower=1> N;  // total number of observations (Conditions x ID)
  int<lower=2> C; // total number of response/rating categories
  int<lower=2> Chalf;
  int<lower=2> nthres;  // number of thresholds (C - 1)
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_disc;  // number of population-level effects
  matrix[N, K_disc] X_disc;  // population-level design matrix
  //int<lower=1> K_nonattmu;  // number of population-level effects
  //matrix[N, K_nonattmu] X_nonattmu;  // population-level design matrix
  matrix[N, nthres] X_gamma;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels / participants
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  int<lower=1> G; // number of group-level effects
  int<lower=1> G_disc; // number of group-level effects
  //int<lower=1> G_nonattmu; // number of group-level effects
  matrix[N, G] Z;
  matrix[N, G_disc] Z_disc;
  //matrix[N, G_nonattmu] Z_nonattmu;
  matrix[N, nthres] Z_gamma;
  int<lower=1> N_2; // number of unique items
  int<lower=1> M_2; // number of coefficients for items (= 1 x parameter, no x by cond)
  int<lower=1> J_2[N]; // which unique item is assigned to which trial
  vector[N] Z_2_mu;
  vector[N] Z_2_disc;
  //vector[N] Z_2_nonattmu;
  //int<lower=1> NC_1;  // number of group-level correlations
  // data for group-level effects of ID 2
  int prior_only;  // should the likelihood be ignored?
  vector[Chalf] prior_alpha_mu;
  vector[Chalf] prior_alpha_scale;
}
parameters {
  vector[K] b;  // population-level effects
    // temporary thresholds for centered predictors
  vector[K_disc] b_disc;  // population-level effects
  //vector[K_nonattmu] b_nonattmu;
  vector[nthres] b_gamma;  // population-level effects
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  vector<lower=0>[M_2] sd_2; // standard deviation for item effects (mean = 0)
  matrix[M_2, N_2] z_2; // standardized group level effects for each item (N_2)
  cholesky_factor_corr[M_2] L_2;  // cholesky factor of correlation matrix
  //ordered[nthres] Intercept[N_1]; // only per-person not per-condition intercept
  //ordered[nthres] mu_cr;
  //real<lower=0> sigma_cr[nthres];
  vector<lower=0>[Chalf] mu_s;
  simplex[Chalf] s_m[N_1];
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  matrix[N_1, G] r_1_mu;
  matrix[N_1, G_disc] r_1_disc;
  //matrix[N_1, G_nonattmu] r_1_nonattmu;
  matrix[N_1, nthres] r_1_gamma;
  matrix[N_2, M_2] r_2;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_2] r_2_mu;
  vector[N_2] r_2_disc;
  //vector[N_2] r_2_nonattmu;
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_mu = r_1[,1:G];
  r_1_disc = r_1[,(G + 1):(G + G_disc)];
 // r_1_nonattmu = r_1[,(G + G_disc + 1):(G + G_disc + G_nonattmu)];
  r_1_gamma = r_1[,(G + G_disc + 1):M_1];
  r_2 = scale_r_cor(z_2, sd_2, L_2);
  r_2_mu = r_2[,1];
  r_2_disc = r_2[,2];
  
 
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
    //vector[N] nonattmu = X_nonattmu * b_nonattmu;
    vector[N] theta;
    matrix[N,nthres] gamma;
    matrix[N,nthres] Intercept;
    
    for(n in 1:N){
      for(j in 1:nthres){
        gamma[n,j] = X_gamma[n,j] * b_gamma[j];
      }
    }
    for(n in 1:N){
      for(j in 1:nthres){
        gamma[n,j]+= r_1_gamma[J_1[n],j] * Z_gamma[n,j];
      }
      Intercept[n,] = 2 * inv_Phi(head(cumulative_sum(softmax(append_row(gamma[n,]', 0))), nthres))';
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += dot_product(r_1_mu[J_1[n],], Z[n,])  + r_2_mu[J_2[n]] * Z_2_mu[n];
      disc[n] += dot_product(r_1_disc[J_1[n],],Z_disc[n,])  + r_2_disc[J_2[n]] * Z_2_disc[n];
     
    }
    //for (n in 1:N) {
      // apply the inverse link function
      // make sure than disc[n] = recollection= 0 for new items
      // by multiplying with isold main effect
      // now: old items = 1 * Phi(disc[n]); new items = 0 * Phi(disc[n])

    //}
    for (n in 1:N){
        disc[n] = X_disc[n,1] * Phi(disc[n]);
       
        target += cumulative_probit_pmf(Y[n], mu[n], disc[n], to_vector(Intercept[n,]), s_m[J_1[n],]);
    }
  }
  // priors including constants

  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - M_1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(to_vector(z_1));
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
  target += lkj_corr_cholesky_lpdf(L_2 | 1);
  target += student_t_lpdf(sd_2 | 3, 0, 2.5)
    - M_2 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(to_vector(z_2));
  target += student_t_lpdf(b[1] | 4, 1, 1.5);
  //target += student_t_lpdf(b[2:K] | 4, 0, 1/sqrt(2));
  target +=  student_t_lpdf(b_disc[1] | 4, -0.3, 0.5);
  //target += normal_lpdf(b_nonattmu[1] | -0.4,0.7);
  //target += normal_lpdf(b_nonattmu[2:K] | 0, 0.2);
  //target += student_t_lpdf(b_disc[2:K] | 4, 0, 0.2);
  s_m ~ dirichlet(mu_s);
  //target += student_t_lpdf(sigma_cr | 3, 0, 2)  - 1 * student_t_lccdf(0 | 3, 0, 2);
  for (i in 1:nthres) {
    target += normal_lpdf(b_gamma[i] | 0,log(100));
  //target += student_t_lpdf(mu_cr[i] | 4,
  //                          -0.75 + (3.0 / (nthres - 1)) * (i - 1),
  //                         0.5);
  //Intercept[,i] ~ normal(mu_cr[i], sigma_cr[i]);
};
}
generated quantities{
  
  matrix[N_1,nthres] crits;
  matrix[N_1,nthres] gammaout;
  for(n in 1:N){
  for(j in 1:nthres){
    gammaout[J_1[n],j] = b_gamma[j] + r_1_gamma[J_1[n],j];
    }
      crits[J_1[n],] = 2 * inv_Phi(head(cumulative_sum(softmax(append_row(gammaout[J_1[n],]', 0))), nthres))';
    }
}
