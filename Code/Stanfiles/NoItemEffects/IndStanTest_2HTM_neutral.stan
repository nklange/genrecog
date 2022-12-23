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
   vector cumulative_probit_pmf(real DO, real DN, real go, real neut, vector s_m, vector a_m_o, vector a_m_n, real old) {
    
     int C = num_elements(a_m_o) + num_elements(a_m_n) + 1;
     real Cold = ceil(C / 2.0);
     vector[C] p;

     if(old){
      for(y in 1:C){
        if (y > Cold){
          p[y] = DO * s_m[C + 1 - y] + (1-DO) * (1-neut)* go * a_m_o[C + 1 - y];
        } else if (y == Cold){
          
          p[y] = (1-DO) * neut;
          
        } else {
          p[y] = (1-DO) * (1-neut) * (1-go) * a_m_n[y];    
        }
        
      }
     } else {
       
       for(y in 1:C){
        if (y < Cold) {
          p[y] = DN * s_m[y] + (1-DN) * (1-neut) *(1-go) * a_m_n[y];
        } else if (y == Cold){
          
          p[y] = (1-DN) * neut;
          
        } else {
          p[y] = (1-DN) * (1-neut) *(go) * a_m_o[C + 1 - y];    
        }
        
      }
     }
  return p;
   }
}
data {
  int<lower=1> N;  // total number of observations (Conditions x ID)
  int<lower=2> C; // total number of response/rating categories
  int<lower=2> Chalf;
  //int<lower=2> nthres;  // number of thresholds (C - 1)
  int Y[N, C];  // response variable
  int<lower=1> K_DO;  // number of population-level effects
  matrix[N, K_DO] X_DO;  // population-level design matrix
  int<lower=1> K_DN;  // number of population-level effects
  matrix[N, K_DN] X_DN;  // population-level design matrix
  int<lower=1> K_go;  // number of population-level effects
  matrix[N, K_go] X_go;  // population-level design matrix
  int<lower=1> K_neut;  // number of population-level effects
  matrix[N, K_neut] X_neut;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels / participants
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  int<lower=1> G_DO; // number of group-level effects
  int<lower=1> G_DN; // number of group-level effects
  int<lower=1> G_go;
  int<lower=1> G_neut;
  matrix[N, G_DO] Z_DO;
  matrix[N, G_DN] Z_DN;
  matrix[N, G_go] Z_go;
   matrix[N, G_neut] Z_neut;
  //int<lower=1> NC_1;  // number of group-level correlations
  // data for group-level effects of ID 2
  int prior_only;  // should the likelihood be ignored?
  vector[Chalf] prior_alpha_mu;
  vector[Chalf] prior_alpha_scale;
  vector[Chalf] prior_alpha_mu_a;
  vector[Chalf] prior_alpha_scale_a;
}
parameters {
  vector[K_DO] b_DO;  // population-level effects
    // temporary thresholds for centered predictors
  vector[K_DN] b_DN;  // population-level effects
  vector[K_go] b_go;  // population-level effects
  vector[K_neut] b_neut;
  //real mu_go;
  //real<lower=0> sigma_go;
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  //ordered[nthres] Intercept[N_1]; // only per-person not per-condition intercept
  //ordered[nthres] mu_cr;
 
  vector<lower=0>[Chalf] mu_s;
  vector<lower=0>[Chalf] mu_a_o;
  vector<lower=0>[Chalf] mu_a_n;
  //real<lower=0> sigma_cr[nthres];
  simplex[Chalf] s_m[N_1];
  simplex[Chalf] a_m_o[N_1];
  simplex[Chalf] a_m_n[N_1];
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  matrix[N_1, G_DO] r_1_DO;
  matrix[N_1, G_DN] r_1_DN;
  matrix[N_1, G_go] r_1_go;
  matrix[N_1, G_neut] r_1_neut;
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_DO = r_1[,1:G_DO];
  r_1_DN = r_1[,(G_DO + 1):(G_DO + G_DN)];
  r_1_go = r_1[,(G_DO + G_DN+1):(G_DO + G_DN + G_neut)];
  r_1_neut = r_1[,(G_DO + G_DN + G_neut +1):(M_1)];
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] DO = X_DO * b_DO;
    // initialize linear predictor term
    // contrast matrix is multiplied over b_disc to make sure that recollection
    // can vary across conditions, i.e., calculate rec cond 1, delta rec cond 2 etc
    // it is coded so that when new: disc = 0
    // problem: as Phi(disc = 0) = 0.5, recollection != 0 for new items
    vector[N] DN = X_DN * b_DN;
    vector[N] go = X_go * b_go;
    vector[N] neut = X_neut * b_neut;
    matrix[N,C] theta;
    for (n in 1:N) {
      // add more terms to the linear predictor
      DO[n] += dot_product(r_1_DO[J_1[n],],Z_DO[n,]);
      DN[n] += dot_product(r_1_DN[J_1[n],],Z_DN[n,]);
      go[n] += dot_product(r_1_go[J_1[n],],Z_go[n,]);
      neut[n] += dot_product(r_1_neut[J_1[n],],Z_neut[n,]); 
    }
    //for (n in 1:N) {
      // apply the inverse link function
      // make sure than disc[n] = recollection= 0 for new items
      // by multiplying with isold main effect
      // now: old items = 1 * Phi(disc[n]); new items = 0 * Phi(disc[n])

    //}
    for (n in 1:N){
        DO[n] =  Phi(DO[n]);
        DN[n] =  Phi(DN[n]);
        go[n] =  Phi(go[n]);
        neut[n] = Phi(neut[n]);
        theta[n,] = to_row_vector(cumulative_probit_pmf(DO[n], DN[n], go[n], neut[n], s_m[J_1[n],], a_m_o[J_1[n],], a_m_n[J_1[n],], X_DO[n,1]));
    }
    for (n in 1:N) {
      target += multinomial_lpmf( to_array_1d(Y[n,]) | to_vector(theta[n,]));
    }
  }
  // priors including constants

  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - M_1 * student_t_lccdf(0 | 3, 0, 2.5);

  target += std_normal_lpdf(to_vector(z_1));
  target += lkj_corr_cholesky_lpdf(L_1 | 4);
  target += normal_lpdf(b_DO[1] | 0,.75);
  target += normal_lpdf(b_DN[1] | 0,.75);
  target += normal_lpdf(b_DO[2:K_DO] | 0,.75);
  target += normal_lpdf(b_DN[2:K_DO] | 0,.75);
  target += normal_lpdf(b_go | 0,0.6);
  target += normal_lpdf(b_neut | -.4,1);
  for (i in 1:Chalf){
    target += student_t_lpdf(mu_s[i] | 4,prior_alpha_mu[i],prior_alpha_scale[i]);

  }
  for (i in 1:Chalf){
    target += student_t_lpdf(mu_a_o[i] | 4,prior_alpha_mu_a[i],prior_alpha_scale_a[i]);
    target += student_t_lpdf(mu_a_n[i] | 4,prior_alpha_mu_a[i],prior_alpha_scale_a[i]);
  }

  s_m ~ dirichlet(mu_s);
  a_m_o ~ dirichlet(mu_a_o);
  a_m_n ~ dirichlet(mu_a_n);

  
  
}
