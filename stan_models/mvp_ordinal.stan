functions{
  // Monotonic Transform
  
  real mo(vector scale, int i) {
    if (i == 0) {
      return 0;
    } else {
      return rows(scale) * sum(scale[1:i]);
    }
  }

  // Numerically stable approximation to Phi()
  // https://github.com/CerulloE1996/dta-ma-mvp-1
  
  real phi_logit_approx(real b) {
    real a;
    a = inv_logit( 1.702*b );
    return(a);
          }
  // Numerically stable approximation to inv_Phi()
  // https://github.com/CerulloE1996/dta-ma-mvp-1
  
  real inv_phi_logit_approx(real b) {
    real a;
    a = (1/1.702)*logit(b);    
    return(a);
   }

}
data {
  int<lower=1> D; // # of response variables
  int<lower=0> N; // # of individuals
  array[N, D] int<lower=1> y; // response vector of ordinal values
  array [N] int X; // ordinal predictor
  int cats; // # of ordinal categories
}
transformed data{
  int n_thr = cats-1; // # of threshold parameters
}
parameters {
  cholesky_factor_corr[D] L; // cholesky factors
  array[N, D] real<lower=0, upper=1> u; // absorbs inequality constraints
  array[D] ordered[n_thr] thr; // threshold parameters
  simplex[3] mono_effects; // monotonic simplex
  array [D] real<lower=0> alpha; // mean parameter
  array [D] real<lower=0> beta; // mean parameter
}
model {
  
  // PRIORS //
  
  alpha ~ normal(0,5);
  beta ~ normal(0,5);
  L ~ lkj_corr_cholesky(4);
  mono_effects ~ dirichlet([1,1,1]);
  for(d in 1:D){
    for(i in 1:n_thr){
      thr[d,i] ~ normal(i+1,1);
    }
  }
  
  // LIKELIHOOD //
  {
    for (n in 1 : N) {
      array [N] vector[D] z;
      array [N] vector[D] mu;
      real prev;
      prev = 0;
      for (d in 1 : D) {
        mu[n,d] = alpha[d] + beta[d]*mo(mono_effects,X[n]);
        // Phi and inv_Phi may overflow and / or be numerically inaccurate
        real t;
        if (y[n,d] == 99) {
          z[n,d] = inv_phi_logit_approx(u[n,d]);
          target += log(1);
        } else if (y[n, d] == 1) {
          real ub = phi_logit_approx((thr[d,1] - (mu[n,d] + prev)) / L[d,d]);
          t = ub * u[n, d];
          z[n,d] = inv_phi_logit_approx(t); // implies utility is positive
          target += log(ub); // Jacobian adjustment
        } else if (y[n,d] == n_thr + 1) {
          real lb = phi_logit_approx((thr[d, n_thr] -(mu[n,d] + prev)) / L[d,d]);
          t = lb + (1 - lb) * u[n,d];
          z[n,d] = inv_phi_logit_approx(t); // implies utility is negative
          target += log1m(lb); // Jacobian adjustment
        } else{
          real lb = phi_logit_approx((thr[d, y[n,d] - 1] -(mu[n,d] + prev)) / L[d,d]);
          real ub = phi_logit_approx((thr[d, y[n,d]    ] -(mu[n,d] + prev)) / L[d,d]);
          t = lb + (ub - lb) * u[n,d];
          z[n,d] = inv_phi_logit_approx(t); // implies utility is negative
          target += log(ub - lb);
        }
        if (d < D) {
          prev = L[d + 1, 1 : d] * head(z[n], d);
        }
        // Jacobian adjustments imply z is truncated standard normal
        // thus utility --- mu + L_Omega * z --- is truncated multivariate normal
      }
    }
  }
}
generated quantities {
  corr_matrix[D] Omega;
  Omega = multiply_lower_tri_self_transpose(L); // correlation matrix
}
