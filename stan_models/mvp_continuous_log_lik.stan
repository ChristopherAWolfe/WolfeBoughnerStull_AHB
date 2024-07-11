functions{
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
  array[N,D] int<lower=1> y; // vector of ordinal response variables
  array [N] real X; // continuous predictor
  int cats; // # of categories
}
transformed data{
  int n_thr = cats-1; // # of threshold parameters
  array[N, D] real<lower=0, upper=1> u; // generate nuisance parameters

  for(n in 1:N){

    for(d in 1:D){

      u[n,d] = uniform_rng(0,1);

    }

  }
}
parameters {
  cholesky_factor_corr[D] L;
  array[D] ordered[n_thr] thr;
  array [D] real<lower=0> beta;
}
generated quantities {
  
  vector[N] lp; // log probability
  
    {
    for (n in 1 : N) {
      array [N] vector[D] z;
      array [N] vector[D] mu;
      vector[D] log_lik;
      real prev;
      prev = 0;
      for (d in 1 : D) {
        mu[n,d] = beta[d]*X[n];
        // Phi and inv_Phi may overflow and / or be numerically inaccurate
        real t;
        if (y[n,d] == 99) {
          z[n,d] = inv_phi_logit_approx(u[n,d]);
          log_lik[d] = log(1);
        } else if (y[n,d] == 1) {
          real ub = phi_logit_approx((thr[d,1] - (mu[n,d] + prev)) / L[d,d]);
          t = ub * u[n, d];
          z[n,d] = inv_phi_logit_approx(t); // implies utility is positive
          log_lik[d] = log(ub); // Jacobian adjustment
        } else if (y[n,d] == n_thr + 1) {
          real lb = phi_logit_approx((thr[d, n_thr] -(mu[n,d] + prev)) / L[d,d]);
          t = lb + (1 - lb) * u[n,d];
          z[n,d] = inv_phi_logit_approx(t); // implies utility is negative
          log_lik [d] = log1m(lb); // Jacobian adjustment
        } else{
          real lb = phi_logit_approx((thr[d, y[n,d] - 1] -(mu[n,d] + prev)) / L[d,d]);
          real ub = phi_logit_approx((thr[d, y[n,d]    ] -(mu[n,d] + prev)) / L[d,d]);
          t = lb + (ub - lb) * u[n,d];
          z[n,d] = inv_phi_logit_approx(t); // implies utility is negative
          log_lik[d] = log(ub - lb);
        }
        if (d < D) {
          prev = L[d + 1, 1 : d] * head(z[n], d);
        }
        // Jacobian adjustments imply z is truncated standard normal
        // thus utility --- mu + L_Omega * z --- is truncated multivariate normal
        
        lp[n] = sum(log_lik);
      }
    }
  }
}
