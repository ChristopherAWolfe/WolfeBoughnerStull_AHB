functions{

  real phi_logit_approx(real b) {
    real a;
    a = inv_logit( 1.702*b );
    return(a);
          }

  real inv_phi_logit_approx(real b) {
    real a;
    a = (1/1.702)*logit(b);    
    return(a);
   }

}
data {
  int<lower=1> D;
  int<lower=0> N;
  array[N, D] int<lower=1> y;
  int cats;
}
transformed data{
  int n_thr = cats-1;
  array[N, D] real<lower=0, upper=1> u; // nuisance that absorbs inequality constraints
  
  for(n in 1:N){

    for(d in 1:D){

      u[n,d] = uniform_rng(0,1);

    }

  }
}
parameters {
  cholesky_factor_corr[D] L;
  array[D] ordered[n_thr] thr;  
}
model {

}
generated quantities {
  vector[N] lp;
  corr_matrix[D] Omega;
  Omega = multiply_lower_tri_self_transpose(L);

    {
    for (n in 1 : N) {
      array [N] vector[D] z;
      vector[D] log_lik;
      real prev;
      prev = 0;
      for (d in 1 : D) {
        // phi_logit_approx and inv_phi_logit_approx_logit_approx may overflow and / or be numerically inaccurate
        real t;
        if (y[n,d] == 99) {
          z[n,d] = inv_phi_logit_approx(u[n,d]);
          log_lik[d] = log(1);
        } else if (y[n,d] == 1) {
          real ub = phi_logit_approx((thr[d,1] - (0 + prev)) / L[d,d]);
          t = ub * u[n, d];
          z[n,d] = inv_phi_logit_approx(t); // implies utility is positive
          log_lik[d] = log(ub); // Jacobian adjustment
        } else if (y[n,d] == n_thr + 1) {
          real lb = phi_logit_approx((thr[d, n_thr] -(0 + prev)) / L[d,d]);
          t = lb + (1 - lb) * u[n,d];
          z[n,d] = inv_phi_logit_approx(t); // implies utility is negative
          log_lik [d] = log1m(lb); // Jacobian adjustment
        } else{
          real lb = phi_logit_approx((thr[d, y[n,d] - 1] -(0 + prev)) / L[d,d]);
          real ub = phi_logit_approx((thr[d, y[n,d]    ] -(0 + prev)) / L[d,d]);
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
