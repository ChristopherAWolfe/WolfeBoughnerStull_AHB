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
}
parameters {
  cholesky_factor_corr[D] L;
  array[N, D] real<lower=0, upper=1> u;
  array[D] ordered[n_thr] thr;
}
model {
  L ~ lkj_corr_cholesky(4);
  for(d in 1:D){
    for(i in 1:n_thr){
      thr[d,i] ~ normal(i+1,1);
    }
  }
  
  // implicit: u is iid standard uniform a priori
  {
    for (n in 1 : N) {
      array [N] vector[D] z;
      real prev;
      prev = 0;
      for (d in 1 : D) {
        // Phi and inv_Phi may overflow and / or be numerically inaccurate
        real t;
        if (y[n,d] == 99) {
          z[n,d] = inv_phi_logit_approx(u[n,d]);
          target += log(1);
        } else if (y[n, d] == 1) {
          real ub = phi_logit_approx((thr[d,1] - (0 + prev)) / L[d,d]);
          t = ub * u[n, d];
          z[n,d] = inv_phi_logit_approx(t); // implies utility is positive
          target += log(ub); // Jacobian adjustment
        } else if (y[n,d] == n_thr + 1) {
          real lb = phi_logit_approx((thr[d, n_thr] -(0 + prev)) / L[d,d]);
          t = lb + (1 - lb) * u[n,d];
          z[n,d] = inv_phi_logit_approx(t); // implies utility is negative
          target += log1m(lb); // Jacobian adjustment
        } else{
          real lb = phi_logit_approx((thr[d, y[n,d] - 1] -(0 + prev)) / L[d,d]);
          real ub = phi_logit_approx((thr[d, y[n,d]    ] -(0 + prev)) / L[d,d]);
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
  Omega = multiply_lower_tri_self_transpose(L);
}
