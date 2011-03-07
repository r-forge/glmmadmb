DATA_SECTION

  init_matrix Zf(1,1250,1,5)
  init_matrix Zg(1,1250,1,5)
  init_matrix X(1,1250,1,2)
  init_vector y(1,1250)
  init_int nobs

PARAMETER_SECTION

  objective_function_value f
  init_vector beta(1,2)
  init_number sigma_f
  init_number sigma_g
  vector eta(1,nobs)
  vector mu(1,nobs)
  random_effects_vector u_f(1,5)
  random_effects_vector u_g(1,5)
TOP_OF_MAIN_SECTION
  arrmblsize = 500000;

  
PROCEDURE_SECTION

  eta = X*beta + sigma_g*(Zg*u_g) + sigma_f*(Zf*u_f);
  mu = exp(eta);
  // binomial log-likelihood (unnormalized)
  int i;
  for (i=1; i<=nobs; i++) {
    f -= log_density_poisson(y(i),mu(i));
  }
  
  f+=0.5*(norm2(u_f)+norm2(u_g));  // log-prior (standard normal)
