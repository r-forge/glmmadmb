DATA_SECTION

  init_int n					// Number of observations
  init_vector y(1,n)				// Observation vector
  init_int p					// Number of fixed effects
  init_matrix X(1,n,1,p)			// Design matrix for fixed effects
  init_int M					// Number of RE blocks (crossed terms)
  init_ivector q(1,M)				// Number of levels of each RE block
  init_ivector m(1,M)				// Number of random effects within each block
  int sum_mq
  init_int ncolZ
  init_matrix Z(1,n,1,ncolZ)			// Design matrix for random effects
  init_imatrix I(1,n,1,ncolZ)			// Index vectors into joint RE vector "u" for each
  init_ivector cor_flag(1,M)			// Indicator for wether each RE block should be correlated
  init_ivector cor_block_start(1,M) 		// Indices for blocks of correlated random effects
  init_ivector cor_block_stop(1,M) 
  init_int numb_cor_params			// Total number of correlation parameters to be estimated
  init_int like_type_flag   			// 0 neg bin 1 poisson
  init_int no_rand_flag   			// 0 have random effects 1 no random effects
  init_int zi_flag				// 1=zi, 0=no zi
  init_int intermediate_maxfn
  init_int has_offset				// 0=no offset, 1=with offset
  init_vector offset(1,n)				// Offset vector
  matrix rr(1,n,1,6)
  matrix phi(1,p,1,p)
 LOC_CALCS
  int i,j;
  phi.initialize();
  for (i=1;i<=p;i++)
  {
    phi(i,i)=1.0;
  }
  dmatrix trr=trans(rr);
  trr(6).fill_seqadd(1,1);
  rr=trans(trr);

  dmatrix TX(1,p,1,n);
  TX=trans(X);
  for (i=1;i<=p;i++)
  {
    double tmp=norm(TX(i));
    TX(i)/=tmp;
    phi(i)/=tmp;
    for (j=i+1;j<=p;j++)
    {
      double a=TX(j)*TX(i);
      TX(j)-=a*TX(i);
      phi(j)-=a*phi(i);
    }
  }
  X=trans(TX);

  sum_mq = 0;
  for (i=1;i<=M;i++)
    sum_mq += m(i)*q(i);

INITIALIZATION_SECTION
 tmpL 1.0
 tmpL1 0.0
 log_alpha 1
 pz .001

PARAMETER_SECTION
 LOC_CALCS
  int alpha_phase = like_type_flag==0 ? 1 : -1;         // Phase 1 if active
  int zi_phase = zi_flag ? 2 : -1;                      // Phase 2 if active
  int rand_phase = no_rand_flag==0 ? 2+zi_flag : -1;    // Right after zi
  int cor_phase = (rand_phase>0) && (sum(cor_flag)>0) ? rand_phase+1 : -1 ; // Right after rand_phase
  ivector ncolS(1,M);
  ncolS = m;                                            // Uncorrelated random effects
  for (int i=1;i<=M;i++)                                // Modifies the correlated ones
    if(cor_flag(i))
      ncolS(i) = m(i)*(m(i)+1)/2;
  int nS = sum(ncolS);             
  cout << "nS=" << nS << endl;
 END_CALCS
  
  init_bounded_number pz(.000001,0.999,zi_phase)
  init_vector beta(1,p,1)     	// Fixed effects
  sdreport_vector real_beta(1,p)     
  init_bounded_vector tmpL(1,ncolZ,-10,10.5,rand_phase)		// Log standard deviations of random effects
  init_bounded_vector tmpL1(1,numb_cor_params,-10,10.5,cor_phase)	// Offdiagonal elements of cholesky-factor of correlation matrix
  init_bounded_number log_alpha(-5.,6.,alpha_phase)	
  sdreport_vector S(1,nS)
  random_effects_vector u(1,sum_mq,rand_phase)    // Pool of random effects 
  objective_function_value g                    	   // Log-likelihood

PRELIMINARY_CALCS_SECTION
  cout << setprecision(4);

PROCEDURE_SECTION
  g=0.0;

  int i;

  if(!no_rand_flag)
    for (i=1;i<=sum_mq;i++)
      n01_prior(u(i));			// u's are N(0,1) distributed

  for(i=1;i<=n;i++)
    log_lik(i,tmpL,tmpL1,u(I(i)),beta,log_alpha,pz);

  if (sd_phase())
  {
    real_beta=beta*phi;

    int i,j,i_m;
    int i1=1, i2=1, ii=1;
    for (i_m=1;i_m<=M;i_m++)
    {
      dvar_matrix L(1,m(i_m),1,m(i_m));
      L.initialize();
      dvar_matrix tmpS(1,m(i_m),1,m(i_m));
      tmpS.initialize();

      if(cor_flag(i_m))
      {
        int ii=1;
        L(1,1)=1;
        for (i=1;i<=m(i_m);i++)
        {
          L(i,i)=1;
          for (int j=1;j<i;j++)
            L(i,j)=tmpL1(i2++);
          L(i)(1,i)/=norm(L(i)(1,i));
        }
        for (i=1;i<=m(i_m);i++)
          L(i)*=exp(tmpL(i1++));
      }
      else
      {
        for (i=1;i<=m(i_m);i++)
          L(i,i)=exp(tmpL(i1++));
      }

      tmpS=L*trans(L);

      for (i=1;i<=m(i_m);i++)
      {
          if(cor_flag(i_m))
            for(j=1;j<i;j++)
              S(ii++) = tmpS(i,j);
          S(ii++) = tmpS(i,i);
      }
    }

  }

SEPARABLE_FUNCTION void n01_prior(const prevariable&  u)
 g -= -0.5*log(2.0*3.1415926535) - 0.5*square(u);

SEPARABLE_FUNCTION void log_lik(int _i, const dvar_vector& tmpL,const dvar_vector& tmpL1,const dvar_vector& ui,const dvar_vector& beta,const prevariable& log_alpha,const prevariable& pz)
  
  int i,j, i_m;
  double e1=1.e-20;
  double e2=1.e-20;

  dvariable alpha = e2+exp(log_alpha);

  // Construct random effects vector with proper var-covar structure from u
  dvar_vector b(1,ncolZ);

  int i1=1, i2=1;
  for (i_m=1;i_m<=M;i_m++)
  {
    dvar_matrix L(1,m(i_m),1,m(i_m));
    L.initialize();

    if(cor_flag(i_m))
    {
      L(1,1)=1;
      for (i=1;i<=m(i_m);i++)
      {
        L(i,i)=1;
        for (int j=1;j<i;j++)
          L(i,j)=tmpL1(i2++);
        L(i)(1,i)/=norm(L(i)(1,i));
      }
      for (i=1;i<=m(i_m);i++)
        L(i)*=exp(tmpL(i1++));
    }
    else
    {
      for (i=1;i<=m(i_m);i++)
        L(i,i)=exp(tmpL(i1++));
    }

    int upper = sum(m(1,i_m));
    int lower = upper-m(i_m)+1;

    b(lower,upper) = (L*(ui(lower,upper).shift(1))).shift(lower);	// L*ui

  }

  dvariable eta = X(_i)*beta + Z(_i)*b;

  if(has_offset)
    eta += offset(_i);
  dvariable lambda = e1+mfexp(eta);            	  
  dvariable tau = 1.0+e1+lambda/alpha;
  dvariable tmpl; 				// Log likelihood

  int cph=current_phase();
  switch(like_type_flag)
  {
    case 0:   // neg binomial
      if (cph<2)
        tmpl = -square(log(1.0+y(_i))-log(1.0+lambda));
      else
        tmpl = log_negbinomial_density(y(_i),lambda,tau);
      break;
    case 1:   // Poisson
      tmpl = log_density_poisson(y(_i),lambda);
      break;
    default:
      cerr << "Illegal value for like_type_flag" << endl;
      ad_exit(1);
  }
	
  // Zero inflation part
  if(zi_flag)
    if(y(_i)==0)
      g -= log(e2+pz+(1.0-pz)*exp(tmpl));
    else
      g -= log(e2+(1.0-pz)*exp(tmpl));
  else
    g -= tmpl;

REPORT_SECTION
  report << beta*phi << endl;

TOP_OF_MAIN_SECTION
  arrmblsize=40000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  gradient_structure::set_MAX_NVAR_OFFSET(2000000);

