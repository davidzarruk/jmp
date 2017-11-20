


double Aggregation_error(parameters params, const double *mortsubsidy, const int error_funct){
	
	// VFI parameters
  const int T 		= params.T;

  // Grid for savings: a
  const int na  	= params.na;
  const int nm  	= params.nm;
  const int nh  	= params.nh;
  const int nr    = params.nr;
  const int nl    = params.nl;
  const int nd  	= params.nd;
  const int ny  		= params.ny;


  // Initialize grids
  double agrid[na];
  double dgrid[nd];
  double hgrid[nh];
  double rgrid[nr];
  double lgrid[nl];
  double mgrid[nm];
  double ygrid[ny];
  double P[ny*ny];
  double survival[T];
  double repay_coeff[T];
  double eprocess[T];

  // Fill grids
  grid_initialize(params, agrid, dgrid, hgrid, rgrid, lgrid, mgrid, ygrid, P, survival, repay_coeff, eprocess);

  //----------------------------------------------//
  //--------    VALUES AND POLICIES     ----------//
  //----------------------------------------------//

  // On the host memory
  double *Value, *Value_equiv, *Policyc, *Pricing_guess;
  int *Default, *Renew, *Policya, *Policym, *Policyh, *Policyr, *Policyl;
  
  size_t sizemats     = ny*na*nm*nh*nd*T*sizeof(double);
  size_t sizematsint  = ny*na*nm*nh*nd*T*sizeof(int);
  
  Value         = (double*)malloc(sizemats);
  Value_equiv   = (double*)malloc(sizemats);
  Policyc       = (double*)malloc(sizemats);
  Pricing_guess = (double*)malloc(sizemats);

  Default       = (int*)malloc(sizematsint);
  Renew         = (int*)malloc(sizematsint);
  Policya       = (int*)malloc(sizematsint);
  Policym       = (int*)malloc(sizematsint);
  Policyh       = (int*)malloc(sizematsint);
  Policyr       = (int*)malloc(sizematsint);
  Policyl       = (int*)malloc(sizematsint);
  
  int ind;
  
  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              Value[ind]        = 0.0;
              Value_equiv[ind]  = 0.0;
              Policyc[ind]      = 0.0;
              Pricing_guess[ind] = 0.0;

              Default[ind]      = 0;
              Renew[ind]        = 0;
              Policya[ind]      = 0;
              Policym[ind]      = 0;
              Policyh[ind]      = 0;
              Policyr[ind]      = 0;
              Policyl[ind]      = 0;
            }
          }
        }
      }
    }
  }

  double incomeshock[T];
  for(int i=0; i<T; i++){
    incomeshock[i] = 0.0;
  }

  //--------------------------------------//
  //---   Value function computation  ----//
  //--------------------------------------//

  int print_out = 1;

  // for(int i=0; i<T; i++){
  //   cout << repay_coeff[i] << endl;
  // }

  vfi_compute(params, error_funct, 
              P, survival, hgrid, rgrid, lgrid, mgrid,
              agrid, dgrid, ygrid, eprocess, repay_coeff,
              Value, Value_equiv, Policyc, Pricing_guess, 
              Default, Renew, Policya, Policym, Policyh, Policyr, Policyl, print_out, incomeshock, mortsubsidy);


  //--------------------------------------//
  //------  Probability Paths  -----------//
  //--------------------------------------//

  double *Pcond, *Puncond, *Pcond_debt, *Puncond_debt;

  Pcond   = (double*)malloc(sizemats);
  Puncond = (double*)malloc(sizemats);
  Pcond_debt   = (double*)malloc(sizemats);
  Puncond_debt = (double*)malloc(sizemats);

  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              Pcond[ind]        = 0.0;
              Puncond[ind]      = 0.0;
              Pcond_debt[ind]        = 0.0;
              Puncond_debt[ind]      = 0.0;
            }
          }
        }
      }
    }
  }

  probability_paths(Pcond, Puncond, P, params, survival, Policya, Policym, Policyh);

  probability_paths_debt(Pcond_debt, Puncond_debt, P, params, survival, Policya, Policym, Policyh);


  double error_equilib;
  error_equilib = vfi_error(params, error_funct, 0,
                            P, survival, hgrid, rgrid, lgrid, mgrid, repay_coeff,
                            agrid, dgrid, ygrid, eprocess,
                            Value, Policyc, Pricing_guess, 
                            Default, Renew, Policya, Policym, Policyh, Policyr, Policyl,
                            Pcond, Puncond, Pcond_debt, Puncond_debt);

  //--------------------------------------//
  //--- Bank disimbursements/revenues ----//
  //--------------------------------------//

  double disbu;
  disbu = disbursements(params, P, survival, hgrid, rgrid, mgrid,
                        agrid, dgrid, ygrid, eprocess, repay_coeff, Pricing_guess, 
                        Default, Renew, Policya, Policym, Policyh, Puncond);

  double rep;
  rep = repayments(params, P, survival, hgrid, rgrid, mgrid,
                                  agrid, dgrid, ygrid, eprocess, repay_coeff, 
                                  Pricing_guess, Default, Renew, Policya, Puncond);

  double fore;
  fore = foreclosures(params, P, survival, hgrid, rgrid, mgrid,
                                      agrid, dgrid, ygrid, eprocess, repay_coeff, 
                                      Pricing_guess, Default, Renew, Policya, Puncond);

  double reneg;
  reneg = renegotiated(params, P, survival, hgrid, rgrid, mgrid,
                                      agrid, dgrid, ygrid, eprocess, repay_coeff, 
                                      Pricing_guess, Default, Renew, Policya, Puncond);

  cout << endl;
  cout << "Disbursements: " << disbu << "\t\tRevenues: " << rep + fore + reneg << " (" << rep << " + " << fore << " + " << reneg << ")" << endl;
  cout << endl;


	return(error_equilib);
}