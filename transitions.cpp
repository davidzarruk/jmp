

//====================================================================//
//         This selects the correct parameters each period            //
//====================================================================//
parameters parameter_selection(parameters params, double *mortgage_subsidy, double *incomeshock, double *repay_coeff, int *colorojo, int *tipo, const int ittrans, const int Ttrans,
                          const double rshock, const double Pashock, const double ppsishock, const double oomegashock,
                          const double sunk_shock, const double interm_shock, const double ltax_shock, const int periods_shock, const int periods_tax,
                          const double *qpath, const double *Ppath, const double *ddeltapath, 
                          const double *lumpsumpath, const double *mortsubsidy, const double *incshock, const double *survival,
                          const double q0, const double Ph0, const double ddeltabar0, const double r0, const double ppsi0, 
                          const double oomega0, const double sunk0, const double interm0, const double Pa0, const double lumpsum0, const double ltax0,
                          const int print_out){


  const int T        = params.T;
  const int na       = params.na;
  const int nm       = params.nm;
  const int nh       = params.nh;
  // const int nr       = params.nr;
  // const int nl       = params.nl;
  // const int nd       = params.nd;
  const int ny       = params.ny;
  const double rrho  = params.rrho;

  int ind;

  parameters params_out = params;

  // Last period is steady state
    if(ittrans == Ttrans+1){
      *tipo = 0;

      params_out.q                   = q0;
      params_out.Ph_today            = Ph0;
      params_out.Ph_tomorrow         = Ph0;
      params_out.ddeltabar_today     = ddeltabar0;
      params_out.ddeltabar_tomorrow  = ddeltabar0;
      params_out.r                   = r0;
      params_out.ppsi                = ppsi0;
      params_out.oomega              = oomega0;
      params_out.sunk                = sunk0;
      params_out.interm              = interm0;
      params_out.Pa                  = Pa0;
      params_out.lumpsum             = lumpsum0;
      params_out.ltax                = ltax0;

      if(print_out == 0){
        cout << "Estado estacionario final";
      } else{
        cout << "Estado estacionario final" << endl;
      }

    } else if(ittrans == 0){
      *tipo = 0;

      params_out.q                   = q0;
      params_out.Ph_today            = Ph0;
      params_out.Ph_tomorrow         = Ph0;
      params_out.ddeltabar_today     = ddeltabar0;
      params_out.ddeltabar_tomorrow  = ddeltabar0;
      params_out.r                   = r0;
      params_out.ppsi                = ppsi0;
      params_out.oomega              = oomega0;
      params_out.sunk                = sunk0;
      params_out.interm              = interm0;
      params_out.Pa                  = Pa0;
      params_out.lumpsum             = lumpsum0;
      params_out.ltax                = 0.0;

      if(print_out == 0){
        cout << "Estado estacionario final";
      } else{
        cout << "Estado estacionario final" << endl;
      }

    } else{
      *tipo = 1;

      // Guess of prices    
      params_out.q                  = qpath[ittrans-1];
      params_out.Ph_today           = Ppath[ittrans-1];
      params_out.ddeltabar_today    = ddeltapath[ittrans-1];
      params_out.lumpsum            = lumpsumpath[ittrans-1];

      if(ittrans == Ttrans){
        params_out.Ph_tomorrow        = Ph0;
        params_out.ddeltabar_tomorrow = ddeltabar0;
      } else{
        params_out.Ph_tomorrow        = Ppath[ittrans];
        params_out.ddeltabar_tomorrow = ddeltapath[ittrans];
      }

     

      if(print_out == 0){
        cout << "Iteracion: " << ittrans << " " << endl;
      } else{
        cout << "Iteracion: " << ittrans << endl;
      }
    }

    // Shock on interest rate
    if(ittrans <= periods_shock-1 && ittrans > 0){
      *colorojo = 1;
      params_out.r    = rshock;
      params_out.interm = interm_shock;
    } else{
      // *colorojo = 0;
      params_out.r  = r0;
      params_out.interm = interm0;
    }

      *colorojo = 1;

    // Other shocks
    if(ittrans <= periods_shock-1 && ittrans > 0){

      params_out.Pa     = Pashock;
      params_out.ppsi   = ppsishock;
      params_out.oomega = oomegashock;
      params_out.sunk   = sunk_shock;

    } else{

      params_out.Pa     = Pa0;
      params_out.ppsi   = ppsi0;
      params_out.oomega = oomega0;
      params_out.sunk   = sunk0;

    }

    if(ittrans <= periods_tax && ittrans > 0){

      params_out.ltax   = ltax_shock;

      // Subsidies only on first period
      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

                mortgage_subsidy[ind] = mortsubsidy[ind];
              }
            }
          }
        }
      }

    } else{

      // params_out.ltax   = ltax0;

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

                mortgage_subsidy[ind] = 0.0;
              }
            }
          }
        }
      }

    }

    
    if(ittrans <= periods_shock-1 && ittrans > 0){
      for(int i=0; i<T; i++){
        incomeshock[i] = incshock[i];
      }
    } else{
      for(int i=0; i<T; i++){
        incomeshock[i] = 0.0;
      }
    }

    repay_coefficient_shock(T, r0 + interm0, rshock + interm_shock, periods_shock - ittrans, rrho, survival, repay_coeff);

    // cout << params_out.q                     << endl;
    //   cout << params_out.Ph_today            << endl;
    //   cout << params_out.Ph_tomorrow         << endl;
    //   cout << params_out.ddeltabar_today     << endl;
    //   cout << params_out.ddeltabar_tomorrow  << endl;
    //   cout << params_out.r                   << endl;
    //   cout << params_out.ppsi                << endl;
    //   cout << params_out.oomega              << endl;
    //   cout << params_out.sunk                << endl;
    //   cout << params_out.interm              << endl;
    //   cout << params_out.Pa                  << endl;
    //   cout << params_out.lumpsum             << endl;
    //   cout << params_out.ltax                << endl;

    //   cout << params_out.Pa     << endl;
    //   cout << params_out.ppsi   << endl;
    //   cout << params_out.oomega << endl;
    //   cout << params_out.sunk   << endl;

    //   cout << params_out.r  << endl;
    //   cout << params_out.interm << endl;

    //   cout << params_out.q                  << endl;
    //   cout << params_out.Ph_today           << endl;
    //   cout << params_out.ddeltabar_today    << endl;
    //   cout << params_out.lumpsum            << endl;

    //   cout << params_out.Ph_tomorrow        << endl;
    //   cout << params_out.ddeltabar_tomorrow << endl;

    return(params_out);
}





//======================================
//         Transitional Dynamics:
//         Transitory shock
//======================================


double transition_error(parameters params, const int error_funct,
                      const int Ttrans, const double rshock, const double Pashock, const double ppsishock, const double oomegashock,
                      const double sunk_shock, const double interm_shock, const double ltax_shock, const int periods_shock, const int periods_tax,
                      const double *qpath, const double *Ppath, const double *ddeltapath, 
                      const double *lumpsumpath, const double *mortsubsidy, const int *subs_eligible, const int *subs_target, 
                      const double prob_mistake, const double *incshock, const int baseline,
                      const std::string tipo_run){
	
  int colorojo = 0;

	// VFI parameters
  const int T 		   = params.T;
  const int na  	   = params.na;
  const int nm  	   = params.nm;
  const int nh  	   = params.nh;
  const int nr       = params.nr;
  const int nl       = params.nl;
  const int nd  	   = params.nd;
  const int ny       = params.ny;
  const double rrho  = params.rrho;

  double q0           = params.q;
  double Ph0          = params.Ph_today;
  double ddeltabar0   = params.ddeltabar_today;
  double r0           = params.r;
  double ppsi0        = params.ppsi;
  double oomega0      = params.oomega;
  double sunk0        = params.sunk;
  double interm0      = params.interm;
  double Pa0          = params.Pa;
  double lumpsum0     = params.lumpsum;
  double ltax0        = params.ltax;

  int ind;
  size_t sizemats = ny*na*nm*nh*T*sizeof(double);
  double *mortgage_subsidy;
  mortgage_subsidy = (double*)malloc(sizemats);

  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int ih=0; ih<nh; ih++){
        for(int im=0; im<nm; im++){
          for(int ia=0; ia<na; ia++){
            ind = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

            mortgage_subsidy[ind] = 0.0;
          }
        }
      }
    }
  }


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

//==================================================================================================================================//
//======================                INITIALIZE      -       INITIZALIZE       -       INITIALIZE       =========================//
//==================================================================================================================================//

  // On the host memory
  double *Value, *Value_equiv, *Value_tomorrow, *Value_equiv_tomorrow, *Policyc, *Pricing_guess, *Pcond, *Puncond, *Pcond_yest, *Puncond_yest, *Pricing_tomorrow, *Puncond_tomorrow;
  double *Value_without, *Value_equiv_without, *Policyc_without, *Pricing_guess_without;

  int *Default, *Renew, *Policya, *Policym, *Policyh, *Policyr, *Policyl;
  int *Default_without, *Renew_without, *Policya_without, *Policym_without, *Policyh_without, *Policyr_without, *Policyl_without;

  double *Value_all, *Value_equiv_all, *Policyc_all, *Pricing_guess_all, *Pcond_all, *Puncond_all;
  double *Value_without_all, *Value_equiv_without_all, *Policyc_without_all, *Pricing_guess_without_all;
  int *Default_all, *Renew_all, *Policya_all, *Policym_all, *Policyh_all, *Policyr_all, *Policyl_all;
  int *Default_without_all, *Renew_without_all, *Policya_without_all, *Policym_without_all, *Policyh_without_all, *Policyr_without_all, *Policyl_without_all;

  int *elig_manana, *target_manana;
  
  sizemats                        = ny*na*nm*nh*nd*T*sizeof(double);
  size_t sizematsint              = ny*na*nm*nh*nd*T*sizeof(int);
  size_t sizemats_all             = ny*na*nm*nh*nd*T*(Ttrans+2)*sizeof(double);
  size_t sizemats_without_all     = ny*na*nm*nh*nd*T*periods_tax*sizeof(double);
  size_t sizematsint_all          = ny*na*nm*nh*nd*T*(Ttrans+2)*sizeof(int);
  size_t sizematsint_without_all  = ny*na*nm*nh*nd*T*periods_tax*sizeof(int);
  

  Value             = (double*)malloc(sizemats);
  Value_equiv       = (double*)malloc(sizemats);
  Value_tomorrow    = (double*)malloc(sizemats);
  Value_equiv_tomorrow    = (double*)malloc(sizemats);
  Policyc           = (double*)malloc(sizemats);
  Pricing_guess     = (double*)malloc(sizemats);

  Value_without             = (double*)malloc(sizemats);
  Value_equiv_without       = (double*)malloc(sizemats);
  Policyc_without           = (double*)malloc(sizemats);
  Pricing_guess_without     = (double*)malloc(sizemats);

  Pcond             = (double*)malloc(sizemats);
  Puncond           = (double*)malloc(sizemats);
  Pcond_yest        = (double*)malloc(sizemats);
  Puncond_yest      = (double*)malloc(sizemats);
  Puncond_tomorrow  = (double*)malloc(sizemats);
  Pricing_tomorrow  = (double*)malloc(sizemats);


  Default        = (int*)malloc(sizematsint);
  Renew          = (int*)malloc(sizematsint);
  Policya        = (int*)malloc(sizematsint);
  Policym        = (int*)malloc(sizematsint);
  Policyh        = (int*)malloc(sizematsint);
  Policyr        = (int*)malloc(sizematsint);
  Policyl        = (int*)malloc(sizematsint);

  elig_manana        = (int*)malloc(sizematsint);
  target_manana      = (int*)malloc(sizematsint);

  Default_without        = (int*)malloc(sizematsint);
  Renew_without          = (int*)malloc(sizematsint);
  Policya_without        = (int*)malloc(sizematsint);
  Policym_without        = (int*)malloc(sizematsint);
  Policyh_without        = (int*)malloc(sizematsint);
  Policyr_without        = (int*)malloc(sizematsint);
  Policyl_without        = (int*)malloc(sizematsint);

  Value_without_all         = (double*)malloc(sizemats_without_all);
  Value_equiv_without_all   = (double*)malloc(sizemats_without_all);
  Policyc_without_all       = (double*)malloc(sizemats_without_all);
  Pricing_guess_without_all = (double*)malloc(sizemats_without_all);

  Default_without_all        = (int*)malloc(sizematsint_without_all);
  Renew_without_all          = (int*)malloc(sizematsint_without_all);
  Policya_without_all        = (int*)malloc(sizematsint_without_all);
  Policym_without_all        = (int*)malloc(sizematsint_without_all);
  Policyh_without_all        = (int*)malloc(sizematsint_without_all);
  Policyr_without_all        = (int*)malloc(sizematsint_without_all);
  Policyl_without_all        = (int*)malloc(sizematsint_without_all);


  // Matrices with everything during transition
  Value_all         = (double*)malloc(sizemats_all);
  Value_equiv_all   = (double*)malloc(sizemats_all);
  Policyc_all       = (double*)malloc(sizemats_all);
  Pricing_guess_all = (double*)malloc(sizemats_all);
  Pcond_all         = (double*)malloc(sizemats_all);
  Puncond_all       = (double*)malloc(sizemats_all);
  
  Default_all       = (int*)malloc(sizematsint_all);
  Renew_all         = (int*)malloc(sizematsint_all);
  Policya_all       = (int*)malloc(sizematsint_all);
  Policym_all       = (int*)malloc(sizematsint_all);
  Policyh_all       = (int*)malloc(sizematsint_all);
  Policyr_all       = (int*)malloc(sizematsint_all);
  Policyl_all       = (int*)malloc(sizematsint_all);
  
  int ind2;


  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              Value[ind]          = 0.0;
              Value_equiv[ind]    = 0.0;
              Value_tomorrow[ind] = 0.0;
              Value_equiv_tomorrow[ind] = 0.0;
              Policyc[ind]        = 0.0;
              Pricing_guess[ind]            = 1.0;

              Value_without[ind]          = 0.0;
              Value_equiv_without[ind]          = 0.0;
              Policyc_without[ind]        = 0.0;
              Pricing_guess_without[ind]  = 0.0;

              Pcond[ind]          = 0.0;
              Puncond[ind]        = 0.0;
              Pcond_yest[ind]     = 0.0;
              Puncond_yest[ind]   = 0.0;
              Puncond_tomorrow[ind]   = 0.0;
              Pricing_tomorrow[ind]   = 0.0;

              Default[ind]        = 0;
              Renew[ind]          = 0;
              Policya[ind]        = 0;
              Policym[ind]        = 0;
              Policyh[ind]        = 0;
              Policyr[ind]        = 0;
              Policyl[ind]        = 0;

              elig_manana[ind]        = 0;
              target_manana[ind]      = 0;

              Default_without[ind]        = 0;
              Renew_without[ind]          = 0;
              Policya_without[ind]        = 0;
              Policym_without[ind]        = 0;
              Policyh_without[ind]        = 0;
              Policyr_without[ind]        = 0;
              Policyl_without[ind]        = 0;
            }
          }
        }
      }
    }
  }


  for(int ittrans=0; ittrans <= Ttrans+1; ittrans++){
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){

                ind = ittrans*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                ind2 = (ittrans-1)*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                
                Value_all[ind]         = 0.0;
                Value_equiv_all[ind]         = 0.0;
                Policyc_all[ind]       = 0.0;
                Pricing_guess_all[ind] = 1.0;
                Pcond_all[ind]         = 0.0;
                Puncond_all[ind]       = 0.0;

                Default_all[ind]       = 0;
                Renew_all[ind]         = 0;
                Policya_all[ind]       = 0;
                Policym_all[ind]       = 0;
                Policyh_all[ind]       = 0;
                Policyr_all[ind]       = 0;
                Policyl_all[ind]       = 0;

                if(ittrans <= periods_tax){

                  Value_without_all[ind2]         = 0.0;
                  Value_equiv_without_all[ind2]         = 0.0;
                  Policyc_without_all[ind2]       = 0.0;
                  Pricing_guess_without_all[ind2] = 1.0;

                  Default_without_all[ind2]       = 0;
                  Renew_without_all[ind2]         = 0;
                  Policya_without_all[ind2]       = 0;
                  Policym_without_all[ind2]       = 0;
                  Policyh_without_all[ind2]       = 0;
                  Policyr_without_all[ind2]       = 0;
                  Policyl_without_all[ind2]       = 0;
                }
              }
            }
          }
        }
      }
    }
  }

  double incomeshock[T];



//==================================================================================================================================//
//======================                COMPUTING THE OPTIMAL POLICIES - BACKWARDS INDUCTION               =========================//
//==================================================================================================================================//

  //---------------------//
  //---   Backwards  ----//
  //---------------------//

  // Empiezo en el estado estacionario final e itero hacia atras. 
  // Cada periodo, encuentro value y policies given next periods'

  int tipo;
  int print_out = 0;

  // Display el bank balance sheets
  int bank_balance_slow = 0;

  cout << endl;

  for(int ittrans=Ttrans+1; ittrans>=0; ittrans--){

    // Select correct parameters for transition period
    params = parameter_selection(params, mortgage_subsidy, incomeshock, repay_coeff, &colorojo, &tipo, ittrans, Ttrans,
                                rshock, Pashock, ppsishock, oomegashock, sunk_shock, interm_shock, ltax_shock, periods_shock, periods_tax, qpath, Ppath, ddeltapath, 
                                lumpsumpath, mortsubsidy, incshock, survival, q0, Ph0, ddeltabar0, r0, ppsi0, oomega0, sunk0, interm0, Pa0, lumpsum0, ltax0, print_out);


    
    if(tipo == 0){

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = ittrans*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  ind2 = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

                  Value[ind2]= 0.0;
                  Value_equiv[ind2]= 0.0;
                  Policyc[ind2]= 0;
                  Pricing_guess[ind2]= 0.0;

                  Default[ind2]= 0;
                  Renew[ind2]= 0;
                  Policya[ind2]= 0;
                  Policym[ind2]= 0;
                  Policyh[ind2]= 0;
                  Policyr[ind2]= 0;
                  Policyl[ind2]= 0;

                  Pricing_guess[ind2]     = 0.0;
                }
              }
            }
          }
        }
      }

      vfi_compute(params, error_funct,
                  P, survival, hgrid, rgrid, lgrid, mgrid,
                  agrid, dgrid, ygrid, eprocess, repay_coeff,
                  Value, Value_equiv, Policyc, Pricing_guess, 
                  Default, Renew, Policya, Policym, Policyh, Policyr, Policyl, print_out, incomeshock, mortgage_subsidy);
      cout << "VFI steady state " << endl;

    } else{

      vfi_transition_compute(params, error_funct, Value_tomorrow, Value_equiv_tomorrow, 
                            P, survival, hgrid, rgrid, lgrid, mgrid,
                            agrid, dgrid, ygrid, eprocess, repay_coeff,
                            Value, Value_equiv, Policyc, Pricing_guess, Pricing_tomorrow, 
                            Default, Renew, Policya, Policym, Policyh, Policyr, Policyl, print_out, incomeshock, mortgage_subsidy);
    }


    // Computo las policies without subsidies
    if(ittrans <= periods_tax && baseline == 0 && ittrans > 0){

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

                mortgage_subsidy[ind] = 0.0;
              }
            }
          }
        }
      }

      cout << " Without subsidies";

      vfi_transition_compute(params, error_funct, Value_tomorrow, Value_equiv_tomorrow, 
                            P, survival, hgrid, rgrid, lgrid, mgrid,
                            agrid, dgrid, ygrid, eprocess, repay_coeff,
                            Value_without, Value_equiv_without, Policyc_without, Pricing_guess_without, Pricing_tomorrow, 
                            Default_without, Renew_without, Policya_without, Policym_without, Policyh_without, 
                            Policyr_without, Policyl_without, print_out, incomeshock, mortgage_subsidy);

    }


    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = ittrans*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                ind2 = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

                Value_all[ind]          = Value[ind2];
                Value_equiv_all[ind]    = Value_equiv[ind2];
                Policyc_all[ind]        = Policyc[ind2];
                Pricing_guess_all[ind]  = Pricing_guess[ind2];

                Default_all[ind]        = Default[ind2];
                Renew_all[ind]          = Renew[ind2];
                Policya_all[ind]        = Policya[ind2];
                Policym_all[ind]        = Policym[ind2];
                Policyh_all[ind]        = Policyh[ind2];
                Policyr_all[ind]        = Policyr[ind2];
                Policyl_all[ind]        = Policyl[ind2];

                Value_tomorrow[ind2]    = Value[ind2];
                Value_equiv_tomorrow[ind2]    = Value_equiv[ind2];
                Pricing_tomorrow[ind2]  = Pricing_guess[ind2];
                Pricing_guess[ind2]     = 0.0;
              }
            }
          }
        }
      }
    }

    cout << endl;
  }

  
  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int ih=0; ih<nh; ih++){
        for(int im=0; im<nm; im++){
          for(int ia=0; ia<na; ia++){
            ind = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

            mortgage_subsidy[ind] = 0.0;
          }
        }
      }
    }
  }


  //---------------------//
  //---    Forward   ----//
  //---------------------//


  // Given las policies, encuentro distribucion y agregados en cada periodo

  double steady_distr[ny];

  steady_distribution(P, params, steady_distr);
  
  //--------------------------------------//
  //----   Distribution computation  -----//
  //--------------------------------------//

  double error_accum = 0.0;
  double *d_error_accum;
  d_error_accum = &error_accum;

  cout << endl;

  // params.q                   = q0;
  // params.Ph_today            = Ph0;
  // params.Ph_tomorrow         = Ph0;
  // params.ddeltabar_today     = ddeltabar0;
  // params.ddeltabar_tomorrow  = ddeltabar0;
  // params.ppsi                = ppsi0;
  // params.oomega              = oomega0;
  // params.sunk                = sunk0;
  // params.lumpsum             = lumpsum0;
  // params.ltax                = 0.0;

  // params.r                   = r0;
  // params.interm              = interm0;
  // params.Pa     = Pa0;

  // for(int i=0; i<T; i++){
  //   incomeshock[i] = 0.0;
  // }



  // repay_coefficient_shock(T, r0 + interm0, rshock + interm_shock, periods_shock - 1, rrho, survival, repay_coeff);

  // for(int i=0; i<T; i++){
  //   cout << repay_coeff[i] << endl;
  // }

  // Initial distribution from steady state
  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              Pcond[ind]        = 0.0;
              Puncond[ind]      = 0.0;
            }
          }
        }
      }
    }
  }

  probability_paths(Pcond, Puncond, P, params, 
                    survival, Policya, Policym, Policyh);

  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

              Pcond_all[ind]        = Pcond[ind];
              Puncond_all[ind]      = Puncond[ind];
            }
          }
        }
      }
    }
  }


  cout << "----------------------------------------------------------------------------------------------------" << endl;


  error_sequence(params, error_funct, 0, P, survival, hgrid, rgrid, mgrid, agrid, dgrid, ygrid, eprocess, lgrid,
                  Value, Value_equiv, Policyc, Pricing_guess, Default, Renew, Policya, Policyl, Policym, Policyh, Policyr, 
                  Value, Value_equiv, Policyc, Pricing_guess, Default, Renew, Policya, Policyl, Policym, Policyh, Policyr, 
                  Pcond, Puncond, mortgage_subsidy, repay_coeff, 
                  Ttrans, rshock, Pashock, ppsishock, oomegashock, incomeshock, sunk_shock, interm_shock, ltax_shock, periods_shock, qpath, Ppath,
                  d_error_accum, colorojo, baseline, subs_eligible, subs_target, prob_mistake, tipo_run);    


  // vfi_error(params, error_funct, 0,
  //           P, survival, hgrid, rgrid, mgrid,
  //           agrid, dgrid, ygrid, eprocess,
  //           Value, Policyc, Pricing_guess, 
  //           Default, Renew, Policya, Policym, Policyh, Policyr,
  //           Pcond, Puncond);

  
  //--------------------------------------//
  //--- Bank disimbursements/revenues ----//
  //--------------------------------------//
  
    //   cout << "ppsi: " << params.ppsi << endl;
    // cout << "r: " << params.r << endl;
    // cout << "oomega: " << params.oomega << endl;
    // cout << "sunk: " << params.sunk << endl;
    // cout << "Pa: " << params.Pa << endl;
    // cout << "Ph_today: " << params.Ph_today << endl;
    // cout << "Ph_tomorrow: " << params.Ph_tomorrow << endl;
    // cout << "dd_today: " << params.ddeltabar_today << endl;
    // cout << "dd_tomorrow: " << params.ddeltabar_tomorrow << endl;
    // cout << "q: " << params.q << endl;
    // cout << "interm: " << params.interm << endl;
    // cout << "lumpsum: " << params.lumpsum << endl;

    // for(int it = 0; it<T; it++){
    //   cout << repay_coeff[it] << endl;
    // }


  if(bank_balance_slow == 1){
    
    balance_sheet_trans(params, P, survival, hgrid, rgrid, mgrid,
                        agrid, dgrid, ygrid, eprocess, Puncond, Value_all, Pricing_guess_all, Default_all, Renew_all, 
                        Policya_all, Policym_all, Policyh_all, Policyr_all, 0, Ttrans,
                        rshock, Pashock, ppsishock, interm_shock, sunk_shock, periods_shock, qpath, Ppath, ddeltapath, r0, Ph0, ddeltabar0, interm0, sunk0);
  
    cout << "----------------------------------------------------------------------------------------------------" << endl;
  }

  pricing_test_all(params, P, survival, hgrid, rgrid, mgrid, agrid, dgrid, ygrid, eprocess, repay_coeff, 
                    Puncond, Pcond, Value_all, Pricing_guess_all, Default_all, Renew_all, Policya_all, 
                    Policym_all, Policyh_all, Policyr_all, 0, Ttrans, qpath, Ppath, ddeltapath, rshock, Ppath[0], ddeltapath[0],
                    prob_mistake, mortsubsidy, subs_eligible, subs_target);


  //--------------------------------------//
  //---         Fast Forward          ----//
  //--------------------------------------//

  // Ahora si itero hacia adelante
  for(int ittrans=1; ittrans<=Ttrans; ittrans++){ 

    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                ind2 = ittrans*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

                Pcond_all[ind2]        = Pcond[ind];
                Puncond_all[ind2]      = Puncond[ind];
              }
            }
          }
        }
      }
    }

    // Select correct parameters for transition period
    params = parameter_selection(params, mortgage_subsidy, incomeshock, repay_coeff, &colorojo, &tipo, ittrans, Ttrans,
                                rshock, Pashock, ppsishock, oomegashock, sunk_shock, interm_shock, ltax_shock, periods_shock, periods_tax, qpath, Ppath, ddeltapath, 
                                lumpsumpath, mortsubsidy, incshock, survival, q0, Ph0, ddeltabar0, r0, ppsi0, oomega0, sunk0, interm0, Pa0, lumpsum0, ltax0, print_out);


    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                ind2 = ittrans*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

                Value[ind]             = Value_all[ind2];
                Value_equiv[ind]       = Value_equiv_all[ind2];
                Policyc[ind]           = Policyc_all[ind2];
                Pricing_guess[ind]     = Pricing_guess_all[ind2];

                Default[ind]           = Default_all[ind2];
                Renew[ind]             = Renew_all[ind2];
                Policya[ind]           = Policya_all[ind2];
                Policym[ind]           = Policym_all[ind2];
                Policyh[ind]           = Policyh_all[ind2];
                Policyr[ind]           = Policyr_all[ind2];
                Policyl[ind]           = Policyl_all[ind2];
              }
            }
          }
        }
      }
    }


  // Error en EE inicial
  if(baseline == 0 && ittrans <= periods_tax){

    // for(int it=0; it<T; it++){
    //   for(int iy=0; iy<ny; iy++){
    //     for(int id=0; id<nd; id++){
    //       for(int ih=0; ih<nh; ih++){
    //         for(int im=0; im<nm; im++){
    //           for(int ia=0; ia<na; ia++){
    //             ind = (ittrans-1)*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
    //             ind2 = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

    //             Value_without[ind2]         = Value_without_all[ind];
    //             Policyc_without[ind2]       = Policyc_without_all[ind];
    //             Pricing_guess_without[ind2] = Pricing_guess_without_all[ind];
                
    //             Default_without[ind2]      = Default_without_all[ind];
    //             Renew_without[ind2]        = Renew_without_all[ind];
    //             Policya_without[ind2]      = Policya_without_all[ind];
    //             Policym_without[ind2]      = Policym_without_all[ind];
    //             Policyh_without[ind2]      = Policyh_without_all[ind];
    //             Policyr_without[ind2]      = Policyr_without_all[ind];
    //             Policyl_without[ind2]      = Policyl_without_all[ind];
    //           }
    //         }
    //       }
    //     }
    //   }
    // }


    if(ittrans == 1){
      cout << "Ver default de estos parceros subsidiados: " << endl;
  
      int ihp;
      int imp;
      int iap;
  
      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
  
                  ihp = Policyh[ind];
                  imp = Policym[ind];
                  iap = Policya[ind];
                  
                  for(int idp=0; idp<nd; idp++){
                    for(int iyp=0; iyp<ny; iyp++){
                      ind2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
  
                      if(subs_eligible[ind] == 1){
                        elig_manana[ind2] == 1;
                      }
                      if(subs_eligible[ind] == 1 && subs_target[ind] == 1){
                        target_manana[ind2] == 1;
                      }
                    }
                  }
  
                }
              }
            }
          }
        }
      }
  
    }

    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                ind2 = ittrans*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

                Value[ind]             = Value_all[ind2];
                Value_equiv[ind]       = Value_equiv_all[ind2];
                Policyc[ind]           = Policyc_all[ind2];
                Pricing_guess[ind]     = Pricing_guess_all[ind2];

                Default[ind]           = Default_all[ind2];
                Renew[ind]             = Renew_all[ind2];
                Policya[ind]           = Policya_all[ind2];
                Policym[ind]           = Policym_all[ind2];
                Policyh[ind]           = Policyh_all[ind2];
                Policyr[ind]           = Policyr_all[ind2];
                Policyl[ind]           = Policyl_all[ind2];
              }
            }
          }
        }
      }
    }


    error_sequence(params, error_funct, ittrans, P, survival, hgrid, rgrid, mgrid, agrid, dgrid, ygrid, eprocess, lgrid,
                    Value, Value_equiv, Policyc, Pricing_guess, Default, Renew, Policya, Policyl, Policym, Policyh, Policyr, 
                    Value_without, Value_equiv_without, Policyc_without, Pricing_guess_without, Default_without, Renew_without, Policya_without, Policyl_without, Policym_without, Policyh_without, Policyr_without, 
                    Pcond, Puncond, mortgage_subsidy, repay_coeff, 
                    Ttrans, rshock, Pashock, ppsishock, oomegashock, incomeshock, sunk_shock, interm_shock, ltax_shock, periods_shock, qpath, Ppath,
                    d_error_accum, colorojo, baseline, subs_eligible, subs_target, prob_mistake, tipo_run);
  } else{
    error_sequence(params, error_funct, ittrans, P, survival, hgrid, rgrid, mgrid, agrid, dgrid, ygrid, eprocess, lgrid,
                    Value, Value_equiv, Policyc, Pricing_guess, Default, Renew, Policya, Policyl, Policym, Policyh, Policyr, 
                    Value, Value_equiv, Policyc, Pricing_guess, Default, Renew, Policya, Policyl, Policym, Policyh, Policyr, 
                    Pcond, Puncond, mortgage_subsidy, repay_coeff, 
                    Ttrans, rshock, Pashock, ppsishock, oomegashock, incomeshock, sunk_shock, interm_shock, ltax_shock, periods_shock, qpath, Ppath,
                    d_error_accum, colorojo, baseline, subs_eligible, subs_target, prob_mistake, tipo_run);
  }

    double defo = 0.0;
    double mean_m = 0.0;
    double pipolseria = 0.0;

    if(ittrans == 2){

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  
                  if(im>0 && elig_manana[ind] == 1){
                    defo   = defo + Puncond[ind]*(double)Default[ind];
                    mean_m = mean_m + Puncond[ind]*mgrid[im];
                    pipolseria = pipolseria + Puncond[ind];
                  }
                }
              }
            }
          }
        }
      }

      cout << "Y el default de los parceros manana es = " << mean_m << endl;
    

    }



    // for(int it = 0; it<T; it++){
    //   cout << repay_coeff[it] << endl;
    // }
    //     cout << "ppsi: " << params.ppsi << endl;
    // cout << "r: " << params.r << endl;
    // cout << "oomega: " << params.oomega << endl;
    // cout << "sunk: " << params.sunk << endl;
    // cout << "Pa: " << params.Pa << endl;
    // cout << "Ph_today: " << params.Ph_today << endl;
    // cout << "Ph_tomorrow: " << params.Ph_tomorrow << endl;
    // cout << "dd_today: " << params.ddeltabar_today << endl;
    // cout << "dd_tomorrow: " << params.ddeltabar_tomorrow << endl;
    // cout << "q: " << params.q << endl;
    // cout << "interm: " << params.interm << endl;
    // cout << "lumpsum: " << params.lumpsum << endl;

    if(bank_balance_slow == 1){
      balance_sheet_trans(params, P, survival, hgrid, rgrid, mgrid,
                          agrid, dgrid, ygrid, eprocess, Puncond, Value_all, Pricing_guess_all, Default_all, Renew_all, 
                          Policya_all, Policym_all, Policyh_all, Policyr_all, ittrans, Ttrans,
                          rshock, Pashock, ppsishock, interm_shock, sunk_shock, periods_shock, qpath, Ppath, ddeltapath, r0, Ph0, ddeltabar0, interm0, sunk0);
  
      cout << "----------------------------------------------------------------------------------------------------" << endl;
    }

    pricing_test_all(params, P, survival, hgrid, rgrid, mgrid, agrid, dgrid, ygrid, eprocess, repay_coeff, 
                      Puncond, Pcond, Value_all, Pricing_guess_all, Default_all, Renew_all, Policya_all, 
                      Policym_all, Policyh_all, Policyr_all, ittrans, Ttrans, qpath, Ppath, ddeltapath, params.r, params.Ph_tomorrow, params.ddeltabar_tomorrow,
                      prob_mistake, mortsubsidy, subs_eligible, subs_target);


    // Evolution of distribution
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

                Pcond_yest[ind]   = Pcond[ind];
                Puncond_yest[ind] = Puncond[ind];

                Pcond[ind]   = 0.0;
                Puncond[ind] = 0.0;
              }
            }
          }
        }
      }
    }
    
    // Update distribution of today, given yesterday's and policies
    if(baseline == 0 && ittrans <= periods_tax){
      probability_paths_trans(Pcond_yest, Puncond_yest, Pcond, Puncond, P, params, 
                              survival, Policya, Policym, Policyh, Policya_without, 
                              Policym_without, Policyh_without, prob_mistake, mortsubsidy, subs_eligible, subs_target);
    } else{
      probability_paths_trans(Pcond_yest, Puncond_yest, Pcond, Puncond, P, params, 
                              survival, Policya, Policym, Policyh, Policya, 
                              Policym, Policyh, prob_mistake, mortsubsidy, subs_eligible, subs_target);
    }


  }

  Color::Modifier green(Color::FG_GREEN);
  Color::Modifier def(Color::FG_DEFAULT);

  cout << endl;
  cout << green << "\e[1m Error: " << error_accum << "\e[0m" << def << endl;

  return(error_accum);
}
