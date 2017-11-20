

//======================================
//         Value function
//======================================

void vfi_compute(parameters params, const int error_funct, 
                  const double *P, const double *survival, const double *hgrid, const double *rgrid, const double *lgrid, const double *mgrid,
                  const double *agrid, const double *dgrid, const double *ygrid, const double *eprocess, const double *repay_coeff,
                  double *Value, double *Value_equiv, double *Policyc, double *Pricing_guess, 
                  int *Default, int *Renew, int *Policya, int *Policym, int *Policyh, int *Policyr, int *Policyl,
                  const int print_out, const double *incshock, const double *mortsubsidy){

  // VFI parameters
  const int maxiter         = params.maxiter;
  const int uti             = params.uti;
  const double tol          = params.tol;
  const double convergence  = params.convergence;
  const int T               = params.T;
  const int Tretirement     = params.Tretirement;

  // Grid sizes
  const int na      = params.na;
  const int nm      = params.nm;
  const int nh      = params.nh;
  const int nr      = params.nr;
  const int nl      = params.nl;
  const int nd      = params.nd;
  const int ny      = params.ny;

  // Preferences
  const double ssigma   = params.ssigma;
  const double rrho     = params.rrho;
  const double ppsi     = params.ppsi;
  const double bbeta    = params.bbeta;
  const double kkappa   = params.kkappa;

  const double tthetalab = params.tthetalab;
  const double eetalab  = params.eetalab;

  const double sunk       = params.sunk;
  const double interm     = params.interm;
  const double rec_probab = params.rec_probab;
  const double oomega     = params.oomega;

  // Equilibrium objects
  const double ddeltabar = params.ddeltabar_today;
  const double r         = params.r;
  const double Ph        = params.Ph_today;
  const double q         = params.q;
  const double Pa        = params.Pa;
  const double fcost     = params.fcost;
  const double refcost   = params.refcost;
  const double pension   = params.pension;
  const double sstax     = params.sstax;
  const double ltax      = params.ltax;
  const double lumpsum   = params.lumpsum;

  const int equivalent   = params.compute_equivalent;

  const double multiplier = params.multiplier;

  int ind;

  ind = 21*ny*nh*nm*na + 1*nh*nm*na + 1*nm*na + 1*na + 10;
  cout << "mortsubsidy = " << mortsubsidy[ind] << endl;
 
//----------------------------------------------//
//--------    VALUES AND POLICIES     ----------//
//----------------------------------------------//
  
  double *Pricing, *Pricing_guess2;
  
  size_t sizemats = ny*na*nm*nh*nd*T*sizeof(double);
  size_t sizematssubs = ny*na*nm*nh*T*sizeof(double);

  Pricing         = (double*)malloc(sizemats);
  Pricing_guess2  = (double*)malloc(sizemats);

  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              Pricing[ind]        = 0.0;
              Pricing_guess2[ind] = 0.0;
            }
          }
        }
      }
    }
  }

  // for(int i=0; i<T; i++){
  //   cout << incshock[i] << endl;
  // }

  size_t sizematsint  = ny*na*nm*nh*nd*T*sizeof(int);
  size_t sizeincshock = T*sizeof(double);

  // On the device memory
  double *dValue, *dValue_equiv, *dPolicyc, *dPricing_guess, *dPricing, *dPricing_guess2, *dincshock, *dmortsubsidy;
  int *dDefault, *dRenew, *dPolicya, *dPolicym, *dPolicyh, *dPolicyr, *dPolicyl;

  cudaMalloc((void**)&dValue,          sizemats);
  cudaMalloc((void**)&dValue_equiv,    sizemats);
  cudaMalloc((void**)&dPolicyc,        sizemats);
  cudaMalloc((void**)&dPricing_guess,  sizemats);
  cudaMalloc((void**)&dPricing,        sizemats);
  cudaMalloc((void**)&dPricing_guess2, sizemats);

  cudaMalloc((void**)&dmortsubsidy,    sizematssubs);
  cudaMalloc((void**)&dincshock,       sizeincshock);
  
  cudaMalloc((void**)&dDefault,        sizematsint);
  cudaMalloc((void**)&dRenew,          sizematsint);
  cudaMalloc((void**)&dPolicya,        sizematsint);
  cudaMalloc((void**)&dPolicym,        sizematsint);
  cudaMalloc((void**)&dPolicyh,        sizematsint);
  cudaMalloc((void**)&dPolicyr,        sizematsint);
  cudaMalloc((void**)&dPolicyl,        sizematsint);

  cudaMemcpy(dValue,        Value,        sizemats,    cudaMemcpyHostToDevice);
  cudaMemcpy(dValue_equiv,  Value_equiv,  sizemats,    cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicyc,      Policyc,      sizemats,    cudaMemcpyHostToDevice);

  cudaMemcpy(dmortsubsidy,  mortsubsidy,  sizematssubs,    cudaMemcpyHostToDevice);
  cudaMemcpy(dincshock,     incshock,     sizeincshock,cudaMemcpyHostToDevice);

  cudaMemcpy(dDefault,      Default,      sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dRenew,        Renew,        sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicya,      Policya,      sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicym,      Policym,      sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicyh,      Policyh,      sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicyr,      Policyr,      sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicyl,      Policyl,      sizematsint, cudaMemcpyHostToDevice);


  // Grids on device
  double *dagrid, *dmgrid, *dhgrid, *drgrid, *dlgrid, *ddgrid, *dygrid, *dP, *deprocess, *dsurvival, *drepay_coeff;

  size_t sizeAgrid      = na*sizeof(double);
  size_t sizeMgrid      = nm*sizeof(double);
  size_t sizeHgrid      = nh*sizeof(double);
  size_t sizeRgrid      = nr*sizeof(double);
  size_t sizeLgrid      = nl*sizeof(double);
  size_t sizeDgrid      = nd*sizeof(double);
  size_t sizeYgrid      = ny*sizeof(double);
  size_t sizeP          = ny*ny*sizeof(double);
  size_t sizeEprocess   = T*sizeof(double);
  size_t sizeSurvival   = T*sizeof(double);
  size_t sizeRepay      = T*sizeof(double);

  cudaMalloc((void**)&dagrid,       sizeAgrid);
  cudaMalloc((void**)&dmgrid,       sizeMgrid);
  cudaMalloc((void**)&dhgrid,       sizeHgrid);
  cudaMalloc((void**)&drgrid,       sizeRgrid);
  cudaMalloc((void**)&dlgrid,       sizeLgrid);
  cudaMalloc((void**)&ddgrid,       sizeDgrid);
  cudaMalloc((void**)&dygrid,       sizeYgrid);
  cudaMalloc((void**)&dP,           sizeP);
  cudaMalloc((void**)&deprocess,    sizeEprocess);
  cudaMalloc((void**)&dsurvival,    sizeSurvival);
  cudaMalloc((void**)&drepay_coeff, sizeRepay);

  // Transfer
  cudaMemcpy(dagrid,   agrid,   sizeAgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(dmgrid,   mgrid,   sizeMgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(dhgrid,   hgrid,   sizeHgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(drgrid,   rgrid,   sizeRgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(dlgrid,   lgrid,   sizeLgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(ddgrid,   dgrid,   sizeDgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(dygrid,   ygrid,   sizeYgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(dP,       P,       sizeP, cudaMemcpyHostToDevice);
  cudaMemcpy(deprocess,eprocess,   sizeEprocess, cudaMemcpyHostToDevice);
  cudaMemcpy(dsurvival,survival,   sizeSurvival, cudaMemcpyHostToDevice);
  cudaMemcpy(drepay_coeff,repay_coeff,   sizeRepay, cudaMemcpyHostToDevice);


//----------------------------------------------//
//--------   DIMENSIONS OF PARALLEL   ----------//
//----------------------------------------------//

  int xpar = nd;
  int ypar = nh;
  int zpar = ny;

  if(xpar*ypar*zpar > 512){
    cout << "-----------------------------" << endl;
    cout << "--- DIMENSIONES ESTAN MAL ---" << endl;
    cout << "-----------------------------" << endl;
  }

  // Hacia derecha va a ser H, hacia abajo va a ser M
  dim3 dimBlock(nd, nh, ny);
  dim3 dimGrid(nm, na);


//----------------------------------------------//
//--------  VALUE FUNCTION ITERATION  ----------//
//----------------------------------------------//

  // VFI starts
  double err   = 1;
  int itnum   = 1;

  double errprice = 1;
  int itnumprice = 1;

  clock_t t;
  clock_t t0;
  t0  = clock();
  t   = t0;

//----------------------------------------------//
//--------   COMPUTATION OF PRICES    ----------//
//----------------------------------------------//

  // Initial Pricing function
  while(errprice > 0.000001){

    cudaMemcpy(dPricing, Pricing, sizemats, cudaMemcpyHostToDevice);
    cudaMemcpy(dPricing_guess, Pricing_guess, sizemats, cudaMemcpyHostToDevice);

    for(int it=T-1; it>=0; it--){
      Pmort<<<dimGrid,dimBlock>>>(T, na, nm, nh, nd, ny, rrho, r, Ph, ddeltabar, sunk, interm, rec_probab, dP, ddgrid, 
                                  dmgrid, dhgrid, drgrid, dsurvival, drepay_coeff, it, dDefault, dRenew, dPolicya, dPricing, dPricing_guess);
      cudaDeviceSynchronize();      
    }

    cudaMemcpy(Pricing,   dPricing,   sizemats, cudaMemcpyDeviceToHost);
    cudaMemcpy(Pricing_guess, dPricing_guess, sizemats, cudaMemcpyDeviceToHost);

    errprice = 0.0;

    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                errprice = errprice + absol(Pricing_guess[ind] - Pricing[ind]);
                Pricing_guess[ind] = Pricing[ind];
                Pricing[ind] = 0.0;
              }
            }
          }
        }
      }
    }
    
    // if (itnumprice % 10 == 0){
    //   cout << "     Pricing Function - It. #: " << itnumprice << ". Error " << errprice << endl;
    // }
    itnumprice++;
  }


  // cout << "     -------------------------------------------" << endl;
  // cout << "     Pricing Function - It. #: " << itnumprice << ". Error " << errprice << endl;
  // cout << "     -------------------------------------------" << endl;

  double conv_prov;

  while((err > tol) & (itnum < maxiter)){
   
  //--------------------------------------//
  //--------  Policy functions  ----------//
  //--------------------------------------//
 
    // Copy back from host (CPU) to device (GPU)
    cudaMemcpy(dPricing_guess, Pricing_guess, sizemats, cudaMemcpyHostToDevice);

    // State variables
    for(int it=T-1; it>=0; it--){
      if(print_out == 1){
        if(it == (T-1)){
        cout << it+1 << ", " << flush;
        } else{
        cout << it+1 << ", " << flush;
        }
      } else{
        if(it % 10 == 0){
          cout << "." << flush;
        }
      }
     
      vfi<<<dimGrid,dimBlock>>>(T, Tretirement, na, nm, nh, nr, nl, nd, ny, uti, rrho, bbeta, Ph, q, Pa, ddeltabar, ssigma, ppsi, kkappa, 
                                tthetalab, eetalab, fcost, refcost, pension, sstax, ltax, lumpsum, oomega, rec_probab, sunk, dincshock, dmortsubsidy,
                                dagrid, dmgrid, dhgrid, drgrid, dlgrid, ddgrid, dygrid, dP, deprocess, dsurvival, drepay_coeff, 
                                it, equivalent, multiplier, dValue, dValue_equiv, dDefault, dRenew, dPolicya, dPolicym, dPolicyh, dPolicyr, dPolicyl, dPolicyc, dPricing_guess); 

      cudaDeviceSynchronize();
    }
    

  //--------------------------------------//
  //--------  Pricing function  ----------//
  //--------------------------------------//

    errprice = 1.0;
    itnumprice = 1;

    while(errprice > 0.000001){

      cudaMemcpy(dPricing, Pricing, sizemats, cudaMemcpyHostToDevice);
      cudaMemcpy(dPricing_guess2, Pricing_guess2, sizemats, cudaMemcpyHostToDevice);

      for(int it=T-1; it>=0; it--){
        Pmort<<<dimGrid,dimBlock>>>(T, na, nm, nh, nd, ny, rrho, r, Ph, ddeltabar, sunk, interm, rec_probab, dP, 
                                    ddgrid, dmgrid, dhgrid, drgrid, dsurvival, drepay_coeff, it, dDefault, dRenew, dPolicya, dPricing, dPricing_guess2);
        cudaDeviceSynchronize();      
      }

      cudaMemcpy(Pricing,   dPricing,   sizemats, cudaMemcpyDeviceToHost);
      cudaMemcpy(Pricing_guess2, dPricing_guess2, sizemats, cudaMemcpyDeviceToHost);

      errprice = 0.0;
      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  errprice = errprice + absol(Pricing_guess2[ind] - Pricing[ind]);
                  Pricing_guess2[ind] = Pricing[ind];
                  Pricing[ind] = 0.0;
                }
              }
            }
          }
        }
      }
      
      if (itnumprice % 10 == 0){
        // cout << "     Pricing Function - It. #: " << itnumprice << ". Error " << errprice << endl;
      }
      itnumprice++;
    }


    // cout << "     -------------------------------------------" << endl;
    // cout << "     Pricing Function - It. #: " << itnumprice << ". Error " << errprice << endl;
    // cout << "     -------------------------------------------" << endl;

    if(convergence != 0.0){ 
      conv_prov = 0.0;
      if (err <= 0.01 & err > 0.001){
        conv_prov = 0.5;
      } else if(err <= 0.001 & err > 0.0002){
        conv_prov = 0.7;
      } else if(err <= 0.002 & err > 0.0001){
        conv_prov = 0.9;
      } else if(err <= 0.0001 & err > 0.000075){
        conv_prov = 0.9;
      } else if(err <= 0.000075 & err > 0.00005){
        conv_prov = 0.9;
      } else if(err <= 0.00005){
        conv_prov = 0.9;
      }
    }
    
    // Compare new and old pricing functions 
    err = 0.0;
    for(int it=T-1; it>=0; it--){
      for(int iy=0; iy<ny; iy++){
        for(int ia=0; ia<na; ia++){
          for(int im=0; im<nm; im++){
            for(int ih=0; ih<nh; ih++){
              for(int id=0; id<nd; id++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                
                err = err + absol(Pricing_guess2[ind] - Pricing_guess[ind])/sizemats;

              }
            }
          }
        }
      }
    }

    
    t = clock() - t0;

    if (itnum % 1 == 0 && print_out == 1){
      // cout << "     -------------------------------------------" << endl;
      // cout << "FUNCTION: Iteration number: " << itnum << ". Time: " << 1000000*((double)t)/CLOCKS_PER_SEC << " seconds. Error " << err << endl;
      cout << "  It: " << itnum << ". Time: " << ((double)t/60)/CLOCKS_PER_SEC << " minutes. Error " << err << ", conv= " << conv_prov << endl;
      // cout << "     -------------------------------------------" << endl;
    }

    itnum = itnum + 1;

    for(int it=T-1; it>=0; it--){
      for(int iy=0; iy<ny; iy++){
        for(int ia=0; ia<na; ia++){
          for(int im=0; im<nm; im++){
            for(int ih=0; ih<nh; ih++){
              for(int id=0; id<nd; id++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

                if((err > tol) & (itnum < maxiter)){
                  Pricing_guess[ind] = conv_prov*Pricing_guess[ind] + (1-conv_prov)*Pricing_guess2[ind];
                  Pricing_guess2[ind] = 0.0;
                } 
              }
            }
          }
        }
      }
    }
  }

  // Copy back from device (GPU) to host (CPU)
  cudaMemcpy(Value,         dValue,         sizemats,    cudaMemcpyDeviceToHost);
  cudaMemcpy(Value_equiv,   dValue_equiv,         sizemats,    cudaMemcpyDeviceToHost);
  cudaMemcpy(Policyc,       dPolicyc,       sizemats,    cudaMemcpyDeviceToHost);
  cudaMemcpy(Pricing_guess, dPricing_guess, sizemats,    cudaMemcpyDeviceToHost);
  cudaMemcpy(Default,       dDefault,       sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Renew,         dRenew,         sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Policya,       dPolicya,       sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Policym,       dPolicym,       sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Policyh,       dPolicyh,       sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Policyr,       dPolicyr,       sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Policyl,       dPolicyl,       sizematsint, cudaMemcpyDeviceToHost);

  // Free variables in device memory
  cudaFree(dValue);
  cudaFree(dValue_equiv);
  cudaFree(dPolicyc);
  cudaFree(dPricing_guess);
  cudaFree(dDefault);
  cudaFree(dRenew);
  cudaFree(dPolicya);
  cudaFree(dPolicym);
  cudaFree(dPolicyh);
  cudaFree(dPolicyr);
  cudaFree(dPolicyl);
  cudaFree(dPricing);
  cudaFree(dPricing_guess2);

  cudaFree(dincshock);
  cudaFree(dmortsubsidy);

  cudaFree(dagrid);
  cudaFree(dmgrid);
  cudaFree(dhgrid);
  cudaFree(drgrid);
  cudaFree(dlgrid);
  cudaFree(ddgrid);
  cudaFree(dygrid);
  cudaFree(dP);
  cudaFree(deprocess);
  cudaFree(dsurvival);
  cudaFree(drepay_coeff);

}



//======================================
//         Value function
//======================================

void vfi_transition_compute(parameters params, const int error_funct, const double *Value_future, const double *Value_equiv_future, 
                            const double *P, const double *survival, const double *hgrid, const double *rgrid, const double *lgrid, const double *mgrid,
                            const double *agrid, const double *dgrid, const double *ygrid, const double *eprocess, const double *repay_coeff,
                            double *Value, double *Value_equiv, double *Policyc, double *Pricing_guess, const double *Pricing_tomorrow, 
                            int *Default, int *Renew, int *Policya, int *Policym, int *Policyh, int *Policyr, int *Policyl,
                            const int print_out, const double *incshock, const double *mortsubsidy){

  // VFI parameters
  const int maxiter         = params.maxiter;
  const int uti             = params.uti;
  const double tol          = params.tol;
  const double convergence  = params.convergence;
  const int T               = params.T;
  const int Tretirement     = params.Tretirement;

  // Grid sizes
  const int na = params.na;
  const int nm = params.nm;
  const int nh = params.nh;
  const int nr = params.nr;
  const int nl = params.nl;
  const int nd = params.nd;
  const int ny = params.ny;

  // Preferences
  const double ssigma   = params.ssigma;
  const double rrho     = params.rrho;
  const double ppsi     = params.ppsi;
  const double bbeta    = params.bbeta;
  const double kkappa   = params.kkappa;
  const double tthetalab = params.tthetalab;
  const double eetalab  = params.eetalab;

  const double sunk       = params.sunk;
  const double interm     = params.interm;
  const double rec_probab = params.rec_probab;
  const double oomega     = params.oomega;

  // Equilibrium objects
  const double ddeltabar_today     = params.ddeltabar_today;
  const double ddeltabar_tomorrow  = params.ddeltabar_tomorrow;
  const double r                   = params.r;
  const double Ph_today            = params.Ph_today;
  const double Ph_tomorrow         = params.Ph_tomorrow;
  const double q                   = params.q;
  const double Pa                  = params.Pa;
  const double fcost               = params.fcost;
  const double refcost             = params.refcost;
  const double pension             = params.pension;
  const double sstax               = params.sstax;
  const double ltax                = params.ltax;
  const double lumpsum             = params.lumpsum;

  const int equivalent             = params.compute_equivalent;

  const double multiplier          = params.multiplier;

  int ind;

  ind = 21*ny*nh*nm*na + 1*nh*nm*na + 1*nm*na + 1*na + 10;
  cout << "mortsubsidy = " << mortsubsidy[ind] << endl;
  
//----------------------------------------------//
//--------    VALUES AND POLICIES     ----------//
//----------------------------------------------//
  
  double *Pricing_guess2;
  
  size_t sizemats     = ny*na*nm*nh*nd*T*sizeof(double);
  size_t sizematssubs = ny*na*nm*nh*T*sizeof(double);

  Pricing_guess2      = (double*)malloc(sizemats);

  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              Pricing_guess2[ind]       = 0.0;
            }
          }
        }
      }
    }
  }


  size_t sizematsint  = ny*na*nm*nh*nd*T*sizeof(int);
  size_t sizeincshock = T*sizeof(double);

  // On the device memory
  double *dValue, *dValue_equiv, *dValue_future, *dValue_equiv_future, *dPolicyc, *dPricing_tomorrow, *dPricing_guess, *dPricing_guess2, *dincshock, *dmortsubsidy;
  int *dDefault, *dRenew, *dPolicya, *dPolicym, *dPolicyh, *dPolicyr, *dPolicyl;

  cudaMalloc((void**)&dValue,            sizemats);
  cudaMalloc((void**)&dValue_equiv,      sizemats);
  cudaMalloc((void**)&dValue_future,     sizemats);
  cudaMalloc((void**)&dValue_equiv_future,     sizemats);
  cudaMalloc((void**)&dPolicyc,          sizemats);
  cudaMalloc((void**)&dPricing_tomorrow, sizemats);
  cudaMalloc((void**)&dPricing_guess,    sizemats);
  cudaMalloc((void**)&dPricing_guess2,   sizemats);

  cudaMalloc((void**)&dmortsubsidy,      sizematssubs);
  cudaMalloc((void**)&dincshock,         sizeincshock);

  cudaMalloc((void**)&dDefault,          sizematsint);
  cudaMalloc((void**)&dRenew,            sizematsint);
  cudaMalloc((void**)&dPolicya,          sizematsint);
  cudaMalloc((void**)&dPolicym,          sizematsint);
  cudaMalloc((void**)&dPolicyh,          sizematsint);
  cudaMalloc((void**)&dPolicyr,          sizematsint);
  cudaMalloc((void**)&dPolicyl,          sizematsint);

  cudaMemcpy(dValue,        Value,        sizemats,    cudaMemcpyHostToDevice);
  cudaMemcpy(dValue_equiv,  Value_equiv,  sizemats,    cudaMemcpyHostToDevice);
  cudaMemcpy(dValue_future, Value_future, sizemats,    cudaMemcpyHostToDevice);
  cudaMemcpy(dValue_equiv_future, Value_equiv_future, sizemats,    cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicyc,      Policyc,      sizemats,    cudaMemcpyHostToDevice);

  cudaMemcpy(dmortsubsidy,  mortsubsidy,  sizematssubs,    cudaMemcpyHostToDevice);
  cudaMemcpy(dincshock,     incshock,     sizeincshock,cudaMemcpyHostToDevice);

  cudaMemcpy(dDefault,      Default,      sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dRenew,        Renew,        sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicya,      Policya,      sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicym,      Policym,      sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicyh,      Policyh,      sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicyr,      Policyr,      sizematsint, cudaMemcpyHostToDevice);
  cudaMemcpy(dPolicyl,      Policyl,      sizematsint, cudaMemcpyHostToDevice);


  // Grids on device
  double *dagrid, *dmgrid, *dhgrid, *drgrid, *dlgrid, *ddgrid, *dygrid, *dP, *deprocess, *dsurvival, *drepay_coeff;

  size_t sizeAgrid      = na*sizeof(double);
  size_t sizeMgrid      = nm*sizeof(double);
  size_t sizeHgrid      = nh*sizeof(double);
  size_t sizeRgrid      = nr*sizeof(double);
  size_t sizeLgrid      = nl*sizeof(double);
  size_t sizeDgrid      = nd*sizeof(double);
  size_t sizeYgrid      = ny*sizeof(double);
  size_t sizeP          = ny*ny*sizeof(double);
  size_t sizeEprocess   = T*sizeof(double);
  size_t sizeSurvival   = T*sizeof(double);
  size_t sizeRepay      = T*sizeof(double);

  cudaMalloc((void**)&dagrid, sizeAgrid);
  cudaMalloc((void**)&dmgrid, sizeMgrid);
  cudaMalloc((void**)&dhgrid, sizeHgrid);
  cudaMalloc((void**)&drgrid, sizeRgrid);
  cudaMalloc((void**)&dlgrid, sizeLgrid);
  cudaMalloc((void**)&ddgrid, sizeDgrid);
  cudaMalloc((void**)&dygrid, sizeYgrid);
  cudaMalloc((void**)&dP, sizeP);
  cudaMalloc((void**)&deprocess, sizeEprocess);
  cudaMalloc((void**)&dsurvival, sizeSurvival);
  cudaMalloc((void**)&drepay_coeff, sizeRepay);

  // Transfer
  cudaMemcpy(dagrid,   agrid,   sizeAgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(dmgrid,   mgrid,   sizeMgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(dhgrid,   hgrid,   sizeHgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(drgrid,   rgrid,   sizeRgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(dlgrid,   lgrid,   sizeLgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(ddgrid,   dgrid,   sizeDgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(dygrid,   ygrid,   sizeYgrid, cudaMemcpyHostToDevice);
  cudaMemcpy(dP,       P,       sizeP, cudaMemcpyHostToDevice);
  cudaMemcpy(deprocess,eprocess,   sizeEprocess, cudaMemcpyHostToDevice);
  cudaMemcpy(dsurvival,survival,   sizeSurvival, cudaMemcpyHostToDevice);
  cudaMemcpy(drepay_coeff,repay_coeff,   sizeRepay, cudaMemcpyHostToDevice);


//----------------------------------------------//
//--------   DIMENSIONS OF PARALLEL   ----------//
//----------------------------------------------//

  int xpar = nd;
  int ypar = nh;
  int zpar = ny;

  if(xpar*ypar*zpar > 512){
    cout << "-----------------------------" << endl;
    cout << "--- DIMENSIONES ESTAN MAL ---" << endl;
    cout << "-----------------------------" << endl;
  }

  // Hacia derecha va a ser H, hacia abajo va a ser M
  dim3 dimBlock(nd, nh, ny);
  dim3 dimGrid(nm, na);


//----------------------------------------------//
//--------  VALUE FUNCTION ITERATION  ----------//
//----------------------------------------------//

  // VFI starts
  double err   = 1;
  int itnum   = 1;

  clock_t t;
  clock_t t0;
  t0  = clock();
  t   = t0;

  //----------------------------------------------//
  //--------   COMPUTATION OF PRICES    ----------//
  //----------------------------------------------//

  cudaMemcpy(dPricing_tomorrow, Pricing_tomorrow, sizemats, cudaMemcpyHostToDevice);
  cudaMemcpy(dPricing_guess, Pricing_guess, sizemats, cudaMemcpyHostToDevice);

  for(int it=T-1; it>=0; it--){
    Pmort<<<dimGrid,dimBlock>>>(T, na, nm, nh, nd, ny, rrho, r, Ph_tomorrow, ddeltabar_tomorrow, sunk, interm, rec_probab, dP, 
                                ddgrid, dmgrid, dhgrid, drgrid, dsurvival, drepay_coeff, it, dDefault, dRenew, dPolicya, dPricing_guess, dPricing_tomorrow);
    cudaDeviceSynchronize();      
  }

  cudaMemcpy(Pricing_guess,   dPricing_guess,   sizemats, cudaMemcpyDeviceToHost);


  double conv_prov;

  while((err > tol) & (itnum < maxiter)){
   
  //--------------------------------------//
  //--------  Policy functions  ----------//
  //--------------------------------------//
 
    // Copy back from host (CPU) to device (GPU)
    cudaMemcpy(dPricing_guess, Pricing_guess, sizemats, cudaMemcpyHostToDevice);

    // State variables
    for(int it=T-1; it>=0; it--){
      if(print_out == 1){
        if(it == (T-1)){
        cout << it+1 << ", " << flush;
        } else{
        cout << it+1 << ", " << flush;
        }
      } else{
        if(it % 10 == 0){
          cout << "." << flush;
        }
      }
     
      vfi_continuation<<<dimGrid,dimBlock>>>(T, Tretirement, na, nm, nh, nr, nl, nd, ny, uti, rrho, bbeta, Ph_today, q, Pa, ddeltabar_today, ssigma, ppsi, kkappa, 
                                            tthetalab, eetalab, fcost, refcost, pension, sstax, ltax, lumpsum, oomega, rec_probab, sunk, dincshock, dmortsubsidy,
                                            dagrid, dmgrid, dhgrid, drgrid, dlgrid, ddgrid, dygrid, dP, deprocess, dsurvival, drepay_coeff, 
                                            it, equivalent, multiplier, dValue_future, dValue_equiv_future, dValue, dValue_equiv, dDefault, dRenew, dPolicya, dPolicym, dPolicyh, dPolicyr, dPolicyl, dPolicyc, dPricing_guess); 

      cudaDeviceSynchronize();
    }
    

  //--------------------------------------//
  //--------  Pricing function  ----------//
  //--------------------------------------//


    cudaMemcpy(dPricing_tomorrow, Pricing_tomorrow, sizemats, cudaMemcpyHostToDevice);
    cudaMemcpy(dPricing_guess2, Pricing_guess2, sizemats, cudaMemcpyHostToDevice);

    for(int it=T-1; it>=0; it--){
      Pmort<<<dimGrid,dimBlock>>>(T, na, nm, nh, nd, ny, rrho, r, Ph_tomorrow, ddeltabar_tomorrow, sunk, interm, rec_probab, dP, 
                                  ddgrid, dmgrid, dhgrid, drgrid, dsurvival, drepay_coeff, it, dDefault, dRenew, dPolicya, dPricing_guess2, dPricing_tomorrow);
      cudaDeviceSynchronize();      
    }

    cudaMemcpy(Pricing_guess2, dPricing_guess2, sizemats, cudaMemcpyDeviceToHost);


    if(convergence != 0.0){ 
        conv_prov = 0.0;
        if (err <= 0.01 & err > 0.001){
          conv_prov = 0.5;
        } else if(err <= 0.001 & err > 0.0002){
          conv_prov = 0.7;
        } else if(err <= 0.002 & err > 0.0001){
          conv_prov = 0.9;
        } else if(err <= 0.0001 & err > 0.000075){
          conv_prov = 0.9;
        } else if(err <= 0.000075 & err > 0.00005){
          conv_prov = 0.9;
        } else if(err <= 0.00005){
          conv_prov = 0.9;
        }
    }
    
    // Compare new and old pricing functions 
    err = 0.0;
    for(int it=T-1; it>=0; it--){
      for(int iy=0; iy<ny; iy++){
        for(int ia=0; ia<na; ia++){
          for(int im=0; im<nm; im++){
            for(int ih=0; ih<nh; ih++){
              for(int id=0; id<nd; id++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                
                err = err + absol(Pricing_guess2[ind] - Pricing_guess[ind])/sizemats;

              }
            }
          }
        }
      }
    }
    
    t = clock() - t0;

    if (itnum % 1 == 0 && print_out == 1){
      // cout << "     -------------------------------------------" << endl;
      // cout << "FUNCTION: Iteration number: " << itnum << ". Time: " << 1000000*((double)t)/CLOCKS_PER_SEC << " seconds. Error " << err << endl;
      cout << "  It: " << itnum << ". Time: " << ((double)t/60)/CLOCKS_PER_SEC << " minutes. Error " << err << ", conv= " << conv_prov << endl;
      // cout << "     -------------------------------------------" << endl;
    }

    itnum = itnum + 1;

    for(int it=T-1; it>=0; it--){
      for(int iy=0; iy<ny; iy++){
        for(int ia=0; ia<na; ia++){
          for(int im=0; im<nm; im++){
            for(int ih=0; ih<nh; ih++){
              for(int id=0; id<nd; id++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

                if((err > tol) & (itnum < maxiter)){
                  Pricing_guess[ind] = conv_prov*Pricing_guess[ind] + (1-conv_prov)*Pricing_guess2[ind];
                  Pricing_guess2[ind] = 0.0;
                } 
              }
            }
          }
        }
      }
    }

  }

  // Copy back from device (GPU) to host (CPU)
  cudaMemcpy(Value,         dValue,         sizemats,    cudaMemcpyDeviceToHost);
  cudaMemcpy(Value_equiv,   dValue_equiv,   sizemats,    cudaMemcpyDeviceToHost);
  cudaMemcpy(Policyc,       dPolicyc,       sizemats,    cudaMemcpyDeviceToHost);
  cudaMemcpy(Pricing_guess, dPricing_guess, sizemats,    cudaMemcpyDeviceToHost);
  cudaMemcpy(Default,       dDefault,       sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Renew,         dRenew,         sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Policya,       dPolicya,       sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Policym,       dPolicym,       sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Policyh,       dPolicyh,       sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Policyr,       dPolicyr,       sizematsint, cudaMemcpyDeviceToHost);
  cudaMemcpy(Policyl,       dPolicyl,       sizematsint, cudaMemcpyDeviceToHost);

  // Free variables in device memory
  cudaFree(dValue);
  cudaFree(dValue_equiv);
  cudaFree(dValue_future);
  cudaFree(dValue_equiv_future);
  cudaFree(dPolicyc);
  cudaFree(dDefault);
  cudaFree(dRenew);
  cudaFree(dPolicya);
  cudaFree(dPolicym);
  cudaFree(dPolicyh);
  cudaFree(dPolicyr);
  cudaFree(dPolicyl);

  cudaFree(dmortsubsidy);
  cudaFree(dincshock);

  cudaFree(dPricing_tomorrow);
  cudaFree(dPricing_guess);
  cudaFree(dPricing_guess2);

  cudaFree(dagrid);
  cudaFree(dmgrid);
  cudaFree(dhgrid);
  cudaFree(drgrid);
  cudaFree(dlgrid);
  cudaFree(ddgrid);
  cudaFree(dygrid);
  cudaFree(dP);
  cudaFree(deprocess);
  cudaFree(dsurvival);
  cudaFree(drepay_coeff);
}



//======================================
//         Aggregate paths
//======================================


// This is the distribution in steady state
void probability_paths(double* Pcond, double* Puncond,
                       const double* P,
                       const parameters params,
                        const double *survival,
                        const int* Policya,
                        const int* Policym,
                        const int* Policyh){

  const int T  = params.T;
  const int na = params.na;
  const int nm = params.nm;
  const int nh = params.nh;
  const int nd = params.nd;
  const int ny = params.ny;

  const int Tretirement  = params.Tretirement;

  const double rrho = params.rrho;

  // Steady state distribution for y
  double steady_distr[ny];

  steady_distribution(P, params, steady_distr);
  
  //--------------------------------------//
  //----   Distribution computation  -----//
  //--------------------------------------//

  int ind;
  int indp;  
  int indp2;  


  for(int it=0; it<T-1; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              if(it == 0 && ih == 0 && im == 0 && ia == 0){
                Pcond[ind]        = steady_distr[iy]*(1/(double)nd);
                Puncond[ind]      = steady_distr[iy]*(1/(double)nd);
              }
              
              if(Pcond[ind] > 0){
                for(int iyp=0; iyp<ny; iyp++){
                  for(int idp=0; idp<nd; idp++){
                    for(int ihp=0; ihp<nh; ihp++){
                      for(int imp=0; imp<nm; imp++){
                        for(int iap=0; iap<na; iap++){
                          
                          // Before retirement there is income uncertainty
                          if(it < Tretirement){

                            indp  = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                            indp2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;

                            if(imp == Policym[ind] && ihp == Policyh[ind] && iap == Policya[ind]){
                              Pcond[indp]   = Pcond[indp]   + rrho*P[iy*ny+iyp]*(1/(double)nd)*Pcond[ind];
                              Puncond[indp] = Puncond[indp] + rrho*P[iy*ny+iyp]*(1/(double)nd)*Puncond[ind]*survival[it];

                              Pcond[indp2]   = Pcond[indp2]   + (1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Pcond[ind];
                              Puncond[indp2] = Puncond[indp2] + (1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Puncond[ind]*survival[it];
                            }

                          // After retirement individuals keep the last productivity shock and receive benefits according to that
                          } else{

                            indp  = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                            indp2 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;

                            if(imp == Policym[ind] && ihp == Policyh[ind] && iap == Policya[ind] && iyp == iy){
                              Pcond[indp]   = Pcond[indp]   + rrho*(1/(double)nd)*Pcond[ind];
                              Puncond[indp] = Puncond[indp] + rrho*(1/(double)nd)*Puncond[ind]*survival[it];

                              Pcond[indp2]   = Pcond[indp2]   + (1-rrho)*(1/(double)nd)*Pcond[ind];
                              Puncond[indp2] = Puncond[indp2] + (1-rrho)*(1/(double)nd)*Puncond[ind]*survival[it];
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
        }
      }
    }
  }
}



// This is the distribution during a transition
void probability_paths_trans(const double* Pcond_initial, const double* Puncond_initial,
                             double* Pcond_final, double* Puncond_final,
                             const double* P,
                             const parameters params,
                             const double *survival,
                             const int* Policya_with,
                             const int* Policym_with,
                             const int* Policyh_with,
                             const int* Policya_without,
                             const int* Policym_without,
                             const int* Policyh_without,
                             const double prob_mistake,
                             const double *mortsubsidy,
                             const int *subs_eligible,
                             const int *subs_target){

  const int T  = params.T;
  const int na = params.na;
  const int nm = params.nm;
  const int nh = params.nh;
  const int nd = params.nd;
  const int ny = params.ny;

  const int Tretirement  = params.Tretirement;

  const double rrho = params.rrho;

  // Steady state distribution for y
  double steady_distr[ny];

  steady_distribution(P, params, steady_distr);
  
  //--------------------------------------//
  //----   Distribution computation  -----//
  //--------------------------------------//

  int ind;
  int indp;  
  int indp2;  

  for(int it=0; it<T-1; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              if(it == 0 && ih == 0 && im == 0 && ia == 0){
                Pcond_final[ind]        = steady_distr[iy]*(1/(double)nd);
                Puncond_final[ind]      = steady_distr[iy]*(1/(double)nd);
              }
              
              if(Pcond_initial[ind] > 0){
                for(int iyp=0; iyp<ny; iyp++){
                  for(int idp=0; idp<nd; idp++){
                    for(int ihp=0; ihp<nh; ihp++){
                      for(int imp=0; imp<nm; imp++){
                        for(int iap=0; iap<na; iap++){
                          
                          // Errores son: 1. No darle a defaulters, 2. Darle a non-defaulters => hay 4 casos

                          if(subs_eligible[ind] == 0){ 

                            // Non-eligible households - no mistakes made

                            indp  = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                            indp2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;

                            // Before retirement there is income uncertainty
                            if(it < Tretirement){

                              if(imp == Policym_without[ind] && ihp == Policyh_without[ind] && iap == Policya_without[ind]){
                                Pcond_final[indp]   = Pcond_final[indp]   + rrho*P[iy*ny+iyp]*(1/(double)nd)*Pcond_initial[ind];
                                Puncond_final[indp] = Puncond_final[indp] + rrho*P[iy*ny+iyp]*(1/(double)nd)*Puncond_initial[ind]*survival[it];

                                Pcond_final[indp2]   = Pcond_final[indp2]   + (1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Pcond_initial[ind];
                                Puncond_final[indp2] = Puncond_final[indp2] + (1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Puncond_initial[ind]*survival[it];
                              }


                            // After retirement individuals keep the last productivity shock and receive benefits according to that
                            } else{

                              if(imp == Policym_without[ind] && ihp == Policyh_without[ind] && iap == Policya_without[ind] && iyp == iy){
                                Pcond_final[indp]   = Pcond_final[indp]   + rrho*(1/(double)nd)*Pcond_initial[ind];
                                Puncond_final[indp] = Puncond_final[indp] + rrho*(1/(double)nd)*Puncond_initial[ind]*survival[it];

                                Pcond_final[indp2]   = Pcond_final[indp2]   + (1-rrho)*(1/(double)nd)*Pcond_initial[ind];
                                Puncond_final[indp2] = Puncond_final[indp2] + (1-rrho)*(1/(double)nd)*Puncond_initial[ind]*survival[it];
                              }

                            }

                          } else{

                            // Eligible households - mistakes made on selection process

                            // Before retirement there is income uncertainty
                            if(it < Tretirement){

                              indp  = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                              indp2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;

                              if(subs_target[ind] == 1){ // If we should give him subsidy

                                // No errors
                                if(imp == Policym_with[ind] && ihp == Policyh_with[ind] && iap == Policya_with[ind]){
                                  Pcond_final[indp]   = Pcond_final[indp]   + (1-prob_mistake)*rrho*P[iy*ny+iyp]*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp] = Puncond_final[indp] + (1-prob_mistake)*rrho*P[iy*ny+iyp]*(1/(double)nd)*Puncond_initial[ind]*survival[it];

                                  Pcond_final[indp2]   = Pcond_final[indp2]   + (1-prob_mistake)*(1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp2] = Puncond_final[indp2] + (1-prob_mistake)*(1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Puncond_initial[ind]*survival[it];
                                }

                                // With errors
                                if(imp == Policym_without[ind] && ihp == Policyh_without[ind] && iap == Policya_without[ind]){
                                  Pcond_final[indp]   = Pcond_final[indp]   + prob_mistake*rrho*P[iy*ny+iyp]*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp] = Puncond_final[indp] + prob_mistake*rrho*P[iy*ny+iyp]*(1/(double)nd)*Puncond_initial[ind]*survival[it];

                                  Pcond_final[indp2]   = Pcond_final[indp2]   + prob_mistake*(1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp2] = Puncond_final[indp2] + prob_mistake*(1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Puncond_initial[ind]*survival[it];
                                }

                              } else{ // If we should not give him subsidy

                                // No errors
                                if(imp == Policym_without[ind] && ihp == Policyh_without[ind] && iap == Policya_without[ind]){
                                  Pcond_final[indp]   = Pcond_final[indp]   + (1-prob_mistake)*rrho*P[iy*ny+iyp]*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp] = Puncond_final[indp] + (1-prob_mistake)*rrho*P[iy*ny+iyp]*(1/(double)nd)*Puncond_initial[ind]*survival[it];

                                  Pcond_final[indp2]   = Pcond_final[indp2]   + (1-prob_mistake)*(1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp2] = Puncond_final[indp2] + (1-prob_mistake)*(1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Puncond_initial[ind]*survival[it];
                                }

                                // With errors
                                if(imp == Policym_with[ind] && ihp == Policyh_with[ind] && iap == Policya_with[ind]){
                                  Pcond_final[indp]   = Pcond_final[indp]   + prob_mistake*rrho*P[iy*ny+iyp]*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp] = Puncond_final[indp] + prob_mistake*rrho*P[iy*ny+iyp]*(1/(double)nd)*Puncond_initial[ind]*survival[it];

                                  Pcond_final[indp2]   = Pcond_final[indp2]   + prob_mistake*(1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp2] = Puncond_final[indp2] + prob_mistake*(1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Puncond_initial[ind]*survival[it];
                                }
                              }

                            // After retirement individuals keep the last productivity shock and receive benefits according to that
                            } else{

                              indp  = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                              indp2 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;

                              if(subs_target[ind] == 1){ // If we should give him subsidy

                                if(imp == Policym_with[ind] && ihp == Policyh_with[ind] && iap == Policya_with[ind] && iyp == iy){
                                  Pcond_final[indp]   = Pcond_final[indp]   + (1-prob_mistake)*rrho*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp] = Puncond_final[indp] + (1-prob_mistake)*rrho*(1/(double)nd)*Puncond_initial[ind]*survival[it];

                                  Pcond_final[indp2]   = Pcond_final[indp2]   + (1-prob_mistake)*(1-rrho)*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp2] = Puncond_final[indp2] + (1-prob_mistake)*(1-rrho)*(1/(double)nd)*Puncond_initial[ind]*survival[it];
                                }

                                if(imp == Policym_without[ind] && ihp == Policyh_without[ind] && iap == Policya_without[ind] && iyp == iy){
                                  Pcond_final[indp]   = Pcond_final[indp]   + prob_mistake*rrho*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp] = Puncond_final[indp] + prob_mistake*rrho*(1/(double)nd)*Puncond_initial[ind]*survival[it];

                                  Pcond_final[indp2]   = Pcond_final[indp2]   + prob_mistake*(1-rrho)*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp2] = Puncond_final[indp2] + prob_mistake*(1-rrho)*(1/(double)nd)*Puncond_initial[ind]*survival[it];
                                }

                              } else{

                                if(imp == Policym_without[ind] && ihp == Policyh_without[ind] && iap == Policya_without[ind] && iyp == iy){
                                  Pcond_final[indp]   = Pcond_final[indp]   + (1-prob_mistake)*rrho*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp] = Puncond_final[indp] + (1-prob_mistake)*rrho*(1/(double)nd)*Puncond_initial[ind]*survival[it];

                                  Pcond_final[indp2]   = Pcond_final[indp2]   + (1-prob_mistake)*(1-rrho)*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp2] = Puncond_final[indp2] + (1-prob_mistake)*(1-rrho)*(1/(double)nd)*Puncond_initial[ind]*survival[it];
                                }

                                if(imp == Policym_with[ind] && ihp == Policyh_with[ind] && iap == Policya_with[ind] && iyp == iy){
                                  Pcond_final[indp]   = Pcond_final[indp]   + prob_mistake*rrho*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp] = Puncond_final[indp] + prob_mistake*rrho*(1/(double)nd)*Puncond_initial[ind]*survival[it];

                                  Pcond_final[indp2]   = Pcond_final[indp2]   + prob_mistake*(1-rrho)*(1/(double)nd)*Pcond_initial[ind];
                                  Puncond_final[indp2] = Puncond_final[indp2] + prob_mistake*(1-rrho)*(1/(double)nd)*Puncond_initial[ind]*survival[it];
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
            }
          }
        }
      }
    }
  }
}




//============================================================================================================//
//         Aggregate de los deudores sin que se les borre la deuda, a ver cual es la base de deudores         //
//         para calcular debt-to-income ratio con la verdadera base                                           //
//============================================================================================================//

void probability_paths_debt(double* Pcond_debt, double* Puncond_debt,
                            const double* P,
                            const parameters params,
                            const double *survival,
                            const int* Policya,
                            const int* Policym,
                            const int* Policyh){

  const int T  = params.T;
  const int na = params.na;
  const int nm = params.nm;
  const int nh = params.nh;
  const int nd = params.nd;
  const int ny = params.ny;

  // Steady state distribution for y
  double steady_distr[ny];

  steady_distribution(P, params, steady_distr);
  
  //--------------------------------------//
  //----   Distribution computation  -----//
  //--------------------------------------//

  int ind;
  int indp;  


  for(int it=0; it<T-1; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              if(it == 0 && ih == 0 && im == 0 && ia == 0){
                Pcond_debt[ind]        = steady_distr[iy]*(1/(double)nd);
                Puncond_debt[ind]      = steady_distr[iy]*(1/(double)nd);
              }
              
              if(Pcond_debt[ind] > 0){
                for(int iyp=0; iyp<ny; iyp++){
                  for(int idp=0; idp<nd; idp++){
                    for(int ihp=0; ihp<nh; ihp++){
                      for(int imp=0; imp<nm; imp++){
                        for(int iap=0; iap<na; iap++){
                          indp  = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                          
                          if(imp == Policym[ind] && ihp == Policyh[ind] && iap == Policya[ind]){
                            Pcond_debt[indp]   = Pcond_debt[indp]   + P[iy*ny+iyp]*(1/(double)nd)*Pcond_debt[ind];
                            Puncond_debt[indp] = Puncond_debt[indp] + P[iy*ny+iyp]*(1/(double)nd)*Puncond_debt[ind]*survival[it];
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
      }
    }
  }
}