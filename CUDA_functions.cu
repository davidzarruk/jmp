

//======================================
//         Utility function
//======================================

__device__ double u(const double c, const double h, const double l, 
                    const double ssigma, const double ppsi, const int uti, 
                    const double kkappa, const double tthetalab, const double eetalab){
  
  double utility = 0.0; 
  
  if(uti == 1){
    // CES
    utility = (powf(powf(ppsi*powf(c, kkappa) + (1-ppsi)*powf(h, kkappa), (1/kkappa)), 1-ssigma) / (1-ssigma)) - (tthetalab*pow(l, 1 + eetalab)/(1 + eetalab)) ;
  } else if(uti == 2){
    // Utility function 2 
    utility = powf(powf(c, ppsi)*powf(h, 1-ppsi), 1-ssigma) / (1-ssigma);
  }
  
  if(c <= 0 || h <= 0){
    utility = powf(-10, 15);
  }

  return(utility);
}


__device__ double mortg_function(const double m, const double Pm, const double oomega, const double h, const double Ph, const double repay_coeff){
  
  double mortgage = 0.0; 
  
  if(m*(1+repay_coeff) <= (1-oomega)*Ph*h){
  // if(m*Pm <= (1-oomega)*Ph*h){

    mortgage = m*Pm;

  } else{

    mortgage = -10000.0;

  }

  return(mortgage);
}

__device__ double maximumab(const double a, const double b){
  
  double ans = a; 
  
  if(b >= a){
    ans = b;
  }

  return(ans);
}


//======================================
//         Pricing function
//======================================


__global__ void Pmort(const int T, const int na, const int nm, const int nh, const int nd, const int ny,
                      const double rrho, const double r, const double Ph, const double ddeltabar, const double sunk, const double interm, const double rec_probab,
                      const double *P, 
                      const double *dgrid, 
                      const double *mgrid, 
                      const double *hgrid, 
                      const double *rgrid, 
                      const double *survival,
                      const double *repay_coeff,
                      const int it,
                      const int* Def,
                      const int* Renew,
                      const int* Policya,
                      double* pricing,
                      double* pricing_guess){

  const int id  = threadIdx.x;
  const int ih  = threadIdx.y;
  const int iy  = threadIdx.z;

  const int im  = blockIdx.x;
  const int ia  = blockIdx.y;

  // If mortgage is equal to zero, the price is not relevant.
  if(im > 0){
    int ind;
    int ind2;
    int ind3;
    
    int iap = 0;
  
    // For every state variable, I compute the pricing function
    ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
  
    if(it == T-1){
  
      pricing[ind] = 0.0;
  
    } else{
  
      // Expected value is sum over tomorrow's possible shocks times probabilities: P[iy][iyp]*(1/nd)
      for(int iyp=0; iyp<ny; iyp++){
        for(int idp=0; idp<nd; idp++){
  
          ind2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + ia;
  
          iap = Policya[ind2];
          ind3 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + iap;
          
          // pricing[ind] = pricing[ind] + ((survival[it]*rrho/(mgrid[im]*(1+r+interm)))*P[iy*ny+iyp]*(1/(double)nd)*(Def[ind2]*(1-Renew[ind2])*Ph*(1-dgrid[idp] - ddeltabar)*hgrid[ih]*(1-sunk) + // If he defaults 
          pricing[ind] = pricing[ind] + ((survival[it]*rrho/(mgrid[im]*(1+r+interm)))*P[iy*ny+iyp]*(1/(double)nd)*(Def[ind2]*(1-Renew[ind2])*(Ph*(1-dgrid[idp] - ddeltabar)*hgrid[ih]*(1-sunk) - Ph*dgrid[idp]*hgrid[ih]) + // If he defaults 
                                                                                                          (1-Def[ind2])*((1-Renew[ind2])*(mgrid[im] + pricing_guess[ind3]*mgrid[im]) + // If he pays and continues with mortgage 
                                                                                                          (1-Def[ind2])*Renew[ind2]*(mgrid[im] + repay_coeff[it+1]*mgrid[im]))));
        }
      }
    }
  }
}


//======================================
//     Value Function Computation
//======================================

__global__ void vfi(const int T, const int Tretirement, const int na, const int nm, const int nh, const int nr, const int nl,
                    const int nd, const int ny, const int uti, const double rrho, 
                    const double bbeta, const double Ph, const double q, 
                    const double Pa, const double ddeltabar, const double ssigma, 
                    const double ppsi, const double kkappa, const double tthetalab, const double eetalab,
                    const double fcost, const double refcost, const double pension, const double sstax, const double ltax,
                    const double lumpsum, const double oomega, const double rec_probab, const double sunk,
                    const double *incshock, const double *mortsubsidy,
                    const double *agrid, const double *mgrid, const double *hgrid, const double *rgrid, const double *lgrid,
                    const double *dgrid, const double *ygrid, const double *P, 
                    const double *eprocess, const double *survival, const double *repay_coeff,
                    const int it,
                    const int equivalent,
                    const double multiplier,
                    double* Value,
                    double* Value_equiv,
                    int* Default,
                    int* Renew,
                    int* Policya,
                    int* Policym,
                    int* Policyh,
                    int* Policyr,
                    int* Policyl,
                    double* Policyc,
                    double* Pricing_guess){

  int ind;
  int ind1;
  int ind2;
  int indsubs;
  
  double VV     = 0.0;
  double VV_eq  = 0.0;

  // Value normal
  double VVk    = powf(-10,11);  double VVn    = powf(-10,11);  double VVd    = powf(-10,11);
  double Vexk   = 0.0;           double Vexn   = 0.0;           double Vexd   = 0.0;

  // Value de consumption equivalent
  double VVk_eq    = powf(-10,11);  double VVn_eq    = powf(-10,11);  double VVd_eq    = powf(-10,11);
  double Vexk_eq   = 0.0;           double Vexn_eq   = 0.0;           double Vexd_eq   = 0.0;

  double cck    = 0;             double ccn    = 0;             double ccd    = 0;    // Consumption
  
  int hhk       = 0;             int hhn       = 0;             int hhd       = 0;    // Home ownership
  int hrk       = 0;             int hrn       = 0;             int hrd       = 0;    // Home renting
  int mmk       = 0;             int mmn       = 0;             int mmd       = 0;    // Mortgage
  int aak       = 0;             int aan       = 0;             int aad       = 0;    // Savings
  int llk       = 0;             int lln       = 0;             int lld       = 0;    // Labor
    
  double yy;
  double aa;
  double ll;
  double aaprime;
  double mm;
  double mmprime;
  double hh;
  double hhprime;
  double hhrent;
  double ddelta;
  double cons;
  double pprice;
  double mort_received;
  double mortgage_subsidy;
  double refinance_cost;

  // State variables that are parallelized
  // const int im  = blockIdx.x * blockDim.x + threadIdx.x;
  const int im  = blockIdx.x;
  const int ia  = blockIdx.y;
  const int id  = threadIdx.x;
  const int ih  = threadIdx.y;
  const int iy  = threadIdx.z;
  
  aa      = agrid[ia];
  mm      = mgrid[im];
  hh      = hgrid[ih];
  ddelta  = dgrid[id];

  ind     = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
  indsubs = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;
  
  mortgage_subsidy = mortsubsidy[indsubs];

  // Control variables
  for(int il=0; il<nl; il++){
    for(int iap=0; iap<na; iap++){
      for(int ihre=0; ihre<nr; ihre++){

        ll      = lgrid[il];
        aaprime = agrid[iap];
        hhrent  = rgrid[ihre];
        
        if(it < Tretirement){
          yy = ygrid[iy]*eprocess[it]*ll*(1-sstax-ltax);
        } else{
          yy = ygrid[iy]*pension;
        }

        // Keeping the same mortgage
        cons = aa + q*hh + yy*(1-incshock[it]) - mm - q*hhrent - Pa*aaprime - Ph*(ddelta + ddeltabar)*hh - lumpsum;
        
        Vexk       = 0.0;
        Vexk_eq    = 0.0;
        if(it < T-1){
          for(int idp=0; idp<nd; idp++){

            if(it < Tretirement){  // Income uncertainty before retirement
              for(int iyp=0; iyp<ny; iyp++){
                
                ind1 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + iap;
                ind2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + 0*na + iap;
                
                Vexk = Vexk + P[iy*ny+iyp]*(1/(double)nd)*(rrho*Value[ind1] +  // Keeps mortg
                                                        (1-rrho)*Value[ind2]);    // Mortg disappears

                Vexk_eq = Vexk_eq + P[iy*ny+iyp]*(1/(double)nd)*(rrho*Value_equiv[ind1] +  // Keeps mortg
                                                               (1-rrho)*Value_equiv[ind2]);    // Mortg disappears
              }
            } else{   // Certainty after retirement
              
              ind1 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + iap;
              ind2 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + 0*na + iap;
              
              Vexk = Vexk + (1/(double)nd)*(rrho*Value[ind1] +  // Keeps mortg
                                            (1-rrho)*Value[ind2]);    // Mortg disappears

              Vexk_eq = Vexk_eq + (1/(double)nd)*(rrho*Value_equiv[ind1] +  // Keeps mortg
                                                  (1-rrho)*Value_equiv[ind2]);    // Mortg disappears
            }
          }
        }
        
        VV    = u(cons, hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexk;

        if(equivalent == 1){
          VV_eq = u(cons*(1+multiplier), hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexk_eq;
        }

        if(VV > VVk){
          VVk    = VV;
          VVk_eq = VV_eq;
          hhk    = ih;
          hrk    = ihre;
          mmk    = im;
          aak    = iap;
          cck    = cons;
          llk    = il;
        }
        
        // Defaulting => Household loses savings
        cons = maximumab(aa - rec_probab*((1+repay_coeff[it])*mm - Ph*(1-ddelta - ddeltabar)*hh*(1-sunk)), 0) + yy*(1-incshock[it]) - q*hhrent - Pa*aaprime - lumpsum;

        Vexd       = 0.0;
        Vexd_eq    = 0.0;
        if(it < T-1){
          for(int idp=0; idp<nd; idp++){

            if(it < Tretirement){  // Income uncertainty before retirement
              for(int iyp=0; iyp<ny; iyp++){
                ind1 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + 0*nm*na + 0*na + iap;

                Vexd = Vexd + P[iy*ny+iyp]*(1/(double)nd)*Value[ind1];

                Vexd_eq = Vexd_eq + P[iy*ny+iyp]*(1/(double)nd)*Value_equiv[ind1];
              }
            } else{   // Certainty after retirement
              ind1 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + 0*nm*na + 0*na + iap;
              
              Vexd = Vexd + (1/(double)nd)*Value[ind1];

              Vexd_eq = Vexd_eq + (1/(double)nd)*Value_equiv[ind1];
            }
          }
        }
        
        VV    = u(cons, hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexd;

        if(equivalent == 1){
          VV_eq = u(cons*(1+multiplier), hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexd_eq;
        }

        if(VV > VVd){
          VVd    = VV;
          VVd_eq = VV_eq;
          hhd    = 0;
          hrd    = ihre;
          mmd    = 0;
          aad    = iap;
          ccd    = cons;
          lld    = il;
        }
          
        // New mortgage
        for(int imp=0; imp<nm; imp++){
          for(int ihp=0; ihp<nh; ihp++){
            
            if(im == 0){
              if(imp > 0){
                refinance_cost = fcost;      // Issuing new mortgage
              } else{
                refinance_cost = 0.0;
              }
            } else{
              if(imp > 0){
                refinance_cost = refcost;  // Refinancing mortgage
              } else{
                refinance_cost = 0.0;      // Paying total debt
              }
            }

            mmprime      = mgrid[imp];
            hhprime      = hgrid[ihp];
            
            ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ihp*nm*na + imp*na + iap;
            pprice = Pricing_guess[ind];

            mort_received = mortg_function(mmprime, pprice, oomega, hhprime, Ph, repay_coeff[it]);
            
            cons = aa + Ph*(1-ddelta - ddeltabar)*hh + q*hhprime + yy*(1-incshock[it]) + mort_received - refinance_cost*(1+repay_coeff[it])*mmprime + mortgage_subsidy - (1+repay_coeff[it])*mm - q*hhrent - Ph*hhprime - Pa*aaprime - lumpsum;
            
            Vexn       = 0.0;
            Vexn_eq    = 0.0;
            if(it < T-1){
              for(int idp=0; idp<nd; idp++){

                if(it < Tretirement){  // Income uncertainty before retirement
                  for(int iyp=0; iyp<ny; iyp++){
                    
                    ind1 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                    ind2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;
                    
                    Vexn = Vexn + P[iy*ny+iyp]*(1/(double)nd)*((rrho * Value[ind1]) + 
                                                            ((1-rrho) * Value[ind2]));

                    Vexn_eq = Vexn_eq + P[iy*ny+iyp]*(1/(double)nd)*((rrho * Value_equiv[ind1]) + 
                                                                     ((1-rrho) * Value_equiv[ind2]));
                  }
                } else{   // Certainty after retirement
                  ind1 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                  ind2 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;
                  
                  Vexn = Vexn + (1/(double)nd)*((rrho * Value[ind1]) + 
                                                ((1-rrho) * Value[ind2]));

                  Vexn_eq = Vexn_eq + (1/(double)nd)*((rrho * Value_equiv[ind1]) + 
                                                      ((1-rrho) * Value_equiv[ind2]));
                }
              }
            }
            
            VV    = u(cons, hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexn;

            if(equivalent == 1){
              VV_eq = u(cons*(1+multiplier), hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexn_eq;
            }

            if(VV > VVn){
              VVn    = VV;
              VVn_eq = VV_eq;
              hhn    = ihp;
              hrn    = ihre;
              mmn    = imp;
              aan    = iap;
              ccn    = cons;
              lln    = il;
            }
          }
        }
      }
    }
  }

  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
    
  if((VVk >= VVd) & (VVk >= VVn)){
    Value[ind]      = VVk;
    Value_equiv[ind]= VVk_eq;
    Policya[ind]    = aak;
    Policyh[ind]    = hhk;
    Policyr[ind]    = hrk;
    Policyl[ind]    = llk;
    Policym[ind]    = mmk;
    Policyc[ind]    = cck;
    Default[ind]    = 0;
    Renew[ind]      = 0;
  } else if((VVd > VVk) & (VVd > VVn)){
    Value[ind]      = VVd;
    Value_equiv[ind]= VVd_eq;
    Policya[ind]    = aad;
    Policyh[ind]    = hhd;
    Policyr[ind]    = hrd;
    Policyl[ind]    = lld;
    Policym[ind]    = mmd;
    Policyc[ind]    = ccd;
    Default[ind]    = 1;
    Renew[ind]      = 0;
  } else{
    Value[ind]      = VVn;
    Value_equiv[ind]= VVn_eq;
    Policya[ind]    = aan;
    Policyh[ind]    = hhn;
    Policyr[ind]    = hrn;
    Policyl[ind]    = lln;
    Policym[ind]    = mmn;
    Policyc[ind]    = ccn;
    Default[ind]    = 0;
    Renew[ind]      = 1;
  }

}



//==================================================================
//     Value Function Computation with different continuation
//==================================================================

__global__ void vfi_continuation(const int T, const int Tretirement, const int na, const int nm, const int nh, const int nr, const int nl,
                                  const int nd,    const int ny, const int uti,
                                  const double rrho,    const double bbeta, const double Ph_today, const double q,  const double Pa, const double ddeltabar_today,
                                  const double ssigma,  const double ppsi,  const double kkappa, const double tthetalab, const double eetalab,
                                  const double fcost, const double refcost, const double pension, const double sstax, const double ltax,
                                  const double lumpsum, const double oomega, const double rec_probab, const double sunk,
                                  const double *incshock, const double *mortsubsidy,
                                  const double *agrid, 
                                  const double *mgrid, 
                                  const double *hgrid, 
                                  const double *rgrid, 
                                  const double *lgrid, 
                                  const double *dgrid, 
                                  const double *ygrid, 
                                  const double *P, 
                                  const double *eprocess, 
                                  const double *survival,
                                  const double *repay_coeff,
                                  const int it,
                                  const int equivalent,
                                  const double multiplier,
                                  const double* Value_future,
                                  const double* Value_equiv_future,
                                  double* Value,
                                  double* Value_equiv,
                                  int* Default,
                                  int* Renew,
                                  int* Policya,
                                  int* Policym,
                                  int* Policyh,
                                  int* Policyr,
                                  int* Policyl,
                                  double* Policyc,
                                  double* Pricing_guess){
  
  int ind;
  int ind1;
  int ind2;
  int indsubs;
  
  double VV     = 0.0;
  double VV_eq  = 0.0;

  double VVk    = powf(-10,11);  double VVn    = powf(-10,11);  double VVd    = powf(-10,11);
  double Vexk   = 0.0;           double Vexn   = 0.0;           double Vexd   = 0.0;

  double VVk_eq    = powf(-10,11);  double VVn_eq    = powf(-10,11);  double VVd_eq    = powf(-10,11);
  double Vexk_eq   = 0.0;           double Vexn_eq   = 0.0;           double Vexd_eq   = 0.0;

  double cck    = 0;             double ccn    = 0;             double ccd    = 0;    // Consumption
  
  int hhk    = 0;  int hhn    = 0;  int hhd    = 0;    // Home ownership
  int hrk    = 0;  int hrn    = 0;  int hrd    = 0;    // Home renting
  int mmk    = 0;  int mmn    = 0;  int mmd    = 0;    // Mortgage
  int aak    = 0;  int aan    = 0;  int aad    = 0;    // Savings
  int llk    = 0;  int lln    = 0;  int lld    = 0;    // Labor

  double yy;
  double aa;
  double ll;
  double aaprime;
  double mm;
  double mmprime;
  double hh;
  double hhprime;
  double hhrent;
  double ddelta;
  double cons;
  double pprice;
  double mort_received;
  double mortgage_subsidy;
  double refinance_cost;

  // State variables that are parallelized
  // const int im  = blockIdx.x * blockDim.x + threadIdx.x;
  const int im  = blockIdx.x;
  const int ia  = blockIdx.y;
  const int id  = threadIdx.x;
  const int ih  = threadIdx.y;
  const int iy  = threadIdx.z;

  aa      = agrid[ia];
  mm      = mgrid[im];
  hh      = hgrid[ih];
  ddelta  = dgrid[id];
 
  ind     = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
  indsubs = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;
  
  mortgage_subsidy = mortsubsidy[indsubs];


  // Control variables
  for(int il=0; il<nl; il++){
    for(int iap=0; iap<na; iap++){
      for(int ihre=0; ihre<nr; ihre++){

        ll      = lgrid[il];
        aaprime = agrid[iap];
        hhrent  = rgrid[ihre];

        if(it < Tretirement){
          yy = ygrid[iy]*eprocess[it]*ll*(1-sstax-ltax);
        } else{
          yy = ygrid[iy]*pension;
        }

        // Keeping the same mortgage
        cons = aa + q*hh + yy*(1-incshock[it]) - mm - q*hhrent - Pa*aaprime - Ph_today*(ddelta + ddeltabar_today)*hh - lumpsum;
        
        Vexk       = 0.0;
        Vexk_eq       = 0.0;
        if(it < T-1){
          for(int idp=0; idp<nd; idp++){

            if(it < Tretirement){  // Income uncertainty before retirement
              for(int iyp=0; iyp<ny; iyp++){
                
                ind1 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + iap;
                ind2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + 0*na + iap;
                
                Vexk = Vexk + P[iy*ny+iyp]*(1/(double)nd)*(rrho*Value_future[ind1] +  // Keeps mortg
                                                        (1-rrho)*Value_future[ind2]);    // Mortg disappears

                Vexk_eq = Vexk_eq + P[iy*ny+iyp]*(1/(double)nd)*(rrho*Value_equiv_future[ind1] +  // Keeps mortg
                                                                (1-rrho)*Value_equiv_future[ind2]);    // Mortg disappears
              }
            } else{   // Certainty after retirement
              
              ind1 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + iap;
              ind2 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + 0*na + iap;
              
              Vexk = Vexk + (1/(double)nd)*(rrho*Value_future[ind1] +  // Keeps mortg
                                            (1-rrho)*Value_future[ind2]);    // Mortg disappears

              Vexk_eq = Vexk_eq + (1/(double)nd)*(rrho*Value_equiv_future[ind1] +  // Keeps mortg
                                                 (1-rrho)*Value_equiv_future[ind2]);    // Mortg disappears
            }
          }
        }
        
        VV    = u(cons, hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexk;
        
        if(equivalent == 1){
          VV_eq = u(cons*(1+multiplier), hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexk_eq;
        }
        
        if(VV > VVk){
          VVk    = VV;
          VVk_eq = VV_eq;
          hhk    = ih;
          hrk    = ihre;
          mmk    = im;
          aak    = iap;
          cck    = cons;
          llk    = il;
        }
        
        // Defaulting => Household loses savings
        cons = maximumab(aa - rec_probab*((1+repay_coeff[it])*mm - Ph_today*(1-ddelta - ddeltabar_today)*hh*(1-sunk)), 0) + yy*(1-incshock[it]) - q*hhrent - Pa*aaprime - lumpsum;
        
        Vexd       = 0.0;
        Vexd_eq       = 0.0;
        if(it < T-1){
          for(int idp=0; idp<nd; idp++){

            if(it < Tretirement){  // Income uncertainty before retirement
              for(int iyp=0; iyp<ny; iyp++){
                ind1 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + iap;
  
                Vexd = Vexd + P[iy*ny+iyp]*(1/(double)nd)*Value_future[ind1];
                Vexd_eq = Vexd_eq + P[iy*ny+iyp]*(1/(double)nd)*Value_equiv_future[ind1];
              }
            } else{   // Certainty after retirement
              ind1 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + iap;
              
              Vexd = Vexd + (1/(double)nd)*Value_future[ind1];
              Vexd_eq = Vexd_eq + (1/(double)nd)*Value_equiv_future[ind1];
            }
          }
        }
        
        VV    = u(cons, hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexd;

        if(equivalent == 1){
          VV_eq = u(cons*(1+multiplier), hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexd_eq;
        }

        if(VV > VVd){
          VVd    = VV;
          VVd_eq = VV_eq;
          hhd    = 0;
          hrd    = ihre;
          mmd    = 0;
          aad    = iap;
          ccd    = cons;
          lld    = il;
        }
          
        // New mortgage
        for(int imp=0; imp<nm; imp++){
          for(int ihp=0; ihp<nh; ihp++){
            
            if(im == 0){
              if(imp > 0){
                refinance_cost = fcost;      // Issuing new mortgage
              } else{
                refinance_cost = 0.0;
              }
            } else{
              if(imp > 0){
                refinance_cost = refcost;  // Refinancing mortgage
              } else{
                refinance_cost = 0.0;      // Paying total debt
              }
            }
  
            mmprime      = mgrid[imp];
            hhprime      = hgrid[ihp];
            
            ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ihp*nm*na + imp*na + iap;
            pprice = Pricing_guess[ind];
            
            mort_received = mortg_function(mmprime, pprice, oomega, hhprime, Ph_today, repay_coeff[it]);
  
            cons = aa + Ph_today*(1-ddelta - ddeltabar_today)*hh + q*hhprime + yy*(1-incshock[it]) - refinance_cost*(1+repay_coeff[it])*mmprime + mort_received + mortgage_subsidy - (1+repay_coeff[it])*mm - q*hhrent - Ph_today*hhprime - Pa*aaprime - lumpsum;
            
            Vexn       = 0.0;
            Vexn_eq    = 0.0;
            if(it < T-1){
              for(int idp=0; idp<nd; idp++){

                if(it < Tretirement){  // Income uncertainty before retirement
                  for(int iyp=0; iyp<ny; iyp++){
                    
                    ind1 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                    ind2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;
                    
                    Vexn = Vexn + P[iy*ny+iyp]*(1/(double)nd)*((rrho * Value_future[ind1]) + 
                                                            ((1-rrho) * Value_future[ind2]));

                    Vexn_eq = Vexn_eq + P[iy*ny+iyp]*(1/(double)nd)*((rrho * Value_equiv_future[ind1]) + 
                                                                    ((1-rrho) * Value_equiv_future[ind2]));
                  }
                } else{   // Certainty after retirement
                  ind1 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                  ind2 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;
                  
                  Vexn = Vexn + (1/(double)nd)*((rrho * Value_future[ind1]) + 
                                                ((1-rrho) * Value_future[ind2]));

                  Vexn_eq = Vexn_eq + (1/(double)nd)*((rrho * Value_equiv_future[ind1]) + 
                                                      ((1-rrho) * Value_equiv_future[ind2]));
                }
              }
            }
            
            VV    = u(cons, hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexn;

            if(equivalent == 1){
              VV_eq = u(cons*(1+multiplier), hhrent, ll, ssigma, ppsi, uti, kkappa, tthetalab, eetalab) + bbeta*survival[it]*Vexn_eq;
            }
            
            if(VV > VVn){
              VVn    = VV;
              VVn_eq = VV_eq;
              hhn    = ihp;
              hrn    = ihre;
              mmn    = imp;
              aan    = iap;
              ccn    = cons;
              lln    = il;
            }
          }
        }
      }
    }
  }

  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
    
  if((VVk >= VVd) & (VVk >= VVn)){
    Value[ind]      = VVk;
    Value_equiv[ind]= VVk_eq;
    Policya[ind]    = aak;
    Policyh[ind]    = hhk;
    Policyr[ind]    = hrk;
    Policyl[ind]    = llk;
    Policym[ind]    = mmk;
    Policyc[ind]    = cck;
    Default[ind]    = 0;
    Renew[ind]      = 0;
  } else if((VVd > VVk) & (VVd >= VVn)){
    Value[ind]      = VVd;
    Value_equiv[ind]= VVd_eq;
    Policya[ind]    = aad;
    Policyh[ind]    = hhd;
    Policyr[ind]    = hrd;
    Policyl[ind]    = lld;
    Policym[ind]    = mmd;
    Policyc[ind]    = ccd;
    Default[ind]    = 1;
    Renew[ind]      = 0;
  } else{
    Value[ind]      = VVn;
    Value_equiv[ind]= VVn_eq;
    Policya[ind]    = aan;
    Policyh[ind]    = hhn;
    Policyr[ind]    = hrn;
    Policyl[ind]    = lln;
    Policym[ind]    = mmn;
    Policyc[ind]    = ccn;
    Default[ind]    = 0;
    Renew[ind]      = 1;
  }

}