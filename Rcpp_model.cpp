//
//  Created by David Zarruk Valencia on June, 2016.
//  Copyright (c) 2016 David Zarruk Valencia. All rights reserved.
//

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <unistd.h>
#include <omp.h>
using std::vector;
using namespace std;


//****************************************************//
//            1. EXOGENOUS PROCESSES                  //
//****************************************************//

//======================================
//         Utility function
//======================================

// [[Rcpp::export]]
double u(double c, double h, double ssigma, double ppsi, int uti){
  
  double utility = 0.0; 
  double kappa = -0.5;
  
  if(uti == 1){
    // CES
    utility = pow(pow(ppsi*pow(c, kappa) + (1-ppsi)*pow(h, kappa), (1/kappa)), 1-ssigma) / (1-ssigma);
  } else if(uti == 2){
    // Utility function 2 
    utility = pow(pow(c, ppsi)*pow(h, 1-ppsi), 1-ssigma) / (1-ssigma);
  }
  
  if(c <= 0 || h < 0){
    utility = pow(-10, 11);
  }

  return(utility);
}


//======================================
//         Grids
//======================================

// [[Rcpp::export]]
vector<double> grida(int na, double amin, double amax){
  
  vector<double> gridaa;
  gridaa.resize(na);
  
  double size = na;
  double astep = (amax - amin) /(size - 1);
  double it = 0;
  
  for(int i = 0; i < na; i++){
    gridaa[i] = amin + it*astep;
    it++;
  }
  
  return(gridaa);
}

// [[Rcpp::export]]
vector<double> gridd(int nd, double dmin, double dmax){
  
  vector<double> gridaa;
  gridaa.resize(nd);
  
  double size = nd;
  double astep = (dmax - dmin) /(size - 1);
  double it = 0;
  
  for(int i = 0; i < nd; i++){
    gridaa[i] = dmin + it*astep;
    it++;
  }
  
  return(gridaa);
}

// [[Rcpp::export]]
vector<double> gridm(int nm, double mmin, double mmax){
  
  vector<double> gridmm;
  gridmm.resize(nm);
  
  double size = nm;
  double mstep = (mmax - mmin) /(size - 1);
  double it = 0;
  
  for(int i = 0; i < nm; i++){
    gridmm[i] = mmin + it*mstep;
    it++;
  }
  
  return(gridmm);
}

// [[Rcpp::export]]
vector<double> gridh(int nh, double hmin, double hmax){
  
  vector<double> gridhh;
  gridhh.resize(nh);
  
  double size = nh;
  double hstep = (hmax - hmin) /(size - 1);
  double it = 0;
  
  for(int i = 0; i < nh; i++){
    gridhh[i] = hmin + it*hstep;
    it++;
  }
  
  return(gridhh);
}

// [[Rcpp::export]]
vector<double> gridy(int ny, double ssigma_y, double llambda_y, double m){
  
  // This grid is made with Tauchen (1986)
  vector<double> gridyy;
  gridyy.resize(ny);
  
  double size = ny;
  double ssigma_yy = sqrt(pow(ssigma_y, 2) / (1 - pow(llambda_y, 2)));
  double ystep = 2*ssigma_yy*m / (size-1);
  
  double it = 0;
  
  for(int i = 0; i < ny; i++){
    gridyy[i] = -m*sqrt(pow(ssigma_y, 2) / (1 - pow(llambda_y, 2))) + it*ystep;
    it = it + 1;
  }
  
  return(gridyy);
}


// [[Rcpp::export]]
double normCDF(double value){
  return 0.5 * erfc(-value * M_SQRT1_2);
}


// [[Rcpp::export]]
vector<vector<double> > yprob(int ny, double ssigma_y, double llambda_y, double m){
  
  // This grid is made with Tauchen (1986)
  vector<double> gridyy;

  gridyy = gridy(ny, ssigma_y, llambda_y, m);
  
  double w = gridyy[1] - gridyy[0];
  
  vector<vector<double> > P;
  P.resize(ny);
  for(int i=0; i<ny; i++){
    P[i].resize(ny);
  }
  
  for(int j = 0; j < ny; j++){
    for(int k = 0; k < ny; k++){
      if(k == 0){
        P[j][k] = normCDF((gridyy[k] - llambda_y*gridyy[j] + (w/2))/ssigma_y);
      } else if(k == ny-1){
        P[j][k] = 1 - normCDF((gridyy[k] - llambda_y*gridyy[j] - (w/2))/ssigma_y);
      } else{
        P[j][k] = normCDF((gridyy[k] - llambda_y*gridyy[j] + (w/2))/ssigma_y) - normCDF((gridyy[k] - llambda_y*gridyy[j] - (w/2))/ssigma_y);
      }
    }
  }
  
  return(P);
}


//****************************************************//
//         3. ENDOGENOUS THINGS                       //
//****************************************************//

//======================================
//         Production
//======================================

// [[Rcpp::export]]
double Y(double K, double L, double Z, double aalpha){
  
  double prod = Z * pow(K, aalpha) * pow(L, 1-aalpha);
  
  return(prod);
}

//======================================
//         Prices
//======================================

// Before tax prices:

// [[Rcpp::export]]
double rfun(double K, double L, double Z, double aalpha, double ddeltak){
  
  double interest = Z * aalpha * pow(L/K, 1-aalpha) - ddeltak;
  
  return(interest);
}

// [[Rcpp::export]]
double q(double K, double L, double Z, double aalpha, double ddeltak, double tauP, double ddeltah){
  
  double rr = rfun(K, L, Z, aalpha, ddeltak);
  
  return(rr + tauP + ddeltah);
}


double absol(double a){
  if(a < 0){
    a = -a;
  }
  return a;
}


void Pmort(double rrho,
             double r, int ny, int na, int nm, int nh, int nd, vector<vector<double> > P,
             vector<double> dgrid, vector<double> mgrid, vector<double> hgrid,
             double* Def,
             double* Renew,
             double Ph,
             double* pricing_guess){
  
  // Initial guess for pricing function: price is equal to 1
  double *pricing;
  size_t sizemats = ny*na*nm*nh*nd*sizeof(double);
  
  pricing = (double*)malloc(sizemats);
  
  for(int iy=0; iy<ny; iy++){
    for(int id=0; id<nd; id++){
      for(int ih=0; ih<nh; ih++){
        for(int im=0; im<nm; im++){
          for(int ia=0; ia<na; ia++){
            pricing[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0.0;
          }
        }
      }
    }
  }
  
  double tol = 0.0001;
  int maxiter = 100;
  double err=1.0;
  int ind;
  int ind2;

  int itnum = 1;
  
  while(err > tol & itnum < maxiter){
    
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              
              // For every state variable, I compute the pricing function
              ind = iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              // Expected value is sum over tomorrow's possible shocks times probabilities: P[iy][iyp]*(1/nd)
              for(int iyp=0; iyp<ny; iyp++){
                for(int idp=0; idp<nd; idp++){
                  
                  ind2 = iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + ia;
                  
                  pricing[ind] = pricing[ind] + (rrho/(mgrid[im]*(1+r)))*P[iy][iyp]*(1/(double)nd)*(Def[ind2]*(1-dgrid[idp])*Ph*hgrid[ih] + // If he defaults 
                                                                                              (1-Def[ind2])*((1-Renew[ind2])*(mgrid[im] + pricing_guess[ind2]*mgrid[im]) + // If he pays and continues with mortgage 
                                                                                                              Renew[ind2]*((1+r)/(1+r-rrho))*mgrid[im]));
                  
                }
              }
            }
          }
        }
      }
    }
    
    err = 0.0;
    
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              err = err + absol(pricing_guess[ind] - pricing[ind]);
              pricing_guess[ind] = pricing[ind];
              pricing[ind] = 0.0;
            }
          }
        }
      }
    }
    if (itnum % 5 == 0){
      // cout << "     Pricing Function - It. #: " << itnum << ". Error " << err << endl;
    }
    itnum++;
  }
  
  cout << "     -------------------------------------------" << endl;
  cout << "     Pricing Function - It. #: " << itnum << ". Error " << err << endl;
  cout << "     -------------------------------------------" << endl;
}

// Supply of housing - linear supply function
double housing_supply(double slope, double hPrice){
  
  double sup = (hPrice - (1-0.69*slope))/slope;
  return(sup);
  
}

//======================================
//         Value Function Iteration
//======================================

// [[Rcpp::export]]
vector<vector<vector<vector<vector<double> > > > > vfi(double tol, int maxiter, int ut,
            int na, double amin, double amax,
            int nm, double mmin, double mmax,
            int nh, double hmin, double hmax, 
            int ny, double ssigma_y, double llambda_y, double m, 
            int nd, double dmin, double dmax,
            double ssigma, double ppsi, double rrho, double bbeta,
            double r, double Ph, double q, double Pa, double ddeltabar) {

  // Initialization  
  vector<vector<vector<vector<vector<double> > > > > Val;
  Val.resize(ny);
  
  for(int i=0; i<ny; i++){
    Val[i].resize(na);
    for(int j=0; j<na; j++){
      Val[i][j].resize(nm);
      for(int l=0; l<nm; l++){
        Val[i][j][l].resize(nh);
        for(int k=0; k<nh; k++){
          Val[i][j][l][k].resize(nd);
          for(int d=0; d<nd; d++){
            Val[i][j][l][k][d]   = 0.0;
          }
        }
      }
    }
  }

  double *Value_guess, *Value, *Default, *Renew, *Policya, *Policym, *Policyh, *Policyc, *Pricing_guess;
  size_t sizemats = ny*na*nm*nh*nd*sizeof(double);
  
  Value_guess = (double*)malloc(sizemats);
  Value       = (double*)malloc(sizemats);
  Default     = (double*)malloc(sizemats);
  Renew       = (double*)malloc(sizemats);
  Policya     = (double*)malloc(sizemats);
  Policym     = (double*)malloc(sizemats);
  Policyh     = (double*)malloc(sizemats);
  Policyc     = (double*)malloc(sizemats);
  Pricing_guess     = (double*)malloc(sizemats);
  
  for(int iy=0; iy<ny; iy++){
    for(int id=0; id<nd; id++){
      for(int ih=0; ih<nh; ih++){
        for(int im=0; im<nm; im++){
          for(int ia=0; ia<na; ia++){
            Value_guess[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = -14.25;
            Value[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0.0;
            Default[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0.0;
            Renew[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0.0;
            Policya[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0.0;
            Policym[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0.0;
            Policyh[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0.0;
            Policyc[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0.0;
            Pricing_guess[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 1.0;
          }
        }
      }
    }
  }


  // Preamble
  vector<double> agrid = grida(na, amin, amax);
  vector<double> dgrid = gridd(nd, dmin, dmax);
  vector<double> mgrid = gridm(nm, mmin, mmax);
  vector<double> hgrid = gridh(nh, hmin, hmax);
  vector<double> ygrid = gridy(ny, ssigma_y, llambda_y, m);

  for(int i = 0; i<ny; i++){
    ygrid[i] = exp(ygrid[i]);
  }
  
  vector<vector<double> > P = yprob(ny, ssigma_y, llambda_y, m);
  
  
  Pmort(rrho, r, ny, na, nm, nh, nd, P, dgrid, mgrid, hgrid, Default, Renew, Ph, Pricing_guess);
  
  // VFI starts
  double err   = 1;
  int itnum = 1;
  int ind;
  double t0;
  double t;

  double VV     = 0.0;
  double VVk    = pow(-10,11);
  double Vexk   = 0.0;
  double VVn    = pow(-10,11);
  double Vexn   = 0.0;
  double VVd    = pow(-10,11);
  double Vexd   = 0.0;
  
  double yy;
  double aa;
  double aaprime;
  double mm;
  double mmprime;
  double hh;
  double hhprime;
  double hhrent;
  double ddelta;
  double cons;
  double pprice;
  
  t0  = omp_get_wtime();
  t   = t0;

  while(err > tol & itnum < maxiter){
    
    // Update the Pricing_guess with the policy functions up to this iteration
    Pmort(rrho, r, ny, na, nm, nh, nd, P, dgrid, mgrid, hgrid, Default, Renew, Ph, Pricing_guess);
    
    // State variables
      
// #pragma omp parallel for shared(Val, Pricing_guess, Value, Value_guess, Default, Renew, P, agrid, hgrid, mgrid, dgrid, ygrid, nd, ny, na, nh, nm, rrho, bbeta, Pa, Ph, ddeltabar, q, r, err, t, t0, ssigma, ppsi, ut) private(yy, aa, mm, hh, ddelta, aaprime, hhrent, cons, Vexk, VVk, Vexd, VVd, Vexn, VVn, VV, mmprime, hhprime, pprice)
    for(int iy=0; iy<ny; iy++){
      for(int ia=0; ia<na; ia++){
        
        for(int im=0; im<nm; im++){
          for(int ih=0; ih<nh; ih++){
            for(int id=0; id<nd; id++){
              
              VVk        = pow(-10,11);
              VVn        = pow(-10,11);
              VVd        = pow(-10,11);
              
              yy      = ygrid[iy];
              aa      = agrid[ia];
              mm      = mgrid[im];
              hh      = hgrid[ih];
              ddelta  = dgrid[ih];
              
              // Control variables
              for(int iap=0; iap<na; iap++){
                for(int ihre=0; ihre<nh; ihre++){
                  aaprime = agrid[iap];
                  hhrent  = hgrid[ihre];
                  
                  // Keeping the same mortgage
                  cons = aa + q*hh + yy - mm - q*hhrent - Pa*aaprime;
                  
                  Vexk       = 0.0;
                  for(int iyp=0; iyp<ny; iyp++){
                    for(int idp=0; idp<nd; idp++){
                      Vexk = Vexk + P[iy][iyp]*(1/(double)nd)*(rrho*Value_guess[iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + iap] +  // Keeps mortg
                                                              (1-rrho)*Value_guess[iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + 0*na + iap]);    // Mortg disappears
                    }
                  }
                    
                  VV = u(cons, hhrent, ssigma, ppsi, ut) + bbeta*Vexk;
                  
                  if(VV > VVk){
                    VVk = VV;
                  }
                  
//                   // Defaulting
//                   cons = aa + yy - q*hhrent - Pa*aaprime;
//                   
//                   Vexd       = 0.0;
//                   for(int iyp=0; iyp<ny; iyp++){
//                     for(int idp=0; idp<nd; idp++){
//                       Vexd = Vexd + P[iy][iyp]*(1/(double)nd)*Value_guess[iyp*nd*nh*nm*na + idp*nh*nm*na + iap];
//                     }
//                   }
//                   
//                   VV = u(cons, hhrent, ssigma, ppsi, ut) + bbeta*Vexd;
//                   
//                   if(iy == 5 & ia == 5 & im == 5 & ih == 5 & cons >0 & hhrent>0){
//                     cout << u(cons, hhrent, ssigma, ppsi, ut) << endl;
//                   }
//                   
//                   if(VV > VVd){
//                     VVd = VV;
//                   }
//                     
                  // New mortgage
//                   for(int imp=0; imp<nm; imp++){
//                     for(int ihp=0; ihp<nh; ihp++){
//                       mmprime      = mgrid[imp];
//                       hhprime      = hgrid[ihp];
//                       
//                       pprice = Pricing_guess[iy*nd*nh*nm*na + id*nh*nm*na + ihp*nm*na + imp*na + iap];
//                       
//                       cons = aa + Ph*(1-ddelta*ddeltabar)*hh + q*hhprime + yy + mmprime*pprice - ((1+r)/(1+r-rrho))*mm - q*hh - Ph*hhprime - Pa*aa;
//                       
//                       Vexn       = 0.0;
//                       for(int iyp=0; iyp<ny; iyp++){
//                         for(int idp=0; idp<nd; idp++){
//                           Vexn = Vexn + P[iy][iyp]*(1/(double)nd)*((rrho * Value_guess[iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap]) + 
//                                                                   ((1-rrho) * Value_guess[iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap]));
//                         }
//                       }
//                       
//                       VV = u(cons, hhrent, ssigma, ppsi, ut) + bbeta*Vexn;
// 
//                       if(iy == 5 & ia == 5 & im == 5 & ih == 5 & cons >0 & hhrent>0){
//                         // cout << u(cons, hhrent, ssigma, ppsi, ut) << endl;
//                       }
//                       
//                       if(VV > VVn){
//                         VVn = VV;
//                       }
//                     }
//                   }
                }
              }
              
              // if((VVk > VVd) & (VVk > VVn)){
                Value[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = VVk;
                Default[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0;
                Renew[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0;
//               } else if((VVd > VVk) & (VVd > VVn)){
//                 Value[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = VVd;
//                 Default[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 1;
//                 Renew[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0;
//               } else{
//                 Value[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = VVn;
//                 Default[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0;
//                 Renew[iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 1;
//               }
//               if(iy == 5 & ia == 5 & im == 5 & ih == 5){
//                 cout << VVk << " " << VVn << " " << VVd << endl;
//               }
            }
          }
        }
      }
    }
    
#pragma omp barrier
    
    err = 0.0;
    
    for(int iy=0; iy<ny; iy++){
      for(int ia=0; ia<na; ia++){
        for(int im=0; im<nm; im++){
          for(int ih=0; ih<nh; ih++){
            for(int id=0; id<nd; id++){
              ind = iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              err = err + absol(Value[ind] - Value_guess[ind]);
              Value_guess[ind] = Value[ind];
              Value[ind] = 0.0;
            }
          }
        }
      }
    }
    
    t = omp_get_wtime() - t0;
    if (itnum % 1 == 0){
      cout << "     -------------------------------------------" << endl;
      cout << "FUNCTION: Iteration number: " << itnum << ". Time: " << 1000000*((float)t)/CLOCKS_PER_SEC << " seconds. Error " << err << endl;
      cout << "     -------------------------------------------" << endl;
    }
    itnum = itnum + 1;
  }

  for(int iy=0; iy<ny; iy++){
    for(int ia=0; ia<na; ia++){
      for(int im=0; im<nm; im++){
        for(int ih=0; ih<nh; ih++){
          for(int id=0; id<nd; id++){
            ind = iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
            
            Val[iy][ia][im][ih][id] = Value_guess[ind];
          }
        }
      }
    }
  }
  
  return(Val);
}
