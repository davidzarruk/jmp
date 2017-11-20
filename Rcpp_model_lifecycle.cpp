//
//  Created by David Zarruk Valencia on June, 2016.
//  Copyright (c) 2016 David Zarruk Valencia. All rights reserved.
//

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
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
  
  if(c <= 0 || h <= 0){
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



void surv(double* survival){
  
  std::ifstream theFile ("survival.csv");
  
  std::string line;
  std::vector<std::vector<std::string> > values;
  int it=0;
  std::string line_value;
  std::vector<std::string> line_values;
  std::stringstream ss;
  while(std::getline(theFile, line)){
    ss<<line;
    
    while(std::getline(ss, line_value, ',')){
      
      line_values.push_back(line_value);
      survival[it] = (double) 1 - ::atof(line_value.c_str());
    }
    values.push_back(line_values);
    it=it+1;
    line_value.clear();
    line_values.clear();
    ss.clear();
  }
}

void eproc(double* eprocess, const double T){
  
  int dirk = 0;
  
  if(dirk == 1){
    // Using data from Conesa, Kitao, Krueger that claim to use Hansen (1993)
    string line;
    vector<vector<string> > values;
    int it = 0;
    ifstream theFile2 ("incprofile.csv");
    
    string line_value;
    vector<string> line_values;
    stringstream ss;
    while(getline(theFile2, line)){
      ss<<line;
      
      while(getline(ss, line_value, ',')){
        
        line_values.push_back(line_value);
        eprocess[it] = (double) ::atof(line_value.c_str());
        
      }
      values.push_back(line_values);
      it=it+1;
      line_value.clear();
      line_values.clear();
      ss.clear();
    }
    
    for(int i = 45; i<T; i++){
      eprocess[i] = 0;
    }
  } else{
    // Using paper of Hansen (1993)
    for(int i = 0; i < T; i++){
      if(i<4){
        eprocess[i] = 0.78 + (i+1)*(1.14-0.78)/5;
      } 
      else if(i>=4 && i<14){
        eprocess[i] = 1.14 + (i-4)*(1.37-1.14)/10;
      } 
      else if(i>=14 && i<24){
        eprocess[i] = 1.37 + (i-14)*(1.39-1.37)/10;
      } 
      else if(i>=24 && i<34){
        eprocess[i] = 1.39 + (i-24)*(1.33-1.39)/10;
      } 
      else if(i>=34 && i<44){
        eprocess[i] = 1.33 + (i-34)*(0.89-1.33)/10;
      } 
      else if(i>=44 && i<49){
        // eprocess[i] = 0.89 + (i-44)*(0-0.89)/10;
        eprocess[i] = 0;
      } 
      else{
        eprocess[i] = 0;
      }
      // cout << eprocess[i] << endl;
    }
  }
}


void export_arrays(const int T,
                   const int ny,
                   const int nd,
                   const int na,
                   const int nh,
                   const int nm,
                   const double* Value,
                   const int* Default,
                   const int* Renew,
                   const int* Policym,
                   const int* Policya,
                   const int* Policyh,
                   const int* Policyr,
                   const double* Policyc,
                   const double* Pricing_guess){
  
  int ind; 
  
  ofstream Valuefile ("Value.txt");
  if (Valuefile.is_open())
  {
    Valuefile << T << "\n";
    Valuefile << na << "\n";
    Valuefile << ny << "\n";
    Valuefile << nd << "\n";
    Valuefile << nh << "\n";
    Valuefile << nm << "\n";
    
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Valuefile << Value[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Valuefile.close();
  }
  else cout << "Unable to open file";
  
  ofstream Defaultfile ("Default.txt");
  if (Defaultfile.is_open())
  {
    Defaultfile << T << "\n";
    Defaultfile << na << "\n";
    Defaultfile << ny << "\n";
    Defaultfile << nd << "\n";
    Defaultfile << nh << "\n";
    Defaultfile << nm << "\n";
    
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Defaultfile << Default[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Defaultfile.close();
  }
  else cout << "Unable to open file";
  
  ofstream Renewfile ("Renew.txt");
  if (Renewfile.is_open())
  {
    Renewfile << T << "\n";
    Renewfile << na << "\n";
    Renewfile << ny << "\n";
    Renewfile << nd << "\n";
    Renewfile << nh << "\n";
    Renewfile << nm << "\n";
    
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Renewfile << Renew[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Renewfile.close();
  }
  else cout << "Unable to open file";
  
  ofstream Policyhfile ("Policyh.txt");
  if (Policyhfile.is_open())
  {
    Policyhfile << T << "\n";
    Policyhfile << na << "\n";
    Policyhfile << ny << "\n";
    Policyhfile << nd << "\n";
    Policyhfile << nh << "\n";
    Policyhfile << nm << "\n";
    
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Policyhfile << Policyh[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Policyhfile.close();
  }
  else cout << "Unable to open file";
  
  ofstream Policymfile ("Policym.txt");
  if (Policymfile.is_open())
  {
    Policymfile << T << "\n";
    Policymfile << na << "\n";
    Policymfile << ny << "\n";
    Policymfile << nd << "\n";
    Policymfile << nh << "\n";
    Policymfile << nm << "\n";
    
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Policymfile << Policym[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Policymfile.close();
  }
  else cout << "Unable to open file";
  
  ofstream Policyafile ("Policya.txt");
  if (Policyafile.is_open())
  {
    Policyafile << T << "\n";
    Policyafile << na << "\n";
    Policyafile << ny << "\n";
    Policyafile << nd << "\n";
    Policyafile << nh << "\n";
    Policyafile << nm << "\n";
    
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Policyafile << Policya[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Policyafile.close();
  }
  else cout << "Unable to open file";
  
  ofstream Policyrfile ("Policyr.txt");
  if (Policyrfile.is_open())
  {
    Policyrfile << T << "\n";
    Policyrfile << na << "\n";
    Policyrfile << ny << "\n";
    Policyrfile << nd << "\n";
    Policyrfile << nh << "\n";
    Policyrfile << nm << "\n";
    
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Policyrfile << Policyr[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Policyrfile.close();
  }
  else cout << "Unable to open file";
  
  ofstream Policycfile ("Policyc.txt");
  if (Policycfile.is_open())
  {
    Policycfile << T << "\n";
    Policycfile << na << "\n";
    Policycfile << ny << "\n";
    Policycfile << nd << "\n";
    Policycfile << nh << "\n";
    Policycfile << nm << "\n";
    
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Policycfile << Policyc[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Policycfile.close();
  }
  else cout << "Unable to open file";
  
  ofstream Pricing_guessfile ("Pricing_guess.txt");
  if (Pricing_guessfile.is_open())
  {
    Pricing_guessfile << T << "\n";
    Pricing_guessfile << na << "\n";
    Pricing_guessfile << ny << "\n";
    Pricing_guessfile << nd << "\n";
    Pricing_guessfile << nh << "\n";
    Pricing_guessfile << nm << "\n";
    
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Pricing_guessfile << Pricing_guess[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Pricing_guessfile.close();
  }
  else cout << "Unable to open file";
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
             double r, int ny, int na, int nm, int nh, int nd, int T, vector<vector<double> > P,
             vector<double> dgrid, vector<double> mgrid, vector<double> hgrid, double* survival,
             int* Def,
             int* Renew,
             int* Policya,
             double Ph,
             double* pricing_guess){
  
  double survival_accum[T];
  survival_accum[0] = survival[0];
  for(int it=1; it<T; it++){
    survival_accum[it] = survival_accum[it-1]*survival[it];
  }

  double repay_coeff[T];
  for(int it=0; it<T; it++){
    repay_coeff[it] = 0.0;
    for(int itp=it; itp<T; itp++){
      repay_coeff[it] = repay_coeff[it] + survival_accum[it]*pow(rrho/(1+r),itp-it);
    }
  }

  // Initial guess for pricing function: price is equal to 1
  double *pricing;
  size_t sizemats = ny*na*nm*nh*nd*T*sizeof(double);
  
  pricing = (double*)malloc(sizemats);
  
  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              pricing[it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia] = 0.0;
            }
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
  int ind3;
  
  int itnum = 1;
  int iap = 0;
  
  // while(err > tol & itnum < maxiter){
  while(err > tol){
  
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                
                // For every state variable, I compute the pricing function
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                
                // Expected value is sum over tomorrow's possible shocks times probabilities: P[iy][iyp]*(1/nd)
                for(int iyp=0; iyp<ny; iyp++){
                  for(int idp=0; idp<nd; idp++){
                    
                    ind2 = it*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + ia;

                    iap = Policya[ind2];
                    ind3 = it*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + iap;
                    
                    pricing[ind] = pricing[ind] + (survival[it]*rrho/(mgrid[im]*(1+r)))*P[iy][iyp]*(1/(double)nd)*(Def[ind2]*(1-dgrid[idp])*Ph*hgrid[ih] + // If he defaults 
                                                                                                                  (1-Def[ind2])*((1-Renew[ind2])*(mgrid[im] + pricing_guess[ind3]*mgrid[im]) + // If he pays and continues with mortgage 
                                                                                                                                  Renew[ind2]*repay_coeff[it]*mgrid[im]));
                    
//                     if(iy == 2 && ia == 4 && im == 1 && ih == 2 && id == 1){
//                       // cout << pricing[ind] << "-:" << Def[ind2] << "-" << (1-dgrid[idp])*Ph*hgrid[ih] << "-" << Renew[ind2] << "-" << mgrid[im] + pricing_guess[ind3]*mgrid[im] << "-" << repay_coeff[it]*mgrid[im] << endl;
//                       cout << pricing[ind] << "-:" << Def[ind2] << "-" << (1-dgrid[idp])*Ph*hgrid[ih] << "-" << dgrid[idp] << "-" << Ph << "-" << hgrid[ih] << endl;
//                     }
                  }
                }
//                 if(iy == 0 && ia == 4 && im == 0 && ih == 1 && id == 0){
//                   cout << "-----" << endl;
//                 }
              }
            }
          }
        }
      }
    }
    // cout << "-----------" << endl;
    err = 0.0;
    
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                err = err + absol(pricing_guess[ind] - pricing[ind]);
                pricing_guess[ind] = pricing[ind];
                pricing[ind] = 0.0;
                
              }
            }
          }
        }
      }
    }
    
     if (itnum % 1 == 0){
       cout << "     Pricing Function - It. #: " << itnum << ". Error " << err << endl;
     }
    itnum++;
  }
  
  // cout << "     -------------------------------------------" << endl;
  // cout << "     Pricing Function - It. #: " << itnum << ". Error " << err << endl;
  // cout << "     -------------------------------------------" << endl;
}

// Supply of housing - linear supply function
double housing_supply(double slope, double hPrice){
  
  double sup = (hPrice - (1-0.69*slope))/slope;
  return(sup);
  
}

void setgrids(int T, double* survival, double* eprocess){
  
  // Survival probabilities
  double survival_long[71];
  surv(survival_long);
  
  for(int it=0; it<T; it++){
    survival[it] = 1;
    for(int i=0; i<71; i++){
      if((i>=it*5) & (i<(it+1)*5)){
        survival[it] = survival[it]*survival_long[i];
      }
    }
  }
  survival[T-1] = 0;
  
  // Income process over life cycle
  double eprocess_long[71];
  eproc(eprocess_long, 71);
  
  for(int it=0; it<T; it++){
    eprocess[it] = 0;
    for(int i=0; i<71; i++){
      if((i>=it*5) & (i<(it+1)*5)){
        eprocess[it] = eprocess[it] + eprocess_long[i]/5;
      }
    }
    if(it>T-4){
      eprocess[it] = eprocess[T-4]*0.5;
    }
  }
  
}

//======================================
//         Value Function Iteration
//======================================

// [[Rcpp::export]]
vector<vector<vector<vector<vector<vector<vector<double> > > > > > > vfi(double tol, int maxiter, int ut, int T,
                                                                          int na, double amin, double amax,
                                                                          int nm, double mmin, double mmax,
                                                                          int nh, double hmin, double hmax, 
                                                                          int ny, double ssigma_y, double llambda_y, double m, 
                                                                          int nd, double dmin, double dmax,
                                                                          double ssigma, double ppsi, double rrho, double bbeta,
                                                                          double r, double Ph, double q, double Pa, double ddeltabar, double fcost) {

  // Initialization of output element: All - contains policies, value, pricing, etc
  vector<vector<vector<vector<vector<vector<vector<double> > > > > > > All;
  All.resize(T);
  
  for(int it=0; it<T; it++){
    All[it].resize(ny);
      for(int iy=0; iy<ny; iy++){
        All[it][iy].resize(na);
      for(int ia=0; ia<na; ia++){
        All[it][iy][ia].resize(nm);
        for(int im=0; im<nm; im++){
          All[it][iy][ia][im].resize(nh);
          for(int ih=0; ih<nh; ih++){
            All[it][iy][ia][im][ih].resize(nd);
            for(int id=0; id<nd; id++){
              All[it][iy][ia][im][ih][id].resize(35);
            }
          }
        }
      }
    }
  }
  
  double *Value, *Policyc, *Pricing_guess, *Pricing2;
  int *Default, *Renew, *Policya, *Policym, *Policyh, *Policyr;

  size_t sizemats = ny*na*nm*nh*nd*T*sizeof(double);
  size_t sizematsint = ny*na*nm*nh*nd*T*sizeof(int);
  
  Value         = (double*)malloc(sizemats);
  Policyc       = (double*)malloc(sizemats);
  Pricing_guess = (double*)malloc(sizemats);
  Pricing2      = (double*)malloc(sizemats);
  Default       = (int*)malloc(sizematsint);
  Renew         = (int*)malloc(sizematsint);
  Policya       = (int*)malloc(sizematsint);
  Policym       = (int*)malloc(sizematsint);
  Policyh       = (int*)malloc(sizematsint);
  Policyr       = (int*)malloc(sizematsint);
  
  int ind;
  int ind1;
  int ind2;
  
  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              Value[ind]        = 0.0;
              Policyc[ind]      = 0.0;
              Pricing_guess[ind] = 1.0;
              Pricing2[ind]     = 1.0;
              Default[ind]      = 0;
              Renew[ind]        = 0;
              Policya[ind]      = 0;
              Policym[ind]      = 0;
              Policyh[ind]      = 0;
              Policyr[ind]      = 0;
            }
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
  
  // Survival probabilities
  double survival[T];
  double eprocess[T];
  setgrids(T, survival, eprocess);

  double survival_accum[T];
  survival_accum[0] = survival[0];
  for(int it=1; it<T; it++){
    survival_accum[it] = survival_accum[it-1]*survival[it];
  }
  
  double repay_coeff[T];
  for(int it=0; it<T; it++){
    repay_coeff[it] = 0.0;
    for(int itp=it; itp<T; itp++){
      repay_coeff[it] = repay_coeff[it] + survival_accum[it]*pow(rrho/(1+r),itp-it);
    }
  }
  
  
  // VFI starts
  double err   = 1;
  int itnum = 1;
  double t0;
  double t;

  double VV     = 0.0;

  double VVk    = pow(-10,11);  double VVn    = pow(-10,11);  double VVd    = pow(-10,11);
  
  double hhk    = 0;  double hhn    = 0;  double hhd    = 0;    // Home ownership
  double hrk    = 0;  double hrn    = 0;  double hrd    = 0;    // Home renting
  double mmk    = 0;  double mmn    = 0;  double mmd    = 0;    // Mortgage
  double aak    = 0;  double aan    = 0;  double aad    = 0;    // Savings
  double cck    = 0;  double ccn    = 0;  double ccd    = 0;    // Consumption
  
  double Vexk   = 0.0;  double Vexn   = 0.0;  double Vexd   = 0.0;
  
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

  Pmort(rrho, r, ny, na, nm, nh, nd, T, P, dgrid, mgrid, hgrid, survival, Default, Renew, Policya, Ph, Pricing_guess);

  while(err > tol & itnum < maxiter){
    
    // State variables
    for(int it=T-1; it>=0; it--){
      
      // if (it % 1 == 0){
      //   cout << "Age is:" << it+1 << "... " << endl;
      // }
      cout << it+1 << ", " << flush;
      
      // Parallelization  
      // #pragma omp parallel for shared(it, All, Pricing_guess, Pricing2, Value, Default, Renew, Policya, Policyh, Policyr, Policym, Policyc, P, agrid, hgrid, mgrid, dgrid, ygrid, nd, ny, na, nh, nm, rrho, bbeta, Pa, Ph, ddeltabar, q, r, err, tol, itnum, maxiter, t, t0, ssigma, ppsi, ut, T, fcost, survival, repay_coeff, survival_accum) private(ind, ind1, ind2, yy, aa, mm, hh, ddelta, aaprime, hhrent, cons, Vexk, Vexd, Vexn, VV, mmprime, hhprime, pprice, VVk, hhk, hrk, mmk, aak, cck, VVd, hhd, hrd, mmd, aad, ccd, VVn, hhn, hrn, mmn, aan, ccn)
      for(int iy=0; iy<ny; iy++){
        for(int ia=0; ia<na; ia++){
          
          for(int im=0; im<nm; im++){
            for(int ih=0; ih<nh; ih++){
              for(int id=0; id<nd; id++){
                
                VVk        = pow(-10,11);
                VVn        = pow(-10,11);
                VVd        = pow(-10,11);
                
                yy      = ygrid[iy]*eprocess[it];
                aa      = agrid[ia];
                mm      = mgrid[im];
                hh      = hgrid[ih];
                ddelta  = dgrid[id];
                
                // Control variables
                for(int iap=0; iap<na; iap++){
                  for(int ihre=0; ihre<nh; ihre++){
                    aaprime = agrid[iap];
                    hhrent  = hgrid[ihre];
                    
                    // Keeping the same mortgage
                    cons = aa + q*hh + yy - mm - q*hhrent - Pa*aaprime - Ph*ddelta*ddeltabar*hh;
                    
                    Vexk       = 0.0;
                    if(it < T-1){
                      for(int idp=0; idp<nd; idp++){
                        if(it<T-3){  // Income uncertainty before retirement
                          for(int iyp=0; iyp<ny; iyp++){
                            
                            ind1 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + iap;
                            ind2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + 0*na + iap;
                            
                            Vexk = Vexk + P[iy][iyp]*(1/(double)nd)*(rrho*Value[ind1] +  // Keeps mortg
                                                                    (1-rrho)*Value[ind2]);    // Mortg disappears
                          }
                        } else{   // Certainty after retirement
                          
                          ind1 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + im*na + iap;
                          ind2 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ih*nm*na + 0*na + iap;
                          
                          Vexk = Vexk + (1/(double)nd)*(rrho*Value[ind1] +  // Keeps mortg
                                                        (1-rrho)*Value[ind2]);    // Mortg disappears
                        }
                      }
                    }
                    
                    VV = u(cons, hhrent, ssigma, ppsi, ut) + bbeta*survival[it]*Vexk;
                    
                    if(VV > VVk){
                      VVk = VV;
                      hhk = ih;
                      hrk = ihre;
                      mmk = im;
                      aak = iap;
                      cck = cons;
                    }
                    
                    // Defaulting
                    cons = aa + yy - q*hhrent - Pa*aaprime;
                    
                    Vexd       = 0.0;
                    if(it < T-1){
                      for(int idp=0; idp<nd; idp++){
                        if(it<T-3){  // Income uncertainty before retirement
                          for(int iyp=0; iyp<ny; iyp++){
                            ind1 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + iap;
  
                            Vexd = Vexd + P[iy][iyp]*(1/(double)nd)*Value[ind1];
                          }
                        } else{   // Certainty after retirement
                          ind1 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + iap;
                          
                          Vexd = Vexd + (1/(double)nd)*Value[ind1];
                        }
                      }
                    }
                    
                    VV = u(cons, hhrent, ssigma, ppsi, ut) + bbeta*survival[it]*Vexd;
                    
                    if(VV > VVd){
                      VVd = VV;
                      hhd = 0;
                      hrd = ihre;
                      mmd = 0;
                      aad = iap;
                      ccd = cons;
                    }
                      
                    // New mortgage
                    for(int imp=0; imp<nm; imp++){
                      for(int ihp=0; ihp<nh; ihp++){
                        
                        mmprime      = mgrid[imp];
                        hhprime      = hgrid[ihp];
                        
                        ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ihp*nm*na + imp*na + iap;
                        pprice = Pricing_guess[ind];
                        
                        cons = aa + Ph*(1-ddelta*ddeltabar)*hh + q*hhprime + yy - fcost + mmprime*pprice - repay_coeff[it]*mm - q*hhrent - Ph*hhprime - Pa*aaprime;
                        // cons = aa + Ph*(1-ddelta*ddeltabar)*hh + q*hhprime + yy - fcost*pow((hhprime-hh), 2) + mmprime*pprice - repay_coeff[it]*mm - q*hhrent - Ph*hhprime - Pa*aaprime;
                        
                        Vexn       = 0.0;
                        if(it < T-1){
                          for(int idp=0; idp<nd; idp++){
                            if(it<T-3){  // Income uncertainty before retirement
                              for(int iyp=0; iyp<ny; iyp++){
                                
                                ind1 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                                ind2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;
                                
                                Vexn = Vexn + P[iy][iyp]*(1/(double)nd)*((rrho * Value[ind1]) + 
                                                                        ((1-rrho) * Value[ind2]));
                              }
                            } else{   // Certainty after retirement
                              ind1 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                              ind2 = (it+1)*ny*nd*nh*nm*na + iy*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;
                              
                              Vexn = Vexn + (1/(double)nd)*((rrho * Value[ind1]) + 
                                                            ((1-rrho) * Value[ind2]));
                            }
                          }
                        }
                        
                        VV = u(cons, hhrent, ssigma, ppsi, ut) + bbeta*survival[it]*Vexn;
  
                        if(VV > VVn){
                          VVn = VV;
                          hhn = ihp;
                          hrn = ihre;
                          mmn = imp;
                          aan = iap;
                          ccn = cons;
                        }
                      }
                    }
                  }
                }
                
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  
                if((VVk > VVd) & (VVk > VVn)){
                  Value[ind]    = VVk;
                  Policya[ind]    = aak;
                  Policyh[ind]    = hhk;
                  Policyr[ind]    = hrk;
                  Policym[ind]    = mmk;
                  Policyc[ind]    = cck;
                  Default[ind]  = 0;
                  Renew[ind]    = 0;
                } else if((VVd > VVk) & (VVd > VVn)){
                  Value[ind] = VVd;
                  Policya[ind]    = aad;
                  Policyh[ind]    = hhd;
                  Policyr[ind]    = hrd;
                  Policym[ind]    = mmd;
                  Policyc[ind]    = ccd;
                  Default[ind] = 1;
                  Renew[ind] = 0;
                } else{
                  Value[ind] = VVn;
                  Policya[ind]    = aan;
                  Policyh[ind]    = hhn;
                  Policyr[ind]    = hrn;
                  Policym[ind]    = mmn;
                  Policyc[ind]    = ccn;
                  Default[ind] = 0;
                  Renew[ind] = 1;
                }
              }
            }
          }
        }
      }
    }
    
    // #pragma omp barrier
    
    Pmort(rrho, r, ny, na, nm, nh, nd, T, P, dgrid, mgrid, hgrid, survival, Default, Renew, Policya, Ph, Pricing2);
    
    err = 0.0;
    
    for(int it=T-1; it>=0; it--){
      for(int iy=0; iy<ny; iy++){
        for(int ia=0; ia<na; ia++){
          for(int im=0; im<nm; im++){
            for(int ih=0; ih<nh; ih++){
              for(int id=0; id<nd; id++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                
                err = err + absol(Pricing2[ind] - Pricing_guess[ind])/sizemats;
                // Pricing_guess[ind] = (Pricing2[ind] + Pricing_guess[ind])/2;
                Pricing_guess[ind] = Pricing2[ind];
                Pricing2[ind] = 0.0;
              }
            }
          }
        }
      }
    }
    
    t = omp_get_wtime() - t0;
    // if (itnum % 1 == 0){
    //   cout << "     -------------------------------------------" << endl;
    //   cout << "FUNCTION: Iteration number: " << itnum << ". Time: " << 1000000*((double)t)/CLOCKS_PER_SEC << " seconds. Error " << err << endl;
    //   cout << "     -------------------------------------------" << endl;
    // }
    if (itnum % 1 == 0){
      // cout << "     -------------------------------------------" << endl;
      // cout << "FUNCTION: Iteration number: " << itnum << ". Time: " << 1000000*((double)t)/CLOCKS_PER_SEC << " seconds. Error " << err << endl;
      cout << "  It: " << itnum << ". Time: " << ((double)t/60)/CLOCKS_PER_SEC << " minutes. Error " << err << endl;
      // cout << "     -------------------------------------------" << endl;
    }
    itnum = itnum + 1;
  }
  
  // I export the relevant matrices: Value, policies, pricing, etc.
  for(int it=T-1; it>=0; it--){
    for(int iy=0; iy<ny; iy++){
      for(int ia=0; ia<na; ia++){
        for(int im=0; im<nm; im++){
          for(int ih=0; ih<nh; ih++){
            for(int id=0; id<nd; id++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              All[it][iy][ia][im][ih][id][0] = Pricing_guess[ind];    // Pricing function
              All[it][iy][ia][im][ih][id][1] = Value[ind];            // Value function
              All[it][iy][ia][im][ih][id][2] = (double)Policya[ind];            // Value function
              All[it][iy][ia][im][ih][id][3] = (double)Policyh[ind];            // Value function
              All[it][iy][ia][im][ih][id][4] = (double)Policyr[ind];            // Value function
              All[it][iy][ia][im][ih][id][5] = (double)Policym[ind];            // Value function
              All[it][iy][ia][im][ih][id][6] = Policyc[ind];            // Value function
              All[it][iy][ia][im][ih][id][7] = (double)Default[ind];          // Default Policy
              All[it][iy][ia][im][ih][id][8] = (double)Renew[ind];            // Renew Policy
            }
          }
        }
      }
    }
  }
  
  //--------------------------------------//
  //--------  Exporting arrays  ----------//
  //--------------------------------------//
  

  export_arrays(T, ny, nd, na, nh, nm,
                Value, Default, Renew, 
                Policym, Policya, 
                Policyh, Policyr, Policyc, Pricing_guess);
  
//-----------------------------------------------------------//
//-----------------  Probability Paths  ---------------------//
//-----------------------------------------------------------//

  // Matrices
  double *Pcond, *Puncond;

  Pcond       = (double*)malloc(sizemats);
  Puncond         = (double*)malloc(sizemats);

  int indp;
  int indp2;
  
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
  

  // Steady state distribution for y
  vector<double> steady_distr;
  vector<vector<double> > Paux;
  vector<vector<double> > Paux2;
  
  steady_distr.resize(ny);
  Paux.resize(ny);
  Paux2.resize(ny);
  for(int iy = 0; iy<ny; iy++){
    Paux[iy].resize(ny);
    Paux2[iy].resize(ny);
    for(int iyp = 0; iyp<ny; iyp++){
      Paux[iy][iyp]=P[iy][iyp];
      
      Paux2[iy][iyp]=0.0;
    }
  }

  for(int itt = 0; itt<10000; itt++){

    for(int iy=0; iy<ny; iy++){
      for(int iyp=0; iyp<ny; iyp++){
        
        Paux2[iy][iyp]=0.0;
        
        for(int is=0; is<ny; is++){
          Paux2[iy][iyp] = Paux2[iy][iyp] + P[iy][is]*Paux[is][iyp];
        }
      }
    }

    for(int iy=0; iy<ny; iy++){
      for(int iyp=0; iyp<ny; iyp++){
        Paux[iy][iyp] = Paux2[iy][iyp];
      }
    }
  }

  for(int iy=0; iy<ny; iy++){
    steady_distr[iy] = Paux[0][iy];
  }

  //--------------------------------------//
  //----   Distribution computation  -----//
  //--------------------------------------//
  

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
                          indp = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
                          indp2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;
                          
                          if(imp == Policym[ind] && ihp == Policyh[ind] && iap == Policya[ind]){
                            Pcond[indp]   = Pcond[indp]   + rrho*P[iy][iyp]*(1/(double)nd)*Pcond[ind];
                            Puncond[indp] = Puncond[indp] + rrho*P[iy][iyp]*(1/(double)nd)*Puncond[ind]*survival[it];

                            Pcond[indp2]   = Pcond[indp2]   + (1-rrho)*P[iy][iyp]*(1/(double)nd)*Pcond[ind];
                            Puncond[indp2] = Puncond[indp2] + (1-rrho)*P[iy][iyp]*(1/(double)nd)*Puncond[ind]*survival[it];
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

  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              All[it][iy][ia][im][ih][id][9] = Pcond[ind];
              All[it][iy][ia][im][ih][id][10] = Puncond[ind];
            }
          }
        }
      }
    }
  }

  
  //--------------------------------------//
  //----  Aggregate paths in economy  ----//
  //--------------------------------------//

  // Lifecycle Average Paths
  double C[T];
  double H[T];
  double Hp[T];
  double M[T];
  double Mp[T];
  double R[T];
  double NR[T];         // % of net renters
  double NetRent[T];    // Net rent
  double OwnerOccupied[T];
  double A[T];
  double PP[T];
  double RR[T];
  double Def[T];
  double Vivos[T];
  
  // Aggregate Statistics
  double people             = 0.0;
  double owners_w_mortgage  = 0.0;
  double mortgage_debt      = 0.0;
  double housing_value      = 0.0;
  double default_rate       = 0.0;
  double housing_total      = 0.0;
  double rent_total         = 0.0;
  
  for(int it=0; it<T; it++){
    C[it]   = 0.0;
    H[it]   = 0.0;
    Hp[it]  = 0.0;
    M[it]   = 0.0;
    Mp[it]  = 0.0;
    R[it]   = 0.0;
    NR[it]  = 0.0;    // Net renter
    NetRent[it]   = 0.0;
    OwnerOccupied[it]   = 0.0;
    A[it]   = 0.0;
    PP[it]  = 0.0;
    RR[it]  = 0.0;
    Def[it] = 0.0;
    
    Vivos[it] = 0.0;
  }
  
  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              // Means:
              C[it] = C[it] + Pcond[ind]*Policyc[ind];
              H[it] = H[it] + Pcond[ind]*hgrid[Policyh[ind]];
              M[it] = M[it] + Pcond[ind]*mgrid[Policym[ind]];
              R[it] = R[it] + Pcond[ind]*hgrid[Policyr[ind]];
              A[it] = A[it] + Pcond[ind]*agrid[Policya[ind]];
              PP[it] = PP[it] + Pcond[ind]*Pricing_guess[ind];
              
              // If rents more than the size owned is a "net renter"
              if(Policyr[ind] > Policyh[ind]){
                NR[it] = NR[it] + Pcond[ind];
                NetRent[it] = NetRent[it] + Pcond[ind]*(hgrid[Policyr[ind]] - hgrid[Policyh[ind]]);
                OwnerOccupied[it] = OwnerOccupied[it] + Pcond[ind]*(hgrid[Policyh[ind]]);
              } else{
                OwnerOccupied[it] = OwnerOccupied[it] + Pcond[ind]*(hgrid[Policyr[ind]]);
              }

              if(ih>0){
                Hp[it] = Hp[it] + Pcond[ind];
              }
              
              if(im > 0){
                RR[it] = RR[it] + Pcond[ind]*(double)Renew[ind];
                Def[it] = Def[it] + Pcond[ind]*(double)Default[ind];
                Mp[it] = Mp[it] + Pcond[ind];
                default_rate = default_rate + Puncond[ind]*(double)Default[ind];
              }

              people = people + Puncond[ind];
              if(Policym[ind] > 0 && Policyh[ind] > 0){
                owners_w_mortgage = owners_w_mortgage + Puncond[ind];
              }
              
              ind2 = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + Policyh[ind]*nm*na + Policym[ind]*na + Policya[ind];

              mortgage_debt = mortgage_debt + Puncond[ind]*mgrid[Policym[ind]]*Pricing_guess[ind2];
              housing_value = housing_value + Puncond[ind]*hgrid[Policyh[ind]]*Ph;
              
              housing_total = housing_total + Puncond[ind]*hgrid[Policyh[ind]];
              rent_total = rent_total + Puncond[ind]*hgrid[Policyr[ind]];
              
              Vivos[it] = Vivos[it] + Puncond[ind];
            }
          }
        }
      }
    }
    
    // Paths
    All[it][0][0][0][0][0][11] = C[it];
    All[it][0][0][0][0][0][12] = H[it];
    All[it][0][0][0][0][0][13] = M[it];
    All[it][0][0][0][0][0][14] = R[it];
    All[it][0][0][0][0][0][15] = A[it];
    All[it][0][0][0][0][0][16] = Vivos[it];
    All[it][0][0][0][0][0][17] = Def[it];
    All[it][0][0][0][0][0][18] = RR[it];
    All[it][0][0][0][0][0][19] = PP[it];
    All[it][0][0][0][0][0][20] = Mp[it];
    All[it][0][0][0][0][0][21] = Hp[it];
    All[it][0][0][0][0][0][22] = NR[it];
    All[it][0][0][0][0][0][23] = NetRent[it];
    All[it][0][0][0][0][0][24] = OwnerOccupied[it];
    
  }
  
  // Aggregate Statistics
  All[0][0][0][0][0][0][34] = owners_w_mortgage;
  All[1][0][0][0][0][0][34] = mortgage_debt;
  All[2][0][0][0][0][0][34] = housing_value;
  All[3][0][0][0][0][0][34] = people;
  All[4][0][0][0][0][0][34] = default_rate;
  All[5][0][0][0][0][0][34] = housing_total;
  All[6][0][0][0][0][0][34] = rent_total;
  
  return(All);
}

