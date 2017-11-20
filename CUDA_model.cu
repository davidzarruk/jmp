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
#include <nlopt.hpp>
using std::vector;
using namespace std;

//****************************************************//
//            1. Parameters                  //
//****************************************************//


class parameters{
  public:
	int maxiter;
  int uti;
  double tol;
  double convergence;
  int T;
  int Tretirement;
  int yearspp;

  // Grid for savings: a
  int na;
  double amin;
  double amax;

  // Grid for mortgages: m
  int nm;
  double mmin;
  double mmax;

  // Grid for housing: h
  int nh;
  double hmin;
  double hmax;

  // Grid for renting: r
  int nr;
  double rmin;
  double rmax;

  // Grid for labor: l
  int nl;
  double lmax;

  // Grid for deoreciation: ddelta
  int nd;
  double dmin;
  double dmax;

  // Grid for income shocks: y
  int ny;
  double ssigma_y;
  double llambda_y;
  double m_y;

  // Preferences
  double ssigma;
  double rrho;
  double ppsi;
  double bbeta;
  double kkappa;

  double tthetalab;
  double eetalab;

  double oomega;
  double sunk;
  double interm;
  double rec_probab;

  // Equilibrium objects
  double ddeltabar_today;
  double ddeltabar_tomorrow;
  double ddeltaf;
  double r;
  double Ph_today;
  double Ph_tomorrow;
  double q;
  double Pa;
  double housing_supply;
  double fcost;
  double refcost;
  double pension;
  double sstax;
  double ltax;
  double lumpsum;
  double Atech;

  int compute_equivalent;

  double multiplier;

  double *d_rental;
  double *d_housing;

	void load(const char*);
};


//****************************************************//
//            1. Import functions                  //
//****************************************************//


#include "Colormod.h" // namespace Color
#include "grid_initialization.cpp"
#include "export_arrays.cpp"
#include "CUDA_functions.cu"
#include "main_functions.cpp"
#include "Aggregation_functions.cpp"
#include "bank_aggregates.cpp"
#include "transitions.cpp"
#include "steady_state.cpp"
#include "maximizaciones.cpp"
#include "transition_max.cpp"
#include "import_arrays.cpp"



//======================================
//         Value Function Iteration
//======================================


double mindosnum(const double a, const double b){
  double res = a;
  if(b<a){
    res = b;
  }
  return(res);
}


double maxdosnum(const double a, const double b){
  double res = a;
  if(b>a){
    res = b;
  }
  return(res);
}


int main(int argc, char *argv[])
{

	cout.precision(6);
	cout.setf(std::ios::fixed);


  //----------------------------------------------//
  //---------    SETTING PARAMETERS    -----------//
  //----------------------------------------------//

  // -------- STEADY STATE PARAMETERS   ----------//

  // VFI parameters
  double tol         = pow(10,-4.0);
  const int uti            = 1;
  int maxiter        = 10;
  const double convergence = 0.7; // Entre mas alto, mas lenta es la convergencia

  // Years per period
  double yearspp = 5;

  // Demographics
  int T;
  int Tretirement;

  T           = 30;
  Tretirement = 23;

  // Grid for savings: a
  int na       = 35;  // 25
  const double amin  = 0;
  double amax  = 1.5;

  // Grid for mortgages: m
  int nm       = 5;  // 11
  const double mmin  = 0.0;
  double mmax  = 0.5;

  // Grid for housing: h
  int nh       = 3;
  const double hmin  = 0.0;
  double hmax  = 3.5;

  // Grid for renting: r
  const int nr       = 9;
  const double rmin  = 0.0001;
  double rmax  = 3.5;

  // Grid for labor: l
  int nl        = 3;
  double lmax   = 0.4;

  // Grid for depreciation: ddelta
  int nd             = 3;
  double ssigma_eps  = 0.05;

  // Grid for income shocks: y
  int ny                  = 3;
  const double ssigma_y   = 0.1*pow(yearspp, 0.5); // 0.44
  const double llambda_y  = pow(0.98, yearspp); // pow(0.98, yearspp) ; pow(0.95, yearspp)
  double m_y              = 1; // 2

  // Preferences
  const double ssigma  = 2;
  double rrho          = 0.8;
  double ppsi          = 0.65;
  double bbeta         = pow(0.965445931, yearspp);
  const double kkappa  = -0.1;

  double Atech         = 2.0;

  double tthetalab     = 30.0;    // So steady state labor is around 0.36
  double eetalab       = 2.0;     // Rogerson - macro Frisch elasticity

  double oomega        = 0.2;
  double sunk          = 0.0;
  double interm        = 0.0;

  // Equilibrium objects
  double r              = pow(1.03, yearspp)-1; // Average of 10yr Treasury bond for 2003-2007: 2.07%
  double Pa             = 1/(1+r);

  // Optimizer initial values
  double Ph             = 1.0;
  double q              = 0.236237;
  double dmin           = -0.050488;
  double dmax           = 0.150493;
  double ddeltabar      = 0.038;
  double ddeltaf        = 1.0;
  double fcost          = 0.153209;
  double refcost        = 0.05;
  double housing_supply = 21.5;
  double pension        = 0.6;
  double sstax          = 0.1;
  double ltax           = 0.0;
  double rec_probab     = 0.0;

  double lumpsum        = 0.0;

  // Transitional dynamics' parameters
  int     Ttrans          = 0;
  int     periods_shock   = 0;
  int     periods_taxes   = periods_shock;

  double  rshock          = r;
  double  Pashock         = Pa;
  double  ppsishock       = ppsi;
  double  oomegashock     = oomega;
  double  sunk_shock      = sunk;
  double  interm_shock    = interm;
  double  ltaxshock       = ltax;
  double  Phinf;
  double  ddeltabarinf;
  double  qinf;
  double  interminf;

  double prob_mistake = 0.0;


  int permanent       = 0;
  int baseline        = 1;

  double multiplier       = 0.008;
  int compute_equivalent  = 0;

  double *Phseq, *Phnoextseq, *qseq, *ddeltabarseq, *lumpsumseq, *mort_subsidy, *qnoextseq, *incshock;
  int *subs_eligible, *subs_target;
  size_t sizemats, sizematssubs;
  int ind;

  // Inputs to see what to compute: initial ss, transition dynamics, etc.
  std::string stage = argv[1];
  std::string tipo  = "tr";


  int maxim = 3;

  const int experimento = 35;


  tol         = pow(10,-4.0);
  // maxiter = 20;

  // Every period: 2 years
  yearspp         = 2;
  T               = 30;
  Tretirement     = 22.5;

  Atech           = 0.8;

  // Grids
  ny              = 3;
  nd              = 3;
  na              = 11;
  nl              = 3;
  lmax            = 0.4;
  // nh = 7;
  // dmax = 0.5;

  mmax            = 0.4;
  hmax            = 3;
  rmax            = hmax;
  amax            = 1.0;

  Ph              = 1.0;
  r               = pow(1.0207, yearspp)-1; // Average of 10yr Treasury bond for 2003-2007: 2.07%
  Pa              = 1/(1+r);

  ppsi            = 0.84;
  bbeta           = pow(0.966499, yearspp);
  rrho            = 1 - yearspp/25;

  ssigma_eps      = 0.04;

  q               = 0.043131;
  dmin            = -0.349340;
  dmax            = 0.459476;
  fcost           = 0.015;
  refcost         = 0.025;

  // fcost           = 0.0;
  // refcost         = 0.0;

  ddeltabar       = 0.00;

  oomega          = -0.25;    // Esto es como Campbell y Cocco

  sunk            = 0.22;  // Mitman (pg18) de Pennington and Cross: 22%, Chatterjee tiene 17%

  interm          = pow(1.000, yearspp)-1;
  rec_probab      = 0.0;

  tthetalab       = 5.0;
  eetalab         = 0.5;

  sstax           = 0.1;
  ltax            = 0.0;

  // ltax            = 0.0035;

  pension         = 0.104;

  lumpsum         = -0.036;

  maxim = 9;

  q               = 0.047873;
  dmin            = -0.345419;
  dmax            = 0.464803;
  bbeta           = pow(0.96, yearspp);
  // r               = pow(1.023, yearspp)-1;



  // -------- TRANSITIONAL DYNAMICS PARAMETERS   ----------//

  housing_supply  = 40.65;
  ddeltabar       = 0.0;


  // TRANSITIONAL DYNAMICS - SHOCKS

  // Final steady state - same as initial
  interminf       = interm;
  ddeltabarinf    = ddeltabar;
  Phinf           = Ph;
  qinf            = q;


  // Shocking interest rate
  Ttrans          = 10;
  // rshock          = pow(1.05, yearspp)-1;
  rshock          = r;
  Pashock         = Pa;
  ppsishock       = ppsi;

  oomegashock     = 0.25;
  // oomegashock     = oomega;

  sunk_shock      = sunk;
  interm_shock    = pow(1.000, yearspp)-1;
  Pashock         = Pa;

  ltaxshock       = 0.00;
  // ltaxshock       = 0.000;

  periods_shock   = 4;
  periods_taxes   = 1;
  permanent       = 0;
  
  // Transitional dynamics objects
  sizemats        = Ttrans*sizeof(double);
  Phseq           = (double*)malloc(sizemats);
  Phnoextseq      = (double*)malloc(sizemats);
  qseq            = (double*)malloc(sizemats);
  qnoextseq       = (double*)malloc(sizemats);
  ddeltabarseq    = (double*)malloc(sizemats);
  lumpsumseq      = (double*)malloc(sizemats);

  // Subsidies only in the first period after shock
  sizematssubs    = ny*na*nm*nh*T*sizeof(double);
  mort_subsidy    = (double*)malloc(sizematssubs);

  sizematssubs    = ny*na*nm*nh*nd*T*sizeof(int);
  subs_eligible   = (int*)malloc(sizematssubs);
  subs_target     = (int*)malloc(sizematssubs);

  sizemats        = T*sizeof(double);
  incshock        = (double*)malloc(sizemats);

  // Ph
  for(int it=0; it<Ttrans; it++){
    Phseq[it]         = Ph;
    ddeltabarseq[it]  = ddeltabar;
    qseq[it]          = qinf;
    qnoextseq[it]          = qinf;
    lumpsumseq[it]    = lumpsum;
  }

  lumpsumseq[0] = -0.043;

  for(int it=0; it<T; it++){
    if(it >= 0 & it < 10/yearspp){
      incshock[it] = 0.128;
    } else if(it >= 10/yearspp && it < 20/yearspp){
      incshock[it] = 0.111;
    } else if(it >= 20/yearspp && it < 30/yearspp){
      incshock[it] = 0.088;
    } else if(it >= 30/yearspp && it < 40/yearspp){
      incshock[it] = 0.096;
    } else if(it >= 40/yearspp && it < 45/yearspp){
      incshock[it] = 0.044;
    } else{
      incshock[it] = 0.0;
    }

    incshock[it] = incshock[it]*1.0;
    // incshock[it] = incshock[it]*0.0;
  }

  Phseq[0] = 0.765;
  Phseq[1] = 0.815886;
  Phseq[2] = 0.876203;
  Phseq[3] = 0.993659;
  Phseq[4] = 0.992694;
  Phseq[5] = 0.992827;
  Phseq[6] = 0.995083;
  Phseq[7] = 0.996250;
  Phseq[8] = 0.999546;
  Phseq[9] = 0.999950;

  qnoextseq[0]  = q - 0.0057;
  qnoextseq[1]  = q - 0.0045;
  qnoextseq[2]  = q - 0.0028;
  qnoextseq[3]  = q - 0.00015;
  qnoextseq[4]  = q - 0.0001;
  qnoextseq[5]  = q - 0.00005;
  qnoextseq[6]  = q - 0.000;
  qnoextseq[7]  = q - 0.000;
  qnoextseq[8]  = q - 0.000;
  qnoextseq[9]  = q - 0.000;


  // Mortgage subsidy
  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int ih=0; ih<nh; ih++){
        for(int im=0; im<nm; im++){
          for(int ia=0; ia<na; ia++){
            ind = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

            mort_subsidy[ind]   = 0.0;

            for(int id=0; id<nd; id++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

              subs_eligible[ind]  = 0;
              subs_target[ind]    = 0;
            }
          }
        }
      }
    }
  }


  double *Value, *Value_future, *Policyc, *Pricing_guess, *Pcond, *Puncond;
  int *Default, *Default0, *Renew, *Policya, *Policym, *Policyh, *Policyr, *Policyl, *Changer;

  sizemats     = ny*na*nm*nh*nd*T*sizeof(double);
  size_t sizematsint  = ny*na*nm*nh*nd*T*sizeof(int);

  size_t sizematsintinfo  = ny*na*nm*nh*T*sizeof(int);
  
  Value         = (double*)malloc(sizemats);
  Value_future  = (double*)malloc(sizemats);
  Policyc       = (double*)malloc(sizemats);
  Pricing_guess = (double*)malloc(sizemats);
  Puncond       = (double*)malloc(sizemats);
  Pcond         = (double*)malloc(sizemats);

  Default0      = (int*)malloc(sizematsint);
  Default       = (int*)malloc(sizematsint);
  Renew         = (int*)malloc(sizematsint);
  Policya       = (int*)malloc(sizematsint);
  Policym       = (int*)malloc(sizematsint);
  Policyh       = (int*)malloc(sizematsint);
  Policyr       = (int*)malloc(sizematsint);
  Policyl       = (int*)malloc(sizematsint);

  Changer       = (int*)malloc(sizematsintinfo);

  int indd;

  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int ih=0; ih<nh; ih++){
        for(int im=0; im<nm; im++){
          for(int ia=0; ia<na; ia++){
            indd = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

            Changer[indd]      = 0;

            for(int id=0; id<nd; id++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              
              Value[ind]         = 0.0;
              Value_future[ind]  = 0.0;
              Policyc[ind]       = 0.0;
              Pricing_guess[ind] = 0.0;
              Puncond[ind]       = 0.0;
              Pcond[ind]         = 0.0;

              Default[ind]       = 0;
              Default0[ind]      = 0;
              Renew[ind]         = 0;
              Policya[ind]       = 0;
              Policym[ind]       = 0;
              Policyh[ind]       = 0;
              Policyr[ind]       = 0;
              Policyl[ind]       = 0;
            }
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

  int find_eq = 1;



	// //----------------------------------------------//
	// //---------   INITIAL STEADY STATE   -----------//
	// //----------------------------------------------//

  if (stage == "ssinitial"){

    clock_t t_start;
    clock_t *d_t_start;
    t_start = clock();
    d_t_start = &t_start;
  
    int iteraciones = 1;
    int *d_iteraciones;
    d_iteraciones =  &iteraciones;

    double min_upto = 1000.0;
    double *d_min_upto;
    d_min_upto =  &min_upto;

    double rental = 0.0;
    double *d_rental;
    d_rental = &rental;
  
    double housing = 0.0;
    double *d_housing;
    d_housing = &housing;

    // Minimum up to this point
    double qmin = q;
    double *d_q_upto;
    d_q_upto =  &qmin;

    double dmaxmin = dmax;
    double *d_dmax_upto;
    d_dmax_upto =  &dmaxmin;

    double dminmin = dmin;
    double *d_dmin_upto;
    d_dmin_upto =  &dminmin;

    double ddeltabarmin = ddeltabar;
    double *d_ddeltabar_upto;
    d_ddeltabar_upto =  &ddeltabarmin;

    double m_ymin = m_y;
    double *d_m_y_upto;
    d_m_y_upto =  &m_ymin;

    double bbetamin = bbeta;
    double *d_bbeta_upto;
    d_bbeta_upto =  &bbetamin;

    double fcostmin = fcost;
    double *d_fcost_upto;
    d_fcost_upto =  &fcostmin;
  
    double refcostmin = refcost;
    double *d_refcost_upto;
    d_refcost_upto =  &refcostmin;
  

    //Loading the structure
    pricesolver_eq_24 paramstructura  = {maxiter, uti, tol, convergence, T, Tretirement, yearspp, na, amin, amax, 
                                        nm, mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, nl, lmax, nd,
                                        ny, ssigma_y, llambda_y, m_y, ssigma, rrho, ppsi, bbeta,
                                        kkappa, tthetalab, eetalab, oomega, sunk, interm, rec_probab, ddeltabar, 
                                        ddeltaf, r, Ph, Pa, housing_supply, pension, sstax, ltax, fcost, refcost, Atech, compute_equivalent, multiplier, lumpsum,
                                        d_iteraciones, d_min_upto, d_t_start, 
                                        d_q_upto, d_dmin_upto, d_dmax_upto, d_ddeltabar_upto, d_bbeta_upto, d_fcost_upto, d_refcost_upto,
                                        d_rental, d_housing};

    //Set up the optimization algorrithm
    nlopt_opt opt;

    opt = nlopt_create(NLOPT_LN_NELDERMEAD, 3);// Dimension 2. Algoritthm cobyla    
    // opt = nlopt_create(NLOPT_GN_CRS2_LM, 4);// Dimension 2. Algoritthm cobyla    

    nlopt_set_min_objective(opt, price_zero_eq_24, &paramstructura);
    nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance

    //const double tolerance=1.0e-5;   q, fcost, dvariance
    double LB[3]={0.001, -0.35, 0.2};
    double UB[3]={0.9, -0.1, 0.6};

    nlopt_set_lower_bounds(opt, LB);
    nlopt_set_upper_bounds(opt, UB);

    nlopt_set_maxeval(opt, 5000);

    //  double init[6sion]={0.001,0.01,0.01,0.01,0.000001,0.1};
    double init[3]={0.002, 0.1, 0.1};
    nlopt_set_initial_step(opt, init);

    // Si es suficientemente pequenho el error, lo pongo en cero y paro
    nlopt_set_stopval(opt, 0.01);

    // Initialize at:
    double xtest[3] = {};
    xtest[0] = q;      // r
    xtest[1] = dmin;      // r
    xtest[2] = dmax;      // r

    //Starting the optimization algorithm
    double minf;
    nlopt_optimize(opt, xtest, &minf);

    // Resultados
    vector<double> Res;
    Res.resize(4+1);
    for(int i=0; i<4; i++){
      Res[i]=xtest[i];
    }
    Res[4]=minf;


  } else{

    double rental = 0.0;
    double *d_rental;
    d_rental = &rental;
  
    double housing = 0.0;
    double *d_housing;
    d_housing = &housing;

    parameters params = {maxiter, uti, tol, convergence, T, Tretirement, yearspp, na, amin, amax, 
                        nm, mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, nl, lmax, nd, dmin, dmin+dmax,
                        ny, ssigma_y, llambda_y, m_y, ssigma, rrho, ppsi, bbeta, kkappa, tthetalab, eetalab, 
                        oomega, sunk, interm, rec_probab, ddeltabar, ddeltabar, ddeltaf, r, Ph, Ph, q, Pa, housing_supply, 
                        fcost, refcost, pension, sstax, ltax, lumpsum, Atech, compute_equivalent, multiplier, d_rental, d_housing};

    // Fill grids
    grid_initialize(params, agrid, dgrid, hgrid, rgrid, lgrid, mgrid, ygrid, P, survival, repay_coeff, eprocess);

    // I read Policy function of first period after transition
    import_basic(Default0, Default, Renew, Policym, Policya, Policyh, Policyr, Policyl, Puncond, Pcond);


    double pipol      = 0.0;

    double costo      = 0.0;
    double targetpop  = 0.0;
    double elig       = 0.0;

    // In first period, households receive a lump-sum transfer
    for(int it=0; it<Ttrans; it++){
      lumpsumseq[it]    = lumpsum;
    }

    lumpsumseq[0]   = -0.043;

    double targetDTI;

    prob_mistake = 0.016; // Such that in equilibrium, strategic default is 10%

    if(stage == "noext"){

      // Este se corre ANTES DE BASELINE para guardar 

      clock_t t_start;
      clock_t *d_t_start;
      t_start = clock();
      d_t_start = &t_start;
    
      int iteraciones = 1;
      int *d_iteraciones;
      d_iteraciones =  &iteraciones;
    
      double min_upto = 1000.0;
      double *d_min_upto;
      d_min_upto =  &min_upto;
    

      //Loading the structure
      transitions_qs paramstructura = {maxiter, uti, tol, convergence, T, Tretirement, yearspp, na, amin, amax, nm, mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, nl, lmax,
                                      nd, dmin, dmin+dmax, ny, ssigma_y, llambda_y, m_y, 
                                      ssigma, rrho, ppsi, bbeta, kkappa, tthetalab, eetalab, oomega, sunk, interm, rec_probab, ddeltabar, ddeltaf, 
                                      r, Ph, q, Pa, housing_supply, fcost, refcost, pension, sstax, ltax, Atech, compute_equivalent, multiplier, lumpsum, Ttrans, 
                                      rshock, Pashock, ppsishock, oomegashock, sunk_shock, interm_shock, ltaxshock,
                                      Phinf, qinf, ddeltabarinf,
                                      periods_shock, periods_taxes, experimento, permanent, baseline, tipo,
                                      qnoextseq, lumpsumseq, mort_subsidy, subs_eligible, subs_target, prob_mistake, incshock,
                                      d_iteraciones, d_min_upto, d_t_start, d_rental, d_housing};


      //Set up the optimization algorrithm
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 10);// Dimension 2. Algoritthm cobyla    
    
      nlopt_set_min_objective(opt, transition_eq_noext, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance
    
      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[10]={0.60, 0.7, 0.7, 0.9, 0.9, 0.95, 0.95, 0.95, 0.95, 0.95};
      double UB[10]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    
      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);
    
      nlopt_set_maxeval(opt, 300);
    
      double init[10]={0.03, 0.03, 0.01, 0.01, 0.002, 0.002, 0.001, 0.001, 0.0001, 0.0001};
      nlopt_set_initial_step(opt, init);
    
      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);
    
      // Initialize at:
      double xtest[10] = {};
      // Ph
      xtest[0] = Phseq[0];
      xtest[1] = Phseq[1];
      xtest[2] = Phseq[2];
      xtest[3] = Phseq[3];
      xtest[4] = Phseq[4];
      xtest[5] = Phseq[5];
      xtest[6] = Phseq[6];
      xtest[7] = Phseq[7];
      xtest[8] = Phseq[8];
      xtest[9] = Phseq[9];

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);


    } else if(stage == "base_policy"){

      // Baseline policy: TARP (45%) and HAMP (55%) are implemented

      targetDTI  = 0.275;   // In the data is 31% - this level ensures 45% in TARP and 55% in HAMP

      ltaxshock  = 0.0;
      ltax       = 0.0035;

      Phseq[0] = 0.774241;
      Phseq[1] = 0.820033;
      Phseq[2] = 0.888033;
      Phseq[3] = 0.999999;
      Phseq[4] = 0.997753;
      Phseq[5] = 0.995629;
      Phseq[6] = 0.996954;
      Phseq[7] = 0.997434;
      Phseq[8] = 0.999592;
      Phseq[9] = 0.999937;

      qnoextseq[0]  = q - 0.0056;
      qnoextseq[1]  = q - 0.0043;
      qnoextseq[2]  = q - 0.0027;
      qnoextseq[3]  = q - 0.00012;
      qnoextseq[4]  = q - 0.0001;
      qnoextseq[5]  = q - 0.00005;
      qnoextseq[6]  = q - 0.000;
      qnoextseq[7]  = q - 0.000;
      qnoextseq[8]  = q - 0.000;
      qnoextseq[9]  = q - 0.000;

      double def_tar = 0.0;

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  indd = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

                  subs_eligible[ind] = 0;
                  subs_target[ind]   = 0;
                  mort_subsidy[indd] = 0.0;

                  if(im>0 && ih > 0){
                    pipol = pipol + Puncond[ind];
                  }
                  
                  tipo = "base_policy";

                  // Subsidio a los que tienen >- targetDTI para bajarlos a ese punto
                  if(im>0 && ih>0 && Policyl[ind]>0 && mgrid[im]/(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih]) > targetDTI){

                    subs_eligible[ind]  = 1;  // Elegibles los que tienen PTI > targetPTI
                    mort_subsidy[indd] = maxdosnum(0, 1-(targetDTI*(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih])/mgrid[im]));


                    if(Default[ind] == 1){
                      subs_target[ind] = 1;   // Target son los que planean hacer default

                      def_tar = def_tar + Puncond[ind]*(double)Default[ind];
                      targetpop = targetpop + Puncond[ind];
                    }

                    costo = costo + mgrid[im]*mort_subsidy[ind]*Puncond[ind];
                    elig = elig + Puncond[ind];
                  }

                  if(it == 21 && iy == 1 && ia == 10 && im == 1 && ih == 1 && id == 2){
                    cout << "El subsidio essss = " << mort_subsidy[indd] << endl;
                  }


                }
              }
            }
          }
        }
      }

      cout << "Default = " << def_tar/targetpop << endl;
      cout << endl;
      cout << "Costo = " << costo << ", Population affected = " << targetpop/pipol << ", target = " << elig/pipol << endl;
      cout << endl;

      compute_equivalent = 1;
      multiplier = 0.01;

      baseline = 0;


    } else if(stage == "subsidy_only"){

      // Expanding HAMP and eliminating TARP

      targetDTI  = 0.23;

      ltaxshock  = 0.0;
      ltax       = 0.0025;

      Phseq[0] = 0.776;
      Phseq[1] = 0.821227;
      Phseq[2] = 0.874368;
      Phseq[3] = 0.999999;
      Phseq[4] = 0.997753;
      Phseq[5] = 0.995629;
      Phseq[6] = 0.996954;
      Phseq[7] = 0.997434;
      Phseq[8] = 0.999592;
      Phseq[9] = 0.999937;

      qnoextseq[0]  = q - 0.0056;
      qnoextseq[1]  = q - 0.0043;
      qnoextseq[2]  = q - 0.0027;
      qnoextseq[3]  = q - 0.00012;
      qnoextseq[4]  = q - 0.0001;
      qnoextseq[5]  = q - 0.00005;
      qnoextseq[6]  = q - 0.000;
      qnoextseq[7]  = q - 0.000;
      qnoextseq[8]  = q - 0.000;
      qnoextseq[9]  = q - 0.000;

      double def_tar = 0.0;

      double mean_age_subsidized = 0.0;
      double mean_age_all = 0.0;

      double pip_subsidized = 0.0;
      double pip_all = 0.0;

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  indd = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

                  subs_eligible[ind] = 0;
                  subs_target[ind]   = 0;
                  mort_subsidy[indd] = 0.0;

                  if(im>0 && ih > 0){
                    pipol = pipol + Puncond[ind];
                  }

                  
                  tipo = "subsidy_only";

                  // Some summary statistics
                  if(im>0 && ih>0 && Policyl[ind]>0){
                    mean_age_all = mean_age_all + agrid[ia]*Puncond[ind];
                    pip_all = pip_all + Puncond[ind];
                    if(mgrid[im]/(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih]) > targetDTI){
                      mean_age_subsidized = mean_age_subsidized + agrid[ia]*Puncond[ind];
                      pip_subsidized = pip_subsidized + Puncond[ind];
                    }

                  }


                  // Subsidio a los que tienen >- targetDTI para bajarlos a ese punto
                  if(im>0 && ih>0 && Policyl[ind]>0 && mgrid[im]/(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih]) > targetDTI){

                    subs_eligible[ind]  = 1;
                    mort_subsidy[indd] = maxdosnum(0, 1-(targetDTI*(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih])/mgrid[im]));

                    if(Default[ind]==1){
                      subs_target[ind] = 1;

                      def_tar = def_tar + Puncond[ind]*(double)Default[ind];
                      targetpop = targetpop + Puncond[ind];
                    }

                    costo = costo + mgrid[im]*mort_subsidy[ind]*Puncond[ind];
                    elig = elig + Puncond[ind];
                  }

                  if(it == 21 && iy == 1 && ia == 10 && im == 1 && ih == 1 && id == 2){
                    cout << "El subsidio essss = " << mort_subsidy[indd] << endl;
                  }

                }
              }
            }
          }
        }
      }


      cout << "age subsidies = " << mean_age_subsidized/pip_subsidized << endl;
      cout << "age all = " << mean_age_all/pip_all << endl;

      cout << "Default = " << def_tar/targetpop << endl;
      cout << endl;
      cout << "Costo = " << costo << ", Population affected = " << targetpop/pipol << ", target = " << elig/pipol << endl;
      cout << endl;


      
      compute_equivalent = 1;
      multiplier = 0.004;

      baseline = 0;


    } else if(stage == "bailout_only"){

      ltaxshock  = 0.0;
      ltax       = 0.0047;

      Phseq[0] = 0.765;
      Phseq[1] = 0.815886;
      Phseq[2] = 0.876203;
      Phseq[3] = 0.993659;
      Phseq[4] = 0.992694;
      Phseq[5] = 0.992827;
      Phseq[6] = 0.995083;
      Phseq[7] = 0.996250;
      Phseq[8] = 0.999546;
      Phseq[9] = 0.999950;

      qnoextseq[0]  = q - 0.0059;
      qnoextseq[1]  = q - 0.0047;
      qnoextseq[2]  = q - 0.0028;
      qnoextseq[3]  = q - 0.00015;
      qnoextseq[4]  = q - 0.0001;
      qnoextseq[5]  = q - 0.00005;
      qnoextseq[6]  = q - 0.000;
      qnoextseq[7]  = q - 0.000;
      qnoextseq[8]  = q - 0.000;
      qnoextseq[9]  = q - 0.000;

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  indd = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

                  tipo = "bailout_only";

                  subs_eligible[ind]  = 1;
                  subs_target[ind] = 1;

                  mort_subsidy[indd]   = 0.0;

                }
              }
            }
          }
        }
      }

      compute_equivalent = 1;
      multiplier = 0.005;

      cout << endl;
      cout << "Costo = " << costo << ", Population affected = " << targetpop/pipol << ", target = " << elig/pipol << endl;
      cout << endl;

      baseline = 0;

    } else if(stage == "first_best"){

      // Mismo subsidio de HAMP pero con un componente adicional por edad
      // Elegibilidad cambia: solo se da subsidio a hogares que cambiaron de default decision after shock

      targetDTI = 0.275;

      ltaxshock       = 0.0;
      ltax       = 0.0025;

      Phseq[0] = 0.79;
      Phseq[1] = 0.84;
      Phseq[2] = 0.895;
      Phseq[3] = 0.9983;
      Phseq[4] = 0.997175;
      Phseq[5] = 0.995046;
      Phseq[6] = 0.996910;
      Phseq[7] = 0.997412;
      Phseq[8] = 0.999589;
      Phseq[9] = 0.999941;

      qnoextseq[0]  = q - 0.0055;
      qnoextseq[1]  = q - 0.0039;
      qnoextseq[2]  = q - 0.0023;
      qnoextseq[3]  = q - 0.00012;
      qnoextseq[4]  = q - 0.0001;
      qnoextseq[5]  = q - 0.00005;
      qnoextseq[6]  = q - 0.000;
      qnoextseq[7]  = q - 0.000;
      qnoextseq[8]  = q - 0.000;
      qnoextseq[9]  = q - 0.000;


      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  indd = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

                  if(im>0 && ih>0 && Default[ind]==1 && Default0[ind]==0){
                    Changer[indd]      = 1;
                  }
                }
              }
            }
          }
        }
      }


      // Export eligible households
      ostringstream ss;
      ss << "matrices/Eligibles.txt";
      ofstream Eligibles (ss.str().c_str());
      if (Eligibles.is_open())
      {
        for(int it=0; it<T; it++){
          for(int iy=0; iy<ny; iy++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;
                  
                  if(im>0 && ih>0 && Changer[ind] == 1){
                    Eligibles << 1 << "\n";
                  } else{
                    Eligibles << 0 << "\n";
                  }
                }
              }
            }
          }
        }
        Eligibles.close();
      }
      else cout << "Unable to open file";

      // export_eligible(T, ny, nd, na, nh, nm, Changer);

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  indd = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

                  mort_subsidy[indd] = 0.0;

                  if(im>0 && ih > 0){
                    pipol = pipol + Puncond[ind];
                  }


                  tipo = "first_best";

                  // First best - perfect information

                  if(im>0 && ih>0 && Changer[indd] == 1){
                    subs_eligible[ind]  = 1;
                    mort_subsidy[indd] = 1.0*mindosnum(1.0, 1-(targetDTI*(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih])/mgrid[im]))  + 0.013*pow(T-it,1.0);

                    if(im>0 && ih>0 && Default[ind]==1 && Default0[ind]==0){
                      subs_target[ind] = 1;
                      costo = costo + mgrid[im]*mort_subsidy[ind]*Puncond[ind];
                    }
                  }
                }
              }
            }
          }
        }
      }

      compute_equivalent = 1;
      multiplier = 0.002;

      cout << endl;
      cout << "Costo = " << costo << ", Population affected = " << targetpop/pipol << ", target = " << elig/pipol << endl;
      cout << endl;

      baseline = 0;


    } else if(stage == "second_best"){

      // Modificar tamanho del subsidio de HAMP - elegibilidad es la misma

      targetDTI  = 0.275;

      ltaxshock       = 0.0;
      ltax       = 0.0038;

      Phseq[0] = 0.77;
      Phseq[1] = 0.818;
      Phseq[2] = 0.876;
      Phseq[3] = 0.998;
      Phseq[4] = 0.997175;
      Phseq[5] = 0.995046;
      Phseq[6] = 0.996910;
      Phseq[7] = 0.997412;
      Phseq[8] = 0.999589;
      Phseq[9] = 0.999941;

      qnoextseq[0]  = q - 0.006;
      qnoextseq[1]  = q - 0.0043;
      qnoextseq[2]  = q - 0.0028;
      qnoextseq[3]  = q - 0.00015;
      qnoextseq[4]  = q - 0.0001;
      qnoextseq[5]  = q - 0.00005;
      qnoextseq[6]  = q - 0.000;
      qnoextseq[7]  = q - 0.000;
      qnoextseq[8]  = q - 0.000;
      qnoextseq[9]  = q - 0.000;

      double def_tar = 0.0;

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  indd = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

                  mort_subsidy[indd] = 0.0;

                  if(im>0 && ih > 0){
                    pipol = pipol + Puncond[ind];
                  }

                  
                  tipo = "second_best";

                  // Subsidio a los que tienen >- targetDTI para bajarlos a ese punto
                  if(im>0 && ih>0 && Policyl[ind]>0 && mgrid[im]/(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih]) > targetDTI){

                    subs_eligible[ind]  = 1;

                    if(it < 25){
                      mort_subsidy[indd] = 1.0*mindosnum(1.0, 1-(targetDTI*(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih])/mgrid[im])) + 0.05*pow(it,1.0);
                    } else{
                      mort_subsidy[indd] = 0.0;
                    }

                    if(Default[ind]==1){
                      subs_target[ind] = 1;

                      def_tar = def_tar + Puncond[ind]*(double)Default[ind];
                      targetpop = targetpop + Puncond[ind];
                    }

                    costo = costo + mgrid[im]*mort_subsidy[ind]*Puncond[ind];
                    elig = elig + Puncond[ind];
                  }

                }
              }
            }
          }
        }
      }

      cout << "Default = " << def_tar/targetpop << endl;
      cout << endl;
      cout << "Costo = " << costo << ", Population affected = " << targetpop/pipol << ", target = " << elig/pipol << endl;
      cout << endl;


      
      compute_equivalent = 1;
      multiplier = 0.004;

      baseline = 0;

    }
  } 


  if(find_eq == 1){


    clock_t t_start;
    clock_t *d_t_start;
    t_start = clock();
    d_t_start = &t_start;
  
    int iteraciones = 1;
    int *d_iteraciones;
    d_iteraciones =  &iteraciones;
  
    double min_upto = 1000.0;
    double *d_min_upto;
    d_min_upto =  &min_upto;
  
    double rental = 0.0;
    double *d_rental;
    d_rental = &rental;
  
    double housing = 0.0;
    double *d_housing;
    d_housing = &housing;
  
    // for(int it=0; it<T; it++){
    //     cout << incshock[it] << endl;
    //   }

    //Loading the structure
    transitions_qs paramstructura = {maxiter, uti, tol, convergence, T, Tretirement, yearspp, na, amin, amax, nm, mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, nl, lmax,
                                    nd, dmin, dmin+dmax, ny, ssigma_y, llambda_y, m_y, 
                                    ssigma, rrho, ppsi, bbeta, kkappa, tthetalab, eetalab, oomega, sunk, interm, rec_probab, ddeltabar, ddeltaf, 
                                    r, Ph, q, Pa, housing_supply, fcost, refcost, pension, sstax, ltax, Atech, compute_equivalent, multiplier, lumpsum, Ttrans, 
                                    rshock, Pashock, ppsishock, oomegashock, sunk_shock, interm_shock, ltaxshock,
                                    Phinf, qinf, ddeltabarinf,
                                    periods_shock, periods_taxes, experimento, permanent, baseline, tipo,
                                    qnoextseq, lumpsumseq, mort_subsidy, subs_eligible, subs_target, prob_mistake, incshock,
                                    d_iteraciones, d_min_upto, d_t_start, d_rental, d_housing};


    //Set up the optimization algorrithm
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LN_NELDERMEAD, 10);// Dimension 2. Algoritthm cobyla    
    // opt = nlopt_create(NLOPT_GN_DIRECT_L, 5);// Dimension 2. Algoritthm cobyla
    // opt = nlopt_create(NLOPT_GN_ESCH, 5);// Dimension 2. Algoritthm cobyla
    //  opt = nlopt_create(NLOPT_LN_BOBYQA, 6);// Dimension 2. Algoritthm cobyla
  
    nlopt_set_min_objective(opt, transition_eq_noext, &paramstructura);
    nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance
  
    //const double tolerance=1.0e-5;   q, fcost, dvariance
    double LB[10]={0.60, 0.7, 0.7, 0.9, 0.9, 0.95, 0.95, 0.95, 0.95, 0.95};
    double UB[10]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    // double LB[6]={0.18,0.18,0.75, 0.85, 1, 1};
    // double UB[6]={0.32,0.32,1, 1, 1.2, 1.15};
  
    nlopt_set_lower_bounds(opt, LB);
    nlopt_set_upper_bounds(opt, UB);
  
    nlopt_set_maxeval(opt, 300);
  
    //  double init[6]={0.001,0.01,0.01,0.01,0.000001,0.1};
    double init[10]={0.03, 0.03, 0.01, 0.001, 0.0005, 0.0005, 0.00005, 0.00005, 0.000005, 0.000005};
    nlopt_set_initial_step(opt, init);
  
    // Si es suficientemente pequenho el error, lo pongo en cero y paro
    nlopt_set_stopval(opt, 0.01);
  
    // Initialize at:
    double xtest[10] = {};
    // Ph
    xtest[0] = Phseq[0];
    xtest[1] = Phseq[1];
    xtest[2] = Phseq[2];
    xtest[3] = Phseq[3];
    xtest[4] = Phseq[4];
    xtest[5] = Phseq[5];
    xtest[6] = Phseq[6];
    xtest[7] = Phseq[7];
    xtest[8] = Phseq[8];
    xtest[9] = Phseq[9];

    //Starting the optimization algorithm
    double minf;
    nlopt_optimize(opt, xtest, &minf);
  }


}