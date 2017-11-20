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



//======================================
//         Value Function Iteration
//======================================



int main(int argc, char *argv[])
{

	cout.precision(6);
	cout.setf(std::ios::fixed);

  // VFI parameters
  const double tol         = pow(10,-5.0);
  const int uti            = 1;
  const int maxiter        = 10;
  const double convergence = 0.7; // Entre mas alto, mas lenta es la convergencia

  // Demographics
  const int T        = 12;

  // Grid for savings: a
  const int na       = 35;  // 25
  const double amin  = 0;
  const double amax  = 1.2;

  // Grid for mortgages: m
  const int nm       = 5;  // 11
  const double mmin  = 0.0;
  const double mmax  = 1.5;

  // Grid for housing: h
  const int nh       = 3;
  const double hmin  = 0.0;
  const double hmax  = 5;

  // Grid for renting: r
  const int nr       = 13;
  const double rmin  = 0.0001;
  const double rmax  = 4;

  // Grid for depreciation: ddelta
  const int nd       = 3;

  // Grid for income shocks: y
  const int ny            = 5;
  const double ssigma_y   = 0.1*pow(5,0.5); // 0.44
  const double llambda_y  = pow(0.98,5.0); // pow(0.98,5.0) ; pow(0.95,5.0)
  const double m_y        = 1.5; // 2

  // Preferences
  const double ssigma  = 2;
  const double rrho    = 0.8;
  const double ppsi    = 0.65;
  const double bbeta   = pow(0.964181319,5.0);
  const double kkappa  = -0.1;

  // Equilibrium objects
  double r              = pow(1.02, 5.0)-1;
  double Pa             = 1/(1+r);
  double housing_supply = 21.5;

  // Optimizer initial values
  const double Ph        = 1.0;
  const double q         = 0.235224;
  const double dmin      = -0.067345;
  const double dmax      = 0.156301;
  double ddeltabar = 0.038;
  double ddeltaf   = 1.0;
  const double fcost     = 0.169573;

  // Transitional dynamics' parameters
  int Ttrans;
  double rshock;
  double Pashock;
  double ppsishock;
  int periods_shock;

  const int experimento   = 2;

  if(experimento == 0){

    // Shocking interest rate
    Ttrans        = 3;
    rshock        = r;
    Pashock       = Pa;
    ppsishock     = ppsi;
    periods_shock = 0;

  } else if(experimento == 1){

    // Shocking interest rate
    Ttrans        = 3;
    rshock        = pow(1.05, 5.0)-1;
    Pashock       = Pa;
    ppsishock     = ppsi;
    periods_shock = 3;

  } else if(experimento == 2){

    // Shocking ppsi
    Ttrans        	= 6;
    rshock        	= r;
    Pashock       	= Pa;
    ppsishock     	= 0.7;
    periods_shock 	= 1;
    ddeltaf 		= 0.25;
    ddeltabar 		= 0.019;

  } else if(experimento == 3){

    // Shocking ppsi
    Ttrans        = 3;
    rshock        = pow(1.02, 5.0)-1;
    Pashock       = Pa;
    ppsishock     = 0.7;
    periods_shock = 1;

  }

  // COn ddeltaf mayor
  // const double q = 0.23;
  // const double dmin =  -0.2;
  // const double dmax =  0.703288;
  // const double ddeltabar = 1.085;
  // const double ddeltaf  = 1.0;
  // const double fcost = 0.158946;


  std::string stage = argv[1];

	// //----------------------------------------------//
	// //---------   INITIAL STEADY STATE   -----------//
	// //----------------------------------------------//

  if (stage == "initial"){

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

    int maxim = 3;
    
    if(maxim == 1){
      //Loading the structure
      pricesolver_eq_15 paramstructura  = {maxiter, uti, tol, convergence, T, na, amin, amax, 
                                          nm, mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, nd, 
                                          ny, ssigma_y, llambda_y, ssigma, rrho, ppsi, 
                                          kkappa, ddeltabar, ddeltaf, r, Ph, Pa, housing_supply,
                                          d_iteraciones, d_min_upto, d_t_start, 
                                          d_q_upto, d_dmin_upto, d_dmax_upto, d_ddeltabar_upto, d_bbeta_upto, d_fcost_upto, d_m_y_upto,
                                          d_rental, d_housing};

      //Set up the optimization algorrithm
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 6);// Dimension 2. Algoritthm cobyla    
      // opt = nlopt_create(NLOPT_GN_CRS2_LM, 6);// Dimension 2. Algoritthm cobyla    
      // opt = nlopt_create(NLOPT_GN_DIRECT_L, 6);// Dimension 2. Algoritthm cobyla
      // opt = nlopt_create(NLOPT_GN_ESCH, 6);// Dimension 2. Algoritthm cobyla
      //  opt = nlopt_create(NLOPT_LN_BOBYQA, 6);// Dimension 2. Algoritthm cobyla

      nlopt_set_min_objective(opt, price_zero_eq_15, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance

      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[6]={0.22, 0.15, 0.59, 0.04, -0.2,  0.5};
      double UB[6]={0.4, 0.35, 1.05, 0.3,  -0.05, 2.0};

      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);

      nlopt_set_maxeval(opt, 5000);

      //  double init[6sion]={0.001,0.01,0.01,0.01,0.000001,0.1};
      double init[6]={0.02, 0.1, 0.1, 0.05, 0.1, 0.4};
      nlopt_set_initial_step(opt, init);

      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);

      // Initialize at:
      double xtest[6] = {};
      xtest[0] = q;      // r
      xtest[1] = dmax;      // r
      xtest[2] = bbeta;      // r
      xtest[3] = fcost;      // r
      xtest[4] = dmin;      // r
      xtest[5] = m_y;      // r

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);

      // Resultados
      vector<double> Res;
      Res.resize(6+1);
      for(int i=0; i<6; i++){
        Res[i]=xtest[i];
      }
      Res[6]=minf;

      cout << "Equilibrio encontrado: " << endl;
      cout << "Rental prices q = " << Res[0]<< endl;
      cout << "Dmax = " << Res[1]<< endl;
      cout << "Ddeltabar = " << ddeltabar<< endl;
      cout << "Bbeta = " << Res[2] << endl;
      cout << "Fcost = " << Res[3] << endl;
      cout << "dmin = " << Res[4] << endl;
      cout << "Minimum = " << Res[5] << endl;

    } else if(maxim == 2){
      //Loading the structure
      pricesolver_eq_16 paramstructura  = {maxiter, uti, tol, convergence, T, na, amin, amax, nm, 
                                            mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, nd, dmin, dmax, ny, ssigma_y, llambda_y, 
                                            ssigma, rrho, ppsi, bbeta, kkappa, ddeltabar, ddeltaf, r, Ph, q, Pa, 
                                            housing_supply, fcost, 
                                            d_iteraciones, d_min_upto, d_t_start, 
                                            d_q_upto, d_dmin_upto, d_dmax_upto, d_ddeltabar_upto, d_bbeta_upto, d_fcost_upto, d_m_y_upto,
                                            d_rental, d_housing};

      //Set up the optimization algorrithm
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 1);// Dimension 2. Algoritthm cobyla    
      // opt = nlopt_create(NLOPT_GN_CRS2_LM, 6);// Dimension 2. Algoritthm cobyla    
      // opt = nlopt_create(NLOPT_GN_DIRECT_L, 6);// Dimension 2. Algoritthm cobyla
      // opt = nlopt_create(NLOPT_GN_ESCH, 6);// Dimension 2. Algoritthm cobyla
      //  opt = nlopt_create(NLOPT_LN_BOBYQA, 6);// Dimension 2. Algoritthm cobyla

      nlopt_set_min_objective(opt, price_zero_eq_16, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance

      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[1]={0.5};
      double UB[1]={2.0};

      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);

      nlopt_set_maxeval(opt, 5000);

      //  double init[6sion]={0.001,0.01,0.01,0.01,0.000001,0.1};
      double init[1]={0.4};
      nlopt_set_initial_step(opt, init);

      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);

      // Initialize at:
      double xtest[1] = {};
      xtest[0] = m_y;      // r

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);

      // Resultados
      vector<double> Res;
      Res.resize(1);
      for(int i=0; i<1; i++){
        Res[i]=xtest[i];
      }
      Res[1]=minf;


    } else if(maxim == 3){
      //Loading the structure
      pricesolver_eq_14 paramstructura  = {maxiter, uti, tol, convergence, T, na, amin, amax, 
                                          nm, mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, nd, 
                                          ny, ssigma_y, llambda_y, m_y, ssigma, rrho, ppsi, 
                                          kkappa, ddeltabar, ddeltaf, r, Ph, Pa, housing_supply,
                                          d_iteraciones, d_min_upto, d_t_start, 
                                          d_q_upto, d_dmin_upto, d_dmax_upto, d_ddeltabar_upto, d_bbeta_upto, d_fcost_upto, 
                                          d_rental, d_housing};

      //Set up the optimization algorrithm
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 5);// Dimension 2. Algoritthm cobyla    

      nlopt_set_min_objective(opt, price_zero_eq_14, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance

      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[5]={0.13, 0.15, 0.59, 0.02, -0.2};
      double UB[5]={0.35, 0.35, 1.05, 0.3,  -0.05};

      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);

      nlopt_set_maxeval(opt, 5000);

      //  double init[6sion]={0.001,0.01,0.01,0.01,0.000001,0.1};
      double init[5]={0.01, 0.05, 0.01, 0.02, 0.05};
      nlopt_set_initial_step(opt, init);

      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);

      // Initialize at:
      double xtest[5] = {};
      xtest[0] = q;      // r
      xtest[1] = dmax;      // r
      xtest[2] = bbeta;      // r
      xtest[3] = fcost;      // r
      xtest[4] = dmin;      // r

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);

      // Resultados
      vector<double> Res;
      Res.resize(5+1);
      for(int i=0; i<5; i++){
        Res[i]=xtest[i];
      }
      Res[5]=minf;


    } else if(maxim == 4){
      // Encuentro 1 de equilibrio
      GenEqParameters_eq paramstructura  = {maxiter, uti, tol, convergence, T,
                                          na, amin, amax, nm, mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, nd, dmin, dmax, ny, ssigma_y, llambda_y, m_y, ssigma,
                                          rrho, ppsi, bbeta, kkappa, ddeltabar, ddeltaf, r, Ph, Pa, housing_supply, fcost, 
                                          d_iteraciones, d_min_upto, d_t_start, 
                                            d_q_upto, d_dmin_upto, d_dmax_upto, d_ddeltabar_upto, d_bbeta_upto, d_fcost_upto, d_m_y_upto,
                                            d_rental, d_housing};

      //Set up the optimization algorrithm
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 1);// Dimension 2. Algoritthm cobyla    

      nlopt_set_min_objective(opt, price_zero_eq, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance

      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[1]={0.22};
      double UB[1]={0.3};

      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);

      nlopt_set_maxeval(opt, 5000);

      //  double init[6sion]={0.001,0.01,0.01,0.01,0.000001,0.1};
      double init[1]={0.01};
      nlopt_set_initial_step(opt, init);

      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);

      // Initialize at:
      double xtest[1] = {};
      xtest[0] = q;      // r

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);

    } 


  } else if(stage == "final"){

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
  
    //Loading the structure
    transitions_qs paramstructura = {maxiter, uti, tol, convergence, T, na, amin, amax, nm, mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, 
                                    nd, dmin, dmin+dmax, ny, ssigma_y, llambda_y, m_y, 
                                    ssigma, rrho, ppsi, bbeta, kkappa, ddeltabar, ddeltaf, 
                                    r, Ph, q, Pa, housing_supply, fcost, Ttrans, rshock, Pashock, ppsishock, periods_shock,
                                    d_iteraciones, d_min_upto, d_t_start, d_rental, d_housing};


    if(Ttrans == 3){
      //Set up the optimization algorrithm
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 6);// Dimension 2. Algoritthm cobyla    
  
      nlopt_set_min_objective(opt, transition_eq_qs, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance
  
      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[6]={0.65, 0.8, 0.85, 1, 1, 1};
      double UB[6]={1.05, 1.05, 1.05, 1.2, 1.2, 1.2};
  
      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);
  
      nlopt_set_maxeval(opt, 300);
  
      double init[6]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
      nlopt_set_initial_step(opt, init);
  
      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);
  
      // Initialize at:
      double xtest[6] = {};
      // Ph
      xtest[0] = 0.8552047;
      xtest[1] = 0.98105;
      xtest[2] = 0.997210;
      // xtest[0] = 1.0;
      // xtest[1] = 1.0;
      // xtest[2] = 1.0;
      // ddeltabar
      xtest[3] = 1.191;
      xtest[4] = 1.15;
      xtest[5] = 1.11;
      // xtest[3] = ddeltabar;
      // xtest[4] = ddeltabar;
      // xtest[5] = ddeltabar;

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);

    } else if(Ttrans == 5){
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 10);// Dimension 2. Algoritthm cobyla    
    
      nlopt_set_min_objective(opt, transition_eq_qs, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance
    
      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[10]={0.75, 0.8, 0.85, 0.8, 0.85, 1, 1, 1, 1, 1};
      double UB[10]={1.05, 1.05, 1.05, 1.05, 1.05, 1.2, 1.2, 1.2, 1.2, 1.1};
    
      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);
    
      nlopt_set_maxeval(opt, 1000);
    
      //  double init[6]={0.001,0.01,0.01,0.01,0.000001,0.1};
      double init[10]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
      nlopt_set_initial_step(opt, init);
    
      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);
    
      // Initialize at:
      double xtest[10] = {};
      // Ph
      xtest[0] = 0.84;
      xtest[1] = 0.981461;
      xtest[2] = 0.995;
      xtest[3] = 0.999;
      xtest[4] = 1.002579;
      // xtest[0] = 1.0;
      // xtest[1] = 1.0;
      // xtest[2] = 1.0;
      // xtest[3] = 1.0;
      // xtest[4] = 1.0;
      // ddeltabar
      xtest[5] = 1.165;
      xtest[6] = 1.084;
      xtest[7] = 1.06;
      xtest[8] = 1.045;
      xtest[9] = 1.077;
      // xtest[5] = ddeltabar;
      // xtest[6] = ddeltabar;
      // xtest[7] = ddeltabar;
      // xtest[8] = ddeltabar;
      // xtest[9] = ddeltabar;

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);

    } else if(Ttrans == 6){
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 12);// Dimension 2. Algoritthm cobyla    
    
      nlopt_set_min_objective(opt, transition_eq_qs, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance
    
      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[12]={0.75, 0.8, 0.85, 0.8, 0.85, 0.8, 0, 0, 0, 0, 0, 0};
      double UB[12]={1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 0.4, 0.2, 0.2, 0.2, 0.1, 0.1};
    
      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);
    
      nlopt_set_maxeval(opt, 1000);
    
      //  double init[6]={0.001,0.01,0.01,0.01,0.000001,0.1};
      double init[12]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
      nlopt_set_initial_step(opt, init);
    
      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);
    
      // Initialize at:
      double xtest[12] = {};
      // Ph
      xtest[0] = 0.897913;
      xtest[1] = 0.979154;
      xtest[2] = 0.997792;
      xtest[3] = 1.010844;
      xtest[4] = 1.004875;
      xtest[5] = 1.000967;
      // xtest[0] = 1.0;
      // xtest[1] = 1.0;
      // xtest[2] = 1.0;
      // xtest[3] = 1.0;
      // xtest[4] = 1.0;
      // xtest[5] = 1.0;
      // ddeltabar
      xtest[6] = 0.083;
      xtest[7] = 0.051;
      xtest[8] = 0.016;
      xtest[9] = 0.010;
      xtest[10] = 0.018;
      xtest[11] = 0.021;
      // xtest[6] = ddeltabar;
      // xtest[7] = ddeltabar;
      // xtest[8] = ddeltabar;
      // xtest[9] = ddeltabar;
      // xtest[10] = ddeltabar;
      // xtest[11] = ddeltabar;

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);

    } else if(Ttrans == 8){
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 16);// Dimension 2. Algoritthm cobyla    
    
      nlopt_set_min_objective(opt, transition_eq_qs, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance
    
      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[16]={0.75, 0.8, 0.85, 0.8, 0.85, 0.85, 0.8, 0.85, 1, 1, 1, 1, 1, 1, 1, 1};
      double UB[16]={1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.2, 1.15, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1};
    
      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);
    
      nlopt_set_maxeval(opt, 300);
    
      //  double init[6]={0.001,0.01,0.01,0.01,0.000001,0.1};
      double init[16]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
      nlopt_set_initial_step(opt, init);
    
      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);
    
      // Initialize at:
      double xtest[16] = {};
      // Ph
      xtest[0] = 0.78;
      xtest[1] = 0.98;
      xtest[2] = 0.995;
      xtest[3] = 1.0;
      xtest[4] = 1.0;
      xtest[5] = 1.0;
      xtest[6] = 1.0;
      xtest[7] = 1.0;
      // ddeltabar
      xtest[8] = 1.11;
      xtest[9] = 1.07;
      xtest[10] = 1.02;
      xtest[11] = 1.02;
      xtest[12] = 1.04;
      xtest[13] = 1.05;
      xtest[14] = 1.06;
      xtest[15] = 1.076;

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);
    }

  } else if(stage == "noext"){

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
  
    //Loading the structure
    transitions_qs paramstructura = {maxiter, uti, tol, convergence, T, na, amin, amax, nm, mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, 
                                    nd, dmin, dmin+dmax, ny, ssigma_y, llambda_y, m_y, 
                                    ssigma, rrho, ppsi, bbeta, kkappa, ddeltabar, ddeltaf, 
                                    r, Ph, q, Pa, housing_supply, fcost, Ttrans, rshock, Pashock, ppsishock, periods_shock,
                                    d_iteraciones, d_min_upto, d_t_start, d_rental, d_housing};

    if(Ttrans == 3){

      //Set up the optimization algorrithm
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 3);// Dimension 2. Algoritthm cobyla    
      // opt = nlopt_create(NLOPT_GN_DIRECT_L, 5);// Dimension 2. Algoritthm cobyla
      // opt = nlopt_create(NLOPT_GN_ESCH, 5);// Dimension 2. Algoritthm cobyla
      //  opt = nlopt_create(NLOPT_LN_BOBYQA, 6);// Dimension 2. Algoritthm cobyla
    
      nlopt_set_min_objective(opt, transition_eq_noext, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance
    
      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[3]={0.75, 0.85, 0.95};
      double UB[3]={1.05, 1.05, 1.05};
    
      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);
    
      nlopt_set_maxeval(opt, 300);
    
      //  double init[6]={0.001,0.01,0.01,0.01,0.000001,0.1};
      double init[3]={0.01, 0.01, 0.01};
      nlopt_set_initial_step(opt, init);
    
      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);
    
      // Initialize at:
      double xtest[3] = {};
      // Ph
      xtest[0] = 0.825;
      xtest[1] = 0.985;
      xtest[2] = 1.0;

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);

    }  else if(Ttrans == 6){
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 6);// Dimension 2. Algoritthm cobyla    
    
      nlopt_set_min_objective(opt, transition_eq_noext, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance
    
      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[6]={0.75, 0.8, 0.85, 0.8, 0.85, 0.8};
      double UB[6]={1.05, 1.05, 1.05, 1.05, 1.05, 1.05};
    
      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);
    
      nlopt_set_maxeval(opt, 1000);
    
      //  double init[6]={0.001,0.01,0.01,0.01,0.000001,0.1};
      double init[6]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
      nlopt_set_initial_step(opt, init);
    
      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);
    
      // Initialize at:
      double xtest[6] = {};
      // Ph
      xtest[0] = 0.935;
      xtest[1] = 0.985;
      xtest[2] = 0.997792;
      xtest[3] = 1.0;
      xtest[4] = 1.0;
      xtest[5] = 1.0;
      // xtest[0] = 1.0;
      // xtest[1] = 1.0;
      // xtest[2] = 1.0;
      // xtest[3] = 1.0;
      // xtest[4] = 1.0;
      // ddeltabar

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);

    } else if(Ttrans == 5){

      //Set up the optimization algorrithm
      nlopt_opt opt;
      opt = nlopt_create(NLOPT_LN_NELDERMEAD, 5);// Dimension 2. Algoritthm cobyla    
      // opt = nlopt_create(NLOPT_GN_DIRECT_L, 5);// Dimension 2. Algoritthm cobyla
      // opt = nlopt_create(NLOPT_GN_ESCH, 5);// Dimension 2. Algoritthm cobyla
      //  opt = nlopt_create(NLOPT_LN_BOBYQA, 6);// Dimension 2. Algoritthm cobyla
    
      nlopt_set_min_objective(opt, transition_eq_noext, &paramstructura);
      nlopt_set_xtol_rel(opt, 1.0e-5); //Tolerance
    
      //const double tolerance=1.0e-5;   q, fcost, dvariance
      double LB[5]={0.75, 0.85, 0.95, 0.95, 0.95};
      double UB[5]={1.05, 1.05, 1.05, 1.05, 1.05};
      // double LB[6]={0.18,0.18,0.75, 0.85, 1, 1};
      // double UB[6]={0.32,0.32,1, 1, 1.2, 1.15};
    
      nlopt_set_lower_bounds(opt, LB);
      nlopt_set_upper_bounds(opt, UB);
    
      nlopt_set_maxeval(opt, 300);
    
      //  double init[6]={0.001,0.01,0.01,0.01,0.000001,0.1};
      double init[5]={0.01, 0.01, 0.01, 0.01, 0.01};
      nlopt_set_initial_step(opt, init);
    
      // Si es suficientemente pequenho el error, lo pongo en cero y paro
      nlopt_set_stopval(opt, 0.01);
    
      // Initialize at:
      double xtest[5] = {};
      // Ph
      xtest[0] = 0.82;
      xtest[1] = 0.97;
      xtest[2] = 0.997;
      xtest[3] = 0.998;
      xtest[4] = 1.0;

      //Starting the optimization algorithm
      double minf;
      nlopt_optimize(opt, xtest, &minf);

    }
  }

}

