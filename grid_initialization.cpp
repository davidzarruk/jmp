
//======================================
//         Grids
//======================================

void grida(const int na, const double amin, const double amax, double* gridaa){
  
  double size = na;
  double astep = (amax - amin) /(size - 1);
  double it = 0;
  
  for(int i = 0; i < na; i++){
    gridaa[i] = amin + it*astep;
    it++;
    // cout << gridaa[i] << endl;
  }
}

void gridd(const int nd, const double dmin, const double dmax, double* gridaa){
    
  double size = nd;
  double astep = (dmax - dmin) /(size - 1);
  double it = 0;
  
  for(int i = 0; i < nd; i++){
    gridaa[i] = dmin + it*astep;
    it++;
  }
}

void gridm(const int nm, const double mmin, const double mmax, double* gridmm){
  
  double size = nm;
  double mstep = (mmax - mmin) /(size - 1);
  double it = 0;
  
  for(int i = 0; i < nm; i++){
    gridmm[i] = mmin + it*mstep;
    it++;
  }
}

void gridh(const int nh, const double hmin, const double hmax, double* gridhh){
  
  double size = nh;
  double hstep = (hmax - hmin) /(size - 1);
  double it = 0;
  
  for(int i = 0; i < nh; i++){
    gridhh[i] = hmin + it*hstep;
    it++;
  }
}

void gridr(const int nr, const double rmin, const double rmax, double* gridrr){
  
  double size = nr;
  double rstep = (rmax - rmin) /(size - 1);
  double it = 0;
  
  for(int i = 0; i < nr; i++){
    gridrr[i] = rmin + it*rstep;
    it++;
  }
}

void gridl(const int nl, const double lmax, double* gridll){
  
  double size = nl;
  double lstep = lmax /(size - 1);
  double it = 0;
  
  for(int i = 0; i < nl; i++){
    gridll[i] = it*lstep;
    it++;
  }
}

void gridy(const int ny, const double ssigma_y, const double llambda_y, const double m, double* gridyy){
  
  // This grid is made with Tauchen (1986)
  double size = ny;
  double ssigma_yy = sqrt(pow(ssigma_y, 2) / (1 - pow(llambda_y, 2)));
  double ystep = 2*ssigma_yy*m / (size-1);
  
  double it = 0;
  
  for(int i = 0; i < ny; i++){
    gridyy[i] = -m*sqrt(pow(ssigma_y, 2) / (1 - pow(llambda_y, 2))) + it*ystep;
    it = it + 1;
  }
}

double normCDF(const double value){
  return 0.5 * erfc(-value * M_SQRT1_2);
}


void yprob(const int ny, const double ssigma_y, const double llambda_y, const double m, const double* gridyy, double* P){
  
  // This grid is made with Tauchen (1986)
  const double w = gridyy[1] - gridyy[0];
  
  for(int j = 0; j < ny; j++){
    for(int k = 0; k < ny; k++){
      if(k == 0){
        P[j*ny + k] = normCDF((gridyy[k] - llambda_y*gridyy[j] + (w/2))/ssigma_y);
      } else if(k == ny-1){
        P[j*ny + k] = 1 - normCDF((gridyy[k] - llambda_y*gridyy[j] - (w/2))/ssigma_y);
      } else{
        P[j*ny + k] = normCDF((gridyy[k] - llambda_y*gridyy[j] + (w/2))/ssigma_y) - normCDF((gridyy[k] - llambda_y*gridyy[j] - (w/2))/ssigma_y);
      }
    }
  }
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

void eproc(double* eprocess, const double T, const double Atech){
  
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

    for(int i = 0; i < T; i++){
      eprocess[i] = Atech*eprocess[i];
    }
  }
}

void survival_coefficient(const int T, double* survival, int yearspp){
    
  // Survival probabilities
  double survival_long[71] = {0};
  surv(survival_long);

  // for(int it=0; it<71; it++){
  //   cout << survival_long[it] << endl;
  // }
  
  for(int it=0; it<T; it++){
    survival[it] = 1;
    // cout << survival[it] << endl;
    for(int i=0; i<71; i++){
      if((i>=it*yearspp) & (i<(it+1)*yearspp)){
        survival[it] = survival[it]*survival_long[i];
        // cout << survival_long[i] << endl;
      }
    }
    // cout << survival[it] << endl;
  }
  survival[T] = 0;

}

void repay_coefficient(const int T, const double r, const double rrho, const double* survival, double* repay_coeff){
  
  double survive = 1.0;
  for(int it=0; it<T; it++){
    repay_coeff[it] = 0.0;
    for(int itp=it+1; itp<T; itp++){
      survive = survive*survival[itp-1];
      repay_coeff[it] = repay_coeff[it] + survive*pow(rrho/(1+r),itp-it);

      // cout << "r = " << 1+r << ", rrho = " << rrho << ", power = " << pow(rrho/(1+r),itp-it) << ", surv = " << survive << endl;
    }
    survive = 1.0;
  }
}

void income_process(const int T, const int Tretirement, const int yearspp, double* eprocess, const double Atech){
  
  double eprocess_long[71];
  eproc(eprocess_long, 71, Atech);

  for(int it=0; it<T; it++){
    eprocess[it] = 0;
    for(int i=0; i<71; i++){
      if((i >= it*yearspp) & (i < (it+1)*yearspp)){
        eprocess[it] = eprocess[it] + eprocess_long[i]/yearspp;
      }
    }
    if(it >= Tretirement){
      eprocess[it] = 0.0;
    }

    // cout << eprocess[it] << endl;
  }
}

// When there is a shock, the repay coefficient depends on the interest rates pre and post shock
void repay_coefficient_shock(const int T, const double r0, const double rshock, const int periods, 
                              const double rrho, const double* survival, double* repay_coeff){
  
  double survive = 1.0;
  double discount;
  double r;

  for(int it=0; it<T; it++){
    repay_coeff[it] = 0.0;
    discount = 1.0;
    for(int itp=it+1; itp<T; itp++){
      survive = survive*survival[itp-1];
      if(itp-it <= periods){
        r = rshock;
      } else{
        r = r0;
      }

      discount = discount*rrho/(1+r);

      repay_coeff[it] = repay_coeff[it] + survive*discount;
    }
    survive = 1.0;
  }
}


//======================================
//         Others
//======================================

double absol(double a){
  if(a < 0){
    a = -a;
  }
  return a;
}


void steady_distribution(const double* P,
                         const parameters params, 
                         double *steady_distr){

  const int ny = params.ny;

  // Steady state distribution for y
  double Paux[ny*ny];
  double C[ny*ny];
  
  for(int iy = 0; iy<ny; iy++){
    for(int iyp = 0; iyp<ny; iyp++){
      Paux[iy*ny+iyp]=P[iy*ny+iyp];
      C[iy*ny+iyp]=0.0;
    }
  }

  for(int itt = 0; itt<10000; itt++){
    for(int iy=0; iy<ny; iy++){
      for(int iyp=0; iyp<ny; iyp++){
        
        C[iy*ny+iyp]=0.0;
        
        for(int is=0; is<ny; is++){
          C[iy*ny+iyp] = C[iy*ny+iyp] + P[iy*ny+is]*Paux[is*ny+iyp];
        }
      }
    }

    for(int iy=0; iy<ny; iy++){
      for(int iyp=0; iyp<ny; iyp++){
        Paux[iy*ny+iyp] = C[iy*ny+iyp];
      }
    }
  }

  for(int iy=0; iy<ny; iy++){
    steady_distr[iy] = Paux[0*ny+iy];
    // cout << steady_distr[iy] << endl;
  }

}


// Supply of housing - linear supply function
double housing_supply(double slope, double hPrice){
  
  double sup = (hPrice - (1-0.69*slope))/slope;
  return(sup);
  
}





//======================================
//         Initialize grids
//======================================

void grid_initialize(parameters params, double *agrid, double *dgrid, 
                      double *hgrid, double *rgrid, double *lgrid, double *mgrid, double *ygrid, double *P,
                      double *survival, double *repay_coeff, double *eprocess){
  
  // VFI parameters
  const int T            = params.T;
  const int Tretirement            = params.Tretirement;
  const int yearspp            = params.yearspp;

  // Grid for savings: a
  const int na           = params.na;
  const double amin      = params.amin;
  const double amax      = params.amax;

  // Grid for mortgages: m
  const int nm           = params.nm;
  const double mmin      = params.mmin;
  const double mmax      = params.mmax;

  // Grid for housing: h
  const int nh           = params.nh;
  const double hmin      = params.hmin;
  const double hmax      = params.hmax;

  // Grid for renting: r
  const int nr           = params.nr;
  const double rmin      = params.rmin;
  const double rmax      = params.rmax;

  // Grid for labor: l
  const int nl           = params.nl;
  const double lmax      = params.lmax;

  // Grid for deoreciation: ddelta
  const int nd           = params.nd;
  const double dmin      = params.dmin;

  // Grid for income shocks: y
  const int ny           = params.ny;
  const double ssigma_y  = params.ssigma_y;
  const double llambda_y = params.llambda_y;
  const double m_y       = params.m_y;

  // Preferences
  const double rrho      = params.rrho;
  const double interm    = params.interm;

  const double Atech     = params.Atech;

  // Equilibrium objects
  double r                  = params.r;
  double dmax               = params.dmax;

  //----------------------------------------------//
  //--------   GRIDS AND PROBABILITY    ----------//
  //----------------------------------------------//

  grida(na, amin, amax, agrid);
  gridd(nd, dmin, dmax, dgrid);
  gridm(nm, mmin, mmax, mgrid);
  gridh(nh, hmin, hmax, hgrid);
  gridr(nr, rmin, rmax, rgrid);
  gridl(nl, lmax, lgrid);
  gridy(ny, ssigma_y, llambda_y, m_y, ygrid);
  yprob(ny, ssigma_y, llambda_y, m_y, ygrid, P);

  for(int i=0; i<ny; i++){
    ygrid[i] = exp(ygrid[i]);
  }

  survival_coefficient(T, survival, yearspp);

  repay_coefficient(T, r + interm, rrho, survival, repay_coeff);
  income_process(T, Tretirement, yearspp, eprocess, Atech);


  // cout << dmin << endl;
  // cout << dmax << endl;
  
  for(int it=0; it<T; it++){
     // cout << repay_coeff[it] << endl;
  }

  for(int it=0; it<T; it++){
     // cout << survival[it] << endl;
  }
  
  // double alive = 1;
  // for(int it=0; it<T; it++){
  //    cout << alive << endl;
  //   alive = alive*survival[it];
  // }

  for(int i=0; i<ny; i++){
     // cout << ygrid[i] << endl;
  }

  for(int i=0; i<nm; i++){
     // cout << mgrid[i] << endl;
  }
  
}
