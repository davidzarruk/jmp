

//======================================
//         Error and aggregates
//======================================

double vfi_error(parameters params, const int error_funct, const int ittrans,
                  const double *P, const double *survival, const double *hgrid, const double *rgrid, const double *lgrid, const double *mgrid, const double *repay_coeff,
                  const double *agrid, const double *dgrid, const double *ygrid, const double *eprocess,
                  const double *Value, const double *Policyc, const double *Pricing_guess, 
                  const int *Default, const int *Renew, const int *Policya, 
                  const int *Policym, const int *Policyh, const int *Policyr, const int *Policyl,
                  const double *Pcond, const double *Puncond, const double *Pcond_debt, const double *Puncond_debt){

  // VFI parameters
  const int maxiter         = params.maxiter;
  const int uti             = params.uti;
  const double tol          = params.tol;
  const int T               = params.T;
  const int Tretirement     = params.Tretirement;
  const int yearspp         = params.yearspp;

  // Grid for savings: a
  const int na      = params.na;
  const double amin = params.amin;
  const double amax = params.amax;

  // Grid for mortgages: m
  const int nm      = params.nm;
  const double mmin = params.mmin;
  const double mmax = params.mmax;

  // Grid for housing: h
  const int nh      = params.nh;
  const double hmin = params.hmin;
  const double hmax = params.hmax;

  // Grid for renting: r
  const int nr      = params.nr;
  const double rmin = params.rmin;
  const double rmax = params.rmax;

  // Grid for deoreciation: ddelta
  const int nd      = params.nd;
  const double dmin = params.dmin;
  const double dmax = params.dmax;

  // Grid for income shocks: y
  const int ny            = params.ny;
  const double ssigma_y   = params.ssigma_y;
  const double llambda_y  = params.llambda_y;
  const double m_y        = params.m_y;

  // Preferences
  const double ssigma = params.ssigma;
  const double rrho   = params.rrho;
  const double ppsi   = params.ppsi;
  const double bbeta  = params.bbeta;
  const double kkappa = params.kkappa;

  const double eetalab  = params.eetalab;
  const double tthetalab  = params.tthetalab;


  // Equilibrium objects
  const double ddeltabar      = params.ddeltabar_today;
  const double ddeltaf        = params.ddeltaf;
  const double r              = params.r;
  const double Ph             = params.Ph_today;
  const double q              = params.q;
  const double Pa             = params.Pa;
  const double fcost          = params.fcost;
  const double refcost        = params.refcost;
  const double pension        = params.pension;
  const double sstax          = params.sstax;
  // const double Atech          = params.Atech;
  // const double lumpsum        = params.lumpsum;
  // const double mortsubsidy    = params.mortsubsidy;
  const double housing_supply = params.housing_supply;
  const double interm  = params.interm;
  const double sunk  = params.sunk;
  // const double oomega  = params.oomega;
  const double lumpsum  = params.lumpsum;

  double *d_rental    = params.d_rental;
  double *d_housing   = params.d_housing;

  int ind;
  
  // cout << na << endl;
  // cout << nm << endl;
  // cout << mmax << endl;
  // cout << nh << endl;
  // cout << nr << endl;
  // cout << nd << endl;
  // cout << ny << endl;
  // cout << ssigma << endl;
  // cout << rrho << endl;
  // cout << ppsi << endl;
  // cout << bbeta << endl;
  // cout << kkappa << endl;
  // cout << sunk << endl;
  // cout << interm << endl;
  // cout << oomega << endl;
  // cout << ddeltabar << endl;
  // cout << r << endl;
  // cout << Ph << endl;
  // cout << q << endl;
  // cout << Pa << endl;
  // cout << fcost << endl;

  //--------------------------------------//
  //----  Aggregate paths in economy  ----//
  //--------------------------------------//

  double C[T];
  double L[T];
  double Income[T];
  double H[T];
  double Hp[T];
  double M[T];
  double Mp[T];
  double R[T];
  double NR[T];         // % of net renters
  double NetRent[T];    // Net rent
  double OwnerOccupied[T];
  double A[T];
  double Aup[T];
  double Adn[T];
  double DD[T];
  double RR[T];
  double RRup[T];
  double RRdn[T];
  double RRdenom[T];
  double Prepay[T];
  double PP[T];
  double KK[T];
  double Def[T];
  double Mort_debt[T];
  double equity[T];
  double owners_age[T];
  double Vivos[T];
  double MortHolders[T];

  for(int it=0; it<T; it++){

    C[it]             = 0.0;
    L[it]             = 0.0;
    Income[it]        = 0.0;
    H[it]             = 0.0;
    Hp[it]            = 0.0;
    M[it]             = 0.0;
    Mp[it]            = 0.0;
    R[it]             = 0.0;
    NR[it]            = 0.0;         // % of neit reniters
    NetRent[it]       = 0.0;    // Neit renit
    OwnerOccupied[it] = 0.0;
    A[it]             = 0.0;
    Aup[it]             = 0.0;
    Adn[it]             = 0.0;
    DD[it]            = 0.0;
    RR[it]            = 0.0;
    RRup[it]            = 0.0;
    RRdn[it]            = 0.0;
    RRdenom[it]            = 0.0;
    Prepay[it]        = 0.0;
    PP[it]            = 0.0;
    KK[it]            = 0.0;
    Def[it]           = 0.0;
    Mort_debt[it]     = 0.0;
    equity[it]        = 0.0;
    owners_age[it]    = 0.0;
    Vivos[it]         = 0.0;
    MortHolders[it]   = 0.0;

  }

  int condind;
  double condDD[T*ny];
  double condYpeople[T*ny];
  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      condind = it*ny + iy;
      condDD[condind] = 0.0;
      condYpeople[condind] = 0.0;
    }
  }

  int condinddelta;
  double condDDdelta[T*nd];
  double condDeltapeople[T*nd];
  for(int it=0; it<T; it++){
    for(int id=0; id<nd; id++){
      condinddelta = it*nd + id;
      condDDdelta[condinddelta] = 0.0;
      condDeltapeople[condinddelta] = 0.0;
    }
  }

  int indTxAxY;
  double DcondTxAxY[T*na*ny];
  double DcondTxAxYpeople[T*na*ny];
  for(int it=0; it<T; it++){
    for(int ia=0; ia<na; ia++){
      for(int iy=0; iy<ny; iy++){
        indTxAxY = it*na*ny + ia*ny + iy;

        DcondTxAxY[indTxAxY] = 0.0;
        DcondTxAxYpeople[indTxAxY] = 0.0;
      }
    }
  }

  // Conditional on A and Y
  int indAxY;
  double DcondAxY[na*ny];
  double DcondAxYpeople[na*ny];
  for(int ia=0; ia<na; ia++){
    for(int iy=0; iy<ny; iy++){
      indAxY = ia*ny + iy;

      DcondAxY[indAxY] = 0.0;
      DcondAxYpeople[indAxY] = 0.0;
    }
  }


  // Conditional on A
  double DcondA[na];
  double DcondApeople[na];
  for(int ia=0; ia<na; ia++){
    DcondA[ia] = 0.0;
    DcondApeople[ia] = 0.0;
  }

  // Conditional on H and M
  int indHxM;
  double DcondHxM[nh*nm];
  double DcondHxMpeople[nh*nm];
  for(int ih=0; ih<nh; ih++){
    for(int im=0; im<nm; im++){
      indHxM = ih*nm + im;

      DcondHxM[indHxM] = 0.0;
      DcondHxMpeople[indHxM] = 0.0;
    }
  }


  // Leverage conditional on A and Y
  double LcondAxY[na*ny];
  double LcondAxYpeople[na*ny];
  for(int ia=0; ia<na; ia++){
    for(int iy=0; iy<ny; iy++){
      indAxY = ia*ny + iy;

      LcondAxY[indAxY] = 0.0;
      LcondAxYpeople[indAxY] = 0.0;
    }
  }

  // Housing conditional on A and Y
  double HcondAxY[na*ny];
  double HcondAxYpeople[na*ny];
  for(int ia=0; ia<na; ia++){
    for(int iy=0; iy<ny; iy++){
      indAxY = ia*ny + iy;

      HcondAxY[indAxY] = 0.0;
      HcondAxYpeople[indAxY] = 0.0;
    }
  }

  // Mortgage conditional on A and Y
  double McondAxY[na*ny];
  double McondAxYpeople[na*ny];
  for(int ia=0; ia<na; ia++){
    for(int iy=0; iy<ny; iy++){
      indAxY = ia*ny + iy;

      McondAxY[indAxY] = 0.0;
      McondAxYpeople[indAxY] = 0.0;
    }
  }


  // Aggregate stats
  double statistics[100];
  double people             = 0.0;
  double owners_w_mortgage  = 0.0;
  double owners_w_mortgage2  = 0.0;
  double mortgage_debt      = 0.0;
  double housing_value      = 0.0;
  double default_rate       = 0.0;
  double housing_total      = 0.0;
  double rent_total         = 0.0;
  double ownership          = 0.0;
  double ownership2         = 0.0;
  double rentership         = 0.0;  
  double av_owned           = 0.0;
  double av_rented          = 0.0;
  double earn_owner         = 0.0;
  double earn_renter        = 0.0;
  double exp_rent           = 0.0;
  double exp_cons           = 0.0;
  double net_renters        = 0.0;
  double total_income       = 0.0;
  double under_35           = 0.0;
  double under_45           = 0.0;
  double under_55           = 0.0;
  double under_65           = 0.0;
  double under_75           = 0.0;
  double people_35          = 0.0;
  double people_45          = 0.0;
  double people_55          = 0.0;
  double people_65          = 0.0;
  double people_75          = 0.0;
  double mean_delta         = 0.0;
  double mean_income        = 0.0;
  double mean_delta_unc     = 0.0;
  double mean_income_unc    = 0.0;
  double underwater         = 0.0;
  double def_underwater     = 0.0;
  double overwater          = 0.0;
  double def_overwater      = 0.0;
  double mortgagors         = 0.0;

  double pension_desembols  = 0.0;
  double ss_revenues        = 0.0;

  double def_y1   = 0.0;
  double def_y2   = 0.0;
  double def_y3   = 0.0;
  double def_y4   = 0.0;
  double def_y5   = 0.0;

  double peo_y1   = 0.0;
  double peo_y2   = 0.0;
  double peo_y3   = 0.0;
  double peo_y4   = 0.0;
  double peo_y5   = 0.0;

  double debttt = 0.0;
  double persons = 0.0;

  double conteo = 0.0;
  double conteo2 = 0.0;

  double debt_mort    = 0.0;
  double home_equity  = 0.0;
  double heqne        = 0.0;
  double heq10        = 0.0;
  double heq20        = 0.0;
  double heq25        = 0.0;
  double heq30        = 0.0;
  double heq50        = 0.0;
  double heq75        = 0.0;
  double heq100       = 0.0;
  double hequal100    = 0.0;
  double heqma        = 0.0;
  double avheq = 0.0;

  double av_value     = 0.0;
  double av_price     = 0.0;

  double home_eq[130];
  double home_eq_pdf[130];

  for(int i = 0; i<130; i++){
    home_eq[i] = 0.0;
    home_eq_pdf[i] = 0.0;
  }

  double home_eq_pdf_dist[26];
  double def_dist[26];
  for(int i = 0; i<26; i++){
    home_eq_pdf_dist[i] = 0.0;
    def_dist[i] = 0.0;
  }

  double Atot = 0.0;
  double Mtot = 0.0;
  double assetsToIncome = 0.0;
  double pip = 0.0;
  
  double av_M_payment  = 0.0;
  double av_M_payment_people  = 0.0;

  double DTI[T];
  double DTI_y[T];
  double DTI_m[T];
  double DTI_people[T];

  for(int i = 0; i<T; i++){
    DTI[i] = 0.0;
    DTI_y[i] = 0.0;
    DTI_m[i] = 0.0;
    DTI_people[i] = 0.0;
  }

  int ind2;

  double avLTV = 0.0;
  double reLTV = 0.0;

  double agC = 0.0;
  double agY = 0.0;
  double agH = 0.0;
  double agF = 0.0;
  double agA = 0.0;
  double agA2 = 0.0;

  // Lump sum transfers for: 1. assets of dead, 2. houses of dead, 3. price discount of banks
  double lump_assets = 0.0;
  double lump_houses = 0.0;
  double lump_disco  = 0.0;



  // double DTI = 0.0;
  // int indicator;
  // // double cons1;
  // // double cons2;

  // for(int it=0; it<T; it++){
  //   cout << eprocess[it]/Atech << endl;
  // }

  // double ccnodef = 0.0;
  // double yynodef = 0.0;

  // double ccdef = 0.0;
  // double yydef = 0.0;

  double retornos_h = 0.0;
  double retornos_a = 0.0;

  double labor_supply = 0.0;
  double workers      = 0.0;

  double realtor_profits = 0.0;

  // Correlation coefficients
  double corr_age_equity = 0.0;
  double corr_equity_assets = 0.0;
  
  double mean_age = 0.0;
  double mean_equity = 0.0;
  double mean_assets = 0.0;  
  double tot_homeowners = 0.0;

  double cov_age_equity = 0.0;
  double cov_equity_assets = 0.0;
  
  double std_age = 0.0;
  double std_equity = 0.0;
  double std_assets = 0.0;

  double perc_renew_incshock = 0.0;
  double perc_renew_prishock = 0.0;
  double people_renew        = 0.0;

  double perc_defau_incshock_lowa = 0.0;
  double perc_defau_prishock_lowa = 0.0;
  double people_default_lowa      = 0.0;

  double perc_defau_incshock_higha = 0.0;
  double perc_defau_prishock_higha = 0.0;
  double people_default_higha      = 0.0;

  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

              agA = agA + Puncond[ind]*agrid[Policya[ind]];
              agA2 = agA2 + Puncond[ind]*agrid[ia];
              
              agC = agC + Puncond[ind]*Policyc[ind];
              agH = agH + Puncond[ind]*hgrid[Policyh[ind]];
              if(it < Tretirement){
                agY = agY + Puncond[ind]*ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]];
              } else{
                // agY = agY + Puncond[ind]*ygrid[iy]*pension;
              }

              if(Renew[ind] == 1 && Policym[ind]>0){
                agF = agF + Puncond[ind]*fcost*(1+repay_coeff[it])*mgrid[Policym[ind]];
              }

              // Refinancing
              if(Renew[ind] == 1 && im>0 && Policym[ind]>0 && ia>0){
                if(iy == 1){
                  perc_renew_incshock = perc_renew_incshock + Puncond[ind];
                }
                if(id == nd-1){
                  perc_renew_prishock = perc_renew_prishock + Puncond[ind];
                }
                people_renew = people_renew + Puncond[ind];
              }

              if(Default[ind] == 1 && im>0){
                if(ia == 0){
                  if(iy == 1){
                    perc_defau_incshock_lowa = perc_defau_incshock_lowa + Puncond[ind];
                  }
                  if(id == nd-1){
                    perc_defau_prishock_lowa = perc_defau_prishock_lowa + Puncond[ind];
                  }
                  people_default_lowa = people_default_lowa + Puncond[ind];
                } else{
                  if(iy == 1){
                    perc_defau_incshock_higha = perc_defau_incshock_higha + Puncond[ind];
                  }
                  if(id == nd-1){
                    perc_defau_prishock_higha = perc_defau_prishock_higha + Puncond[ind];
                  }
                  people_default_higha = people_default_higha + Puncond[ind];
                }
              }

              lump_assets = lump_assets + Puncond[ind]*(1-survival[it])*agrid[Policya[ind]];
              lump_houses = lump_houses + Puncond[ind]*(1-survival[it])*hgrid[Policyh[ind]];

              if(im > 0){
                lump_disco  = lump_disco  + Puncond[ind]*Default[ind]*(Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih]*sunk);
              }

              realtor_profits = realtor_profits + Puncond[ind]*Default[ind]*(1-Renew[ind])*Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih]*sunk;




//               if(it == 1 && iy == 1 && im == 1 && ia == 1 && ih == 1 && id == 1){

//                 cout << "PRUEBA:" << endl;
//                 cout << Policyc[ind] << endl;
//                 cout << agrid[ia] + q*hgrid[ih] + eprocess[it]*ygrid[iy]*lgrid[Policyl[ind]]*(1-sstax) -  mgrid[im] - q*rgrid[Policyr[ind]] - Pa*agrid[Policya[ind]] - Ph*(dgrid[id]+ddeltabar)*hgrid[Policyh[ind]] - lumpsum << endl;
//                 cout << "PRUEBA:" << endl;


//               }
//               // if(it < Tretirement && Default[ind] == 0 && Renew[ind] == 0){
//                 ccnodef = ccnodef + Puncond[ind]*Policyc[ind];
//                 // yynodef = yynodef + Puncond[ind]*(agrid[ia] + q*hgrid[ih] + eprocess[it]*ygrid[iy]*lgrid[Policyl[ind]]*(1-sstax) -  mgrid[im] - q*rgrid[Policyr[ind]] - Pa*agrid[Policya[ind]] - Ph*(dgrid[id]+ddeltabar)*hgrid[Policyh[ind]] - lumpsum);
//                 yynodef = yynodef + Puncond[ind]*(eprocess[it]*ygrid[iy]*lgrid[Policyl[ind]]);
//               // }

//               if(it < Tretirement && Default[ind] == 0 && Renew[ind] == 0 && (Policyc[ind] - (agrid[ia] + q*hgrid[ih] + eprocess[it]*ygrid[iy]*lgrid[Policyl[ind]]*(1-sstax) -  mgrid[im] - q*rgrid[Policyr[ind]] - Pa*agrid[Policya[ind]] - Ph*(dgrid[id]+ddeltabar)*hgrid[Policyh[ind]] - lumpsum)) > 0.1){
//                 cout << Policyc[ind] << endl;
//                 cout << agrid[ia] + q*hgrid[ih] + eprocess[it]*ygrid[iy]*lgrid[Policyl[ind]]*(1-sstax) -  mgrid[im] - q*rgrid[Policyr[ind]] - Pa*agrid[Policya[ind]] - Ph*(dgrid[id]+ddeltabar)*hgrid[Policyh[ind]] - lumpsum << endl;
//                 cout << it << " " << iy << " " << ia << " " << im << " " << ih << " " << id << " " << endl;
//               }

              retornos_h = retornos_h + Puncond[ind]*dgrid[id]*Ph*hgrid[ih];
              retornos_a = retornos_a + Puncond[ind]*agrid[ia]*r;


// // cons = maximumab(aa - rec_probab*((1+repay_coeff[it])*mm - Ph*(1-ddelta - ddeltabar)*hh*(1-sunk)), 0) + yy*(1-incshock[it]) - q*hhrent - Pa*aaprime - lumpsum;
              
//               if(it < Tretirement && Default[ind] == 1 && Renew[ind] == 0){
//                 ccdef = ccdef + Puncond[ind]*Policyc[ind];
//                 // yydef = yydef + Puncond[ind]*(agrid[ia] + eprocess[it]*ygrid[iy]*lgrid[Policyl[ind]]*(1-sstax) - q*rgrid[Policyr[ind]] - Pa*agrid[Policya[ind]] - lumpsum);
//                 yydef = yydef + Puncond[ind]*(eprocess[it]*ygrid[iy]*lgrid[Policyl[ind]]);
//               }

//               if(it < Tretirement && Default[ind] == 1 && Renew[ind] == 0 && (Policyc[ind] - (agrid[ia] + eprocess[it]*ygrid[iy]*lgrid[Policyl[ind]]*(1-sstax) - q*rgrid[Policyr[ind]] - Pa*agrid[Policya[ind]] - lumpsum)) > 0.1){
//                 cout << Policyc[ind] << endl;
//                 cout << agrid[ia] + eprocess[it]*ygrid[iy]*lgrid[Policyl[ind]]*(1-sstax) - q*rgrid[Policyr[ind]] - Pa*agrid[Policya[ind]] - lumpsum << endl;
//                 cout << it << " " << iy << " " << ia << " " << im << " " << ih << " " << id << " " << endl;
//               }


              // Means:
              C[it] = C[it] + Pcond[ind]*Policyc[ind];
              L[it] = L[it] + Pcond[ind]*lgrid[Policyl[ind]];
              H[it] = H[it] + Pcond[ind]*hgrid[Policyh[ind]];
              M[it] = M[it] + Pcond[ind]*mgrid[Policym[ind]];
              R[it] = R[it] + Pcond[ind]*rgrid[Policyr[ind]];
              A[it] = A[it] + Pcond[ind]*agrid[Policya[ind]];
              PP[it] = PP[it] + Pcond[ind]*Pricing_guess[ind];


              labor_supply = labor_supply + Puncond[ind]*lgrid[Policyl[ind]];
              if(it < Tretirement){
                workers = workers + Puncond[ind]; 
              }

              Atot = Atot + Puncond[ind]*agrid[Policya[ind]];
              // Mortgage amount disimbursed

              ind2 = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + Policyh[ind]*nm*na + Policym[ind]*na + Policya[ind];

              if(Renew[ind] == 1 && Policym[ind] > 0){
                Mtot = Mtot + Puncond[ind]*mgrid[Policym[ind]]*Pricing_guess[ind2];
              }

              if(it >= 15 && it <= 40){
                assetsToIncome = assetsToIncome + Pcond[ind]*(agrid[Policya[ind]] / (ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]]));
                pip = pip + Pcond[ind];
              }

              av_value = av_value + Puncond[ind]*Value[ind];
              av_price = av_price + Puncond[ind]*Pricing_guess[ind];


              if(it < Tretirement){
                Income[it] = Income[it] + Pcond[ind]*ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]]*(1-sstax);
                ss_revenues = ss_revenues + Puncond[ind]*ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]]*sstax;
              } else{
                Income[it] = Income[it] + Pcond[ind]*ygrid[iy]*pension;
                pension_desembols = pension_desembols + Puncond[ind]*ygrid[iy]*pension;
              }


              ind2 = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + Policyh[ind]*nm*na + Policym[ind]*na + Policya[ind];

              if(Renew[ind] == 1 && Policym[ind] > 0){
                avLTV = avLTV + Puncond[ind]*(mgrid[Policym[ind]]*(1+repay_coeff[it]))/(hgrid[Policyh[ind]]*Ph);
                reLTV = reLTV + Puncond[ind];
              }




              // If rents more than the size owned is a "net renter"
              if(rgrid[Policyr[ind]] > hgrid[Policyh[ind]]){
                NR[it] = NR[it] + Pcond[ind];
                NetRent[it] = NetRent[it] + Pcond[ind]*(rgrid[Policyr[ind]] - hgrid[Policyh[ind]]);
                OwnerOccupied[it] = OwnerOccupied[it] + Pcond[ind]*(hgrid[Policyh[ind]]);
                net_renters = net_renters + Puncond[ind];
              } else{
                OwnerOccupied[it] = OwnerOccupied[it] + Pcond[ind]*(rgrid[Policyr[ind]]);
              }

              if(ih>0){
                Hp[it] = Hp[it] + Pcond[ind];
              }

              if(im > 0){
                if(Policym[ind] > 0){
                  RR[it] = RR[it] + Pcond[ind]*(double)Renew[ind];
                } else{
                  Prepay[it] = Prepay[it] + Pcond[ind]*(double)Renew[ind];
                }

                if(Policym[ind] > 0 && Policym[ind] < im){
                  RRdn[it] = RRdn[it] + Pcond[ind]*(double)Renew[ind];
                  Adn[it] = Adn[it] + Pcond[ind]*agrid[ia];

                } else if(Policym[ind] > 0 && Policym[ind] > im){
                  RRup[it] = RRup[it] + Pcond[ind]*(double)Renew[ind];
                  Aup[it] = Aup[it] + Pcond[ind]*agrid[ia];
                }

                RRdenom[it] = RRdenom[it] + Pcond[ind];
              }
              
              if(im > 0){
                KK[it] = KK[it] + Pcond[ind]*(double)(1-Renew[ind])*(double)(1-Default[ind]);
                Def[it] = Def[it] + Pcond[ind]*(double)Default[ind];
                MortHolders[it] = MortHolders[it] + Pcond[ind];

                debttt = debttt + Puncond[ind]*mgrid[im];
                mean_delta = mean_delta + dgrid[id]*Puncond[ind]*(double)Default[ind];

                Mp[it] = Mp[it] + Pcond[ind];

                if(Default[ind] == 1 && im > 0){
                  default_rate  = default_rate  + Puncond[ind];
                  mean_income   = mean_income   + Puncond[ind]*ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]];
                }

                mean_delta_unc = mean_delta_unc + Puncond[ind]*dgrid[id];
                mean_income_unc = mean_income_unc + Puncond[ind]*ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]];

                // Default conditional on income
                condind = it*ny + iy;
                condDD[condind]       = condDD[condind] + Puncond[ind]*(double)Default[ind];
                condYpeople[condind]  = condYpeople[condind] + Puncond[ind];

                condinddelta = it*nd + id;
                condDDdelta[condinddelta]      = condDDdelta[condinddelta] + Puncond[ind]*(double)Default[ind];
                condDeltapeople[condinddelta]  = condDeltapeople[condinddelta] + Puncond[ind];

                indTxAxY = it*na*ny + ia*ny + iy;
                DcondTxAxY[indTxAxY]        = DcondTxAxY[indTxAxY] + Puncond[ind]*(double)Default[ind];
                DcondTxAxYpeople[indTxAxY]  = DcondTxAxYpeople[indTxAxY] + Puncond[ind];

                indAxY = ia*ny + iy;
                DcondAxY[indAxY] = DcondAxY[indAxY] + Puncond[ind]*(double)Default[ind];
                DcondAxYpeople[indAxY] = DcondAxYpeople[indAxY] + Puncond[ind];

                DcondA[ia] = DcondA[ia] + Puncond[ind]*(double)Default[ind];
                DcondApeople[ia] = DcondApeople[ia] + Puncond[ind];

                if(ih > 0){
                  LcondAxY[indAxY] = LcondAxY[indAxY] + Puncond[ind]*(mgrid[im]*(1+repay_coeff[it]))/(Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih]);
                  LcondAxYpeople[indAxY] = LcondAxYpeople[indAxY] + Puncond[ind];
                }

                HcondAxY[indAxY] = HcondAxY[indAxY] + Puncond[ind]*hgrid[ih];
                HcondAxYpeople[indAxY] = HcondAxYpeople[indAxY] + Puncond[ind];

                McondAxY[indAxY] = McondAxY[indAxY] + Puncond[ind]*mgrid[im];
                McondAxYpeople[indAxY] = McondAxYpeople[indAxY] + Puncond[ind];

                indHxM = ih*nm + im;
                DcondHxM[indHxM] = DcondHxM[indHxM] + Puncond[ind]*(double)Default[ind];
                DcondHxMpeople[indHxM] = DcondHxMpeople[indHxM] + Puncond[ind];
              }



              if(Default[ind] == 1 && ih == 0){
                conteo = conteo + Puncond[ind];
              }

              if(Default[ind] == 1 && ih == 0 && im > 0){
                conteo2 = conteo2 + Puncond[ind];
              }


              if(im > 0){
                mortgagors = mortgagors + Puncond[ind];
              }

              if(im > 0){
                if(iy == 0){
                  peo_y1 = peo_y1 + Puncond[ind];
                  if(iy == 0 && Default[ind] == 1){
                    def_y1 = def_y1 + Puncond[ind];
                  }
                } else if(iy == 1){
                  peo_y2 = peo_y2 + Puncond[ind];
                  if(iy == 1 && Default[ind] == 1){
                    def_y2 = def_y2 + Puncond[ind];
                  }
                } else if(iy == 2){
                  peo_y3 = peo_y3 + Puncond[ind];
                  if(iy == 2 && Default[ind] == 1){
                    def_y3 = def_y3 + Puncond[ind];
                  }
                } else if(iy == 3){
                  peo_y4 = peo_y4 + Puncond[ind];
                  if(iy == 3 && Default[ind] == 1){
                    def_y4 = def_y4 + Puncond[ind];
                  }
                } else if(iy == 4){
                  peo_y5 = peo_y5 + Puncond[ind];
                  if(iy == 4 && Default[ind] == 1){
                    def_y5 = def_y5 + Puncond[ind];
                  }
                }  
              }
              // if(Policym[ind] > 0 & Policyh[ind] == 0 & Puncond[ind] > 0.001){
              //   cout << "alert " << it << " " << mgrid[Policym[ind]] << endl;

              //   indicator = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + Policyh[ind]*nm*na + Policym[ind]*na + Policya[ind];

              //   cons1 = agrid[ia] + Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih] + q*hgrid[Policyh[ind]] + ygrid[iy]*eprocess[it] - fcost + mgrid[Policym[ind]]*Pricing_guess[indicator] - repay_coeff[it]*mgrid[im] - q*hgrid[Policyr[ind]] - Ph*hgrid[Policyh[ind]] - Pa*agrid[Policya[ind]];
                
              //   cons2 = agrid[ia] + q*hgrid[ih] + ygrid[iy]*eprocess[it] - mgrid[im] - q*hgrid[Policyr[ind]] - Pa*agrid[Policya[ind]] - Ph*(dgrid[id]+ddeltabar)*hgrid[ih];

              //   cout << Pricing_guess[indicator] << endl;
              //   cout << Default[ind] << endl;
              //   cout << Renew[ind] << endl;
              //   cout << Policyc[ind] << endl;
              //   cout << cons1 << endl;
              //   cout << cons2 << endl;
              // }

              if(mgrid[im]*(1+repay_coeff[it]) >  Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih] & ih > 0){
                underwater = underwater + Puncond[ind];
                def_underwater = def_underwater + Puncond[ind]*(double)Default[ind];
              }

              if(mgrid[im]*(1+repay_coeff[it]) <=  Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih] & ih > 0){
                overwater = overwater + Puncond[ind];
                def_overwater = def_overwater + Puncond[ind]*(double)Default[ind];
              }

              people = people + Puncond[ind];
              
              if(Policym[ind] > 0 && Policyh[ind] > 0){
                owners_w_mortgage = owners_w_mortgage + Puncond[ind];
              }

              if(im > 0 && ih > 0){
                owners_w_mortgage2 = owners_w_mortgage2 + Puncond[ind];
              }

              if(Policyh[ind] > 0){
                ownership = ownership + Puncond[ind];
                earn_owner = earn_owner + Puncond[ind]*ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]];
              }

              if(ih > 0){
                ownership2 = ownership2 + Puncond[ind];
              }

              if(Policyr[ind] > 0){
                rentership = rentership + Puncond[ind];
                earn_renter = earn_renter + Puncond[ind]*ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]];
              }

              if(it <= 15/yearspp){
                people_35 = people_35 + Puncond[ind];
                if(rgrid[Policyr[ind]] > hgrid[Policyh[ind]]){
                  under_35 = under_35 + Puncond[ind];
                }
              } else if(it > 15/yearspp && it <= 25/yearspp){
                people_45 = people_45 + Puncond[ind];
                if(rgrid[Policyr[ind]] > hgrid[Policyh[ind]]){
                  under_45 = under_45 + Puncond[ind];
                }
              } else if(it > 25/yearspp && it <= 35/yearspp){
                people_55 = people_55 + Puncond[ind];
                if(rgrid[Policyr[ind]] > hgrid[Policyh[ind]]){
                  under_55 = under_55 + Puncond[ind];
                }
              } else if(it > 35/yearspp && it <= 45/yearspp){
                people_65 = people_65 + Puncond[ind];
                if(rgrid[Policyr[ind]] > hgrid[Policyh[ind]]){
                  under_65 = under_65 + Puncond[ind];
                }
              } else if(it > 45/yearspp && it <= 55/yearspp){
                people_75 = people_75 + Puncond[ind];
                if(rgrid[Policyr[ind]] > hgrid[Policyh[ind]]){
                  under_75 = under_75 + Puncond[ind];
                }
              }


              persons = persons + Pcond[ind];

              av_owned = av_owned + Puncond[ind]*hgrid[Policyh[ind]];
              av_rented = av_rented + Puncond[ind]*rgrid[Policyr[ind]];
              
              ind2 = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + Policyh[ind]*nm*na + Policym[ind]*na + Policya[ind];

              mortgage_debt = mortgage_debt + Puncond[ind]*mgrid[Policym[ind]]*Pricing_guess[ind2];

              // home_equity = (hgrid[Policyh[ind]]*Ph - mgrid[Policym[ind]]*Pricing_guess[ind2])/(hgrid[Policyh[ind]]*Ph);
              // debt_mort = (mgrid[Policym[ind]]*Pricing_guess[ind2])/(hgrid[Policyh[ind]]*Ph);

              // for(int i = 0; i<100; i++){
              //   if(Policyh[ind] > 0 && home_equity <= (1+(double)i)/100){
              //     home_eq[i] = home_eq[i] + Puncond[ind];
              //   } 
              // }

              // if(Policyh[ind] > 0){
              //   owners_age[it] = owners_age[it] + Pcond[ind];
              //   Mort_debt[it] = Mort_debt[it] + Pcond[ind]*debt_mort;
              //   equity[it] = equity[it] + Pcond[ind]*home_equity;
              // } 

              if(ih > 0){
                home_equity = (Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih] - mgrid[im]*(1+repay_coeff[it]))/(Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih]);
                debt_mort = (mgrid[Policym[ind]]*Pricing_guess[ind2])/(hgrid[Policyh[ind]]*Ph);

                mean_age = mean_age + Puncond[ind]*it;
                mean_equity = mean_equity + Puncond[ind]*home_equity;
                mean_assets = mean_assets + Puncond[ind]*agrid[ia];

                tot_homeowners = tot_homeowners + Puncond[ind];

              } else{
                home_equity = 0;
              }

              for(int i = -30; i<100; i++){
                if(ih > 0 && home_equity <= (1+(double)i)/100){
                  home_eq[i+30]   = home_eq[i+30] + Puncond[ind];
                } 

                if(ih > 0 && home_equity <= (1+(double)i)/100 && home_equity > ((double)i)/100){
                  home_eq_pdf[i+30]   = home_eq_pdf[i+30] + Puncond[ind];
                } 
              }

              for(int i = -6; i<20; i++){

                if(ih > 0 && home_equity <= 5*(1+(double)i)/100 && home_equity > 5*((double)i)/100){
                  home_eq_pdf_dist[i+6]   = home_eq_pdf_dist[i+6] + Puncond[ind];
                  def_dist[i+6]  = def_dist[i+6] + Puncond[ind]*((double)Default[ind]);
                } 
              }


              if(ih > 0){
                owners_age[it] = owners_age[it] + Pcond[ind];
                Mort_debt[it] = Mort_debt[it] + Pcond[ind]*debt_mort;
                equity[it] = equity[it] + Pcond[ind]*home_equity;
              } 



              if(ih > 0 && home_equity < 0.0){
                heqne = heqne + Puncond[ind];
              } 
              if(ih > 0 && home_equity < 0.1){
                heq10 = heq10 + Puncond[ind];
              } 
              if(ih > 0 && home_equity < 0.2){
                heq20 = heq20 + Puncond[ind];
              } 
              if(ih > 0 && home_equity < 0.25){
                heq25 = heq25 + Puncond[ind];
              } 
              if(ih > 0 && home_equity < 0.3){
                heq30 = heq30 + Puncond[ind];
              } 
              if(ih > 0 && home_equity < 0.5){
                heq50 = heq50 + Puncond[ind];
              } 
              if(ih > 0 && home_equity < 0.75){
                heq75 = heq75 + Puncond[ind];
              } 
              if(ih > 0 && home_equity <= 1){
                heq100 = heq100 + Puncond[ind];
              } 
              if(ih > 0 && home_equity == 1){
                hequal100 = hequal100 + Puncond[ind];
              } 
              if(ih > 0 && home_equity > 1){
                heqma = heqma + Puncond[ind];
              }

              if(ih > 0){
                avheq = avheq + home_equity*Puncond[ind];
              }

              // Debt-to-income ratio
              if(im>0 && ih>0 && Policyl[ind]>0){            
              // if(im>0 && ih>0 && Default[ind] == 0 & Renew[ind] == 0 && Policyl[ind]>0){            
                av_M_payment        = av_M_payment        + Puncond[ind]*mgrid[im]/(ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih]);
                av_M_payment_people = av_M_payment_people + Puncond[ind];

                // DTI[it] = DTI[it] + Puncond[ind]*mgrid[im]/(ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]]);
                // DTI_y[it] = DTI_y[it] + Puncond_debt[ind]*(ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]]);
                // DTI_m[it] = DTI_m[it] + Puncond[ind]*mgrid[im];
                // DTI_people[it] = DTI_people[it] + Puncond_debt[ind];

                // DTI[it] = DTI[it]               + Puncond[ind]*mgrid[im]*(1+repay_coeff[it])/(ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih]);
                // DTI_y[it] = DTI_y[it]           + Puncond[ind]*(ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih]);
                // DTI_m[it] = DTI_m[it]           + Puncond[ind]*mgrid[im]*(1+repay_coeff[it]);
                // DTI_people[it] = DTI_people[it] + Puncond[ind];


                DTI[it] = DTI[it]               + Puncond[ind]*mgrid[im]/(-lumpsum + ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih]);
                DTI_y[it] = DTI_y[it]           + Puncond[ind]*(-lumpsum + ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]] + agrid[ia]*r + q*hgrid[ih]);
                DTI_m[it] = DTI_m[it]           + Puncond[ind]*mgrid[im];
                DTI_people[it] = DTI_people[it] + Puncond[ind];
              }


              housing_value = housing_value + Puncond[ind]*hgrid[Policyh[ind]]*Ph;
              housing_total = housing_total + Puncond[ind]*hgrid[Policyh[ind]];

              rent_total = rent_total + Puncond[ind]*rgrid[Policyr[ind]];

              exp_rent = exp_rent + Puncond[ind]*rgrid[Policyr[ind]]*q;
              exp_cons = exp_cons + Puncond[ind]*Policyc[ind];

              Vivos[it] = Vivos[it] + Puncond[ind];

              total_income = total_income + Puncond[ind]*(ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]]);
            }
          }
        }
      }
    }
  }

  mean_age = mean_age/tot_homeowners;
  mean_equity = mean_equity/tot_homeowners;
  mean_assets = mean_assets/tot_homeowners;

  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

              if(ih > 0){
                home_equity = (Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih] - mgrid[im]*(1+repay_coeff[it]))/(Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih]);

                std_age = std_age + Puncond[ind]*pow(it - mean_age, 2);
                std_equity = std_equity + Puncond[ind]*pow(home_equity - mean_equity, 2);
                std_assets = std_assets + Puncond[ind]*pow(agrid[ia] - mean_assets, 2);

                cov_age_equity = cov_age_equity + Puncond[ind]*(it - mean_age)*(home_equity - mean_equity);
                cov_equity_assets = cov_equity_assets + Puncond[ind]*(home_equity - mean_equity)*(agrid[ia] - mean_assets);

              }
            }
          }
        }
      }
    }
  }


  cout << "renew income shock = " << perc_renew_incshock/people_renew << "renew price shock = " << perc_renew_prishock/people_renew << ", renewers = " << people_renew << endl;
  // cout << "defau income shock = " << perc_defau_incshock/people_default << "defau price shock = " << perc_defau_prishock/people_default << endl;

  // cout << endl;
  // cout << endl;
  // for(int i=0; i<T; i++){
  //   cout << RRup[i]/RRdenom[i] << "  " << endl;
  // }
  // cout << endl;
  // cout << endl;
  // for(int i=0; i<T; i++){
  //   cout << RRdn[i]/RRdenom[i] << "  " << endl;
  // }
  // cout << endl;
  // cout << endl;

  // for(int i=0; i<T; i++){
  //   cout << Aup[i]/RRdenom[i] << "  " << endl;
  // }
  // cout << endl;
  // cout << endl;
  // for(int i=0; i<T; i++){
  //   cout << Adn[i]/RRdenom[i] << "  " << endl;
  // }
  // cout << endl;
  // cout << endl;



  avLTV = avLTV/reLTV;

  assetsToIncome = assetsToIncome / pip;

  mean_delta = mean_delta / default_rate;

  mean_income = mean_income / default_rate;
  
  mean_delta_unc = mean_delta_unc / mortgagors;
  mean_income_unc = mean_income_unc / mortgagors;

  def_underwater = def_underwater/underwater;
  underwater = underwater/mortgagors;

  def_overwater = def_overwater/overwater;
  overwater = overwater/mortgagors;

  av_M_payment = av_M_payment/av_M_payment_people;

  heqne = heqne/ownership;
  heq10 = heq10/ownership;
  heq20 = heq20/ownership;
  heq25 = heq25/ownership;
  heq30 = heq30/ownership;
  heq50 = heq50/ownership;
  heq75 = heq75/ownership;
  heq100 = heq100/ownership;
  hequal100 = hequal100/ownership;
  heqma = heqma/ownership;
  avheq = avheq / ownership;

  lump_assets = lump_assets/people;
  lump_houses = lump_houses/people;
  lump_disco  = lump_disco/people;

  for(int i = -30; i<100; i++){
    home_eq[i+30] = home_eq[i+30]/ownership;
  }

  for(int i = -6; i<20; i++){
    def_dist[i+6] = def_dist[i+6]/home_eq_pdf_dist[i+6];
  }

  for(int i = 0; i<T; i++){
    Mort_debt[i] = Mort_debt[i]/owners_age[i];
    equity[i] = equity[i]/owners_age[i];

    Def[i] = Def[i]/MortHolders[i];
    // cout << Income[i] << endl;
  }

  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      condind = it*ny + iy;

      condDD[condind]       = condDD[condind]/condYpeople[condind];
    }
  }

  for(int it=0; it<T; it++){
    for(int id=0; id<nd; id++){
      condinddelta = it*nd + id;

      condDDdelta[condinddelta]       = condDDdelta[condinddelta]/condDeltapeople[condinddelta];
    }
  }

  for(int it=0; it<T; it++){
    for(int ia=0; ia<na; ia++){
      for(int iy=0; iy<ny; iy++){
        indTxAxY = it*na*ny + ia*ny + iy;

        DcondTxAxY[indTxAxY] = DcondTxAxY[indTxAxY]/DcondTxAxYpeople[indTxAxY];
      }
    }
  }


  // for(int ia=0; ia<na; ia++){
  //   for(int iy=0; iy<ny; iy++){
  //     indAxY = ia*ny + iy;

  //     DcondAxY[indAxY] = DcondAxY[indAxY]/DcondAxYpeople[indAxY];
  //   }
  // }

  for(int ia=0; ia<na; ia++){
    DcondA[ia] = DcondA[ia]/DcondApeople[ia];
  }


  for(int ih=0; ih<nh; ih++){
    for(int im=0; im<nm; im++){
      indHxM = ih*nm + im;

      DcondHxM[indHxM] = DcondHxM[indHxM]/DcondHxMpeople[indHxM];
    }
  }

  for(int ia=0; ia<na; ia++){
    for(int iy=0; iy<ny; iy++){
      indAxY = ia*ny + iy;

      LcondAxY[indAxY] = LcondAxY[indAxY]/LcondAxYpeople[indAxY];
    }
  }


  for(int ia=0; ia<na; ia++){
    for(int iy=0; iy<ny; iy++){
      indAxY = ia*ny + iy;

      HcondAxY[indAxY] = HcondAxY[indAxY]/HcondAxYpeople[indAxY];
    }
  }

  for(int ia=0; ia<na; ia++){
    for(int iy=0; iy<ny; iy++){
      indAxY = ia*ny + iy;

      McondAxY[indAxY] = McondAxY[indAxY]/McondAxYpeople[indAxY];
    }
  }


  for(int i = 0; i<T; i++){
    DTI[i] = DTI[i] / DTI_people[i];
    DTI_y[i] = DTI_y[i] / DTI_people[i];
    DTI_m[i] = DTI_m[i] / DTI_people[i];

    // cout << DTI[i] << " - " << DTI_y[i] << " - " << Income[i] << " - " << DTI_m[i] << endl;
  }

  owners_w_mortgage = owners_w_mortgage/ownership;
  owners_w_mortgage2 = owners_w_mortgage2/ownership2;

  default_rate = default_rate/mortgagors;

  net_renters = net_renters/people;

  av_owned = av_owned/ownership;
  av_rented = av_rented/rentership;

  earn_owner = earn_owner/ownership;
  earn_renter = earn_renter/rentership;

  ownership = ownership/people;
  ownership2 = ownership2/people;
  rentership = rentership/people;

  under_35 = under_35/people_35;
  under_45 = under_45/people_45;
  under_55 = under_55/people_55;
  under_65 = under_65/people_65;
  under_75 = under_75/people_75;

  // Parameters
  statistics[0] = tol;
  statistics[1] = uti;
  statistics[2] = maxiter;
  statistics[3] = T;

  statistics[4] = na;
  statistics[5] = amin;
  statistics[6] = amax;

  statistics[7] = nm;
  statistics[8] = mmin;
  statistics[9] = mmax;

  statistics[10] = nh;
  statistics[11] = hmin;
  statistics[12] = hmax;

  statistics[13] = nr;
  statistics[14] = rmin;
  statistics[15] = rmax;

  statistics[16] = nd;
  statistics[17] = dmin;
  statistics[18] = dmax;

  statistics[19] = ny;
  statistics[20] = ssigma_y;
  statistics[21] = llambda_y;
  statistics[22] = m_y;

  statistics[23] = ssigma;
  statistics[24] = rrho;
  statistics[25] = ppsi;
  statistics[26] = bbeta;
  statistics[27] = kkappa;
  statistics[28] = ddeltabar;
  statistics[29] = r;
  statistics[30] = Ph;
  statistics[31] = q;
  statistics[32] = Pa;
  statistics[33] = fcost;
  statistics[34] = housing_supply;

  statistics[43] = people;
  statistics[44] = owners_w_mortgage;
  statistics[45] = mortgage_debt;
  statistics[46] = housing_value;
  statistics[47] = default_rate;
  statistics[48] = housing_total;

  statistics[49] = yearspp;
  statistics[50] = Tretirement;
  statistics[51] = eetalab;
  statistics[52] = tthetalab;
  statistics[53] = refcost;
  statistics[54] = sunk;


  statistics[60] = heqne;
  statistics[61] = heq10;
  statistics[62] = heq20;
  statistics[63] = heq25;
  statistics[64] = heq30;
  statistics[65] = heq50;
  statistics[66] = heq75;
  statistics[67] = heq100;
  statistics[68] = hequal100;
  statistics[69] = avheq;

  statistics[80] = ownership;
  statistics[81] = av_owned/av_rented;
  statistics[82] = Tretirement;
  statistics[83] = yearspp;
  statistics[84] = av_M_payment;

  statistics[85] = exp_cons/exp_rent;

  statistics[90] = perc_renew_incshock/people_renew;
  statistics[91] = perc_renew_prishock/people_renew;

  statistics[94] = perc_defau_incshock_lowa/people_default_lowa;
  statistics[95] = perc_defau_prishock_lowa/people_default_lowa;
  statistics[96] = people_default_lowa/mortgagors;

  statistics[97] = perc_defau_incshock_higha/people_default_higha;
  statistics[98] = perc_defau_prishock_higha/people_default_higha;
  statistics[99] = people_default_higha/mortgagors;
  

  // for(int i=0; i<54; i++){
  //   cout << statistics[i] << endl;
  // }

  //--------------------------------------//
  //--------  Exporting arrays  ----------//
  //--------------------------------------//

  export_arrays(T, ny, nd, na, nh, nm,
                Value, Default, Renew, 
                Policym, Policya, 
                Policyh, Policyr, Policyl, Policyc, Pricing_guess,
                Pcond, Puncond,
                C, L, Income, H, R, M, A, DD, RR, Prepay, PP, Mp, Hp, KK, RRup, RRdn,
                NR, NetRent, OwnerOccupied, Def, Vivos, home_eq, home_eq_pdf, def_dist, Mort_debt, equity, condDD, 
                condDDdelta, DcondTxAxY, DcondAxY, DcondA, DcondHxM, LcondAxY, LcondAxYpeople, HcondAxY, McondAxY,
                statistics, error_funct);


  //--------------------------------------//
  //--------  Error computation ----------//
  //--------------------------------------//

  double error_equilib = 0.0;

  double ddeltaoutput = ddeltaf*default_rate;

  cout.precision(3);
  cout.setf(std::ios::fixed);
  
  Color::Modifier red(Color::FG_RED);
  Color::Modifier blue(Color::FG_BLUE);
  Color::Modifier green(Color::FG_GREEN);
  Color::Modifier def(Color::FG_DEFAULT);
  Color::Modifier yellow(Color::FG_YELLOW);

    if(error_funct == 1){
      // Etapa 1: equilibrio de estado estacionario inicial
      double def_target = yearspp*0.015;

      error_equilib = error_equilib + pow(housing_total-rent_total, 2.0);

      // error_equilib = error_equilib + 1000*pow(heqne - 0.018, 2.0);
      // error_equilib = error_equilib + 1000*pow(heq10 - 0.07, 2.0);
      // error_equilib = error_equilib + 1000*pow(heq20 - 0.141, 2.0);
      // error_equilib = error_equilib + 1000*pow(heq30 - 0.224, 2.0);

      cout << " " << endl;

      cout << "Estado estacionario inicial: " << endl;
      cout << setprecision(6) << "q = \t\t" << q << setprecision(3) << "\t | \t Rental: \t" << red << rent_total << "\t" << housing_total << def << "\t:Owned \t | \t Under 35:\t" << under_35 << "\t0.663" << endl;
      cout << setprecision(6) << "dmin = \t\t" << dmin << setprecision(3) << "\t | \t Ddelta: \t" << red << ddeltabar << "\t" << ddeltaoutput << def << "\t:Delou \t | \t 35-44:\t\t" << under_45 << "\t0.404" << endl;
      cout << setprecision(6) << "dmax = \t\t" << dmax << setprecision(3) << "\t | \t Own w/mort:\t" << red << owners_w_mortgage << "\t0.66" << def <<"\t:Data \t | \t 45-54:\t\t" << under_55 << "\t0.297" << endl;
      cout << setprecision(6) << "Ph = \t\t" << Ph << setprecision(3) << "\t | \t Mort debt:\t" << (mortgage_debt/housing_value) << "\t0.485\t:Data \t | \t 55-64:\t\t" << under_65 << "\t0.233" << endl;
      cout << setprecision(6) << "Ddelta =\t" << ddeltabar << setprecision(3) << "\t | \t Default:\t" << red << default_rate << "\t" << def_target << def << "\t:Data \t | \t 65-74:\t\t" << under_75 << "\t0.194" << endl;
      cout << setprecision(6) << "fcost =\t\t" << fcost << setprecision(3) << "\t | \t Ownership: \t" << green << ownership << "\t0.66" << def << "\t:Data \t | " << endl;
      cout << setprecision(6) << "bbeta =\t\t" << bbeta << setprecision(3) << "\t | \t Owned/rent: \t" << green << av_owned/av_rented << "\t1.5" << def << "\t:Data \t | \t Heq10:\t" << heq10 << "\tHe(-):\t" << heqne << endl;
      cout << setprecision(6) << "r     =\t\t" << r << setprecision(3) << "\t | \t Earn Ow/re: \t" << earn_owner/earn_renter << "\t2.05\t:Data \t | \t Heq20:\t" << heq20 << "\tHeq50:\t" << heq50 << endl;
      cout << setprecision(6) << "ppsi  =\t\t" << ppsi << setprecision(3) << "\t | \t Cons/rent: \t" << green << exp_cons/exp_rent << "\t1.89" << def << "\t:Data \t | \t Heq25:\t" << heq25 << "\tHeq75:\t" << heq75 << endl;
      cout << setprecision(6) << "ttrans=\t\t" << ittrans << setprecision(3) << "\t" << "\t | \t Net renters: \t" << net_renters << "\t0.35\t:Data \t | \t Heq30:\t" << heq30 << "\tHeq+:\t" << heqma << endl;
      cout << setprecision(6) << "interm=\t\t" << interm << setprecision(3) << "\t | \t Housing exp: \t"  << green << housing_value/total_income << "\t1.80" << def << "\t:Data \t | \t He100:\t" << hequal100 << "\t" << endl;
      cout << setprecision(6) << "yincom=\t\t" << total_income << setprecision(3) << "\t | \t Mean delta: \t" << mean_delta << "\t" << mean_delta_unc << "\t:Data" << endl;
      cout << "    \t" << "\t\t\t | \t Mean incom: \t" << mean_income << "\t" << mean_income_unc << "\t:Data \t | \t Aveq:\t" << avheq << "\t" << endl;
      cout << "    \t\t Underwat:\t" << underwater << "\t" << "\t: Default_under:\t" << def_underwater << "\t" << endl;
      cout << "    \t\t Overwat:\t" << overwater << "\t" << "\t: Default_over:\t" << def_overwater << "\t" << endl;
      cout << "    \t\t Pension:\t" << pension_desembols << "\t" << "\t: ss_revenues:\t" << ss_revenues << "\t"  << endl;
      cout << "    \t\t Av valu:\t" << av_value << endl;
      cout << "    \t\t Av pric:\t" << av_price << endl;
      cout << "    \t\t M/A:    \t" << Mtot/Atot << endl;
      cout << "    \t\t A/Y:    \t" << assetsToIncome << endl;
      cout << "    \t\t Av Mort \t" << av_M_payment << "\t\t\t LTV = \t" << avLTV << endl;
      cout << "    \t\t C:      \t" << agC << "\t\t\t Y =  " << agY << "\t\t\t F =  :" << agF << "\t\t\t A(r) =  " << retornos_a << "\t\t\t H(d) = " << retornos_h << endl;
      cout << "    \t\t A:      \t" << agA << "\t\t\t A2 =  :" << agA2 << endl;
      cout << green << "    \t\t Labor Supply:\t" << labor_supply << " (" << labor_supply/workers << ")" << def << endl;
      cout << "Lump A: \t" << lump_assets << ",\tLump H: \t" << lump_houses << ",\tLump $: " << lump_disco << ",\t LUMP: " << lump_assets+lump_houses+lump_disco << endl;
      cout << "y1 =\t" << def_y1/peo_y1 << endl;
      cout << "y2 =\t" << def_y2/peo_y2 << endl;
      cout << "y3 =\t" << def_y3/peo_y3 << endl;
      cout << "y4 =\t" << def_y4/peo_y4 << endl;
      cout << "y5 =\t" << def_y5/peo_y5 << endl;
      // cout << "sunk =\t" << sunk << endl;
      // cout << "oomega =\t" << oomega << endl;
      cout << " " << endl;

      cout << "ERROR: \t\t" << error_equilib << "       " << endl;
      cout << " " << endl;

    } else if(error_funct == 2){
      // Etapa 2: equilibrio de estado estacionario final
      error_equilib = error_equilib + pow(housing_total-rent_total, 2.0);
      error_equilib = error_equilib + pow(housing_total-housing_supply, 2.0);
      cout << " " << endl;
      cout << "Estado estacionario final: " << endl;
      cout << "Rental:  " << rent_total << ", Owner: " << housing_total << " (q = " << q << ")" << "       " << endl;
      cout << "Housing: " << housing_total << " (Ph = " << Ph << ") - Target: " << housing_supply << "       " << endl;
      cout << "Default rate: " << (default_rate/owners_w_mortgage) << " (ddeltabar = " << ddeltabar << ")" << "       " << endl;
      cout << "ERROR: " << error_equilib << "       " << endl;
      cout << " " << endl;
    } else if(error_funct == 3 | error_funct == 4){
      // Etapa 1: equilibrio de estado estacionario inicial

      double def_target = yearspp*0.015;

      // error_equilib = error_equilib + pow(housing_total-rent_total, 2.0);
      // error_equilib = error_equilib + 1000*pow(owners_w_mortgage - 0.66, 2.0);
      error_equilib = error_equilib + 100000*pow(default_rate - def_target, 2.0);
      // error_equilib = error_equilib + 1000*pow(ownership - 0.66, 2.0);

      // error_equilib = error_equilib + 1000*pow(heqne - 0.018, 2.0);
      // error_equilib = error_equilib + 1000*pow(heq10 - 0.07, 2.0);
      // error_equilib = error_equilib + 1000*pow(heq20 - 0.141, 2.0);
      // error_equilib = error_equilib + 1000*pow(heq30 - 0.224, 2.0);

      cout << " " << endl;

      cout << "Estado estacionario inicial: " << endl;
      cout << setprecision(6) << "q = \t\t" << q << setprecision(3) << "\t | \t Rental: \t" << red << rent_total << "\t" << housing_total << def << "\t:Owned \t | \t Under 35:\t" << under_35 << "\t0.663" << endl;
      cout << setprecision(6) << "dmin = \t\t" << dmin << setprecision(3) << "\t | \t Ddelta: \t" << red << ddeltabar << "\t" << ddeltaoutput << def << "\t:Delou \t | \t 35-44:\t\t" << under_45 << "\t0.404" << endl;
      cout << setprecision(6) << "dmax = \t\t" << dmax << setprecision(3) << "\t | \t Own w/mort:\t" << red << owners_w_mortgage << "\t0.66" << def <<"\t:Data \t | \t 45-54:\t\t" << under_55 << "\t0.297" << endl;
      cout << setprecision(6) << "Ph = \t\t" << Ph << setprecision(3) << "\t | \t Mort debt:\t" << (mortgage_debt/housing_value) << "\t0.485\t:Data \t | \t 55-64:\t\t" << under_65 << "\t0.233" << endl;
      cout << setprecision(6) << "Ddelta =\t" << ddeltabar << setprecision(3) << "\t | \t Default:\t" << red << default_rate << "\t" << def_target << def << "\t:Data \t | \t 65-74:\t\t" << under_75 << "\t0.194" << endl;
      cout << setprecision(6) << "fcost =\t\t" << fcost << setprecision(3) << "\t | \t Ownership: \t" << green << ownership << "\t0.66" << def << "\t:Data \t | " << endl;
      cout << setprecision(6) << "bbeta =\t\t" << bbeta << setprecision(3) << "\t | \t Owned/rent: \t" << green << av_owned/av_rented << "\t1.5" << def << "\t:Data \t | \t Heq10:\t" << heq10 << "\tHe(-):\t" << heqne << endl;
      cout << setprecision(6) << "r     =\t\t" << r << setprecision(3) << "\t | \t Earn Ow/re: \t" << earn_owner/earn_renter << "\t2.05\t:Data \t | \t Heq20:\t" << heq20 << "\tHeq50:\t" << heq50 << endl;
      cout << setprecision(6) << "ppsi  =\t\t" << ppsi << setprecision(3) << "\t | \t Cons/rent: \t" << green << exp_cons/exp_rent << "\t1.89" << def << "\t:Data \t | \t Heq25:\t" << heq25 << "\tHeq75:\t" << heq75 << endl;
      cout << setprecision(6) << "ttrans=\t\t" << ittrans << setprecision(3) << "\t" << "\t | \t Net renters: \t" << net_renters << "\t0.35\t:Data \t | \t Heq30:\t" << heq30 << "\tHeq+:\t" << heqma << endl;
      cout << setprecision(6) << "interm=\t\t" << interm << setprecision(3) << "\t | \t Housing exp: \t"  << green << housing_value/total_income << "\t1.80" << def << "\t:Data \t | \t He100:\t" << hequal100 << "\t" << endl;
      cout << setprecision(6) << "yincom=\t\t" << total_income << setprecision(3) << "\t | \t Mean delta: \t" << mean_delta << "\t" << mean_delta_unc << "\t:Data" << endl;
      cout << "    \t" << "\t\t\t | \t Mean incom: \t" << mean_income << "\t" << mean_income_unc << "\t:Data \t | \t Aveq:\t" << avheq << "\t" << endl;
      cout << "    \t\t Underwat:\t" << underwater << "\t" << "\t: Default_under:\t" << def_underwater << "\t" << endl;
      cout << "    \t\t Overwat:\t" << overwater << "\t" << "\t: Default_over:\t" << def_overwater << "\t" << endl;
      cout << "    \t\t Pension:\t" << pension_desembols << "\t" << "\t: ss_revenues:\t" << ss_revenues << "\t"  << endl;
      cout << "    \t\t Realtor:\t" << realtor_profits/people << endl;
      cout << "    \t\t Av valu:\t" << av_value << endl;
      cout << "    \t\t Av pric:\t" << av_price << endl;
      cout << "    \t\t M/A:    \t" << Mtot/Atot << endl;
      cout << "    \t\t A/Y:    \t" << assetsToIncome << endl;
      cout << "    \t\t Av Mort \t" << av_M_payment << "\t\t\t LTV = \t" << avLTV << endl;
      cout << "    \t\t C:      \t" << agC << "\t\t\t Y =  " << agY << "\t\t\t F =  :" << agF << "\t\t\t A(r) =  " << retornos_a << "\t\t\t H(d) = " << retornos_h << endl;
      cout << "    \t\t A:      \t" << agA << "\t\t\t A2 =  :" << agA2 << endl;
      cout << green << "    \t\t Labor Supply:\t" << labor_supply << " (" << labor_supply/workers << ")" << def << endl;
      cout << "Lump A: \t" << lump_assets << ",\tLump H: \t" << lump_houses << ",\tLump $: " << lump_disco << ",\t LUMP: " << lump_assets+lump_houses+lump_disco << endl;
      cout << "y1 =\t" << def_y1/peo_y1 << endl;
      cout << "y2 =\t" << def_y2/peo_y2 << endl;
      cout << "y3 =\t" << def_y3/peo_y3 << endl;
      cout << "y4 =\t" << def_y4/peo_y4 << endl;
      cout << "y5 =\t" << def_y5/peo_y5 << endl;
      // cout << "sunk =\t" << sunk << endl;
      // cout << "oomega =\t" << oomega << endl;
      cout << " " << endl;

      cout << "ERROR: \t\t" << error_equilib << "       " << endl;
      cout << " " << endl;

    } else if(error_funct == 20){
      // Etapa 1: equilibrio de estado estacionario inicial
      error_equilib = error_equilib + pow(housing_total-rent_total, 2.0);
      error_equilib = error_equilib + pow(housing_total - housing_supply, 2.0);
      cout << " " << endl;

      cout << "Estado estacionario inicial: " << endl;
      cout << setprecision(6) << "q = \t\t" << q << setprecision(3) << "\t | \t Rental: \t" << red << rent_total << "\t" << housing_total << def << "\t:Owned \t | \t Under 35:\t" << under_35 << "\t0.663" << endl;
      cout << setprecision(6) << "dmin = \t\t" << dmin << setprecision(3) << "\t | \t Ddelta: \t" << red << ddeltabar << "\t" << ddeltaoutput << def << "\t:Delou \t | \t 35-44:\t\t" << under_45 << "\t0.404" << endl;
      cout << setprecision(6) << "dmax = \t\t" << dmax << setprecision(3) << "\t | \t Own w/mort:\t" << red << owners_w_mortgage << "\t0.66" << def <<"\t:Data \t | \t 45-54:\t\t" << under_55 << "\t0.297" << endl;
      cout << setprecision(6) << "Ph = \t\t" << Ph << setprecision(3) << "\t | \t Mort debt:\t" << (mortgage_debt/housing_value) << "\t0.485\t:Data \t | \t 55-64:\t\t" << under_65 << "\t0.233" << endl;
      cout << setprecision(6) << "Ddelta =\t" << ddeltabar << setprecision(3) << "\t | \t Default:\t" << red << default_rate << "\t0.076" << def << "\t:Data \t | \t 65-74:\t\t" << under_75 << "\t0.194" << endl;
      cout << setprecision(6) << "fcost =\t\t" << fcost << setprecision(3) << "\t | \t Ownership: \t" << green << ownership << "\t0.66" << def << "\t:Data \t | " << endl;
      cout << setprecision(6) << "bbeta =\t\t" << bbeta << setprecision(3) << "\t | \t Owned/rent: \t" << green << av_owned/av_rented << "\t1.5" << def << "\t:Data \t | \t Heq10:\t" << heq10 << "\t" << endl;
      cout << setprecision(6) << "r     =\t\t" << r << setprecision(3) << "\t | \t Earn Ow/re: \t" << earn_owner/earn_renter << "\t2.05\t:Data \t | \t Heq20:\t" << heq20 << "\tHeq50:\t" << heq50 << endl;
      cout << setprecision(6) << "ppsi  =\t\t" << ppsi << setprecision(3) << "\t | \t Cons/rent: \t" << green << exp_cons/exp_rent << "\t1.89" << def << "\t:Data \t | \t Heq25:\t" << heq25 << "\tHeq75:\t" << heq75 << endl;
      cout << setprecision(6) << "ttrans=\t\t" << ittrans << setprecision(3) << "\t" << "\t | \t Net renters: \t" << net_renters << "\t0.35\t:Data \t | \t Heq30:\t" << heq30 << "\tHeq+:\t" << heqma << endl;
      cout << setprecision(6) << "interm=\t\t" << interm << setprecision(3) << "\t | \t Housing exp: \t" << green << housing_value/total_income << "\t1.80" << def << "\t:Data \t | \t He100:\t" << heq100 << "\t" << endl;
      cout << setprecision(6) << "yincom=\t\t" << total_income << setprecision(3) << "\t | \t Mean delta: \t" << mean_delta << "\t" << mean_delta_unc << "\t:Data" << endl;
      cout << "    \t\t\t" << "\t | \t Mean incom: \t" << mean_income << "\t" << mean_income_unc << "\t:Data \t | \t Aveq:\t" << avheq << "\t" << endl;
      cout << "    \t\t\t Underwat:\t" << underwater << "\t" << "\t: Default_under:\t" << def_underwater << "\t" << endl;
      cout << "y1 =\t" << def_y1/peo_y1 << endl;
      cout << "y2 =\t" << def_y2/peo_y2 << endl;
      cout << "y3 =\t" << def_y3/peo_y3 << endl;
      cout << "y4 =\t" << def_y4/peo_y4 << endl;
      cout << "y5 =\t" << def_y5/peo_y5 << endl;
      // cout << "sunk =\t" << sunk << endl;
      // cout << "oomega =\t" << oomega << endl;      cout << " " << endl;

      cout << "ERROR: \t\t" << error_equilib << "       " << endl;
      cout << " " << endl;

    } else if(error_funct == 5){
      // Etapa 1: equilibrio de estado estacionario inicial
      error_equilib = error_equilib + pow(housing_total-rent_total, 2.0);
      error_equilib = error_equilib + 10000*pow(ddeltabar - ddeltaoutput, 2.0);
      error_equilib = error_equilib + pow(housing_total-housing_supply, 2.0);
      cout << " " << endl;

      cout << "Estado estacionario final: " << endl;
      cout << "Rental:  " << rent_total << ", Owner: " << housing_total << " (q = " << q << ")" << "       " << endl;
      cout << "Housing: " << housing_total << " (Ph = " << Ph << ") - Target: " << housing_supply << "       " << endl;
      cout << "Owners w mortg: " << (owners_w_mortgage/people) << "       " << endl;
      cout << "Mortg debt: " << (mortgage_debt/housing_value) << " (dmin = " << dmin << ", dvar = " << dmax + 0.05 << ") - 0.43" << "       " << endl;
      cout << "Default rate: " << (default_rate/owners_w_mortgage) << " (dmax = " << dmax << ") - 0.0152" << "       " << endl;
      cout << "Delta bar output: " << ddeltaoutput << endl;
      cout << "ERROR: " << error_equilib << "       " << endl;
      cout << " " << endl;
    } else if(error_funct == 16){
      // Etapa 1: equilibrio de estado estacionario inicial
      error_equilib = error_equilib + 1000*pow(owners_w_mortgage - 0.66, 2.0);
      cout << " " << endl;
      // error_equilib = error_equilib + 1000*pow(owners_w_mortgage/people - 0.66, 2.0);

      cout << "Estado estacionario inicial: " << endl;
      cout << setprecision(6) << "q = \t\t" << q << setprecision(3) << "\t | \t Rental: \t" << red << rent_total << "\t" << housing_total << def << "\t:Owned \t | \t Under 35:\t" << under_35 << "\t0.663" << endl;
      cout << setprecision(6) << "dmin = \t\t" << dmin << setprecision(3) << "\t | \t Ddelta: \t" << red << ddeltabar << "\t" << ddeltaoutput << def << "\t:Delou \t | \t 35-44:\t\t" << under_45 << "\t0.404" << endl;
      cout << setprecision(6) << "dmax = \t\t" << dmax << setprecision(3) << "\t | \t Own w/mort:\t" << red << owners_w_mortgage << "\t0.66" << def <<"\t:Data \t | \t 45-54:\t\t" << under_55 << "\t0.297" << endl;
      cout << setprecision(6) << "Ph = \t\t" << Ph << setprecision(3) << "\t | \t Mort debt:\t" << (mortgage_debt/housing_value) << "\t0.485\t:Data \t | \t 55-64:\t\t" << under_65 << "\t0.233" << endl;
      cout << setprecision(6) << "Ddelta =\t" << ddeltabar << setprecision(3) << "\t | \t Default:\t" << red << default_rate << "\t0.076" << def << "\t:Data \t | \t 65-74:\t\t" << under_75 << "\t0.194" << endl;
      cout << setprecision(6) << "fcost =\t\t" << fcost << setprecision(3) << "\t | \t Ownership: \t" << green << ownership << "\t0.66" << def << "\t:Data \t | " << endl;
      cout << setprecision(6) << "bbeta =\t\t" << bbeta << setprecision(3) << "\t | \t Owned/rent: \t" << green << av_owned/av_rented << "\t1.5" << def << "\t:Data \t | \t Underwat:\t" << underwater << "\t" << endl;
      cout << setprecision(6) << "r     =\t\t" << r << setprecision(3) << "\t | \t Earn Ow/re: \t" << earn_owner/earn_renter << "\t2.05\t:Data \t | \t Underdef:\t" << def_underwater << "\t" << endl;
      cout << setprecision(6) << "ppsi  =\t\t" << ppsi << setprecision(3) << "\t | \t Cons/rent: \t" << green << exp_cons/exp_rent << "\t1.89" << def << "\t:Data" << endl;
      cout << setprecision(6) << "ttrans=\t\t" << ittrans << setprecision(3) << "\t" << "\t | \t Net renters: \t" << net_renters << "\t0.35\t:Data" << endl;
      cout << "    \t\t\t" << "\t | \t Housing exp: \t" << green << housing_value/total_income << "\t1.59" << def << "\t:Data" << endl;
      cout << "    \t\t\t" << "\t | \t Mean delta: \t" << mean_delta << "\t \t:Data" << endl;
      cout << "    \t\t\t" << "\t | \t Mean incom: \t" << mean_income << "\t \t:Data" << endl;
      cout << " " << endl;
      cout << "ERROR: \t\t" << error_equilib << "       " << endl;
      cout << " " << endl;

      // cout << "Estado estacionario inicial: " << endl;
      // cout << "Rental:  " << rent_total << ", Owner: " << housing_total << " (q = " << q << ")" << "       " << endl;
      // cout << "Housing: " << housing_total << " (Ph = " << Ph << ")" << "       " << endl;
      // cout << "Owners w mortg: " << (owners_w_mortgage/people) << "  - 0.66     " << endl;
      // cout << "Mortg debt: " << (mortgage_debt/housing_value) << " (dmin = " << dmin << ", dvar = " << dmax + 0.05 << ") - 0.485" << "       " << endl;
      // cout << "Default rate: " << (default_rate/owners_w_mortgage) << " (dmax = " << dmax << ") - 0.0152" << "       " << endl;
    //   cout << "Delta bar output: " << ddeltaoutput << endl;
    //   cout << "Ownership rate: " << ownership << " (0.66)" << endl;
    //   cout << "Owned/rented size: " << av_owned/av_rented << " (1.5)" << endl;
      // cout << "ERROR: " << error_equilib << "       " << endl;
      // cout << " " << endl;
    } 
  
  cout.precision(6);
  cout.setf(std::ios::fixed);


  *d_rental   = rent_total;
  *d_housing  = housing_total;

  return(error_equilib);
}



//======================================
//         Error and aggregates
//======================================

double error_sequence(parameters params, const int error_funct, const int ittrans,
                      const double *P, const double *survival, const double *hgrid, const double *rgrid, const double *mgrid,
                      const double *agrid, const double *dgrid, const double *ygrid, const double *eprocess, const double *lgrid,
                      const double *Value, const double *Value_equiv, const double *Policyc, const double *Pricing_guess, 
                      const int *Default, const int *Renew, const int *Policya, const int *Policyl, 
                      const int *Policym, const int *Policyh, const int *Policyr,
                      const double *Value_without, const double *Value_equiv_without, const double *Policyc_without, const double *Pricing_guess_without, 
                      const int *Default_without, const int *Renew_without, const int *Policya_without, const int *Policyl_without, 
                      const int *Policym_without, const int *Policyh_without, const int *Policyr_without,
                      const double *Pcond, const double *Puncond, const double *mortsubsidy, const double *repay_coeff,
                      const int Ttrans, const double rshock, const double Pashock, const double ppsishock, const double oomegashock, const double *incshock,
                      const double sunk_shock, const double interm_shock, const double ltax_shock, const int periods_shock,
                      const double *qpath, const double *Ppath,
                      double *error_accum, const int colorojo, const int baseline, 
                      const int *subs_eligible, const int *subs_target, const double prob_mistake,
                      const std::string tipo){



  // VFI parameters
  const int maxiter         = params.maxiter;
  const int uti             = params.uti;
  const double tol          = params.tol;
  const int T               = params.T;
  const int Tretirement     = params.Tretirement;

  // Grid for savings: a
  const int na      = params.na;
  const double amin = params.amin;
  const double amax = params.amax;

  // Grid for mortgages: m
  const int nm      = params.nm;
  const double mmin = params.mmin;
  const double mmax = params.mmax;

  // Grid for housing: h
  const int nh      = params.nh;
  const double hmin = params.hmin;
  const double hmax = params.hmax;

  // Grid for renting: r
  const int nr      = params.nr;
  const double rmin = params.rmin;
  const double rmax = params.rmax;

  // Grid for deoreciation: ddelta
  const int nd      = params.nd;
  const double dmin = params.dmin;
  const double dmax = params.dmax;

  // Grid for income shocks: y
  const int ny            = params.ny;
  const double ssigma_y   = params.ssigma_y;
  const double llambda_y  = params.llambda_y;
  const double m_y        = params.m_y;

  // Preferences
  const double ssigma = params.ssigma;
  const double rrho   = params.rrho;
  const double ppsi   = params.ppsi;
  const double bbeta  = params.bbeta;
  const double kkappa = params.kkappa;


  // Equilibrium objects
  const double ddeltabar      = params.ddeltabar_today;
  const double ddeltaf        = params.ddeltaf;
  const double r              = params.r;
  const double Ph             = params.Ph_today;
  const double q              = params.q;
  const double Pa             = params.Pa;
  const double fcost          = params.fcost;
  // const double lumpsum        = params.lumpsum;
  // const double mortsubsidy    = params.mortsubsidy;
  const double housing_supply = params.housing_supply;

  const double oomega     = params.oomega;
  const double sunk       = params.sunk;
  const double interm     = params.interm;
  const double lumpsum    = params.lumpsum;

  const double sstax          = params.sstax;
  const double ltax           = params.ltax;

  double *d_rental    = params.d_rental;
  double *d_housing   = params.d_housing;

  int ind;


  //--------------------------------------//
  //----  Aggregate paths in economy  ----//
  //--------------------------------------//

  // Aggregate stats
  double default_rate       = 0.0;
  double refinance_rate     = 0.0;
  double housing_total      = 0.0;
  double rent_total         = 0.0;
  double ownership          = 0.0;
  double mortgagors         = 0.0;
  double mort_creation      = 0.0;
  double savings            = 0.0;
  double average_Pm         = 0.0;
  double people             = 0.0;
  double subsidy_disburs    = 0.0;
  // double total_income       = 0.0;

  double ss_revenues        = 0.0;
  double l_revenues         = 0.0;

  double DTI[T];
  double DTI_y[T];
  double DTI_m[T];
  double DTI_people[T];

  // double Income[T];

  double pipol[T];

  for(int i = 0; i<T; i++){
    DTI[i] = 0.0;
    DTI_y[i] = 0.0;
    DTI_m[i] = 0.0;
    DTI_people[i] = 0.0;

    // Income[i] = 0.0;

    pipol[i] = 0.0;
  }

  double high_DTI = 0.0;
  double all_DTI  = 0.0;
  double DebtToInc = 0.0;
  double DebtToPeo = 0.0;

  double underwater = 0.0;
  // double mortgagors = 0.0;

  double av_value = 0.0;
  double av_value_eq = 0.0;

  double labor_supply = 0.0;

  int indsubs = 0;

  double incorrect_subs = 0.0;
  double correct_subs   = 0.0;

  double realtor_profits = 0.0;

  double default_rate_def = 0.0;
  double default_rate_ndef = 0.0;

  double av_a = 0.0;
  double av_m = 0.0;
  double av_h = 0.0;
  // double av_y = 0.0;

  double mean_def_with = 0.0;
  double ref_baja      = 0.0;
  double mean_mort     = 0.0;
  double mean_mort_bef = 0.0;
  double mean_mort_out = 0.0;
  double mean_def_out  = 0.0;
  double mean_ref      = 0.0;
  double pip_with      = 0.0;
  double pip_out       = 0.0;

  double home_equity = 0.0;
  double home_eq[130];

  for(int i = 0; i<130; i++){
    home_eq[i] = 0.0;
  }

  double av_value_seq[15];
  double av_value_seq_eq[15];
  double av_pip_seq[15];

  for(int i = 0; i<15; i++){
    av_value_seq[i] = 0.0;
    av_value_seq_eq[i] = 0.0;
    av_pip_seq[i] = 0.0;
  }

  double av_value_age[15];
  double av_value_seq_age[15];
  double av_pip_age[15];

  for(int i = 0; i<15; i++){
    av_value_age[i] = 0.0;
    av_value_seq_age[i] = 0.0;
    av_pip_age[i] = 0.0;
  }

  double wealth_dist[10];
  double av_value_wealth[10];
  double av_value_wealth_eq[10];
  double av_pip_wealth[10];

  for(int i = 0; i<10; i++){
    wealth_dist[i] = 0.0;
    av_value_wealth[i] = 0.0;
    av_value_wealth_eq[i] = 0.0;
    av_pip_wealth[i] = 0.0;
  }

  double av_value_nonowners = 0.0;
  double av_value_nonowners_eq = 0.0;

  double ownership_eq = 0.0;

  double PP[T];
  double DD[T];
  double Eq[T];
  double Pip_Eq[T];
  double Pip_age[T];

  for(int it=0; it<T; it++){
    PP[it]            = 0.0;
    DD[it]            = 0.0;
    Eq[it]            = 0.0;
    Pip_age[it]            = 0.0;
    Pip_Eq[it]            = 0.0;
  }

  double DDa[na];
  double Pip_a[na];

  for(int ia=0; ia<na; ia++){
    DDa[ia]            = 0.0;
    Pip_a[ia]            = 0.0;
  }

  double err = 0.0;

  double av_eq_def = 0.0;
  double av_def    = 0.0;
  double max_eq_def = -10.0;

  double wealth = 0.0;

  int ssigma_welfare = 10;

  double av_value_sigma[11];
  double av_value_sigma_eq[11];
  double pipol_prov[11];

  for(int it=0; it<11; it++){
    av_value_sigma[it]            = 0.0;
    av_value_sigma_eq[it]            = 0.0;
    pipol_prov[it]            = 0.0;
  }

  // for(int it=0; it<T; it++){
  //   for(int iy=0; iy<ny; iy++){
  //     for(int id=0; id<nd; id++){
  //       for(int ih=0; ih<nh; ih++){
  //         for(int im=0; im<nm; im++){
  //           for(int ia=0; ia<na; ia++){
  //             ind     = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

  //             if(iy ==1 && id == 1 && ih == 1 && im==1 && ia==1){
  //               cout << Value[ind] << endl;
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  double Phseq[11];
  Phseq[0] = 1.0;
  Phseq[1] = 0.774241;
  Phseq[2] = 0.820033;
  Phseq[3] = 0.888033;
  Phseq[4] = 0.999999;
  Phseq[5] = 0.997753;
  Phseq[6] = 0.995629;
  Phseq[7] = 0.996954;
  Phseq[8] = 0.997434;
  Phseq[9] = 0.999592;
  Phseq[10] = 0.999937;

  double home_eq_counterfactual = 0.0;



  for(int it=0; it<T; it++){
    for(int iy=0; iy<ny; iy++){
      for(int id=0; id<nd; id++){
        for(int ih=0; ih<nh; ih++){
          for(int im=0; im<nm; im++){
            for(int ia=0; ia<na; ia++){
              ind     = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
              indsubs = it*ny*nh*nm*na + iy*nh*nm*na + ih*nm*na + im*na + ia;

              av_a = av_a + Puncond[ind]*agrid[ia];
              av_m = av_m + Puncond[ind]*mgrid[im];
              av_h = av_h + Puncond[ind]*hgrid[ih];

              PP[it] = PP[it] + Pcond[ind]*Pricing_guess[ind];

              if(im > 0){
                // subsidy_disburs = subsidy_disburs + mortsubsidy[ind]*(Puncond[ind]*(1-Renew[ind])*(1-Default[ind])*mgrid[im] +
                //                                                       Puncond[ind]*Renew[ind]*(1-Default[ind])*mgrid[im]*(1+repay_coeff[it]));

                // Le doy a los correctos (1-prob_mistake) y a los incorrectos prob_mistake
                // subsidy_disburs = subsidy_disburs + Puncond[ind]*(1+repay_coeff[it])*mortsubsidy[indsubs]*mgrid[im]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake) +
                //                                     Puncond[ind]*(1+repay_coeff[it])*mortsubsidy[indsubs]*mgrid[im]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake;
                subsidy_disburs = subsidy_disburs + Puncond[ind]*mortsubsidy[indsubs]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake) +
                                                    Puncond[ind]*mortsubsidy[indsubs]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake;

                correct_subs   = correct_subs + Puncond[ind]*mortsubsidy[indsubs]*mgrid[im]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake);
                incorrect_subs = incorrect_subs + Puncond[ind]*mortsubsidy[indsubs]*mgrid[im]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake;

                realtor_profits = realtor_profits + Puncond[ind]*Default[ind]*(1-Renew[ind])*Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih]*sunk;

                // mortgagors = mortgagors + Puncond[ind];

                // Guys underwater
                if(mgrid[im]*(1+repay_coeff[it]) >  Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih] & ih > 0){
                  underwater = underwater + Puncond[ind];
                }


              }

              labor_supply = labor_supply + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*lgrid[Policyl[ind]] +
                                            Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*lgrid[Policyl_without[ind]] +
                                            Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*lgrid[Policyl[ind]] + 
                                            Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*lgrid[Policyl_without[ind]] + 
                                            Puncond[ind]*(1-subs_eligible[ind])*lgrid[Policyl_without[ind]];


              if(im > 0){
                mortgagors = mortgagors + Puncond[ind];

                Pip_age[it] = Pip_age[it] + Puncond[ind];

                Pip_a[ia] = Pip_a[ia] + Puncond[ind];
                
                default_rate = default_rate + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*(double)Default[ind] + 
                                              Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*(double)Default_without[ind] +
                                              Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*(double)Default[ind] +
                                              Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*(double)Default_without[ind] +
                                              Puncond[ind]*(1-subs_eligible[ind])*(double)Default_without[ind];

                DD[it] = DD[it] + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*(double)Default[ind] + 
                                  Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*(double)Default_without[ind] +
                                  Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*(double)Default[ind] +
                                  Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*(double)Default_without[ind] +
                                  Puncond[ind]*(1-subs_eligible[ind])*(double)Default_without[ind];

                DDa[ia] = DDa[ia] + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*(double)Default[ind] + 
                                  Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*(double)Default_without[ind] +
                                  Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*(double)Default[ind] +
                                  Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*(double)Default_without[ind] +
                                  Puncond[ind]*(1-subs_eligible[ind])*(double)Default_without[ind];
                

                if(subs_eligible[ind] == 1 && subs_target[ind] == 1){
                  mean_def_with = mean_def_with + Puncond[ind]*(double)Default[ind];

                  mean_ref = mean_ref + Puncond[ind]*(double)Renew[ind];

                  if(Policym[ind]<im){
                    ref_baja = ref_baja + Puncond[ind];
                  }

                  mean_mort = mean_mort + Puncond[ind]*mgrid[Policym[ind]];
                  mean_mort_bef = mean_mort_bef + Puncond[ind]*mgrid[im];
                  mean_mort_out = mean_mort_out + Puncond[ind]*mgrid[Policym_without[ind]];

                  mean_def_out = mean_def_out + Puncond[ind]*(double)Default_without[ind];

                  pip_with = pip_with + Puncond[ind];
                  pip_out = pip_out + Puncond[ind];
                }

                if(subs_eligible[ind] == 1 && subs_target[ind] == 1 && Default_without[ind] == 0){
                  err = err + Puncond[ind];
                }

                // mean_def_out = mean_def_out + Puncond[ind]*(subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*(double)Default_without[ind] + 
                //                                             (1-subs_eligible[ind])*(double)Default_without[ind]);


                default_rate_def = default_rate_def + Puncond[ind]*(double)Default[ind];
                default_rate_ndef = default_rate_ndef + Puncond[ind]*(double)Default_without[ind];


                refinance_rate = refinance_rate + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*(double)Renew[ind] + 
                                                  Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*(double)Renew_without[ind] +
                                                  Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*(double)Renew[ind] +
                                                  Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*(double)Renew_without[ind] +
                                                  Puncond[ind]*(1-subs_eligible[ind])*(double)Renew_without[ind];
              }

              if(Policym[ind] > 0 && Renew[ind] == 1 && subs_eligible[ind] == 1 && subs_target[ind] == 1){
                mort_creation = mort_creation + Puncond[ind]*(1-prob_mistake);
              } else if(Policym[ind] > 0 && Renew[ind] == 1 && subs_eligible[ind] == 1 && subs_target[ind] == 0){
                mort_creation = mort_creation + Puncond[ind]*prob_mistake;
              } else if(Policym_without[ind] > 0 && Renew_without[ind] == 1 && subs_eligible[ind] == 1 && subs_target[ind] == 1){
                mort_creation = mort_creation + Puncond[ind]*prob_mistake;
              } else if(Policym_without[ind] > 0 && Renew_without[ind] == 1 && subs_eligible[ind] == 1 && subs_target[ind] == 0){
                mort_creation = mort_creation + Puncond[ind]*(1-prob_mistake);
              } else if(Policym_without[ind] > 0 && Renew_without[ind] == 1 && subs_eligible[ind] == 0){
                mort_creation = mort_creation + Puncond[ind];
              }


              if(Policyh[ind] > 0 && subs_eligible[ind] == 1 && subs_target[ind] == 1){
                ownership = ownership + Puncond[ind]*(1-prob_mistake);
              } else if(Policyh[ind] > 0 && subs_eligible[ind] == 1 && subs_target[ind] == 0){
                ownership = ownership + Puncond[ind]*prob_mistake;
              } else if(Policyh_without[ind] > 0 && subs_eligible[ind] == 1 && subs_target[ind] == 1){
                ownership = ownership + Puncond[ind]*prob_mistake;
              } else if(Policyh_without[ind] > 0 && subs_eligible[ind] == 1 && subs_target[ind] == 0){
                ownership = ownership + Puncond[ind]*(1-prob_mistake);
              } else if(Policyh_without[ind] > 0 && subs_eligible[ind] == 0){
                ownership = ownership + Puncond[ind];
              }



              if(it < Tretirement){
                ss_revenues = ss_revenues + Puncond[ind]*ygrid[iy]*eprocess[it]*(1-incshock[it])*sstax*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*lgrid[Policyl[ind]] + 
                                            Puncond[ind]*ygrid[iy]*eprocess[it]*(1-incshock[it])*sstax*subs_eligible[ind]*subs_target[ind]*prob_mistake*lgrid[Policyl_without[ind]] + 
                                            Puncond[ind]*ygrid[iy]*eprocess[it]*(1-incshock[it])*sstax*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*lgrid[Policyl[ind]] + 
                                            Puncond[ind]*ygrid[iy]*eprocess[it]*(1-incshock[it])*sstax*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*lgrid[Policyl_without[ind]] +
                                            Puncond[ind]*ygrid[iy]*eprocess[it]*(1-incshock[it])*sstax*(1-subs_eligible[ind])*lgrid[Policyl_without[ind]];

                l_revenues  = l_revenues  + Puncond[ind]*ygrid[iy]*eprocess[it]*(1-incshock[it])*ltax*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*lgrid[Policyl[ind]] +
                                            Puncond[ind]*ygrid[iy]*eprocess[it]*(1-incshock[it])*ltax*subs_eligible[ind]*subs_target[ind]*prob_mistake*lgrid[Policyl_without[ind]] +
                                            Puncond[ind]*ygrid[iy]*eprocess[it]*(1-incshock[it])*ltax*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*lgrid[Policyl[ind]] +
                                            Puncond[ind]*ygrid[iy]*eprocess[it]*(1-incshock[it])*ltax*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*lgrid[Policyl_without[ind]] +
                                            Puncond[ind]*ygrid[iy]*eprocess[it]*(1-incshock[it])*ltax*(1-subs_eligible[ind])*lgrid[Policyl_without[ind]];
              }


              // // I only look for payments of already originated mortgages - exlude mortgages originated at same period
              // if(im>0 && ih>0 && Policyl_without[ind]>0 && Renew_without[ind] == 0){            

              //   // DTI[it] = DTI[it]               + Puncond[ind]*mgrid[im]*(1+repay_coeff[it])/(ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]]);
              //   // DTI_y[it] = DTI_y[it]           + Puncond[ind]*(ygrid[iy]*eprocess[it]*lgrid[Policyl[ind]]);
              //   // DTI_m[it] = DTI_m[it]           + Puncond[ind]*mgrid[im]*(1+repay_coeff[it]);
              //   // DTI_people[it] = DTI_people[it] + Puncond[ind];

              //   DTI[it] = DTI[it]               + Puncond[ind]*(1-mortsubsidy[ind])*mgrid[im]/(-lumpsum + ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl_without[ind]] + agrid[ia]*r + q*hgrid[ih]);
              //   DTI_y[it] = DTI_y[it]           + Puncond[ind]*(-lumpsum + ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl_without[ind]] + agrid[ia]*r + q*hgrid[ih]);
              //   DTI_m[it] = DTI_m[it]           + Puncond[ind]*(1-mortsubsidy[ind])*mgrid[im];
              //   DTI_people[it] = DTI_people[it] + Puncond[ind];

              //   DebtToInc = DebtToInc + Puncond[ind]*(1-mortsubsidy[ind])*mgrid[im]/(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl_without[ind]] + agrid[ia]*r + q*hgrid[ih]);
              //   DebtToPeo = DebtToPeo + Puncond[ind];

              //   // Income includes rent and interests
              //   Income[it] = Income[it] + Pcond[ind]*(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl_without[ind]] + agrid[ia]*r + q*hgrid[ih]);

              //   all_DTI = all_DTI + Puncond[ind];
              //   if((1-mortsubsidy[ind])*mgrid[im]/(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl_without[ind]] + agrid[ia]*r + q*hgrid[ih]) > 0.55){
              //     high_DTI = high_DTI + Puncond[ind];
              //   }

              // }

              if(ih > 0){
                ownership_eq = ownership_eq + Puncond[ind];
                home_equity = (Ph*(1-dgrid[id])*hgrid[ih] - mgrid[im]*(1+repay_coeff[it]))/(Ph*(1-dgrid[id])*hgrid[ih]);
                home_eq_counterfactual = (Phseq[ittrans]*(1-dgrid[id])*hgrid[ih] - mgrid[im]*(1+repay_coeff[it]))/(Phseq[ittrans]*(1-dgrid[id])*hgrid[ih]);

                Eq[it] = Eq[it] + Puncond[ind]*home_equity;
                Pip_Eq[it] = Pip_Eq[it] + Puncond[ind];
              } else{
                home_equity = 0.0;
                home_eq_counterfactual = 0.0;
              }

              if(ih > 0 && Default[ind] == 1){
                av_def = av_def + Puncond[ind];
                av_eq_def = av_eq_def + home_equity*Puncond[ind];

                if(home_equity > max_eq_def){
                  max_eq_def = home_equity;
                }
              }

              for(int i = -30; i<100; i++){
                if(ih > 0 && home_equity <= (1+(double)i)/100){
                  home_eq[i+30] = home_eq[i+30] + Puncond[ind];
                } 
              }

              wealth = Pa*agrid[ia] + Phseq[ittrans]*(1-dgrid[id])*hgrid[ih] - mgrid[im]*(1+repay_coeff[it]);

              for(int i = -1; i<8; i++){
                if(wealth <= 0.1*(1+(double)i)+0.1 && wealth > 0.1*((double)i)+0.1){
                  wealth_dist[i+3] = wealth_dist[i+3] + Puncond[ind];

                  av_value_wealth[i+3] = av_value_wealth[i+3] + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Value[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Value_without[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Value[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Value_without[ind] +
                                                          Puncond[ind]*(1-subs_eligible[ind])*Value_without[ind];

                  av_value_wealth_eq[i+3] = av_value_wealth_eq[i+3] + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Value_equiv[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Value_equiv_without[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Value_equiv[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Value_equiv_without[ind] +
                                                                Puncond[ind]*(1-subs_eligible[ind])*Value_equiv_without[ind];
                  
                  av_pip_wealth[i+3] = av_pip_wealth[i+3] + Puncond[ind];
                }
              }

              // for(ssigma_welfare=0; ssigma_welfare<11; ssigma_welfare++){
                
              //   av_value_sigma[ssigma_welfare] = av_value_sigma[ssigma_welfare] + (Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*pow(Value[ind], 1-(double)ssigma_welfare) +
              //                                                         Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*pow(Value_without[ind], 1-(double)ssigma_welfare) +
              //                                                         Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*pow(Value[ind], 1-(double)ssigma_welfare) +
              //                                                         Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*pow(Value_without[ind], 1-(double)ssigma_welfare) +
              //                                                         Puncond[ind]*(1-subs_eligible[ind])*pow(Value_without[ind], 1-(double)ssigma_welfare)) / (1-(double)ssigma_welfare);

              //   av_value_sigma_eq[ssigma_welfare] = av_value_sigma_eq[ssigma_welfare] + (Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*pow(Value_equiv[ind], 1-(double)ssigma_welfare) +
              //                                                         Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*pow(Value_equiv_without[ind], 1-(double)ssigma_welfare) +
              //                                                         Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*pow(Value_equiv[ind], 1-(double)ssigma_welfare) +
              //                                                         Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*pow(Value_equiv_without[ind], 1-(double)ssigma_welfare) +
              //                                                         Puncond[ind]*(1-subs_eligible[ind])*pow(Value_equiv_without[ind], 1-(double)ssigma_welfare)) / (1-(double)ssigma_welfare);
              // }

              // Average total value
              av_value = av_value + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Value[ind] +
                                    Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Value_without[ind] +
                                    Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Value[ind] +
                                    Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Value_without[ind] +
                                    Puncond[ind]*(1-subs_eligible[ind])*Value_without[ind];

              av_value_eq = av_value_eq + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Value_equiv[ind] +
                                          Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Value_equiv_without[ind] +
                                          Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Value_equiv[ind] +
                                          Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Value_equiv_without[ind] +
                                          Puncond[ind]*(1-subs_eligible[ind])*Value_equiv_without[ind];

              // Average value by equity
              for(int i = -2; i<5; i++){
                if(ih > 0 && home_eq_counterfactual > 20*((double)i)/100 && home_eq_counterfactual <= 20*(1+(double)i)/100){
                  av_value_seq[i+2] = av_value_seq[i+2] + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Value[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Value_without[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Value[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Value_without[ind] +
                                                          Puncond[ind]*(1-subs_eligible[ind])*Value_without[ind];

                  av_value_seq_eq[i+2] = av_value_seq_eq[i+2] + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Value_equiv[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Value_equiv_without[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Value_equiv[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Value_equiv_without[ind] +
                                                                Puncond[ind]*(1-subs_eligible[ind])*Value_equiv_without[ind];

                  av_pip_seq[i+2] = av_pip_seq[i+2] + Puncond[ind];
                } 
              }




              if(ih == 0){
                av_value_nonowners = av_value_nonowners + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Value[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Value_without[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Value[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Value_without[ind] +
                                                          Puncond[ind]*(1-subs_eligible[ind])*Value_without[ind];

                av_value_nonowners_eq = av_value_nonowners_eq + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Value_equiv[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Value_equiv_without[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Value_equiv[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Value_equiv_without[ind] +
                                                                Puncond[ind]*(1-subs_eligible[ind])*Value_equiv_without[ind];
              }


              for(int i = 0; i<6; i++){
                if(ih > 0 && it > 5*i && it <= 5*(1+i)){
                  av_value_age[i] = av_value_age[i] + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Value[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Value_without[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Value[ind] +
                                                          Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Value_without[ind] +
                                                          Puncond[ind]*(1-subs_eligible[ind])*Value_without[ind];
                  av_value_seq_age[i] = av_value_seq_age[i] + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Value_equiv[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Value_equiv_without[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Value_equiv[ind] +
                                                                Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Value_equiv_without[ind] +
                                                                Puncond[ind]*(1-subs_eligible[ind])*Value_equiv_without[ind];
                  av_pip_age[i] = av_pip_age[i] + Puncond[ind];
                } 
              }


              people = people + Puncond[ind];

              housing_total = housing_total + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*hgrid[Policyh[ind]] + 
                                              Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*hgrid[Policyh_without[ind]] +
                                              Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*hgrid[Policyh[ind]] +
                                              Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*hgrid[Policyh_without[ind]] +
                                              Puncond[ind]*(1-subs_eligible[ind])*hgrid[Policyh_without[ind]];

              rent_total = rent_total + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*rgrid[Policyr[ind]] +
                                        Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*rgrid[Policyr_without[ind]] +
                                        Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*rgrid[Policyr[ind]] +
                                        Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*rgrid[Policyr_without[ind]] +
                                        Puncond[ind]*(1-subs_eligible[ind])*rgrid[Policyr_without[ind]];

              savings = savings + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*agrid[Policya[ind]] +
                                  Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*agrid[Policya_without[ind]] +
                                  Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*agrid[Policya[ind]] +
                                  Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*agrid[Policya_without[ind]] +
                                  Puncond[ind]*(1-subs_eligible[ind])*agrid[Policya_without[ind]];

              average_Pm = average_Pm + Puncond[ind]*subs_eligible[ind]*subs_target[ind]*(1-prob_mistake)*Pricing_guess[ind] +
                                        Puncond[ind]*subs_eligible[ind]*subs_target[ind]*prob_mistake*Pricing_guess_without[ind] +
                                        Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*prob_mistake*Pricing_guess[ind] +
                                        Puncond[ind]*subs_eligible[ind]*(1-subs_target[ind])*(1-prob_mistake)*Pricing_guess_without[ind] +
                                        Puncond[ind]*(1-subs_eligible[ind])*Pricing_guess_without[ind];

              // total_income = total_income + Puncond[ind]*ygrid[iy]*(1-incshock[it])*eprocess[it];

              pipol[it] = pipol[it] + Puncond[ind];

            }
          }
        }
      }
    }
  }


  for(int i = -30; i<100; i++){
    home_eq[i+30] = home_eq[i+30]/ownership_eq;
  }

  for(int i = 0; i<T; i++){
    DTI[i] = DTI[i] / DTI_people[i];
    DTI_y[i] = DTI_y[i] / DTI_people[i];
    DTI_m[i] = DTI_m[i] / DTI_people[i];

    // cout << pipol[i] << endl;
    // cout << DTI[i] << " - " << DTI_y[i] << " - " << Income[i] << " - " << DTI_m[i] << endl;
  }

  for(int i = 0; i<T; i++){
    DD[i] = DD[i]/Pip_age[i];
    Eq[i] = Eq[i]*Pip_Eq[i];
  }

  for(int i = 0; i<na; i++){
    DDa[i] = DDa[i]/Pip_a[i];
  }

  underwater = underwater/mortgagors;

  high_DTI = high_DTI / all_DTI;
  // cout << "mortgages >= " << default_rate << endl;
  DebtToInc = DebtToInc/DebtToPeo;

  // cout << "media con subsidio = " << mean_def_with/pip_with << ", mortgage mean = " << mean_mort/pip_with << ", mean before = " << mean_mort_bef/pip_with << ", renew = " << mean_ref/pip_with << ", perc. baja = " << ref_baja/pip_with << endl;
  // cout << "media sin subsidio = " << mean_def_out/pip_out << ", mortgage mean = " << mean_mort_out/pip_out << endl;


  // cout << "equity of default = " << max_eq_def << endl;
  // cout << "Err = " << err << endl;

  // cout << "A = " << av_a << endl;
  // cout << "M = " << av_m << endl;
  // cout << "H = " << av_h << endl;

  default_rate = default_rate/mortgagors;
  refinance_rate = refinance_rate/mortgagors;

  // Export arrays
  double statistics[100];
  statistics[0] = tol;
  statistics[1] = uti;
  statistics[2] = maxiter;
  statistics[3] = T;

  statistics[4] = na;
  statistics[5] = amin;
  statistics[6] = amax;

  statistics[7] = nm;
  statistics[8] = mmin;
  statistics[9] = mmax;

  statistics[10] = nh;
  statistics[11] = hmin;
  statistics[12] = hmax;

  statistics[13] = nr;
  statistics[14] = rmin;
  statistics[15] = rmax;

  statistics[16] = nd;
  statistics[17] = dmin;
  statistics[18] = dmax;

  statistics[19] = ny;
  statistics[20] = ssigma_y;
  statistics[21] = llambda_y;
  statistics[22] = m_y;

  statistics[23] = ssigma;
  statistics[24] = rrho;
  statistics[25] = ppsi;
  statistics[26] = bbeta;
  statistics[27] = kkappa;
  statistics[28] = ddeltabar;
  statistics[29] = r;
  statistics[30] = Ph;
  statistics[31] = q;
  statistics[32] = Pa;
  statistics[33] = fcost;
  statistics[34] = housing_supply;

  // for(int it=0; it<35; it++){
  //   cout << statistics[it] << endl;
  // }


  export_basic(T, ny, nd, na, nh, nm, 
                  Value_without, Default, Renew, Policym, Policya, Policyh, Policyr, Policyl, Policyc, Pricing_guess, 
                  Ttrans, rshock, Pashock, ppsishock, oomegashock, sunk_shock, interm_shock, ltax_shock, periods_shock, qpath, Ppath,
                  Pcond, Puncond, home_eq, statistics, PP, DD, DDa, Eq, ittrans,
                  default_rate, tipo);

  //--------------------------------------//
  //--------  Exporting arrays  ----------//
  //--------------------------------------//

  double error_equilib = 0.0;

  double ddeltaoutput = ddeltaf*default_rate;

  cout.precision(3);
  cout.setf(std::ios::fixed);
  
  Color::Modifier red(Color::FG_RED);
  Color::Modifier blue(Color::FG_BLUE);
  Color::Modifier green(Color::FG_GREEN);
  Color::Modifier def(Color::FG_DEFAULT);
  Color::Modifier yellow(Color::FG_YELLOW);


  double *derror_accum = error_accum;

  int summarize = 1;

  // *derror_accum = *derror_accum + pow(housing_total-rent_total, 2.0);
  // *derror_accum = *derror_accum + 10000*pow(ddeltabar - ddeltaoutput, 2.0);

  if(error_funct == 3){

    *derror_accum = *derror_accum + pow(housing_total-rent_total, 2.0);
    *derror_accum = *derror_accum + pow(housing_total - housing_supply, 2.0);
    // *derror_accum = *derror_accum + 1000*pow(ddeltabar - ddeltaoutput, 2.0);

    cout << setprecision(6) << "q = \t\t" << q << setprecision(3) << "\t | \t Rental: \t" << rent_total << "\t" << housing_total << "\tOwned |\t " << green << "\e[1m Error: " << *derror_accum << "\e[0m" << yellow << "\t\t Subsidy: \t" << subsidy_disburs << "   (error = " << incorrect_subs/correct_subs << ")" << def << endl;
    cout << setprecision(6) << "ttrans= \t\t" << ittrans << setprecision(3) << "\t | \t Ddelta: \t" << ddeltabar << "\t" << ddeltaoutput << green << "\t       ----------------------" << def << endl;
    cout << setprecision(6) << "Ph = \t\t" << Ph << setprecision(3) << "\t | \t def:\t" << default_rate << "\tRef:\t" << refinance_rate << "\tMort:\t" << mortgagors << "\tM_cre:\t" << mort_creation << endl;
    cout << setprecision(6) << "r = \t\t" << r << setprecision(3) << "\t | \t sav:\t" << savings << endl;
    cout << setprecision(6) << "Pa  = \t\t" << Pa << setprecision(3) << "\t | \t av_Pm:\t" << average_Pm << endl;
    cout << setprecision(6) << "ppsi = \t\t" << ppsi << setprecision(3) << "\t | \t Ph/q:\t" << red << Ph/q << def << endl;
    cout << setprecision(6) << "tax  = \t\t" << lumpsum << setprecision(3) << endl;
    cout << setprecision(6) << "sunk = \t\t" << sunk << setprecision(3) << "\t | \t inter:\t" << interm << "\t | \t omega:\t" << oomega << "\t | \t lump:\t" << lumpsum << endl;
    cout << "----------------------------------------------------------------------------------------------------" << endl;

  } else if(error_funct == 4){

    if(ittrans < Ttrans-1){
      *derror_accum = *derror_accum + pow(housing_total-rent_total, 2.0);
      *derror_accum = *derror_accum + pow(housing_total - housing_supply, 2.0);
    }

    if(colorojo == 1){

      if(summarize == 0){
        cout << setprecision(6) << red << "q = \t\t" << q << setprecision(3) << "\t | \t Rental: \t" << rent_total << "\t" << housing_total << "\tOwned |\t " << green << "\e[1m Error: " << *derror_accum << "\e[0m" << yellow << "\t\t Subsidy: \t" << subsidy_disburs << "   (error = " << incorrect_subs/correct_subs << ")" << red << endl;
        cout << setprecision(6) << red << "ttrans= \t\t" << ittrans << setprecision(3) << "\t | \t Ddelta: \t" << ddeltabar << "\t" << ddeltaoutput << green << "\t       ----------------------" << red << endl;
        cout << setprecision(6) << red << "Ph = \t\t" << Ph << setprecision(3) << "\t | \t def:\t" << default_rate << "\tRef:\t" << refinance_rate << "\tMort:\t" << mortgagors << "\tM_cre:\t" << mort_creation << endl;
        cout << setprecision(6) << red << "r = \t\t" << r << setprecision(3) << "\t | \t sav:\t" << savings << green << "\t\t taxes:\t" << l_revenues << def << ", collect " << r*subsidy_disburs << red << ", \tDTI = \t" << DebtToInc << ", \tUnd = \t" << underwater << endl;
        cout << setprecision(6) << red << "Pa  = \t\t" << Pa << setprecision(3) << "\t | \t av_Pm:\t" << average_Pm << "\t | \t av_V:\t" << av_value << yellow << ", Veq = " << av_value_eq << red << endl;
        cout << setprecision(6) << red << "ppsi = \t\t" << ppsi << setprecision(3) << "\t | \t L:\t" << labor_supply << endl;
        cout << setprecision(6) << red << "       \t\t"  << setprecision(3) << "\t\t | \t Realtor:\t" << realtor_profits/people << endl;
        cout << setprecision(6) << red << "tax  = \t\t" << lumpsum << setprecision(3) << endl;
        cout << setprecision(6) << red << "sunk = \t\t" << sunk << setprecision(3) << "\t | \t inter:\t" << interm << "\t | \t ltax:\t" << ltax << endl;
        cout << "----------------------------------------------------------------------------------------------------" << endl;

        // cout << yellow;
        // for(int i = -3; i<7; i++){
        //   cout << wealth_dist[i+3] << "  ";
        // }
        // cout << endl;

        // By equity
        cout << green;
        for(int i = -2; i<5; i++){
          cout << av_value_seq[i+2] << "  ";
        }
        cout << yellow << av_value_nonowners;
        cout << endl;

        for(int i = -2; i<5; i++){
          cout << av_value_seq_eq[i+2] << "  ";
        }
        cout << yellow << av_value_nonowners_eq << def << endl;
        cout << "----------------------------------------------------------------------------------------------------" << endl;
        cout << endl;

        // By wealth
        // cout << green;
        // for(int i = -3; i<7; i++){
        //   cout << av_value_wealth[i+3] << "  ";
        // }
        // cout << endl;

        // for(int i = -3; i<7; i++){
        //   cout << av_value_wealth_eq[i+3] << "  ";
        // }
        // cout << endl;
        // cout << endl;

        // By age
        // cout << green;
        // for(int i = 0; i<15; i++){
        //   cout << av_value_age[i] << "  ";
        // }
        // cout << endl;

        // for(int i = 0; i<15; i++){
        //   cout << av_value_seq_age[i] << "  ";
        // }
        // cout << endl;
      } else{
        cout << setprecision(6) << red << "q = \t\t" << q << setprecision(3) << "\t | \t Rental: \t" << rent_total << "\t" << housing_total << "\tOwned " << "\t | \t Ddelta: \t" << ddeltaoutput << green << "\e[1m Error: " << *derror_accum << "\e[0m" << yellow << "\t\t Subsidy: \t" << subsidy_disburs << "   (error = " << incorrect_subs/correct_subs << ")" << red << endl;
        cout << setprecision(6) << red << "Ph = \t\t" << Ph << setprecision(3) << green << "\t\t taxes:\t" << l_revenues << def << ", collect " << r*subsidy_disburs << red << ", \tUnd = \t" << underwater << green << "\t | \t av_V:\t" << av_value << yellow << ", Veq = " << av_value_eq << blue << "\t | \t av_Pm:\t" << average_Pm << red << endl;
        cout << "----------------------------------------------------------------------------------------------------" << endl;

        // for(int i=0; i<11; i++){
        //   cout << av_value_sigma[i] << "  ";
        // }
        // cout << endl;
        // for(int i=0; i<11; i++){
        //   cout << av_value_sigma_eq[i] << "  ";
        // }
        // cout << endl;

        // cout << yellow;
        // for(int i = -3; i<7; i++){
        //   cout << wealth_dist[i+3] << "  ";
        // }
        // cout << endl;

        // By equity
        // cout << green;
        // for(int i = -2; i<5; i++){
        //   cout << av_value_seq[i+2]/av_pip_seq[i+2] << "  ";
        // }
        // cout << yellow << av_value_nonowners;
        // cout << endl;

        // for(int i = -2; i<5; i++){
        //   cout << av_value_seq_eq[i+2]/av_pip_seq[i+2] << "  ";
        // }
        // cout << yellow << av_value_nonowners_eq << def << endl;
        // cout << "----------------------------------------------------------------------------------------------------" << endl;
        // cout << endl;

        // By wealth
        // cout << green;
        // for(int i = -3; i<7; i++){
        //   cout << av_value_wealth[i+3]/av_pip_wealth[i+3] << "  ";
        // }
        // cout << endl;

        // for(int i = -3; i<7; i++){
        //   cout << av_value_wealth_eq[i+3]/av_pip_wealth[i+3] << "  ";
        // }
        // cout << endl;
        // cout << endl;

        // By age
        cout << green;
        for(int i = 0; i<15; i++){
          cout << av_value_age[i]/av_pip_age[i] << "  ";
        }
        cout << endl;

        for(int i = 0; i<15; i++){
          cout << av_value_seq_age[i]/av_pip_age[i] << "  ";
        }
        cout << endl;
      }


    } else{
      if(summarize == 0){
        cout << setprecision(6) << "q = \t\t" << q << setprecision(3) << "\t | \t Rental: \t" << rent_total << "\t" << housing_total << "\tOwned |\t " << green << "\e[1m Error: " << *derror_accum << "\e[0m" << def << endl;
        cout << setprecision(6) << "ttrans= \t\t" << ittrans << setprecision(3) << "\t | \t Ddelta: \t" << ddeltabar << "\t" << ddeltaoutput << green << "\t       ----------------------" << def << endl;
        cout << setprecision(6) << "Ph = \t\t" << Ph << setprecision(3) << "\t | \t def:\t" << default_rate << "\tRef:\t" << refinance_rate << "\tMort:\t" << mortgagors << "\tM_cre:\t" << mort_creation << endl;
        cout << setprecision(6) << "r = \t\t" << r << setprecision(3) << "\t | \t sav:\t" << savings << "\t\t taxes:\t" << l_revenues << ", \tDTI = \t" << DebtToInc << ", \tUnd = \t" << underwater << endl;
        cout << setprecision(6) << "Pa  = \t\t" << Pa << setprecision(3) << "\t | \t av_Pm:\t" << average_Pm << "\t | \t av_V:\t" << av_value << yellow << ", Veq = " << av_value_eq << red <<  endl;
        cout << setprecision(6) << "ppsi = \t\t" << ppsi << setprecision(3) << "\t | \t L:\t" << labor_supply  << endl;
        cout << setprecision(6) << red << "       \t\t"  << setprecision(3) << "\t\t | \t Realtor:\t" << realtor_profits/people << endl;
        cout << setprecision(6) << "tax  = \t\t" << lumpsum << setprecision(3) << endl;
        cout << setprecision(6) << "sunk = \t\t" << sunk << setprecision(3) << "\t | \t inter:\t" << interm << endl;
        cout << "----------------------------------------------------------------------------------------------------" << endl;
      } else if(summarize == 1 && ittrans == 5 || ittrans == 10){
        cout << setprecision(6) << red << "q = \t\t" << q << setprecision(3) << "\t | \t Rental: \t" << rent_total << "\t" << housing_total << "\tOwned " << "\t | \t Ddelta: \t" << ddeltaoutput << green << "\e[1m Error: " << *derror_accum << "\e[0m" << yellow << "\t\t Subsidy: \t" << subsidy_disburs << "   (error = " << incorrect_subs/correct_subs << ")" << red << endl;
        cout << setprecision(6) << red << "Ph = \t\t" << Ph << setprecision(3) << green << "\t\t taxes:\t" << l_revenues << def << ", collect " << r*subsidy_disburs << red << ", \tUnd = \t" << underwater << green << "\t | \t av_V:\t" << av_value << yellow << ", Veq = " << av_value_eq << red << endl;
        cout << "----------------------------------------------------------------------------------------------------" << endl;

      }

    }
  }


  cout.precision(6);
  cout.setf(std::ios::fixed);


  *d_rental   = rent_total;
  *d_housing  = housing_total;

  return(error_equilib);

}
