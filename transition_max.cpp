//------------------------------------------------------------------------//
//            EQUILIBRIO INICIAL: Q, FCOST, DELTAS, DDELTABAR             //
//------------------------------------------------------------------------//

typedef struct transitions_qs{

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
	double ddeltabar;
	double ddeltaf;
	double r;
	double Ph;
	double q;
	double Pa;
	double housing_supply;
	double fcost;
	double refcost;
	double pension;
	double sstax;
	double ltax;

	double Atech;
	int compute_equivalent;
	
	double multiplier;
	double lumpsum;

	// Transitional dynamics' parameters
	int Ttrans;
	double rshock;
	double Pashock;
	double ppsishock;
	double oomegashock;
	double sunk_shock;
	double interm_shock;
	double ltax_shock;

	double Phinf;
	double qinf;
	double ddeltabarinf;

	int periods_shock;
	int periods_tax;
	int experimento;
	int permanente;

	int baseline;
	std::string tipo;

	double *qpath;
    double *lumpsumpath;
    double *mortsubsidypath;
    int *subs_eligible;
    int *subs_target;

    double prob_mistake;

    double *incshock;

	int   *Piteraciones;
	double *Pmin_upto;
	clock_t *Pt_start;

	double *Prental;
	double *Phousing;

}transitions_qs;




// Equilibrium over rental price q
double transition_eq_qs(unsigned n, const double *x,
                              double *grad,
                              void *params){
	//0. Loading parameters
	struct transitions_qs *p=(struct transitions_qs *)params;
  
	const int maxiter 	= p->maxiter;
	const int uti 		= p->uti;
	const double tol 	= p->tol;
	const double convergence 	= p->convergence;
	const int T 		= p->T;
	const int Tretirement 		= p->Tretirement;
	const int yearspp 	= p->yearspp;

	// Grid for savings: a
	const int na  		= p->na;
	const double amin 	= p->amin;
	const double amax 	= p->amax;

	// Grid for mortgages: m
	const int nm  		= p->nm;
	const double mmin 	= p->mmin;
	const double mmax 	= p->mmax;

	// Grid for housing: h
	const int nh  		= p->nh;
	const double hmin 	= p->hmin;
	const double hmax 	= p->hmax;

	// Grid for renting: r
	const int nr  		= p->nr;
	const double rmin 	= p->rmin;
	const double rmax 	= p->rmax;

	// Grid for labor: l
	const int nl  		= p->nl;
	const double lmax 	= p->lmax;

	// Grid for deoreciation: ddelta
	const int nd  		= p->nd;
	const double dmin 	= p->dmin;
	const double dmax 	= p->dmax;

	// Grid for income shocks: y
	const int ny  			= p->ny;
	const double ssigma_y 	= p->ssigma_y;
	const double llambda_y 	= p->llambda_y;
	const double m_y 		= p->m_y;

	// Preferences
	const double ssigma 	= p->ssigma;
	const double rrho 		= p->rrho;
	const double ppsi 		= p->ppsi;
	const double bbeta 		= p->bbeta;
	const double kkappa 	= p->kkappa;
	const double tthetalab 	= p->tthetalab;
	const double eetalab 	= p->eetalab;

	const double oomega 	= p->oomega;
	const double sunk	 	= p->sunk;
	const double interm	 	= p->interm;
	const double rec_probab = p->rec_probab;

	// Equilibrium objects
	const double ddeltabar 		= p->ddeltabar;
	const double ddeltaf 		= p->ddeltaf;
	const double r 				= p->r;
	const double Ph 			= p->Ph;
	double q 					= p->q;
	const double Pa 			= p->Pa;
	const double housing_supply = p->housing_supply;
	const double fcost 			= p->fcost;
	const double refcost		= p->refcost;
	const double pension		= p->pension;
	const double sstax			= p->sstax;
	const double ltax			= p->ltax;
	const double Atech			= p->Atech;
	const double multiplier		= p->multiplier;
	const double lumpsum 		= p->lumpsum;

	// Transitional dynamics' parameters
	const int Ttrans 			= p->Ttrans;
	const double rshock 		= p->rshock;
	const double Pashock 		= p->Pashock;
	const double ppsishock 		= p->ppsishock;
	const double oomegashock	= p->oomegashock;
	const double sunk_shock 	= p->sunk_shock;
	const double interm_shock 	= p->interm_shock;
	const double ltax_shock 	= p->ltax_shock;

	const int periods_shock  	= p->periods_shock;
	const int periods_tax  		= p->periods_tax;
	const int permanente  		= p->permanente;

	const int baseline  		= p->baseline;
	const std::string tipo  		= p->tipo;

	const double prob_mistake 		= p->prob_mistake;

	const double Phinf  		= p->Phinf;
	const double qinf  			= p->qinf;
	const double ddeltabarinf  	= p->ddeltabarinf;

	const int compute_equivalent = p->compute_equivalent;

	// The function:
	clock_t t;
	t = clock();
	clock_t t2;

	int  *d_iteraciones     = p->Piteraciones; 
	double *d_min_upto 		= p->Pmin_upto;
	clock_t *d_t_start  	= p->Pt_start;

	double *d_rental 		= p->Prental;
	double *d_housing 		= p->Phousing;
	
	parameters paramets = {maxiter, uti, tol, convergence, T, Tretirement, yearspp, na, amin, amax, nm, mmin, mmax, 
							nh, hmin, hmax, nr, rmin, rmax, nl, lmax, nd, dmin, dmax,
							ny, ssigma_y, llambda_y, m_y, ssigma, rrho, ppsi, bbeta, kkappa, tthetalab, eetalab, oomega, sunk, interm, rec_probab,
							ddeltabar, ddeltabar, ddeltaf, r, Ph, Ph, q, Pa, housing_supply, 
							fcost, refcost, pension, sstax, ltax, lumpsum, Atech, compute_equivalent, multiplier, d_rental, d_housing};

	double *qpath, *lumpsumpath, *mortsubsidypath, *incshock;
	int *subs_eligible, *subs_target;
	
	qpath 			= p->qpath;
	lumpsumpath 	= p->lumpsumpath;
	mortsubsidypath = p->mortsubsidypath;
	incshock 		= p->incshock;

	subs_eligible 	= p->subs_eligible;
	subs_target 	= p->subs_target;

	double Ppath[Ttrans];
    double ddeltapath[Ttrans];
    
    for(int i = 0; i<Ttrans; i++){
		Ppath[i] = x[i];
		ddeltapath[i] = x[i+Ttrans];
    }

    double residual;
	residual = transition_error(paramets, 3, Ttrans, rshock, Pashock, ppsishock, oomegashock, sunk_shock, interm_shock, ltax_shock,
								periods_shock, periods_tax, qpath, Ppath, ddeltapath, lumpsumpath, mortsubsidypath, subs_eligible, subs_target, prob_mistake, incshock, baseline, tipo);
  
  	if(residual < *d_min_upto){
  		*d_min_upto = residual;
  	}

	t = clock() - t;
	t2 = clock() - *d_t_start;
	float tiempo = ((float)t2)/CLOCKS_PER_SEC;
	tiempo = tiempo / 60.0;

	std::streamsize ss = std::cout.precision();

	std::cout << "(time: " << ((float)t)/CLOCKS_PER_SEC << " seconds. Total: " << setprecision(2) << tiempo << " minutes. It#: " << setprecision(ss) << (*d_iteraciones) << "    -     Error up to: " << *d_min_upto << ")" << std::endl;
	cout << " " << endl;
	*d_iteraciones = *d_iteraciones + 1;

	return(residual);
}



// Equilibrium over rental price q
double transition_eq_noext(unsigned n, const double *x,
                              double *grad,
                              void *params){
	//0. Loading parameters
	struct transitions_qs *p=(struct transitions_qs *)params;
  
	const int maxiter 	= p->maxiter;
	const int uti 		= p->uti;
	const double tol 	= p->tol;
	const double convergence 	= p->convergence;
	const int T 		= p->T;
	const int Tretirement 		= p->Tretirement;
	const int yearspp 	= p->yearspp;

	// Grid for savings: a
	const int na  		= p->na;
	const double amin 	= p->amin;
	const double amax 	= p->amax;

	// Grid for mortgages: m
	const int nm  		= p->nm;
	const double mmin 	= p->mmin;
	const double mmax 	= p->mmax;

	// Grid for housing: h
	const int nh  		= p->nh;
	const double hmin 	= p->hmin;
	const double hmax 	= p->hmax;

	// Grid for renting: r
	const int nr  		= p->nr;
	const double rmin 	= p->rmin;
	const double rmax 	= p->rmax;

	// Grid for labor: l
	const int nl  		= p->nl;
	const double lmax 	= p->lmax;

	// Grid for deoreciation: ddelta
	const int nd  		= p->nd;
	const double dmin 	= p->dmin;
	const double dmax 	= p->dmax;

	// Grid for income shocks: y
	const int ny  			= p->ny;
	const double ssigma_y 	= p->ssigma_y;
	const double llambda_y 	= p->llambda_y;
	const double m_y 		= p->m_y;

	// Preferences
	const double ssigma 	= p->ssigma;
	const double rrho 		= p->rrho;
	const double ppsi 		= p->ppsi;
	const double bbeta 		= p->bbeta;
	const double kkappa 	= p->kkappa;
	const double tthetalab 	= p->tthetalab;
	const double eetalab 	= p->eetalab;

	const double oomega 	= p->oomega;
	const double sunk	 	= p->sunk;
	const double interm	 	= p->interm;
	const double rec_probab = p->rec_probab;

	// Equilibrium objects
	const double ddeltabar 		= p->ddeltabar;
	const double ddeltaf 		= p->ddeltaf;
	const double r 				= p->r;
	const double Ph 			= p->Ph;
	double q 					= p->q;
	const double Pa 			= p->Pa;
	const double housing_supply = p->housing_supply;
	const double fcost 			= p->fcost;
	const double refcost		= p->refcost;
	const double pension		= p->pension;
	const double sstax			= p->sstax;
	const double ltax			= p->ltax;
	const double Atech			= p->Atech;
	const double multiplier		= p->multiplier;
	const double lumpsum 		= p->lumpsum;

	// Transitional dynamics' parameters
	const int Ttrans 			= p->Ttrans;
	const double rshock 		= p->rshock;
	const double Pashock 		= p->Pashock;
	const double ppsishock 		= p->ppsishock;
	const double oomegashock	= p->oomegashock;
	const double sunk_shock		= p->sunk_shock;
	const double interm_shock	= p->interm_shock;
	const double ltax_shock 	= p->ltax_shock;

	const double periods_shock  = p->periods_shock;
	const int periods_tax  		= p->periods_tax;
	const int permanente  		= p->permanente;

	const int baseline  		= p->baseline;
	const std::string tipo  		= p->tipo;

	const double prob_mistake 		= p->prob_mistake;

	const double Phinf  		= p->Phinf;
	const double qinf  			= p->qinf;
	const double ddeltabarinf  	= p->ddeltabarinf;

	const int compute_equivalent = p->compute_equivalent;


	// The function:
	clock_t t;
	t = clock();
	clock_t t2;

	int  *d_iteraciones     = p->Piteraciones; 
	double *d_min_upto 		= p->Pmin_upto;
	clock_t *d_t_start  	= p->Pt_start;

	double *d_rental 		= p->Prental;
	double *d_housing 		= p->Phousing;

	parameters paramets = {maxiter, uti, tol, convergence, T, Tretirement, yearspp, na, amin, amax, nm, mmin, mmax, 
							nh, hmin, hmax, nr, rmin, rmax, nl, lmax, nd, dmin, dmax,
							ny, ssigma_y, llambda_y, m_y, ssigma, rrho, ppsi, bbeta, kkappa, tthetalab, eetalab, oomega, sunk, interm, rec_probab,
							ddeltabar, ddeltabar, ddeltaf, r, Ph, Ph, q, Pa, housing_supply, 
							fcost, refcost, pension, sstax, ltax, lumpsum, Atech, compute_equivalent, multiplier, d_rental, d_housing};


	double *qpath, *lumpsumpath, *mortsubsidypath, *incshock;
	int *subs_eligible, *subs_target;
	
	qpath 			= p->qpath;
	lumpsumpath 	= p->lumpsumpath;
	mortsubsidypath = p->mortsubsidypath;
	incshock 		= p->incshock;

	subs_eligible 	= p->subs_eligible;
	subs_target 	= p->subs_target;

	double Ppath[Ttrans];
    double ddeltapath[Ttrans];
    
    for(int i = 0; i<Ttrans; i++){
		Ppath[i] = x[i];
		ddeltapath[i] = ddeltabar;

    }

    // for(int it=0; it<T; it++){
    //     cout << incshock[it] << endl;
    //   }



    // cout << ddeltapath[0] << ddeltapath[1] << ddeltapath[2] << ddeltapath[3] << ddeltapath[4] << ddeltapath[5] << endl;
    double residual;
	residual = transition_error(paramets, 4, Ttrans, rshock, Pashock, ppsishock, oomegashock, sunk_shock, interm_shock, ltax_shock,
								periods_shock, periods_tax, qpath, Ppath, ddeltapath, lumpsumpath, mortsubsidypath, subs_eligible, subs_target, prob_mistake, incshock, baseline, tipo);


  	if(residual < *d_min_upto){
  		*d_min_upto = residual;
  	}

	t = clock() - t;
	t2 = clock() - *d_t_start;
	float tiempo = ((float)t2)/CLOCKS_PER_SEC;
	tiempo = tiempo / 60.0;

	std::streamsize ss = std::cout.precision();

	std::cout << "(time: " << ((float)t)/CLOCKS_PER_SEC << " seconds. Total: " << setprecision(2) << tiempo << " minutes. It#: " << setprecision(ss) << (*d_iteraciones) << ", error = " << *d_min_upto << ")" << std::endl;
	cout << " " << endl;
	*d_iteraciones = *d_iteraciones + 1;


	return(residual);
}