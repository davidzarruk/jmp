
//------------------------------------------------------------------------//
//            EQUILIBRIO INICIAL: Q, FCOST, DELTAS, DDELTABAR             //
//------------------------------------------------------------------------//

typedef struct GenEqParameters_eq_24{
	
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
	double Pa;
	double housing_supply;

	double pension;
	double sstax;
	double ltax;

	double fcost;
	double refcost;

	double Atech;
	int compute_equivalent;
	double multiplier;
	double lumpsum;

	int   *Piteraciones;
	double *Pmin_upto;
	clock_t *Pt_start;

	double *d_q_upto;
    double *d_dmin_upto;
    double *d_dmax_upto;
    double *d_ddeltabar_upto;
    double *d_bbeta_upto;
    double *d_fcost_upto;
    double *d_refcost_upto;

	double *Prental;
	double *Phousing;

}pricesolver_eq_24;



double price_zero_eq_24(unsigned n, const double *x,
                              double *grad,
                              void *params){
	//0. Loading parameters
	struct GenEqParameters_eq_24 *p=(struct GenEqParameters_eq_24 *)params;
  
	const int maxiter 	= p->maxiter;
	const int uti 		= p->uti;
	const double tol 	= p->tol;
	const double convergence 	= p->convergence;
	const int T 		= p->T;
	const int Tretirement 		= p->Tretirement;
	const int yearspp	= p->yearspp;
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
	const double lmax   = p->lmax;

	// Grid for deoreciation: ddelta
	const int nd  		= p->nd;

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
	const double ddeltabar	= p->ddeltabar;
	const double ddeltaf 	= p->ddeltaf;
	const double r 			= p->r;
	const double Pa 		= p->Pa;
	const double Ph 		= p->Ph;
	const double housing_supply 	= p->housing_supply;
	const double pension	= p->pension;
	const double sstax 		= p->sstax;
	const double ltax 		= p->ltax;
	const double Atech 		= p->Atech;
	const int compute_equivalent = p->compute_equivalent;
	const double multiplier = p->multiplier;
	const double lumpsum 	= p->lumpsum;

	const double fcost 		= p->fcost;
	const double refcost	= p->refcost;
	
	int ind;	
	size_t sizemats = ny*na*nm*nh*nd*T*sizeof(double);
	double *mortsubsidy;
  	mortsubsidy = (double*)malloc(sizemats);
	
	for(int it=0; it<T; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){
							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

							mortsubsidy[ind] = 0.0;
						}
					}
				}
			}
		}
	}

	const int error_funct 	= 3;

	// The function:
	clock_t t;
	t = clock();
	clock_t t2;

	int  *d_iteraciones     = p->Piteraciones; 
	double *d_min_upto 		= p->Pmin_upto;
	clock_t *d_t_start  	= p->Pt_start;

	double *d_q_upto 			= p->d_q_upto;
	double *d_dmin_upto 		= p->d_dmin_upto;
	double *d_dmax_upto 		= p->d_dmax_upto;
	double *d_ddeltabar_upto 	= p->d_ddeltabar_upto;
	double *d_bbeta_upto 		= p->d_bbeta_upto;
	double *d_fcost_upto 		= p->d_fcost_upto;
	double *d_refcost_upto 		= p->d_refcost_upto;

	double *d_rental 		= p->Prental;
	double *d_housing 		= p->Phousing;

	cout << "Parameters: q = " << x[0] << ", dmin = " << x[1] << ", dmax = " << x[1] + x[2] << ", ddeltabar = " << ddeltabar << ", bbeta = " << bbeta << ", fcost = " << fcost << ", refcost = " << refcost <<  endl;

	parameters paramets = {maxiter, uti, tol, convergence, T, Tretirement, yearspp, na, amin, amax, nm, mmin, mmax, nh, hmin, hmax, nr, rmin, rmax, nl, lmax, nd, x[1], (x[1]+x[2]),
						ny, ssigma_y, llambda_y, m_y, ssigma, rrho, ppsi, bbeta, kkappa, tthetalab, eetalab, oomega, sunk, interm, rec_probab, ddeltabar, ddeltabar, ddeltaf, r, Ph, Ph, x[0], Pa, housing_supply, 
						fcost, refcost, pension, sstax, ltax, lumpsum, Atech, compute_equivalent, multiplier, d_rental, d_housing};

  	double residual = Aggregation_error(paramets, mortsubsidy, error_funct);
  	
  	cout << '\a';
  	
  	if(residual < *d_min_upto){
  		*d_min_upto = residual;
  		*d_q_upto = x[0];
  		*d_ddeltabar_upto = ddeltabar;
  		*d_bbeta_upto = bbeta;
  		*d_fcost_upto = fcost;
  		*d_refcost_upto = refcost;
  		*d_dmax_upto = x[1];
  		*d_dmin_upto = x[2];
  	}

	t = clock() - t;
	t2 = clock() - *d_t_start;
	float tiempo = ((float)t2)/CLOCKS_PER_SEC;
	tiempo = tiempo / 60.0;

	std::streamsize ss = std::cout.precision();

	std::cout << "(time: " << ((float)t)/CLOCKS_PER_SEC << " seconds. Total: " << setprecision(2) << tiempo << " minutes. It#: " << setprecision(ss) << (*d_iteraciones) << "    -     Error up to: " << *d_min_upto << ")" << std::endl;
	std::cout << "(q = " << *d_q_upto << ", dmax = " << *d_dmax_upto << ", dmin = " << *d_dmin_upto << ", bbeta = " << pow(*d_bbeta_upto, 1/(double)yearspp) << ", fcost = " << *d_fcost_upto << ", refcost = " << *d_refcost_upto <<")" << std::endl;
	
	cout << " " << endl;
	*d_iteraciones = *d_iteraciones + 1;

	return(residual);
}



