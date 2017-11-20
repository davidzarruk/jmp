

//--------------------------------------//
//----       In steady state       -----//
//--------------------------------------//

// These following functions return disimbursements and bank assets
// in steady state. Every period, total disimbursements are equal to assets.
// Out of steady state, this need not be true. 


//--------------------------------------//
//--  During Transitional dynamics  ----//
//--------------------------------------//

// Total amount disimbursed by banks to households
double balance_sheet_trans(const parameters params, const double *P, 
							const double *survival, const double *hgrid, const double *rgrid, const double *mgrid,
							const double *agrid, const double *dgrid, const double *ygrid, 
							const double *eprocess, 
							const double *Puncond,
							const double *Value_all, const double *Pricing_guess_all, 
							const int *Default_all, const int *Renew_all, const int *Policya_all, 
							const int *Policym_all, const int *Policyh_all, const int *Policyr_all, const int ittrans, const int Ttrans,
							const double rshock, const double Pashock, const double ppsishock, const double interm_shock, const double sunk_shock, const int periods_shock,
							const double *qpath, const double *Ppath, const double *ddeltapath, 
							const double r0, const double Ph0, const double ddeltabar0, const double interm0, const double sunk0){

	const int T  = params.T;
	const int na = params.na;
	const int nm = params.nm;
	const int nh = params.nh;
	const int nd = params.nd;
	const int ny = params.ny;

	double r 			= params.r;
	const double rrho 	= params.rrho;
	
	double interm = 0.0;
	double Ph;
	double ddeltabar;
	double sunk 	= 0.0;

	int ind;
	int ind2;
	int indp;
	int indp2;
	int indprime;

	double mortgage 	= 0.0;
	double disbursed 	= 0.0;
	double price 		= 0.0;
	double repay_total 	= 0.0;
	double housing 		= 0.0;
	double deprec 		= 0.0;
	double foreclosure 	= 0.0;
	double renewing 	= 0.0;

	double repay_coeff[T];

	double steady_distr[ny];
	steady_distribution(P, params, steady_distr);

	int time;

	// On the host memory
	double *Pricing_guess, *Puncond2, *Puncond_yest;
	int *Default, *Renew, *Policya, *Policym, *Policyh;

	size_t sizemats     = ny*na*nm*nh*nd*T*sizeof(double);
	size_t sizematsint  = ny*na*nm*nh*nd*T*sizeof(int);

	Pricing_guess = (double*)malloc(sizemats);
	Puncond_yest  = (double*)malloc(sizemats);
	Puncond2	  = (double*)malloc(sizemats);

	Default       = (int*)malloc(sizematsint);
	Renew         = (int*)malloc(sizematsint);
	Policya       = (int*)malloc(sizematsint);
	Policym       = (int*)malloc(sizematsint);
	Policyh       = (int*)malloc(sizematsint);

	for(int it=0; it<T; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){
							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
							ind2 = ittrans*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;


							Pricing_guess[ind]     = Pricing_guess_all[ind2];
							Default[ind]           = Default_all[ind2];
							Renew[ind]             = Renew_all[ind2];
							Policya[ind]           = Policya_all[ind2];
							Policym[ind]           = Policym_all[ind2];
							Policyh[ind]           = Policyh_all[ind2];

							Puncond_yest[ind] 	= Puncond[ind];
							Puncond2[ind] 		= 0.0;
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

							// Renew solamente, que son los nuevos mortgages
							mortgage 	= mgrid[Policym[ind]]*Renew[ind];

							indprime = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + Policyh[ind]*nm*na + Policym[ind]*na + Policya[ind];

							price = Pricing_guess[indprime];

							if(Policym[ind]>0){
								disbursed = disbursed + Puncond[ind]*mortgage*price;
							}

						}
					}
				}
			}
		}
	}


	for(int it=0; it<T-1; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){

							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

							if(it == 0 && ih == 0 && im == 0 && ia == 0){
								Puncond2[ind] = 0.0;
							}

							if(Puncond_yest[ind] > 0){
								for(int iyp=0; iyp<ny; iyp++){
									for(int idp=0; idp<nd; idp++){
										for(int ihp=0; ihp<nh; ihp++){
											for(int imp=0; imp<nm; imp++){
												for(int iap=0; iap<na; iap++){
													indp  = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
													indp2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;

													if(Renew[ind] == 1){
														if(imp == Policym[ind] && ihp == Policyh[ind] && iap == Policya[ind]){
															Puncond2[indp] = Puncond2[indp] + rrho*P[iy*ny+iyp]*(1/(double)nd)*Puncond_yest[ind]*survival[it];
															Puncond2[indp2] = Puncond2[indp2] + (1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Puncond_yest[ind]*survival[it];
														}
													} else{
														if(imp == Policym[ind] && ihp == Policyh[ind] && iap == Policya[ind]){
															Puncond2[indp2] = Puncond2[indp2] + P[iy*ny+iyp]*(1/(double)nd)*Puncond_yest[ind]*survival[it];
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

	double discount = 1.0;

	for(int itprime = ittrans+1; itprime<T+ittrans; itprime++){
		
		if(itprime<Ttrans+1){
			time = itprime;
			Ph = Ppath[time-1];
			ddeltabar = ddeltapath[time-1];
		} else{
			time = Ttrans+1;
			Ph = Ph0;
			ddeltabar = ddeltabar0;
		}

		// Shock on interest rate
		if(time <= periods_shock){
			r  = rshock;
			interm = interm_shock;
		} else{
			r  = r0;
			interm = interm0;
		}

		if(time > 0 && time <= periods_shock-1){
			sunk   = sunk_shock;
		} else{
			sunk   = sunk0;
		}


		repay_coefficient_shock(T, r0 + interm0, rshock + interm_shock, periods_shock - itprime, rrho, survival, repay_coeff);

		discount = discount*(1/(1+r+interm));

		// cout << interm << endl;
		// cout << "tiempo: " << ittrans << ", time: " << time << ", Ph = " << Ph << ", ddeltabar = " << ddeltabar << ", r = " << r << ", discount = " << discount << ", ittrans = " << ittrans << endl;

		// cout << "Ph = " << Ph << ", ddeltabar = " << ddeltabar << ", r = " << r << endl;


		// Update of policy functions
		for(int it=0; it<T; it++){
			for(int iy=0; iy<ny; iy++){
				for(int id=0; id<nd; id++){
					for(int ih=0; ih<nh; ih++){
						for(int im=0; im<nm; im++){
							for(int ia=0; ia<na; ia++){
								ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
								ind2 = time*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

								Pricing_guess[ind]     = Pricing_guess_all[ind2];
								Default[ind]           = Default_all[ind2];
								Renew[ind]             = Renew_all[ind2];
								Policya[ind]           = Policya_all[ind2];
								Policym[ind]           = Policym_all[ind2];
								Policyh[ind]           = Policyh_all[ind2];
							}
						}
					}
				}
			}
		}


		// Computo los total repayments
		for(int it=0; it<T; it++){
			for(int iy=0; iy<ny; iy++){
				for(int id=0; id<nd; id++){
					for(int ih=0; ih<nh; ih++){
						for(int im=0; im<nm; im++){
							for(int ia=0; ia<na; ia++){

								ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

								mortgage 	= mgrid[im];

								if(im>0){
									repay_total = repay_total + Puncond2[ind]*discount*(1-Renew[ind])*(1-Default[ind])*mortgage;
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

								housing = hgrid[ih];
								deprec 	= dgrid[id];

								if(im>0){
									foreclosure = foreclosure + Puncond2[ind]*discount*(1-Renew[ind])*Default[ind]*housing*(1-deprec - ddeltabar)*Ph*(1-sunk);
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

								mortgage = mgrid[im];

								if(im>0){
									renewing = renewing + Puncond2[ind]*discount*Renew[ind]*(1-Default[ind])*mortgage*(1+repay_coeff[it]);
								}
							}
						}
					}
				}
			}
		}

		
		// cout << "Forecl = " << foreclosure << ", Repay = " << repay_total << ", Renew = " << renewing << endl;

		// Evolution of distribution
		for(int it=0; it<T; it++){
			for(int iy=0; iy<ny; iy++){
				for(int id=0; id<nd; id++){
					for(int ih=0; ih<nh; ih++){
						for(int im=0; im<nm; im++){
							for(int ia=0; ia<na; ia++){
								ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

								Puncond_yest[ind] = Puncond2[ind];
								Puncond2[ind] = 0.0;
							}
						}
					}
				}
			}
		}
	    
		// Computo las distribuciones unicamente con los mortgages iniciales. No acepto nuevos mortgages 
		for(int it=0; it<T-1; it++){
			for(int iy=0; iy<ny; iy++){
				for(int id=0; id<nd; id++){
					for(int ih=0; ih<nh; ih++){
						for(int im=0; im<nm; im++){
							for(int ia=0; ia<na; ia++){

								ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

								if(it == 0 && ih == 0 && im == 0 && ia == 0){
									Puncond2[ind]      = 0.0;
								}

								if(Puncond_yest[ind] > 0){
									for(int iyp=0; iyp<ny; iyp++){
										for(int idp=0; idp<nd; idp++){
											for(int ihp=0; ihp<nh; ihp++){
												for(int imp=0; imp<nm; imp++){
													for(int iap=0; iap<na; iap++){
														indp  = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + imp*na + iap;
														indp2 = (it+1)*ny*nd*nh*nm*na + iyp*nd*nh*nm*na + idp*nh*nm*na + ihp*nm*na + 0*na + iap;

														if(Default[ind] == 0 && Renew[ind] == 0){
															if(imp == Policym[ind] && ihp == Policyh[ind] && iap == Policya[ind]){
																Puncond2[indp] = Puncond2[indp] + rrho*P[iy*ny+iyp]*(1/(double)nd)*Puncond_yest[ind]*survival[it];
																Puncond2[indp2] = Puncond2[indp2] + (1-rrho)*P[iy*ny+iyp]*(1/(double)nd)*Puncond_yest[ind]*survival[it];
															}
														} else{
															if(imp == Policym[ind] && ihp == Policyh[ind] && iap == Policya[ind]){
																Puncond2[indp2] = Puncond2[indp2] + P[iy*ny+iyp]*(1/(double)nd)*Puncond_yest[ind]*survival[it];
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



	cout << "Disbursed: " << disbursed << ", Revenues = " << repay_total+foreclosure+renewing << ", repays: " << repay_total << ", forecl = " << foreclosure << ", renewing = " << renewing << ", Deficit = " << disbursed - (repay_total+foreclosure+renewing) << endl;
	
	return(repay_total);
}








double pricing_test_all(parameters params, const double *P, 
						const double *survival, const double *hgrid, const double *rgrid, const double *mgrid,
						const double *agrid, const double *dgrid, const double *ygrid, const double *eprocess, const double *repay_coeff, 
						const double *Puncond, const double *Pcond,
						const double *Value_all, const double *Pricing_guess_all, 
						const int *Default_all, const int *Renew_all, const int *Policya_all, 
						const int *Policym_all, const int *Policyh_all, const int *Policyr_all, const int ittrans, const int Ttrans,
						const double *qpath, const double *Ppath, const double *ddeltapath,
						double r, const double Ph, const double ddeltabar,
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

	// const double rrho  	 = params.rrho;
	const double sunk  	 = params.sunk;
	const double interm  = params.interm;

	int ind;
	int ind2;
	int ind3;
	int indprime;

	double mortg 		= 0.0;
	double pprice 		= 0.0;
	double expenditure	= 0.0;

	double pricing_prime 	= 0.0;
	double expected 		= 0.0;

	double steady_distr[ny];
	steady_distribution(P, params, steady_distr);

	// On the host memory
	double *Pricing_guess, *Puncond2, *Pcond2, *Puncond_yest, *Pcond_yest;
	int *Default, *Renew, *Policya, *Policym, *Policyh;

	size_t sizemats     = ny*na*nm*nh*nd*T*sizeof(double);
	size_t sizematsint  = ny*na*nm*nh*nd*T*sizeof(int);

	Pricing_guess = (double*)malloc(sizemats);
	Puncond_yest  = (double*)malloc(sizemats);
	Puncond2	  = (double*)malloc(sizemats);
	Pcond_yest 	  = (double*)malloc(sizemats);
	Pcond2	  	  = (double*)malloc(sizemats);

	Default       = (int*)malloc(sizematsint);
	Renew         = (int*)malloc(sizematsint);
	Policya       = (int*)malloc(sizematsint);
	Policym       = (int*)malloc(sizematsint);
	Policyh       = (int*)malloc(sizematsint);

	for(int it=0; it<T; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){
							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
							ind2 = ittrans*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;


							Pricing_guess[ind]     = Pricing_guess_all[ind2];
							Default[ind]           = Default_all[ind2];
							Renew[ind]             = Renew_all[ind2];
							Policya[ind]           = Policya_all[ind2];
							Policym[ind]           = Policym_all[ind2];
							Policyh[ind]           = Policyh_all[ind2];

							Puncond_yest[ind] 	= 0.0;
							Puncond2[ind] 		= 0.0;
							Pcond_yest[ind] 	= 0.0;
							Pcond2[ind] 		= 0.0;

							if(Policym[ind]>0 && ittrans == 0){
								
								Puncond_yest[ind] 	= Puncond[ind];
								Pcond_yest[ind] 	= Puncond[ind];

							} else if(Policym[ind]>0 && ittrans >0 && Renew[ind] == 1){
								
								Puncond_yest[ind] 	= Puncond[ind];
								Pcond_yest[ind] 	= Puncond[ind];
							
							}
						}
					}
				}
			}
		}
	}

	// int indp;
	// int indp2;

	for(int it=0; it<T; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){
							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

							// Disimbursed
							mortg = mgrid[Policym[ind]];

							indprime = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + Policyh[ind]*nm*na + Policym[ind]*na + Policya[ind];

							pprice = Pricing_guess[indprime];

							if(Policym[ind]>0){
								expenditure = expenditure + Puncond_yest[ind]*mortg*pprice;
							}
						}
					}
				}
			}
		}
	}


	probability_paths_trans(Pcond_yest, Puncond_yest, Pcond2, Puncond2,
                             P, params, survival, Policya, Policym, Policyh, Policya, Policym, Policyh, prob_mistake, mortsubsidy, subs_eligible, subs_target);

	
	// Update of policy functions
	for(int it=0; it<T; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){
							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
							ind2 = (ittrans+1)*T*ny*nd*nh*nm*na + it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

							Pricing_guess[ind]     = Pricing_guess_all[ind2];
							Default[ind]           = Default_all[ind2];
							Renew[ind]             = Renew_all[ind2];
							Policya[ind]           = Policya_all[ind2];
							Policym[ind]           = Policym_all[ind2];
							Policyh[ind]           = Policyh_all[ind2];
						}
					}
				}
			}
		}
	}


	// if(ittrans <= periods_shock-1){
	// 	sunk   = sunk_shock;
	// } else{
	// 	sunk   = sunk0;
	// }

	double kk2 = 0.0;
	double rr2 = 0.0;
	double rene = 0.0;

	for(int it=0; it<T; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){
							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

							pricing_prime = 0.0;
							expected = 0.0;

							ind3 = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + Policyh[ind]*nm*na + Policym[ind]*na + Policya[ind];
							pricing_prime = Pricing_guess[ind3];

							if(im > 0){
								kk2 = kk2 + Puncond2[ind]*(1/(1+r+interm))*(1-Renew[ind])*(1-Default[ind])*(mgrid[im] + mgrid[im]*pricing_prime);
								rr2 = rr2 + Puncond2[ind]*(1/(1+r+interm))*(1-Renew[ind])*Default[ind]*(Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih]*(1-sunk) - Ph*dgrid[id]*hgrid[ih]);
	
								rene = rene + Puncond2[ind]*(1/(1+r+interm))*Renew[ind]*(1-Default[ind])*mgrid[im]*(1+repay_coeff[it]);
							}
						}
					}
				}
			}
		}
	}

	Color::Modifier green(Color::FG_GREEN);
	Color::Modifier def(Color::FG_DEFAULT);

	// cout << "Ph = " << Ph << ", ddeltabar = " << ddeltabar << ", r = " << r << endl;

	cout << "Disbursed: " << expenditure << ", Revenues = " << kk2+rr2+rene << ", repays: = " << kk2 << ", forecl = " << rr2 << ", Renewing = " << rene << ", Deficit = " << green << kk2+rr2+rene-expenditure << ", collect = " << (kk2+rr2+rene-expenditure)*r << def << endl;
    cout << "-------------------------------------------------------------------------------" << green << "---------------------" << def << endl;

	return(expenditure - expected);
}










// -------------------------------------------------------//
// 					Estos son individuales				  //
// -------------------------------------------------------//




// Total amount disimbursed by banks to households

double disbursements(parameters params, const double *P, 
					const double *survival, const double *hgrid, const double *rgrid, const double *mgrid,
	                const double *agrid, const double *dgrid, const double *ygrid, const double *eprocess, const double *repay_coeff, 
	                const double *Pricing_guess, const int *Default, const int *Renew, const int *Policya, 
	                const int *Policym, const int *Policyh, const double *Puncond){

	const int T  = params.T;
	const int na = params.na;
	const int nm = params.nm;
	const int nh = params.nh;
	const int nd = params.nd;
	const int ny = params.ny;

	const double r = params.r;

	int ind;
	int indprime;
	double disbursed 	= 0.0;
	double mortgage 	= 0.0;
	double price 		= 0.0;

	for(int it=0; it<T; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){

							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

							// Renew solamente, que son los nuevos mortgages
							mortgage 	= mgrid[Policym[ind]]*Renew[ind];

							indprime = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + Policyh[ind]*nm*na + Policym[ind]*na + Policya[ind];

							price = Pricing_guess[indprime];

							if(Policym[ind]>0){
								disbursed = disbursed + Puncond[ind]*pow((1/(1+r)), it)*mortgage*price;
							}

						}
					}
				}
			}
		}
	}
	return(disbursed);
}


// Total amount repaid by households to bank

double repayments(parameters params, const double *P, 
					const double *survival, const double *hgrid, const double *rgrid, const double *mgrid,
	                const double *agrid, const double *dgrid, const double *ygrid, const double *eprocess, const double *repay_coeff, 
	                const double *Pricing_guess, const int *Default, const int *Renew, const int *Policya, const double *Puncond){

	const int T  = params.T;
	const int na = params.na;
	const int nm = params.nm;
	const int nh = params.nh;
	const int nd = params.nd;
	const int ny = params.ny;

	const double r = params.r;

	int ind;
	double mortgage 	= 0.0;
	double repay_total 	= 0.0;

	for(int it=0; it<T; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){

							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

							mortgage 	= mgrid[im];

							if(im>0){
								repay_total = repay_total + Puncond[ind]*pow((1/(1+r)), it)*(1-Renew[ind])*(1-Default[ind])*mortgage;
							}
						}
					}
				}
			}
		}
	}

	return(repay_total);
}

// Total amount received by banks generated by foreclosures

double foreclosures(parameters params, const double *P, 
					const double *survival, const double *hgrid, const double *rgrid, const double *mgrid,
	                const double *agrid, const double *dgrid, const double *ygrid, const double *eprocess, const double *repay_coeff, 
	                const double *Pricing_guess, const int *Default, const int *Renew, const int *Policya, const double *Puncond){

	const int T  = params.T;
	const int na = params.na;
	const int nm = params.nm;
	const int nh = params.nh;
	const int nd = params.nd;
	const int ny = params.ny;

	const double r = params.r;
	const double Ph         = params.Ph_tomorrow;
	const double ddeltabar  = params.ddeltabar_tomorrow;
	const double sunk		= params.sunk;

	int ind;
	double housing 		= 0.0;
	double deprec 		= 0.0;
	double foreclosure 	= 0.0;

	for(int it=0; it<T; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){

							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

							housing = hgrid[ih];
							deprec 	= dgrid[id];

							if(im>0){
								foreclosure = foreclosure + Puncond[ind]*pow((1/(1+r)), it)*(1-Renew[ind])*Default[ind]*housing*(1-deprec - ddeltabar)*Ph*(1-sunk);
							}
						}
					}
				}
			}
		}
	}

	return(foreclosure);
}



// Total amount received by banks by refinancing

double renegotiated(parameters params, const double *P, 
					const double *survival, const double *hgrid, const double *rgrid, const double *mgrid,
	                const double *agrid, const double *dgrid, const double *ygrid, const double *eprocess, const double *repay_coeff, 
	                const double *Pricing_guess, const int *Default, const int *Renew, const int *Policya, const double *Puncond){

	const int T  = params.T;
	const int na = params.na;
	const int nm = params.nm;
	const int nh = params.nh;
	const int nd = params.nd;
	const int ny = params.ny;

	const double r = params.r;

	int ind;
	double mortgage		= 0.0;
	double renewing 	= 0.0;

	for(int it=0; it<T; it++){
		for(int iy=0; iy<ny; iy++){
			for(int id=0; id<nd; id++){
				for(int ih=0; ih<nh; ih++){
					for(int im=0; im<nm; im++){
						for(int ia=0; ia<na; ia++){

							ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;

							mortgage = mgrid[im];

							if(im>0){
								renewing = renewing + Puncond[ind]*pow((1/(1+r)), it)*Renew[ind]*(1-Default[ind])*mortgage*(1+repay_coeff[it]);
							}
						}
					}
				}
			}
		}
	}

	return(renewing);
}