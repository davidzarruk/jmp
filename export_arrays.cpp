
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
           const int* Policyl,
				   const double* Policyc,
				   const double* Pricing_guess,
				   const double* Pcond,
				   const double* Puncond,
				   const double* C, 
           const double* L, 
           const double* Income, 
				   const double* H, 
				   const double* R, 
				   const double* M, 
				   const double* A, 
				   const double* DD, 
				   const double* RR, 
           const double* Prepay, 
           const double* PP, 
           const double* Mp, 
           const double* Hp, 
           const double* KK, 
           const double* RRup, 
           const double* RRdn, 
           const double* NR, 
           const double* NetRent, 
           const double* OwnerOccupied, 
           const double* Def, 
           const double* Vivos, 
           const double* equity, 
           const double* equity_pdf, 
           const double* def_dist,
           const double* Mort_debt, 
           const double* equity_cycle,
           const double* condDD,
           const double* condDDdelta,
           const double* DcondTxAxY,
           const double* DcondAxY,
           const double* DcondA,
           const double* DcondHxM, 
           const double* LcondAxY,
           const double* LcondAxYpeople,
           const double* HcondAxY,
           const double* McondAxY,
           const double* statistics, 
           const int error_funct){
  
  int ind; 
  ostringstream ss;

  ss << "matrices/Value" << error_funct << ".txt";
  ofstream Valuefile (ss.str().c_str());
  if (Valuefile.is_open())
  {
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

  ss.str("");
  ss.clear();
  ss << "matrices/Default" << error_funct << ".txt";
  ofstream Defaultfile (ss.str().c_str());
  if (Defaultfile.is_open())
  {
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

  ss.str("");
  ss.clear();
  ss << "matrices/Renew" << error_funct << ".txt";
  ofstream Renewfile (ss.str().c_str());
  if (Renewfile.is_open())
  {
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

  ss.str("");
  ss.clear();
  ss << "matrices/Policyh" << error_funct << ".txt";
  ofstream Policyhfile (ss.str().c_str());
  if (Policyhfile.is_open())
  {

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

  ss.str("");
  ss.clear();
  ss << "matrices/Policym" << error_funct << ".txt";
  ofstream Policymfile (ss.str().c_str());
  if (Policymfile.is_open())
  {

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

  ss.str("");
  ss.clear();
  ss << "matrices/Policya" << error_funct << ".txt";
  ofstream Policyafile (ss.str().c_str());
  if (Policyafile.is_open())
  {

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

  ss.str("");
  ss.clear();
  ss << "matrices/Policyr" << error_funct << ".txt";
  ofstream Policyrfile (ss.str().c_str());
  if (Policyrfile.is_open())
  {

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


  ss.str("");
  ss.clear();
  ss << "matrices/Policyl" << error_funct << ".txt";
  ofstream Policylfile (ss.str().c_str());
  if (Policylfile.is_open())
  {

    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Policylfile << Policyl[ind] << "\n";

              }
            }
          }
        }
      }
    }
    Policylfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Policyc" << error_funct << ".txt";
	ofstream Policycfile (ss.str().c_str());
  if (Policycfile.is_open())
  {

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

  ss.str("");
  ss.clear();
  ss << "matrices/Pricing_guess" << error_funct << ".txt";
  ofstream Pricing_guessfile (ss.str().c_str());
  if (Pricing_guessfile.is_open())
  {

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

  ss.str("");
  ss.clear();
  ss << "matrices/Pcond" << error_funct << ".txt";
  ofstream Pcondfile (ss.str().c_str());
  if (Pcondfile.is_open())
  {

    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Pcondfile << Pcond[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Pcondfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Puncond" << error_funct << ".txt";
  ofstream Puncondfile (ss.str().c_str());
  if (Puncondfile.is_open())
  {

    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        for(int id=0; id<nd; id++){
          for(int ih=0; ih<nh; ih++){
            for(int im=0; im<nm; im++){
              for(int ia=0; ia<na; ia++){
                ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                Puncondfile << Puncond[ind] << "\n";
              }
            }
          }
        }
      }
    }
    Puncondfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Cpath" << error_funct << ".txt";
  ofstream Cpathfile (ss.str().c_str());
  if (Cpathfile.is_open())
  {
    for(int it=0; it<T; it++){
	    Cpathfile << C[it] << "\n";
    }
    Cpathfile.close();
  }
  else cout << "Unable to open file";

  


  ss.str("");
  ss.clear();
  ss << "matrices/Lpath" << error_funct << ".txt";
  ofstream Lpathfile (ss.str().c_str());
  if (Lpathfile.is_open())
  {
    for(int it=0; it<T; it++){
      Lpathfile << L[it] << "\n";
    }
    Lpathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Incomepath" << error_funct << ".txt";
  ofstream Incomepathfile (ss.str().c_str());
  if (Incomepathfile.is_open())
  {
    for(int it=0; it<T; it++){
      Incomepathfile << Income[it] << "\n";
    }
    Incomepathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Hpath" << error_funct << ".txt";
  ofstream Hpathfile (ss.str().c_str());
  if (Hpathfile.is_open())
  {
    for(int it=0; it<T; it++){
        Hpathfile << H[it] << "\n";
    }
    Hpathfile.close();
  }
  else cout << "Unable to open file";




  ss.str("");
  ss.clear();
  ss << "matrices/Apath" << error_funct << ".txt";
  ofstream Apathfile (ss.str().c_str());
  if (Apathfile.is_open())
  {

    for(int it=0; it<T; it++){
        Apathfile << A[it] << "\n";
    }
    Apathfile.close();
  }
  else cout << "Unable to open file";



  ss.str("");
  ss.clear();
  ss << "matrices/Mpath" << error_funct << ".txt";
  ofstream Mpathfile (ss.str().c_str());
  if (Mpathfile.is_open())
  {

    for(int it=0; it<T; it++){
	    Mpathfile << M[it] << "\n";
    }
    Mpathfile.close();
  }
  else cout << "Unable to open file";



  ss.str("");
  ss.clear();
  ss << "matrices/Rpath" << error_funct << ".txt";
  ofstream Rpathfile (ss.str().c_str());
  if (Rpathfile.is_open())
  {

    for(int it=0; it<T; it++){
	    Rpathfile << R[it] << "\n";
    }
    Rpathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Defaultpath" << error_funct << ".txt";
  ofstream Defaultpathfile (ss.str().c_str());
  if (Defaultpathfile.is_open())
  {

    for(int it=0; it<T; it++){
	    Defaultpathfile << DD[it] << "\n";
    }
    Defaultpathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Renewpath" << error_funct << ".txt";
  ofstream Renewpathfile (ss.str().c_str());
  if (Renewpathfile.is_open())
  {
    for(int it=0; it<T; it++){
	    Renewpathfile << RR[it] << "\n";
    }
    Renewpathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Prepaypath" << error_funct << ".txt";
  ofstream Prepaypathfile (ss.str().c_str());
  if (Prepaypathfile.is_open())
  {
    for(int it=0; it<T; it++){
      Prepaypathfile << Prepay[it] << "\n";
    }
    Prepaypathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/PPpath" << error_funct << ".txt";
  ofstream PPpathfile (ss.str().c_str());
  if (PPpathfile.is_open())
  {

    for(int it=0; it<T; it++){
      PPpathfile << PP[it] << "\n";
    }
    PPpathfile.close();
  }
  else cout << "Unable to open file";



  ss.str("");
  ss.clear();
  ss << "matrices/Mppath" << error_funct << ".txt";
  ofstream Mppathfile (ss.str().c_str());
  if (Mppathfile.is_open())
  {

    for(int it=0; it<T; it++){
      Mppathfile << Mp[it] << "\n";
    }
    Mppathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Hppath" << error_funct << ".txt";
  ofstream Hppathfile (ss.str().c_str());
  if (Hppathfile.is_open())
  {
    for(int it=0; it<T; it++){
      Hppathfile << Hp[it] << "\n";
    }
    Hppathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/KKpath" << error_funct << ".txt";
  ofstream KKpathfile (ss.str().c_str());
  if (KKpathfile.is_open())
  {

    for(int it=0; it<T; it++){
      KKpathfile << KK[it] << "\n";
    }
    KKpathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/RRuppath" << error_funct << ".txt";
  ofstream RRuppathfile (ss.str().c_str());
  if (RRuppathfile.is_open())
  {

    for(int it=0; it<T; it++){
      RRuppathfile << RRup[it] << "\n";
    }
    RRuppathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/RRdnpath" << error_funct << ".txt";
  ofstream RRdnpathfile (ss.str().c_str());
  if (RRdnpathfile.is_open())
  {

    for(int it=0; it<T; it++){
      RRdnpathfile << RRdn[it] << "\n";
    }
    RRdnpathfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/NRpath" << error_funct << ".txt";
  ofstream NRpathfile (ss.str().c_str());
  if (NRpathfile.is_open())
  {

    for(int it=0; it<T; it++){
      NRpathfile << NR[it] << "\n";
    }
    NRpathfile.close();
  }
  else cout << "Unable to open file";



  ss.str("");
  ss.clear();
  ss << "matrices/NetRentpath" << error_funct << ".txt";
  ofstream NetRentpath (ss.str().c_str());
  if (NetRentpath.is_open())
  {

    for(int it=0; it<T; it++){
      NetRentpath << NetRent[it] << "\n";
    }
    NetRentpath.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/OwnerOccupiedpath" << error_funct << ".txt";
  ofstream OwnerOccupiedpath (ss.str().c_str());
  if (OwnerOccupiedpath.is_open())
  {

    for(int it=0; it<T; it++){
      OwnerOccupiedpath << OwnerOccupied[it] << "\n";
    }
    OwnerOccupiedpath.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Defpath" << error_funct << ".txt";
  ofstream Defpath (ss.str().c_str());
  if (Defpath.is_open())
  {

    for(int it=0; it<T; it++){
      Defpath << Def[it] << "\n";
    }
    Defpath.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Vivospath" << error_funct << ".txt";
  ofstream Vivospath (ss.str().c_str());
  if (Vivospath.is_open())
  {

    for(int it=0; it<T; it++){
      Vivospath << Vivos[it] << "\n";
    }
    Vivospath.close();
  }
  else cout << "Unable to open file";

  

  ss.str("");
  ss.clear();
  ss << "matrices/equitypath" << error_funct << ".txt";
  ofstream equitypath (ss.str().c_str());
  if (equitypath.is_open())
  {

    for(int it=-30; it<100; it++){
      equitypath << equity[it+30] << "\n";
    }
    equitypath.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/equitypath_pdf" << error_funct << ".txt";
  ofstream equitypath_pdf (ss.str().c_str());
  if (equitypath_pdf.is_open())
  {

    for(int it=-30; it<100; it++){
      equitypath_pdf << equity_pdf[it+30] << "\n";
    }
    equitypath_pdf.close();
  }
  else cout << "Unable to open file";



  ss.str("");
  ss.clear();
  ss << "matrices/defaultdist" << error_funct << ".txt";
  ofstream defaultdist (ss.str().c_str());
  if (defaultdist.is_open())
  {

    for(int it=-6; it<20; it++){
      defaultdist << def_dist[it+6] << "\n";
    }
    defaultdist.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/Mort_debt_path" << error_funct << ".txt";
  ofstream Mort_debtfile (ss.str().c_str());
  if (Mort_debtfile.is_open())
  {
    for(int it=0; it<T; it++){
        Mort_debtfile << Mort_debt[it] << "\n";
    }
    Mort_debtfile.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/equity_path_age" << error_funct << ".txt";
  ofstream equityfile (ss.str().c_str());
  if (equityfile.is_open())
  {
    for(int it=0; it<T; it++){
        equityfile << equity_cycle[it] << "\n";
    }
    equityfile.close();
  }
  else cout << "Unable to open file";

  int condind;
  ss.str("");
  ss.clear();
  ss << "matrices/conDD" << error_funct << ".txt";
  ofstream conDD (ss.str().c_str());
  if (conDD.is_open())
  {
    for(int it=0; it<T; it++){
      for(int iy=0; iy<ny; iy++){
        condind = it*ny + iy;
        conDD << condDD[condind] << "\n";
      }
    }
    conDD.close();
  }
  else cout << "Unable to open file";

  int condinddelta;
  ss.str("");
  ss.clear();
  ss << "matrices/conDDdelta" << error_funct << ".txt";
  ofstream conDDdelta (ss.str().c_str());
  if (conDDdelta.is_open())
  {
    for(int it=0; it<T; it++){
      for(int id=0; id<nd; id++){
        condinddelta = it*nd + id;

        conDDdelta << condDDdelta[condinddelta] << "\n";
      }
    }
    conDDdelta.close();
  }
  else cout << "Unable to open file";

  int indTxAxY;
  ss.str("");
  ss.clear();
  ss << "matrices/DcondTAY" << error_funct << ".txt";
  ofstream DcondTAY (ss.str().c_str());
  if (DcondTAY.is_open())
  {
    for(int it=0; it<T; it++){
      for(int ia=0; ia<na; ia++){
        for(int iy=0; iy<ny; iy++){
          indTxAxY = it*na*ny + ia*ny + iy;

          DcondTAY << DcondTxAxY[indTxAxY] << "\n";
        }
      }
    }
    DcondTAY.close();
  }
  else cout << "Unable to open file";

  int indAxY;
  ss.str("");
  ss.clear();
  ss << "matrices/DcondAY" << error_funct << ".txt";
  ofstream DcondAY (ss.str().c_str());
  if (DcondAY.is_open())
  {
    for(int ia=0; ia<na; ia++){
      for(int iy=0; iy<ny; iy++){
        indAxY = ia*ny + iy;

        DcondAY << DcondAxY[indAxY] << "\n";
      }
    }
    DcondAY.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/DcondAA" << error_funct << ".txt";
  ofstream DcondAA (ss.str().c_str());
  if (DcondAA.is_open())
  {
    for(int ia=0; ia<na; ia++){
      DcondAA << DcondA[ia] << "\n";
    }
    DcondAA.close();
  }
  else cout << "Unable to open file";


  int indHxM;
  ss.str("");
  ss.clear();
  ss << "matrices/DcondHM" << error_funct << ".txt";
  ofstream DcondHM (ss.str().c_str());
  if (DcondHM.is_open())
  {
    for(int ih=0; ih<nh; ih++){
      for(int im=0; im<nm; im++){
        indHxM = ih*nm + im;

        DcondHM << DcondHxM[indHxM] << "\n";
      }
    }
    DcondHM.close();
  }
  else cout << "Unable to open file";



  ss.str("");
  ss.clear();
  ss << "matrices/LcondAY" << error_funct << ".txt";
  ofstream LcondAY (ss.str().c_str());
  if (LcondAY.is_open())
  {
    for(int ia=0; ia<na; ia++){
      for(int iy=0; iy<ny; iy++){
        indAxY = ia*ny + iy;

        LcondAY << LcondAxY[indAxY] << "\n";
      }
    }
    LcondAY.close();
  }
  else cout << "Unable to open file";



  ss.str("");
  ss.clear();
  ss << "matrices/PeoplecondAY" << error_funct << ".txt";
  ofstream PeoplecondAY (ss.str().c_str());
  if (PeoplecondAY.is_open())
  {
    for(int ia=0; ia<na; ia++){
      for(int iy=0; iy<ny; iy++){
        indAxY = ia*ny + iy;

        PeoplecondAY << LcondAxYpeople[indAxY] << "\n";
      }
    }
    PeoplecondAY.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/HcondAY" << error_funct << ".txt";
  ofstream HcondAY (ss.str().c_str());
  if (HcondAY.is_open())
  {
    for(int ia=0; ia<na; ia++){
      for(int iy=0; iy<ny; iy++){
        indAxY = ia*ny + iy;

        HcondAY << HcondAxY[indAxY] << "\n";
      }
    }
    HcondAY.close();
  }
  else cout << "Unable to open file";


  ss.str("");
  ss.clear();
  ss << "matrices/McondAY" << error_funct << ".txt";
  ofstream McondAY (ss.str().c_str());
  if (McondAY.is_open())
  {
    for(int ia=0; ia<na; ia++){
      for(int iy=0; iy<ny; iy++){
        indAxY = ia*ny + iy;

        McondAY << McondAxY[indAxY] << "\n";
      }
    }
    McondAY.close();
  }
  else cout << "Unable to open file";



  ss.str("");
  ss.clear();
  ss << "matrices/statisticspath" << error_funct << ".txt";
  ofstream statisticspath (ss.str().c_str());
  if (statisticspath.is_open())
  {

    for(int it=0; it<100; it++){
      statisticspath << statistics[it] << "\n";
    }
    statisticspath.close();
  }
  else cout << "Unable to open file";

  // ---------------------------------------------------------------------------------//
  // I generate latex documents with values I want to directly import on my document  //
  // ---------------------------------------------------------------------------------//

  ss.str("");
  ss.clear();
  ss << "matrices/values/negeqvalue.tex";
  ofstream negeqvalue (ss.str().c_str());
  if (negeqvalue.is_open()){

    negeqvalue <<  setprecision(2) << statistics[60]*100;
    negeqvalue.close();
  }
  else cout << "Unable to open file";  


  ss.str("");
  ss.clear();
  ss << "matrices/values/eq10value.tex";
  ofstream eq10value (ss.str().c_str());
  if (eq10value.is_open()){

    eq10value <<  setprecision(3) << statistics[61]*100;
    eq10value.close();
  }
  else cout << "Unable to open file";  


  ss.str("");
  ss.clear();
  ss << "matrices/values/eq20value.tex";
  ofstream eq20value (ss.str().c_str());
  if (eq20value.is_open()){

    eq20value <<  setprecision(3) << statistics[62]*100;
    eq20value.close();
  }
  else cout << "Unable to open file";  

  ss.str("");
  ss.clear();
  ss << "matrices/values/eq25value.tex";
  ofstream eq25value (ss.str().c_str());
  if (eq25value.is_open()){

    eq25value <<  setprecision(3) << statistics[63]*100;
    eq25value.close();
  }
  else cout << "Unable to open file";  

  ss.str("");
  ss.clear();
  ss << "matrices/values/eq30value.tex";
  ofstream eq30value (ss.str().c_str());
  if (eq30value.is_open()){

    eq30value <<  setprecision(3) << statistics[64]*100;
    eq30value.close();
  }
  else cout << "Unable to open file";  


  ss.str("");
  ss.clear();
  ss << "matrices/values/eq50value.tex";
  ofstream eq50value (ss.str().c_str());
  if (eq50value.is_open()){

    eq50value <<  setprecision(3) << statistics[65]*100;
    eq50value.close();
  }
  else cout << "Unable to open file";  


  ss.str("");
  ss.clear();
  ss << "matrices/values/eq75value.tex";
  ofstream eq75value (ss.str().c_str());
  if (eq75value.is_open()){

    eq75value <<  setprecision(3) << statistics[66]*100;
    eq75value.close();
  }
  else cout << "Unable to open file";  

  ss.str("");
  ss.clear();
  ss << "matrices/values/eq100value.tex";
  ofstream eq100value (ss.str().c_str());
  if (eq100value.is_open()){

    eq100value <<  setprecision(3) << statistics[68]*100;
    eq100value.close();
  }
  else cout << "Unable to open file";  


  ss.str("");
  ss.clear();
  ss << "matrices/values/avheqvalue.tex";
  ofstream avheqvalue (ss.str().c_str());
  if (avheqvalue.is_open()){

    avheqvalue <<  setprecision(3) << statistics[69]*100;
    avheqvalue.close();
  }
  else cout << "Unable to open file";  



  // Targeted moments
  ss.str("");
  ss.clear();
  ss << "matrices/values/ownerswmortvalue.tex";
  ofstream ownerswmort (ss.str().c_str());
  if (ownerswmort.is_open()){

    ownerswmort <<  setprecision(3) << statistics[44]*100;
    ownerswmort.close();
  }
  else cout << "Unable to open file";  


  ss.str("");
  ss.clear();
  ss << "matrices/values/defratevalue.tex";
  ofstream defrate (ss.str().c_str());
  if (defrate.is_open()){

    defrate <<  setprecision(2) << statistics[47]*100;
    defrate.close();
  }
  else cout << "Unable to open file";  


  ss.str("");
  ss.clear();
  ss << "matrices/values/ownershipvalue.tex";
  ofstream ownership (ss.str().c_str());
  if (ownership.is_open()){

    ownership <<  setprecision(3) << statistics[80]*100;
    ownership.close();
  }
  else cout << "Unable to open file";  

  ss.str("");
  ss.clear();
  ss << "matrices/values/relsizevalue.tex";
  ofstream relsize (ss.str().c_str());
  if (relsize.is_open()){

    relsize <<  setprecision(3) << statistics[81];
    relsize.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/renew_incomeshock.tex";
  ofstream renincomeshock (ss.str().c_str());
  if (renincomeshock.is_open()){

    renincomeshock <<  setprecision(4) << 100*statistics[90];
    renincomeshock.close();
  }
  else cout << "Unable to open file"; 

  ss.str("");
  ss.clear();
  ss << "matrices/values/renew_priceshock.tex";
  ofstream renpriceshock (ss.str().c_str());
  if (renpriceshock.is_open()){

    renpriceshock <<  setprecision(4) << 100*statistics[91];
    renpriceshock.close();
  }
  else cout << "Unable to open file"; 

  ss.str("");
  ss.clear();
  ss << "matrices/values/default_incomeshock_lowa.tex";
  ofstream defincomeshock_lowa (ss.str().c_str());
  if (defincomeshock_lowa.is_open()){

    defincomeshock_lowa <<  setprecision(4) << 100*statistics[94];
    defincomeshock_lowa.close();
  }
  else cout << "Unable to open file"; 

  ss.str("");
  ss.clear();
  ss << "matrices/values/default_priceshock_lowa.tex";
  ofstream defpriceshock_lowa (ss.str().c_str());
  if (defpriceshock_lowa.is_open()){

    defpriceshock_lowa <<  setprecision(4) << 100*statistics[95];
    defpriceshock_lowa.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/default_lowa.tex";
  ofstream def_lowa (ss.str().c_str());
  if (def_lowa.is_open()){

    def_lowa <<  setprecision(2) << 100*statistics[96];
    def_lowa.close();
  }
  else cout << "Unable to open file"; 

    ss.str("");
  ss.clear();
  ss << "matrices/values/default_incomeshock_higha.tex";
  ofstream defincomeshock_higha (ss.str().c_str());
  if (defincomeshock_higha.is_open()){

    defincomeshock_higha <<  setprecision(4) << 100*statistics[97];
    defincomeshock_higha.close();
  }
  else cout << "Unable to open file"; 

  ss.str("");
  ss.clear();
  ss << "matrices/values/default_priceshock_higha.tex";
  ofstream defpriceshock_higha (ss.str().c_str());
  if (defpriceshock_higha.is_open()){

    defpriceshock_higha <<  setprecision(4) << 100*statistics[98];
    defpriceshock_higha.close();
  }
  else cout << "Unable to open file"; 

  ss.str("");
  ss.clear();
  ss << "matrices/values/default_higha.tex";
  ofstream def_higha (ss.str().c_str());
  if (def_higha.is_open()){

    def_higha <<  setprecision(2) << 100*statistics[99];
    def_higha.close();
  }
  else cout << "Unable to open file"; 


  // Parameters
  ss.str("");
  ss.clear();
  ss << "matrices/values/T.tex";
  ofstream Tval (ss.str().c_str());
  if (Tval.is_open()){

    Tval << statistics[3];
    Tval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/Tret.tex";
  ofstream Tretval (ss.str().c_str());
  if (Tretval.is_open()){

    Tretval << statistics[82];
    Tretval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/yearspp.tex";
  ofstream yearsppval (ss.str().c_str());
  if (yearsppval.is_open()){

    yearsppval <<  setprecision(1) << statistics[83];
    yearsppval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/bbeta.tex";
  ofstream bbetaval (ss.str().c_str());
  if (bbetaval.is_open()){

    bbetaval <<  setprecision(3) << pow(statistics[26], 1/statistics[83]);
    bbetaval.close();
  }
  else cout << "Unable to open file"; 

  ss.str("");
  ss.clear();
  ss << "matrices/values/dmin.tex";
  ofstream dminval (ss.str().c_str());
  if (dminval.is_open()){

    dminval <<  setprecision(3) << statistics[17];
    dminval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/dmax.tex";
  ofstream dmaxval (ss.str().c_str());
  if (dmaxval.is_open()){

    dmaxval <<  setprecision(3) << statistics[18];
    dmaxval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/rshare.tex";
  ofstream rshare (ss.str().c_str());
  if (rshare.is_open()){

    rshare <<  setprecision(3) << 100/(1+statistics[85]);
    rshare.close();
  }
  else cout << "Unable to open file"; 

  
  ss.str("");
  ss.clear();
  ss << "matrices/values/eeta.tex";
  ofstream eetaval (ss.str().c_str());
  if (eetaval.is_open()){

    eetaval << statistics[51];
    eetaval.close();
  }
  else cout << "Unable to open file"; 

  ss.str("");
  ss.clear();
  ss << "matrices/values/tthetal.tex";
  ofstream tthetalval (ss.str().c_str());
  if (tthetalval.is_open()){

    tthetalval << statistics[52];
    tthetalval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/ssigma.tex";
  ofstream ssigmaval (ss.str().c_str());
  if (ssigmaval.is_open()){

    ssigmaval <<  setprecision(2) << statistics[23];
    ssigmaval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/kkappa.tex";
  ofstream kkappaval (ss.str().c_str());
  if (kkappaval.is_open()){

    kkappaval <<  setprecision(2) << statistics[27];
    kkappaval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/ppsi.tex";
  ofstream ppsival (ss.str().c_str());
  if (ppsival.is_open()){

    ppsival <<  setprecision(2) << statistics[25];
    ppsival.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/llambda_y.tex";
  ofstream llambdaval (ss.str().c_str());
  if (llambdaval.is_open()){

    llambdaval <<  setprecision(2) << pow(statistics[21], 1/statistics[83]);
    llambdaval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/ssigma_y.tex";
  ofstream ssigma_yval (ss.str().c_str());
  if (ssigma_yval.is_open()){

    ssigma_yval <<  setprecision(2) << statistics[20];
    ssigma_yval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/m_y.tex";
  ofstream m_yval (ss.str().c_str());
  if (m_yval.is_open()){

    m_yval << setprecision(2) << statistics[22];
    m_yval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/rrho.tex";
  ofstream rrhoval (ss.str().c_str());
  if (rrhoval.is_open()){

    rrhoval <<  setprecision(2) << statistics[24];
    rrhoval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/r.tex";
  ofstream rval (ss.str().c_str());
  if (rval.is_open()){

    rval <<  setprecision(3) << statistics[29];
    rval.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/DTI.tex";
  ofstream DTI (ss.str().c_str());
  if (DTI.is_open()){

    DTI <<  setprecision(3) << statistics[84]*100;
    DTI.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/fcost.tex";
  ofstream fcost (ss.str().c_str());
  if (fcost.is_open()){

    fcost <<  setprecision(3) << statistics[33];
    fcost.close();
  }
  else cout << "Unable to open file"; 

  ss.str("");
  ss.clear();
  ss << "matrices/values/refcost.tex";
  ofstream refcost (ss.str().c_str());
  if (refcost.is_open()){

    refcost <<  setprecision(3) << statistics[53];
    refcost.close();
  }
  else cout << "Unable to open file"; 


  ss.str("");
  ss.clear();
  ss << "matrices/values/sunk.tex";
  ofstream sunk (ss.str().c_str());
  if (sunk.is_open()){

    sunk <<  setprecision(2) << statistics[54];
    sunk.close();
  }
  else cout << "Unable to open file"; 


  // Standard deviation of idiosyncratic house price risk

  ss.str("");
  ss.clear();
  ss << "matrices/values/stdrisk.tex";
  ofstream stdrisk (ss.str().c_str());
  if (stdrisk.is_open()){

    stdrisk <<  setprecision(3) << 100.0*pow((1.0/12.0), 0.5)*(statistics[18]-statistics[17]);
    stdrisk.close();
  }
  else cout << "Unable to open file"; 
}







void export_basic(const int T,
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
                  const int* Policyl,
                  const double* Policyc,
                  const double* Pricing_guess,
                  const int Ttrans, const double rshock, const double Pashock, const double ppsishock, const double oomegashock,
                  const double sunk_shock, const double interm_shock, const double ltax_shock, const int periods_shock,
                  const double *qpath, const double *Ppath,
                  const double* Pcond,
                  const double* Puncond,
                  const double* home_eq,
                  const double* statistics,
                  const double* PP,
                  const double* DD,
                  const double* DDa,
                  const double* Eq,
                  const int ittrans,
                  const double default_rate,
                  const std::string tipo){


  int ind; 
  ostringstream ss;
  cout << setprecision(9); 
  if(ittrans > -1){
    ss << "matrices/Value" << ittrans << tipo << ".txt";
    ofstream Valuefile (ss.str().c_str());
    if (Valuefile.is_open())
    {
      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  Valuefile << setprecision(9) << Value[ind] << "\n";
                }
              }
            }
          }
        }
      }
      Valuefile.close();
    }
    else cout << "Unable to open file";

    ss.str("");
    ss.clear();
    ss << "matrices/Default" << ittrans << tipo << ".txt";
    ofstream Defaultfile (ss.str().c_str());
    if (Defaultfile.is_open())
    {
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

    ss.str("");
    ss.clear();
    ss << "matrices/Renew" << ittrans << tipo << ".txt";
    ofstream Renewfile (ss.str().c_str());
    if (Renewfile.is_open())
    {
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

    ss.str("");
    ss.clear();
    ss << "matrices/DD_tr" << ittrans << ".txt";
    ofstream Defaultpathfile (ss.str().c_str());
    if (Defaultpathfile.is_open())
    {

      for(int it=0; it<T; it++){
        Defaultpathfile << setprecision(9) << DD[it] << "\n";
      }
      Defaultpathfile.close();
    }
    else cout << "Unable to open file";


    ss.str("");
    ss.clear();
    ss << "matrices/DDa_tr" << ittrans << ".txt";
    ofstream Defaultapathfile (ss.str().c_str());
    if (Defaultapathfile.is_open())
    {

      for(int ia=0; ia<na; ia++){
        Defaultapathfile << setprecision(9) << DDa[ia] << "\n";
      }
      Defaultapathfile.close();
    }
    else cout << "Unable to open file";



    ss.str("");
    ss.clear();
    ss << "matrices/Eq_tr" << ittrans << ".txt";
    ofstream Equitypathfile (ss.str().c_str());
    if (Equitypathfile.is_open())
    {

      for(int it=0; it<T; it++){
        Equitypathfile << setprecision(9) << Eq[it] << "\n";
      }
      Equitypathfile.close();
    }
    else cout << "Unable to open file";



    ss.str("");
    ss.clear();
    ss << "matrices/Policyh" << ittrans << tipo << ".txt";
    ofstream Policyhfile (ss.str().c_str());
    if (Policyhfile.is_open())
    {

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

    ss.str("");
    ss.clear();
    ss << "matrices/Policym" << ittrans << tipo << ".txt";
    ofstream Policymfile (ss.str().c_str());
    if (Policymfile.is_open())
    {

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

    ss.str("");
    ss.clear();
    ss << "matrices/Policya" << ittrans << tipo << ".txt";
    ofstream Policyafile (ss.str().c_str());
    if (Policyafile.is_open())
    {

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

    ss.str("");
    ss.clear();
    ss << "matrices/Policyr" << ittrans << tipo << ".txt";
    ofstream Policyrfile (ss.str().c_str());
    if (Policyrfile.is_open())
    {

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



    ss.str("");
    ss.clear();
    ss << "matrices/Policyl" << ittrans << tipo << ".txt";
    ofstream Policylfile (ss.str().c_str());
    if (Policylfile.is_open())
    {

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  Policylfile << Policyl[ind] << "\n";
                }
              }
            }
          }
        }
      }
      Policylfile.close();
    }
    else cout << "Unable to open file";


    

    ss.str("");
    ss.clear();
    ss << "matrices/Policyc" << ittrans << tipo << ".txt";
    ofstream Policycfile (ss.str().c_str());
    if (Policycfile.is_open())
    {

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  Policycfile << setprecision(9) << Policyc[ind] << "\n";
                }
              }
            }
          }
        }
      }
      Policycfile.close();
    }
    else cout << "Unable to open file";

    ss.str("");
    ss.clear();
    ss << "matrices/Pricing_guess" << ittrans << tipo << ".txt";
    ofstream Pricing_guessfile (ss.str().c_str());
    if (Pricing_guessfile.is_open())
    {

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  Pricing_guessfile << setprecision(9) << Pricing_guess[ind] << "\n";
                }
              }
            }
          }
        }
      }
      Pricing_guessfile.close();
    }
    else cout << "Unable to open file";

    ss.str("");
    ss.clear();
    ss << "matrices/Pcond" << ittrans << tipo << ".txt";
    ofstream Pcondfile (ss.str().c_str());
    if (Pcondfile.is_open())
    {

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  Pcondfile << setprecision(9) << Pcond[ind] << "\n";
                }
              }
            }
          }
        }
      }
      Pcondfile.close();
    }
    else cout << "Unable to open file";


    ss.str("");
    ss.clear();
    ss << "matrices/Puncond" << ittrans << tipo << ".txt";
    ofstream Puncondfile (ss.str().c_str());
    if (Puncondfile.is_open())
    {

      for(int it=0; it<T; it++){
        for(int iy=0; iy<ny; iy++){
          for(int id=0; id<nd; id++){
            for(int ih=0; ih<nh; ih++){
              for(int im=0; im<nm; im++){
                for(int ia=0; ia<na; ia++){
                  ind = it*ny*nd*nh*nm*na + iy*nd*nh*nm*na + id*nh*nm*na + ih*nm*na + im*na + ia;
                  Puncondfile << setprecision(9) << Puncond[ind] << "\n";
                }
              }
            }
          }
        }
      }
      Puncondfile.close();
    }
    else cout << "Unable to open file";


    ss.str("");
    ss.clear();
    ss << "matrices/equitypath_tr" << ittrans << ".txt";
    ofstream equitypath (ss.str().c_str());
    if (equitypath.is_open())
    {

      for(int it=-30; it<100; it++){
        equitypath << setprecision(9) << home_eq[it+30] << "\n";
      }
      equitypath.close();
    }
    else cout << "Unable to open file";


    ss.str("");
    ss.clear();
    ss << "matrices/PPpath_tr" << ittrans << ".txt";
    ofstream PPpathfile (ss.str().c_str());
    if (PPpathfile.is_open())
    {

      for(int it=0; it<T; it++){
        PPpathfile << setprecision(9) << PP[it] << "\n";
      }
      PPpathfile.close();
    }
    else cout << "Unable to open file";




    ss.str("");
    ss.clear();
    ss << "matrices/statisticspath" << ittrans << tipo << ".txt";
    ofstream statisticspath (ss.str().c_str());
    if (statisticspath.is_open())
    {

      for(int it=0; it<100; it++){
        statisticspath << statistics[it] << "\n";
      }
      statisticspath.close();
    }
    else cout << "Unable to open file";    


    // Equilibrium objects
    ss.str("");
    ss.clear();
    ss << "matrices/values/Ppath" << tipo << ".txt";
    ofstream Ppathtr (ss.str().c_str());
    if (Ppathtr.is_open())
    {

      for(int it=0; it<Ttrans; it++){
        Ppathtr << setprecision(9) << Ppath[it] << "\n";
      }
      Ppathtr.close();
    }
    else cout << "Unable to open file";  



    ss.str("");
    ss.clear();
    ss << "matrices/values/qpath" << tipo << ".txt";
    ofstream qpathtr (ss.str().c_str());
    if (qpathtr.is_open())
    {

      for(int it=0; it<Ttrans; it++){
        qpathtr << setprecision(9) << qpath[it] << "\n";
      }
      qpathtr.close();
    }
    else cout << "Unable to open file";  

    ss.str("");
    ss.clear();
    ss << "matrices/values/rshock" << tipo << ".txt";
    ofstream rshocktr (ss.str().c_str());
    if (rshocktr.is_open())
    {
      rshocktr << rshock << "\n";

      rshocktr.close();
    }
    else cout << "Unable to open file";  


    ss.str("");
    ss.clear();
    ss << "matrices/values/Pashock" << tipo << ".txt";
    ofstream Pashocktr (ss.str().c_str());
    if (Pashocktr.is_open())
    {
      Pashocktr << Pashock << "\n";

      Pashocktr.close();
    }
    else cout << "Unable to open file";  


    ss.str("");
    ss.clear();
    ss << "matrices/values/ppsishock" << tipo << ".txt";
    ofstream ppsishocktr (ss.str().c_str());
    if (ppsishocktr.is_open())
    {
      ppsishocktr << ppsishock << "\n";

      ppsishocktr.close();
    }
    else cout << "Unable to open file";  


    ss.str("");
    ss.clear();
    ss << "matrices/values/oomegashock" << tipo << ".txt";
    ofstream oomegashocktr (ss.str().c_str());
    if (oomegashocktr.is_open())
    {
      oomegashocktr << oomegashock << "\n";

      oomegashocktr.close();
    }
    else cout << "Unable to open file";  


    ss.str("");
    ss.clear();
    ss << "matrices/values/sunk_shock" << tipo << ".txt";
    ofstream sunk_shocktr (ss.str().c_str());
    if (sunk_shocktr.is_open())
    {
      sunk_shocktr << sunk_shock << "\n";

      sunk_shocktr.close();
    }
    else cout << "Unable to open file"; 


    ss.str("");
    ss.clear();
    ss << "matrices/values/interm_shock" << tipo << ".txt";
    ofstream interm_shocktr (ss.str().c_str());
    if (interm_shocktr.is_open())
    {
      interm_shocktr << interm_shock << "\n";

      interm_shocktr.close();
    }
    else cout << "Unable to open file"; 


    ss.str("");
    ss.clear();
    ss << "matrices/values/ltax_shock" << tipo << ".txt";
    ofstream ltax_shocktr (ss.str().c_str());
    if (ltax_shocktr.is_open())
    {
      ltax_shocktr << ltax_shock << "\n";

      ltax_shocktr.close();
    }
    else cout << "Unable to open file"; 


    ss.str("");
    ss.clear();
    ss << "matrices/values/periodsshock" << tipo << ".txt";
    ofstream periodsshocktr (ss.str().c_str());
    if (periodsshocktr.is_open())
    {
      periodsshocktr << periods_shock << "\n";

      periodsshocktr.close();
    }
    else cout << "Unable to open file"; 


    ss.str("");
    ss.clear();
    ss << "matrices/values/default" << ittrans << tipo << ".txt";
    ofstream deftr (ss.str().c_str());
    if (deftr.is_open())
    {
      deftr << default_rate << "\n";

      deftr.close();
    }
    else cout << "Unable to open file"; 

  }


}





void export_eligible(const int T,
                      const int ny,
                      const int nd,
                      const int na,
                      const int nh,
                      const int nm,
                      const int* Changer){

  int ind; 
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
              
              // if(im>0 && ih>0 && Changer[ind] == 1){
              //   Eligibles << 1 << "\n";
              // } else{
              //   Eligibles << 0 << "\n";
              // }
            }
          }
        }
      }
    }
    Eligibles.close();
  }
  else cout << "Unable to open file";
}
