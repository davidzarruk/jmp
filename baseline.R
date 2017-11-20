
# Folders and addresses
foldder = "/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/"
CESFold = "/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Data/CEX/fmli/CES_collapsed.dta"
CESlarge = "/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Data/CEX/fmli/CES_nocollapse.dta"
ACSFold = "/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Data/ACS/summary.dta"
SCFFold = "/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Data/SCF/summary.dta"
SCFequity = "/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Data/SCF/Equity distribution/equity_dist.dta"
SCFmortdebt = "/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Data/SCF/Equity distribution/age_equity.dta"
DocumentFold = "Document/"

statistics = as.numeric(readLines("matrices/statisticspath3.txt"));

tol =       statistics[1];
uti =       statistics[2];
maxiter =   statistics[3];
T =         statistics[4];

na =        statistics[5];
amin =      statistics[6];
amax =      statistics[7];

nm =        statistics[8];
mmin =      statistics[9];
mmax =      statistics[10];

nh =        statistics[11];
hmin =      statistics[12];
hmax =      statistics[13];

nr =        statistics[14];
rmin =      statistics[15];
rmax =      statistics[16];

nd =        statistics[17];
dmin =      statistics[18];
dmax =      statistics[19];

ny =        statistics[20];
ssigma_y =  statistics[21];
llambda_y = statistics[22];
m_y =       statistics[23];

ssigma =    statistics[24];
rrho =      statistics[25];
ppsi =      statistics[26];
bbeta =     statistics[27];
kkappa =    statistics[28];
ddeltabar = statistics[29];

eetalab   = statistics[52];
tthetalab = statistics[53];

r =         statistics[30];
Ph =        statistics[31];
q =         statistics[32];
Pa =        statistics[33];
fcost =     statistics[34];
housing_supply = statistics[35];

people            = statistics[44];
owners_w_mortgage = statistics[45];
mortgage_debt     = statistics[46];
housing_value     = statistics[47];
default_rate      = statistics[48];

yearspp           = statistics[50];
Tretirement       = statistics[51];

# Initialize matrices
Value     = array(0, dim = c(T, ny, na, nm, nh, nd))
Policya   = array(0, dim = c(T, ny, na, nm, nh, nd))
Policyh   = array(0, dim = c(T, ny, na, nm, nh, nd))
Policyr   = array(0, dim = c(T, ny, na, nm, nh, nd))
Policym   = array(0, dim = c(T, ny, na, nm, nh, nd))
Policyc   = array(0, dim = c(T, ny, na, nm, nh, nd))
Policyl   = array(0, dim = c(T, ny, na, nm, nh, nd))
Pricing   = array(0, dim = c(T, ny, na, nm, nh, nd))
Default   = array(0, dim = c(T, ny, na, nm, nh, nd))
Renew     = array(0, dim = c(T, ny, na, nm, nh, nd))
Pcond     = array(0, dim = c(T, ny, na, nm, nh, nd))
Puncond   = array(0, dim = c(T, ny, na, nm, nh, nd))

Changer   = array(0, dim = c(T, ny, na, nm, nh))

# Policies
Value_base = array(0, dim = c(T, ny, na, nm, nh, nd))
Value_def  = array(0, dim = c(T, ny, na, nm, nh, nd))
Puncond_base = array(0, dim = c(T, ny, na, nm, nh, nd))
Puncond_def  = array(0, dim = c(T, ny, na, nm, nh, nd))

C       = array(0, dim=c(T))
L       = array(0, dim=c(T))
Income  = array(0, dim=c(T))
H       = array(0, dim=c(T))
Hp      = array(0, dim=c(T))
M       = array(0, dim=c(T))
Mp      = array(0, dim=c(T))
R       = array(0, dim=c(T))
NR      = array(0, dim=c(T))
A       = array(0, dim=c(T))
KK       = array(0, dim=c(T))
Vivos   = array(0, dim=c(T))

D       = array(0, dim=c(T))
D1       = array(0, dim=c(T))
D0       = array(0, dim=c(T))
D9       = array(0, dim=c(T))

Eq       = array(0, dim=c(T))
Eq1       = array(0, dim=c(T))

Da       = array(0, dim=c(na))
Da1       = array(0, dim=c(na))

RR      = array(0, dim=c(T))
Prepay  = array(0, dim=c(T))
PP      = array(0, dim=c(T))
PP1      = array(0, dim=c(T))
PP2      = array(0, dim=c(T))
PP3      = array(0, dim=c(T))
NetRent = array(0, dim=c(T))
OwnerOccupied = array(0, dim=c(T))


# ------------------------------------------ #
#           Steady state policies            #
# ------------------------------------------ #

# Policies and value functions
Val = readLines("matrices/Value3.txt")
Pri = readLines("matrices/Pricing_guess3.txt")
Con = readLines("matrices/Policyc3.txt")
Lab = readLines("matrices/Policyl3.txt")
Sav = readLines("matrices/Policya3.txt")
Mor = readLines("matrices/Policym3.txt")
Ren = readLines("matrices/Policyr3.txt")
Rrenew = readLines("matrices/Renew3.txt")
Ddefault = readLines("matrices/Default3.txt")
Hou = readLines("matrices/Policyh3.txt")
Pco = readLines("matrices/Pcond3.txt")
Pun = readLines("matrices/Puncond3.txt")

# Policies
Val_base = readLines("matrices/Value1base_policy.txt")
Val_def = readLines("matrices/Value1default_policy.txt")
Pun_base = readLines("matrices/Puncond1base_policy.txt")
Pun_def = readLines("matrices/Puncond1default_policy.txt")

#Life-cycle paths: 
C       = as.numeric(readLines("matrices/Cpath3.txt"));
L       = as.numeric(readLines("matrices/Lpath3.txt"));
Income  = as.numeric(readLines("matrices/Incomepath3.txt"));
H       = as.numeric(readLines("matrices/Hpath3.txt"));
M       = as.numeric(readLines("matrices/Mpath3.txt"));
Mp      = as.numeric(readLines("matrices/Mppath3.txt"));
R       = as.numeric(readLines("matrices/Rpath3.txt"));
A       = as.numeric(readLines("matrices/Apath3.txt"));
KK      = as.numeric(readLines("matrices/KKpath3.txt"));

Eq       = as.numeric(readLines("matrices/Eq_tr0.txt"));
Eq1       = as.numeric(readLines("matrices/Eq_tr1.txt"));


D       = as.numeric(readLines("matrices/Defpath3.txt"));
D1       = as.numeric(readLines("matrices/DD_tr1.txt"));
D0       = as.numeric(readLines("matrices/DD_tr0.txt"));
D9       = as.numeric(readLines("matrices/DD_tr9.txt"));
Da       = as.numeric(readLines("matrices/DDa_tr0.txt"));
Da1       = as.numeric(readLines("matrices/DDa_tr1.txt"));

NR      = as.numeric(readLines("matrices/NRpath3.txt"));
RR      = as.numeric(readLines("matrices/Renewpath3.txt"));
RRup      = as.numeric(readLines("matrices/RRuppath3.txt"));
RRdn      = as.numeric(readLines("matrices/RRdnpath3.txt"));
Prepay  = as.numeric(readLines("matrices/Prepaypath3.txt"));

PP      = as.numeric(readLines("matrices/PPpath3.txt"));
PP1      = as.numeric(readLines("matrices/PPpath_tr1.txt"));
PP2      = as.numeric(readLines("matrices/PPpath_tr2.txt"));
PP3      = as.numeric(readLines("matrices/PPpath_tr3.txt"));

OwnerOccupied = as.numeric(readLines("matrices/OwnerOccupiedpath3.txt"));
NetRent = as.numeric(readLines("matrices/NetRentpath3.txt"));
home_equity = as.numeric(readLines("matrices/equitypath3.txt"));
home_equity_pdf = as.numeric(readLines("matrices/equitypath_pdf3.txt"));
default_dist = as.numeric(readLines("matrices/defaultdist3.txt"));

# Distributions over transitions
home_equity1 = as.numeric(readLines("matrices/equitypath_tr1.txt"));
home_equity2 = as.numeric(readLines("matrices/equitypath_tr2.txt"));
home_equity3 = as.numeric(readLines("matrices/equitypath_tr3.txt"));
home_equity4 = as.numeric(readLines("matrices/equitypath_tr4.txt"));
home_equity5 = as.numeric(readLines("matrices/equitypath_tr5.txt"));
home_equity6 = as.numeric(readLines("matrices/equitypath_tr6.txt"));

home_equity7 = as.numeric(readLines("matrices/equitypath_tr7.txt"));
home_equity8 = as.numeric(readLines("matrices/equitypath_tr8.txt"));
home_equity9 = as.numeric(readLines("matrices/equitypath_tr9.txt"));
# 
# home_equity11 = as.numeric(readLines("matrices/equitypath_tr11.txt"));
# home_equity14 = as.numeric(readLines("matrices/equitypath_tr14.txt"));


Chan       = as.numeric(readLines("matrices/Eligibles.txt"));

condDD = as.numeric(readLines("matrices/conDD3.txt"));
condDDdelta = as.numeric(readLines("matrices/conDDdelta3.txt"));

DcondTxAxY = as.numeric(readLines("matrices/DcondTAY3.txt"));
DcondAxY   = as.numeric(readLines("matrices/DcondAY3.txt"));
DcondA     = as.numeric(readLines("matrices/DcondAA3.txt"));
DcondHxM   = as.numeric(readLines("matrices/DcondHM3.txt"));

LcondAxY = as.numeric(readLines("matrices/LcondAY3.txt"));
HcondAxY = as.numeric(readLines("matrices/HcondAY3.txt"));
McondAxY = as.numeric(readLines("matrices/McondAY3.txt"));
PeoplecondAxY = as.numeric(readLines("matrices/PeoplecondAY3.txt"));

equity_cycle = as.numeric(readLines("matrices/equity_path_age3.txt"));
mort_debt = as.numeric(readLines("matrices/Mort_debt_path3.txt"));

DDEFAULT = readLines("matrices/Defaultpath3.txt")

for(it in 1:T){
  for(iy in 1:ny){
    for(ia in 1:na){
      for(im in 1:nm){
        for(ih in 1:nh){
          for(id in 1:nd){
            ind = (it-1)*ny*nd*nh*nm*na + (iy-1)*nd*nh*nm*na + (id-1)*nh*nm*na + (ih-1)*nm*na + (im-1)*na + ia;
            ind2 = (it-1)*ny*nh*nm*na + (iy-1)*nh*nm*na + (ih-1)*nm*na + (im-1)*na + ia;

            Value[it, iy, ia, im, ih, id]   = as.numeric(Val[ind]);
            Pricing[it, iy, ia, im, ih, id] = as.numeric(Pri[ind]);
            Policyc[it, iy, ia, im, ih, id] = as.numeric(Con[ind]);
            Policyl[it, iy, ia, im, ih, id] = as.numeric(Lab[ind]);
            Policya[it, iy, ia, im, ih, id] = as.numeric(Sav[ind]);
            Policym[it, iy, ia, im, ih, id] = as.numeric(Mor[ind]);
            Policyr[it, iy, ia, im, ih, id] = as.numeric(Ren[ind]);
            Policyh[it, iy, ia, im, ih, id] = as.numeric(Hou[ind]);
            Renew[it, iy, ia, im, ih, id] = as.numeric(Rrenew[ind]);
            Default[it, iy, ia, im, ih, id] = as.numeric(Ddefault[ind]);
            Pcond[it, iy, ia, im, ih, id] = as.numeric(Pco[ind]);
            Puncond[it, iy, ia, im, ih, id] = as.numeric(Pun[ind]);
            
            Changer[it, iy, ia, im, ih] = as.numeric(Chan[ind2]);
            
            Value_base[it, iy, ia, im, ih, id]   = as.numeric(Val_base[ind]);
            Value_def[it, iy, ia, im, ih, id]   = as.numeric(Val_def[ind]);
            
            Puncond_base[it, iy, ia, im, ih, id] = as.numeric(Pun_base[ind]);
            Puncond_def[it, iy, ia, im, ih, id] = as.numeric(Pun_def[ind]);
          }
        }
      }
    }
  }
}


eprocess = c(0.710400,0.825600,0.921200,0.958000,0.994800,1.031600,1.068400,1.096800,1.100000,1.103200,1.106400,1.109600,1.109600,
1.100000,1.090400,1.080800,1.071200,1.046400,0.976000,0.905600,0.835200,0.764800,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
0.000000,0.000000)


CCC = 0
YYY = 0

for(it in 1:T){
  for(iy in 1:ny){
    for(ia in 1:na){
      for(im in 1:nm){
        for(ih in 1:nh){
          for(id in 1:nd){
            
            CCC = CCC + Puncond[it, iy, ia, im, ih, id]*Policyc[it, iy, ia, im, ih, id]
          }
        }
      }
    }
  }
}

# 
# lumpsum = 0.036
# nl = 3
# lgrid = seq(0,0.4,length.out = nl)
# it = 2; ia = 2; ih = 2; im = 2; id = 2; iy = 2;
# 
# C2 = 0
# Y2 = 0
# 
# for(it in 1:T){
#   for(iy in 1:ny){
#     for(ia in 1:na){
#       for(im in 1:nm){
#         for(ih in 1:nh){
#           for(id in 1:nd){
#             
#             if(Default[it, iy, ia, im, ih, id] == 0 && Renew[it, iy, ia, im, ih, id] == 0 && it<20){
#               C2 = C2 + Puncond[it, iy, ia, im, ih, id]*Policyc[it, iy, ia, im, ih, id]
#               Y2 = Y2 + Puncond[it, iy, ia, im, ih, id]*(agrid[ia] + q*hgrid[ih] + eprocess[it]*ygrid[iy]*lgrid[Policyl[it, iy, ia, im, ih, id]+1]*(1-0.1) - 
#                         mgrid[im] - q*rgrid[Policyr[it, iy, ia, im, ih, id]+1] - Pa*agrid[Policya[it, iy, ia, im, ih, id]+1] - Ph*(dgrid[id]+ddeltabar)*hgrid[Policyh[it, iy, ia, im, ih, id]+1] + lumpsum)
#             }          
#           }
#         }
#       }
#     }
#   }
# }
# 
# Policyc[it, iy, ia, im, ih, id]
# agrid[ia] + q*hgrid[ih] + eprocess[it]*ygrid[iy]*lgrid[Policyl[it, iy, ia, im, ih, id]+1]*(1-0.1) - 
#   mgrid[im] - q*rgrid[Policyr[it, iy, ia, im, ih, id]+1] - Pa*agrid[Policya[it, iy, ia, im, ih, id]+1] - Ph*(dgrid[id]+ddeltabar)*hgrid[Policyh[it, iy, ia, im, ih, id]+1] + lumpsum


# ------------------------------------------ #
#             Transition policies            #
# ------------------------------------------ #


Ppath_tr = as.numeric(readLines("matrices/values/Ppathtr.txt"));
Ppath_pol = as.numeric(readLines("matrices/values/Ppathbase_policy.txt"));
Ppath_defpol = as.numeric(readLines("matrices/values/Ppathdefault_policy.txt"));

qpath_tr = as.numeric(readLines("matrices/values/qpathtr.txt"));
qpath_pol = as.numeric(readLines("matrices/values/qpathbase_policy.txt"));
qpath_defpol = as.numeric(readLines("matrices/values/qpathdefault_policy.txt"));

Ttrans = 10

defpath_tr = array(0, Ttrans)
defpath_pol = array(0, Ttrans)
defpath_defpol = array(0, Ttrans)
for(i in 1:Ttrans){
  defpath_tr[i] = as.numeric(readLines(paste0("matrices/values/default", i-1, "tr.txt")));
  defpath_pol[i] = as.numeric(readLines(paste0("matrices/values/default", i-1, "base_policy.txt")));
  defpath_defpol[i] = as.numeric(readLines(paste0("matrices/values/default", i-1, "default_policy.txt")));
}

foredataq = data.frame(year = seq(1992,2011.4,by=0.25),
                       foreclosures = c(1.36,1.35,1.33,1.32,1.3,1.3,1.29,1.26,1.25,1.27,1.3,1.32,1.33,1.32,1.31,1.31,1.36,1.37,
                                        1.38,1.38,1.36,1.37,1.38,1.41,1.42,1.44,1.46,1.49,1.48,1.52,1.49,1.46,1.47,1.37,1.47,1.52,
                                        1.55,1.72,1.76,1.79,1.84,1.86,1.83,1.8,1.77,1.63,1.62,1.66,1.72,1.76,1.72,1.73,1.68,1.67,
                                        1.68,1.64,1.63,1.67,1.72,1.84,2.01,2.23,2.55,2.84,3.25,3.79,4.08,4.26,4.61,4.89,5.24,5.37,5.2,4.9,4.8,4.95,4.81,4.65),
                       hpi = c(0.784,0.786,0.778,0.787,0.795,0.807,0.827,0.858,0.869,0.87,0.884,0.896,0.903,0.917,0.945,0.964,0.984,1.007,1.022,1.039,1.065,1.092,1.131,1.181,1.249,1.299,1.353,1.421,1.47,1.528,1.597,1.673,1.765,1.87,1.959,2.054,2.156,2.23,2.308,2.401,2.485,2.591,2.72,2.817,2.919,3.017,3.153,3.317,3.518,3.74,3.943,4.14,4.404,4.654,4.86,5.02,5.169,5.19,5.156,5.184,5.177,5.085,4.941,4.787,4.608,4.426,4.264,4.066,3.823,3.77,3.773,3.745,3.727,3.741,3.681,3.563,3.574,3.606))

datayearly = aggregate(foredataq, by=list(floor(foredataq$year)), FUN=mean, na.rm=TRUE)
datayearly = datayearly[ , !(names(datayearly) %in% c("year"))]
names(datayearly)[names(datayearly) == 'Group.1'] <- 'year'
datayearly$pair = datayearly$year%%2

data_biyear = data.frame(year = array(0,dim(datayearly)[1]/2), 
                         foreclosures = array(0,dim(datayearly)[1]/2), 
                         hpi = array(0,dim(datayearly)[1]/2))

j = 1
for(i in 1:dim(datayearly)[1]){
  data_biyear$foreclosures[j] = data_biyear$foreclosures[j] + datayearly$foreclosures[i]
  data_biyear$hpi[j] = data_biyear$hpi[j] + datayearly$hpi[i]/2
  if(datayearly$pair[i] == 1){
    j = j+1
  } else{
    data_biyear$year[j]         = datayearly$year[i]
  }
}




Value_tr     = array(0, dim = c(T, ny, na, nm, nh, nd))
Policya_tr   = array(0, dim = c(T, ny, na, nm, nh, nd))
Policyh_tr   = array(0, dim = c(T, ny, na, nm, nh, nd))
Policyr_tr   = array(0, dim = c(T, ny, na, nm, nh, nd))
Policym_tr   = array(0, dim = c(T, ny, na, nm, nh, nd))
Policyc_tr   = array(0, dim = c(T, ny, na, nm, nh, nd))
Pricing_tr   = array(0, dim = c(T, ny, na, nm, nh, nd))
Default_tr   = array(0, dim = c(T, ny, na, nm, nh, nd))
Renew_tr     = array(0, dim = c(T, ny, na, nm, nh, nd))
Pcond_tr     = array(0, dim = c(T, ny, na, nm, nh, nd))
Puncond_tr   = array(0, dim = c(T, ny, na, nm, nh, nd))


# Policies and value functions
Val = readLines("matrices/Value1base_policy.txt")
Pri = readLines("matrices/Pricing_guess1base_policy.txt")
Con = readLines("matrices/Policyc1base_policy.txt")
Sav = readLines("matrices/Policya1base_policy.txt")
Mor = readLines("matrices/Policym1base_policy.txt")
Ren = readLines("matrices/Policyr1base_policy.txt")
Rrenew = readLines("matrices/Renew1base_policy.txt")
Ddefault = readLines("matrices/Default1base_policy.txt")
Hou = readLines("matrices/Policyh1base_policy.txt")
Pco = readLines("matrices/Pcond1base_policy.txt")
Pun = readLines("matrices/Puncond1base_policy.txt")


for(it in 1:T){
  for(iy in 1:ny){
    for(ia in 1:na){
      for(im in 1:nm){
        for(ih in 1:nh){
          for(id in 1:nd){
            ind = (it-1)*ny*nd*nh*nm*na + (iy-1)*nd*nh*nm*na + (id-1)*nh*nm*na + (ih-1)*nm*na + (im-1)*na + ia;

            Value_tr[it, iy, ia, im, ih, id]   = as.numeric(Val[ind]);
            Pricing_tr[it, iy, ia, im, ih, id] = as.numeric(Pri[ind]);
            Policyc_tr[it, iy, ia, im, ih, id] = as.numeric(Con[ind]);
            Policya_tr[it, iy, ia, im, ih, id] = as.numeric(Sav[ind]);
            Policym_tr[it, iy, ia, im, ih, id] = as.numeric(Mor[ind]);
            Policyr_tr[it, iy, ia, im, ih, id] = as.numeric(Ren[ind]);
            Policyh_tr[it, iy, ia, im, ih, id] = as.numeric(Hou[ind]);
            Renew_tr[it, iy, ia, im, ih, id] = as.numeric(Rrenew[ind]);
            Default_tr[it, iy, ia, im, ih, id] = as.numeric(Ddefault[ind]);
            Pcond_tr[it, iy, ia, im, ih, id] = as.numeric(Pco[ind]);
            Puncond_tr[it, iy, ia, im, ih, id] = as.numeric(Pun[ind]);

          }
        }
      }
    }
  }
}



defo0 = 0.0
mortg0 = 0.0
price0 = 0.0
def_age0 = array(0, T)
mort_age0 = array(0, T)

defo1 = 0.0
mortg1 = 0.0
price1 = 0.0
def_age1 = array(0, T)
mort_age1 = array(0, T)

for(it in 1:T){
  for(iy in 1:ny){
    for(ia in 1:na){
      for(im in 1:nm){
        for(ih in 1:nh){
          for(id in 1:nd){
            
            price0 = price0 + Puncond[it, iy, ia, im, ih, id]*Pricing[it, iy, ia, im, ih, id]
            price1 = price1 + Puncond_tr[it, iy, ia, im, ih, id]*Pricing_tr[it, iy, ia, im, ih, id]
            
            if(im > 1){
              # Pre-shock
              mortg0 = mortg0 + Puncond[it, iy, ia, im, ih, id]
              defo0 = defo0 + Puncond[it, iy, ia, im, ih, id]*Default[it, iy, ia, im, ih, id];
              
              mort_age0[it] = mort_age0[it] + Puncond[it, iy, ia, im, ih, id]
              def_age0[it] = def_age0[it] + Puncond[it, iy, ia, im, ih, id]*Default[it, iy, ia, im, ih, id];
              
              # Post-shock
              mortg1 = mortg1 + Puncond_tr[it, iy, ia, im, ih, id]
              defo1 = defo1 + Puncond_tr[it, iy, ia, im, ih, id]*Default_tr[it, iy, ia, im, ih, id];
              
              mort_age1[it] = mort_age1[it] + Puncond_tr[it, iy, ia, im, ih, id]
              def_age1[it] = def_age1[it] + Puncond_tr[it, iy, ia, im, ih, id]*Default_tr[it, iy, ia, im, ih, id];
            }
          }
        }
      }
    }
  }
}

defo0 = defo0/mortg0
defo1 = defo1/mortg1





# ------------------------------------------ #
#        Who become more defaulters?         #
# ------------------------------------------ #

y_changers = array(0,ny)
y_holders = array(0,ny)
y_pipol = array(0,ny)

h_changers = array(0,nh)
h_holders = array(0,nh)
h_pipol = array(0,nh)

m_changers = array(0,nm)
m_holders = array(0,nm)
m_pipol = array(0,nm)

t_changers = array(0,T)
t_holders = array(0,T)
t_pipol = array(0,T)

a_changers = array(0,na)
a_holders = array(0,na)
a_pipol = array(0,na)

d_changers = array(0,nd)
d_holders = array(0,nd)
d_pipol = array(0,nd)

probabs       = array(0, dim = c(T, ny, na, nm, nh))
probabs_pipol = array(0, dim = c(T, ny, na, nm, nh))

for(it in 1:T){
  for(iy in 1:ny){
    for(ia in 1:na){
      for(im in 1:nm){
        for(ih in 1:nh){
          for(id in 1:nd){

            if(im > 0 && ih>0){
              if(Changer[it, iy, ia, im, ih] == 1){
                y_changers[iy] = y_changers[iy] + Puncond[it, iy, ia, im, ih, id];
                h_changers[ih] = h_changers[ih] + Puncond[it, iy, ia, im, ih, id];
                m_changers[im] = m_changers[im] + Puncond[it, iy, ia, im, ih, id];
                t_changers[it] = t_changers[it] + Puncond[it, iy, ia, im, ih, id];
                a_changers[ia] = a_changers[ia] + Puncond[it, iy, ia, im, ih, id];
                d_changers[id] = d_changers[id] + Puncond[it, iy, ia, im, ih, id];
              }
              y_pipol[iy] = y_pipol[iy] + Puncond[it, iy, ia, im, ih, id];
              h_pipol[ih] = h_pipol[ih] + Puncond[it, iy, ia, im, ih, id];
              m_pipol[im] = m_pipol[im] + Puncond[it, iy, ia, im, ih, id];
              t_pipol[it] = t_pipol[it] + Puncond[it, iy, ia, im, ih, id];
              a_pipol[ia] = a_pipol[ia] + Puncond[it, iy, ia, im, ih, id];
              d_pipol[id] = d_pipol[id] + Puncond[it, iy, ia, im, ih, id];

              y_holders[iy] = y_holders[iy] + Puncond[it, iy, ia, im, ih, id];
              h_holders[ih] = h_holders[ih] + Puncond[it, iy, ia, im, ih, id];
              m_holders[im] = m_holders[im] + Puncond[it, iy, ia, im, ih, id];
              t_holders[it] = t_holders[it] + Puncond[it, iy, ia, im, ih, id];
              a_holders[ia] = a_holders[ia] + Puncond[it, iy, ia, im, ih, id];
              d_holders[id] = d_holders[id] + Puncond[it, iy, ia, im, ih, id];

              probabs[it, iy, ia, im, ih] = probabs[it, iy, ia, im, ih] + Puncond[it, iy, ia, im, ih, id]
              probabs_pipol[it, iy, ia, im, ih] = probabs_pipol[it, iy, ia, im, ih] + Puncond[it, iy, ia, im, ih, id]
            }
          }
        }
      }
    }
  }
}



# for(it in 1:T){
#   for(iy in 1:ny){
#     for(ia in 1:na){
#       for(im in 1:nm){
#         for(ih in 1:nh){
#           for(id in 1:nd){
# 
#             if(im > 0){
#               y_changers[iy] = y_changers[iy] + Puncond_tr[it, iy, ia, im, ih, id]*(Default_tr[it, iy, ia, im, ih, id]-Default[it, iy, ia, im, ih, id]);
#               y_pipol[iy] = y_pipol[iy] + Puncond_tr[it, iy, ia, im, ih, id];
# 
#               h_changers[ih] = h_changers[ih] + Puncond_tr[it, iy, ia, im, ih, id]*(Default_tr[it, iy, ia, im, ih, id]-Default[it, iy, ia, im, ih, id]);
#               h_pipol[ih] = h_pipol[ih] + Puncond_tr[it, iy, ia, im, ih, id];
# 
#               m_changers[im] = m_changers[im] + Puncond_tr[it, iy, ia, im, ih, id]*(Default_tr[it, iy, ia, im, ih, id]-Default[it, iy, ia, im, ih, id]);
#               m_pipol[im] = m_pipol[im] + Puncond_tr[it, iy, ia, im, ih, id];
# 
#               t_changers[it] = t_changers[it] + Puncond_tr[it, iy, ia, im, ih, id]*(Default_tr[it, iy, ia, im, ih, id]-Default[it, iy, ia, im, ih, id]);
#               t_pipol[it] = t_pipol[it] + Puncond_tr[it, iy, ia, im, ih, id];
# 
#               a_changers[ia] = a_changers[ia] + Puncond_tr[it, iy, ia, im, ih, id]*(Default_tr[it, iy, ia, im, ih, id]-Default[it, iy, ia, im, ih, id]);
#               a_pipol[ia] = a_pipol[ia] + Puncond_tr[it, iy, ia, im, ih, id];
# 
#               d_changers[id] = d_changers[id] + Puncond_tr[it, iy, ia, im, ih, id]*(Default_tr[it, iy, ia, im, ih, id]-Default[it, iy, ia, im, ih, id]);
#               d_pipol[id] = d_pipol[id] + Puncond_tr[it, iy, ia, im, ih, id];
# 
#               y_holders[iy] = y_holders[iy] + Puncond_tr[it, iy, ia, im, ih, id];
#               h_holders[ih] = h_holders[ih] + Puncond_tr[it, iy, ia, im, ih, id];
#               m_holders[im] = m_holders[im] + Puncond_tr[it, iy, ia, im, ih, id];
#               t_holders[it] = t_holders[it] + Puncond_tr[it, iy, ia, im, ih, id];
#               a_holders[ia] = a_holders[ia] + Puncond_tr[it, iy, ia, im, ih, id];
#               d_holders[id] = d_holders[id] + Puncond_tr[it, iy, ia, im, ih, id];
# 
#               probabs[it, iy, ia, im, ih] = probabs[it, iy, ia, im, ih] + Puncond_tr[it, iy, ia, im, ih, id]*(Default_tr[it, iy, ia, im, ih, id]-Default[it, iy, ia, im, ih, id])
#               probabs_pipol[it, iy, ia, im, ih] = probabs_pipol[it, iy, ia, im, ih] + Puncond_tr[it, iy, ia, im, ih, id]
#             }
#           }
#         }
#       }
#     }
#   }
# }

for(iy in 1:ny){
  y_changers[iy] = y_changers[iy]/y_holders[iy]
}

for(ih in 1:nh){
  h_changers[ih] = h_changers[ih]/h_holders[ih]
}

for(im in 1:nm){
  m_changers[im] = m_changers[im]/m_holders[im]
}

for(it in 1:T){
  t_changers[it] = t_changers[it]/t_holders[it]
}

for(ia in 1:na){
  a_changers[ia] = a_changers[ia]/a_holders[ia]
}

for(id in 1:nd){
  d_changers[id] = d_changers[id]/d_holders[id]
}


for(it in 1:T){
  for(iy in 1:ny){
    for(ia in 1:na){
      for(im in 1:nm){
        for(ih in 1:nh){

          if(probabs_pipol[it, iy, ia, im, ih] != 0){
            probabs[it, iy, ia, im, ih] = probabs[it, iy, ia, im, ih]/probabs_pipol[it, iy, ia, im, ih]
          } else{
            probabs[it, iy, ia, im, ih] = 0
          }
        }
      }
    }
  }
}






# 
# # Welfare improvements
# 
# Value_seq     = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# Policyr_seq   = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# Policyc_seq   = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# Policyl_seq   = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# Policya_seq   = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# Policyh_seq   = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# Policym_seq   = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# 
# Default_seq   = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# Renew_seq     = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# Pcond_seq     = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# Puncond_seq   = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# 
# for(ip in 1:11){
#   val = readLines(paste0("matrices/Value", ip-1,"base_policy.txt"));
#   polc = readLines(paste0("matrices/Policyc", ip-1,"base_policy.txt"));
#   poll = readLines(paste0("matrices/Policyl", ip-1,"base_policy.txt"));
#   polr = readLines(paste0("matrices/Policyr", ip-1,"base_policy.txt"));
#   pola = readLines(paste0("matrices/Policya", ip-1,"base_policy.txt"));
#   polh = readLines(paste0("matrices/Policyh", ip-1,"base_policy.txt"));
#   polm = readLines(paste0("matrices/Policym", ip-1,"base_policy.txt"));
#   poldef = readLines(paste0("matrices/Default", ip-1,"base_policy.txt"));
#   polren = readLines(paste0("matrices/Renew", ip-1,"base_policy.txt"));
#   polpuncond = readLines(paste0("matrices/Puncond", ip-1,"base_policy.txt"));
# 
#   for(it in 1:T){
#     for(iy in 1:ny){
#       for(ia in 1:na){
#         for(im in 1:nm){
#           for(ih in 1:nh){
#             for(id in 1:nd){
#               ind = (it-1)*ny*nd*nh*nm*na + (iy-1)*nd*nh*nm*na + (id-1)*nh*nm*na + (ih-1)*nm*na + (im-1)*na + ia;
# 
#               Value_seq[ip, it, iy, ia, im, ih, id]   = as.numeric(val[ind]);
#               Policyc_seq[ip, it, iy, ia, im, ih, id]   = as.numeric(polc[ind]);
#               Policyl_seq[ip, it, iy, ia, im, ih, id]   = as.numeric(poll[ind]);
#               Policyr_seq[ip, it, iy, ia, im, ih, id]   = as.numeric(polr[ind]);
#               Policya_seq[ip, it, iy, ia, im, ih, id]   = as.numeric(pola[ind]);
#               Policyh_seq[ip, it, iy, ia, im, ih, id]   = as.numeric(polh[ind]);
#               Policym_seq[ip, it, iy, ia, im, ih, id]   = as.numeric(polm[ind]);
#               Default_seq[ip, it, iy, ia, im, ih, id]   = as.numeric(poldef[ind]);
#               Renew_seq[ip, it, iy, ia, im, ih, id]   = as.numeric(polren[ind]);
# 
#               Puncond_seq[ip, it, iy, ia, im, ih, id]   = as.numeric(polpuncond[ind]);
#             }
#           }
#         }
#       }
#     }
#   }
# }
# 
# Value_comp     = array(0, dim = c(11, T, ny, na, nm, nh, nd))
# Ptrans = Rtauchen(ny, ssigma_y, llambda_y, 1)
# 
# mV = array(0, dim = c(11))
# 
# for(ip in 2:11){
#   for(it in 1:T){
#     for(iy in 1:ny){
#       for(ia in 1:na){
#         for(im in 1:nm){
#           for(ih in 1:nh){
#             for(id in 1:nd){
# 
#               v_tom = 0;
# 
#               c = Policyc_seq[ip, it, iy, ia, im, ih, id]
#               r = Policyr_seq[ip, it, iy, ia, im, ih, id]
#               l = Policyl_seq[ip, it, iy, ia, im, ih, id]
#               iap = Policya_seq[ip, it, iy, ia, im, ih, id]
#               ihp = Policyh_seq[ip, it, iy, ia, im, ih, id]
#               imp = Policym_seq[ip, it, iy, ia, im, ih, id]
# 
#               if(it<T){
#                 for(iyp in 1:ny){
#                   for(idp in 1:nd){
#                     if(ip == 11){
#                       v_tom = v_tom + (1/nd)*Ptrans[iy,iyp]*(rrho*Value_seq[1, it+1, iyp, iap+1, imp+1, ihp+1, idp] + (1-rrho)*Value_seq[1, it+1, iyp, iap+1, 1, ihp+1, idp])
#                     } else{
#                       v_tom = v_tom + (1/nd)*Ptrans[iy,iyp]*(rrho*Value_comp[ip+1, it+1, iyp, iap+1, imp+1, ihp+1, idp] + (1-rrho)*Value_comp[ip+1, it+1, iyp, iap+1, 1, ihp+1, idp])
#                     }
#                   }
#                 }
#               }
# 
#               Value_comp[ip, it, iy, ia, im, ih, id] = ((ppsi*(c^kkappa) + (1-ppsi)*(r^kkappa))^((1-ssigma)/(kkappa)))/(1-ssigma) - tthetalab*(l^(1+eetalab))/(1+eetalab) + bbeta*surv[it]*v_tom;
# 
#               mV[ip] = mV[ip] + Puncond_seq[ip, it, iy, ia, im, ih, id]*Value_seq[ip, it, iy, ia, im, ih, id]
#             }
#           }
#         }
#       }
#     }
#   }
# }
# 



repay_coeff = c(7.084407,7.041049,6.994780,6.942594,6.883242,6.816373,6.741416,6.657291,6.562933,6.457651,6.340766,6.211675,6.070287,5.916850,5.750890,5.571630,5.377996,5.168729,4.942374,4.696728,4.426358,4.127028,3.796141,3.430416,3.025226,2.573416,2.065566,1.485583,0.808801,0.000000)
dgrid = array(0,nd)
size = nd;
dstep = (dmax - dmin) /(size - 1);
for(i in 1:nd){
  dgrid[i] = dmin + (i-1)*dstep;
}
hgrid = array(0,nh)
size = nh;
hstep = (hmax - hmin) /(size - 1);
for(i in 1:nh){
  hgrid[i] = hmin + (i-1)*hstep;
}
mgrid = array(0,nm)
size = nm;
mstep = (mmax - mmin) /(size - 1);

for(i in 1:nm){
  mgrid[i] = mmin + (i-1)*mstep;
}
# Welfare improvements

Value_bas     = array(0, dim = c(T, ny, na, nm, nh, nd))
Value_sub     = array(0, dim = c(T, ny, na, nm, nh, nd))
Puncond       = array(0, dim = c(T, ny, na, nm, nh, nd))
Punconds       = array(0, dim = c(T, ny, na, nm, nh, nd))


cont = 0
cont2 = 0

for(ip in 2:2){
  val_b = readLines(paste0("matrices/Value", ip-1,"base_policy.txt"));
  val_s = readLines(paste0("matrices/Value", ip-1,"subsidy_only.txt"));
  polpuncond = readLines(paste0("matrices/Puncond", ip-1,"base_policy.txt"));
  polpuncond_s = readLines(paste0("matrices/Puncond", ip-1,"subsidy_only.txt"));
  
  for(it in 1:T){
    for(iy in 1:ny){
      for(ia in 1:na){
        for(im in 1:nm){
          for(ih in 1:nh){
            for(id in 1:nd){
              ind = (it-1)*ny*nd*nh*nm*na + (iy-1)*nd*nh*nm*na + (id-1)*nh*nm*na + (ih-1)*nm*na + (im-1)*na + ia;
              
              if(ih>1){
                home_eq = (0.765*(1-dgrid[id])*hgrid[ih] - mgrid[im]*(1+repay_coeff[it]))/(0.765*(1-dgrid[id])*hgrid[ih])
              }
              Value_bas[it, iy, ia, im, ih, id]   = as.numeric(val_b[ind]);
              Value_sub[it, iy, ia, im, ih, id]   = as.numeric(val_s[ind]);
              
              Puncond[it, iy, ia, im, ih, id]   = as.numeric(polpuncond[ind]);
              Punconds[it, iy, ia, im, ih, id]   = as.numeric(polpuncond_s[ind]);
              
              # if(Punconds[it, iy, ia, im, ih, id] != Puncond[it, iy, ia, im, ih, id]){
              #   print("holalala")
              # }
              # if(Value_bas[it, iy, ia, im, ih, id] != Value_sub[it, iy, ia, im, ih, id]){
              #   print(paste(it, iy, ia, im, ih, id))
              #   # cont = cont + 1
              # }
              # if(home_eq >0 && home_eq < 0.25){
                cont = cont + Puncond[it, iy, ia, im, ih, id]*Value_bas[it, iy, ia, im, ih, id]
                cont2 = cont2 + Punconds[it, iy, ia, im, ih, id]*Value_sub[it, iy, ia, im, ih, id]
              # }
            }
          }
        }
      }
    }
  }
}
