
# -----------------------------------------------------------#
#                  Prelimiaries                              #
# -----------------------------------------------------------#

.libPaths( c( .libPaths(), "~/R/x86_64-pc-linux-gnu-library/3.4/") )

library("Rcpp")
library("MASS")
library("Rtauchen")
library("ggplot2")
library("foreign")
library("scales")
library("gridExtra")
library("fANCOVA")
library("grid")
library("nleqslv")

rm(list=ls())

setwd("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Model_externalities_19/")

set.seed(90)

# Corro un archivo con los parametros
# Sys.setenv("PKG_CXXFLAGS"="-std=c++11 -fopenmp")
Sys.setenv("PKG_CXXFLAGS"=" -fopenmp")
Sys.setenv("OMP_NUM_THREADS"="8")
Sys.setenv("PKG_LIBS"="-lgsl -lgslcblas -lnlopt -fopenmp")


sourceCpp("Rcpp_model_lifecycle.cpp")


# -----------------------------------------------------------#
#                  Parameters                                #
# -----------------------------------------------------------#

source('baseline.R')



hgrid = array(0,nh)
size = nh;
hstep = (hmax - hmin) /(size - 1);
for(i in 1:nh){
  hgrid[i] = hmin + (i-1)*hstep;
}

rgrid = array(0,nr)
size = nr;
rstep = (rmax - rmin) /(size - 1);

for(i in 1:nr){
  rgrid[i] = rmin + (i-1)*rstep;
}

RRR = array(0, dim=c(T))
DDD = array(0, dim=c(T))
KKK = array(0, dim=c(T))
RRent = array(0, dim=c(T))

for(it in 1:T){
  for(iy in 1:ny){
    for(id in 1:nd){
      for(ih in 1:nh){
        for(im in 1:nm){
          for(ia in 1:na){
            
            RRent[it] = RRent[it] + Pcond[it, iy, ia, im, ih, id]*rgrid[Policyr[it, iy, ia, im, ih, id]+1]
            if(im>1){
              RRR[it] = RRR[it] + Pcond[it, iy, ia, im, ih, id]*Renew[it, iy, ia, im, ih, id];
              DDD[it] = DDD[it] + Pcond[it, iy, ia, im, ih, id]*Default[it, iy, ia, im, ih, id];
              KKK[it] = KKK[it] + Pcond[it, iy, ia, im, ih, id]*(1-Renew[it, iy, ia, im, ih, id])*(1-Default[it, iy, ia, im, ih, id]);
            }
          }
        }
      }
    }
  }
}






dist_m = array(0,nm)
dist_a = array(0,na)
dist_h = array(0,nh)
dist_r = array(0,nr)
cross_m = array(0,ny)
people_y = array(0,ny)

mgrid = array(0,nm)
size = nm;
mstep = (mmax - mmin) /(size - 1);

for(i in 1:nm){
  mgrid[i] = mmin + (i-1)*mstep;
}

hgrid = seq(hmin, hmax, length.out=nh)
for(it in 1:T){
  for(iy in 1:ny){
    for(id in 1:nd){
      for(ih in 1:nh){
        for(im in 1:nm){
          for(ia in 1:na){
            dist_m[im] = dist_m[im] + Pcond[it, iy, ia, im, ih, id]/12;
            dist_a[ia] = dist_a[ia] + Pcond[it, iy, ia, im, ih, id]/12;
            
            cross_m[iy] = cross_m[iy] + Puncond[it, iy, ia, im, ih, id]*mgrid[Policym[it, iy, ia, im, ih, id]+1];
            people_y[iy] = people_y[iy] + Puncond[it, iy, ia, im, ih, id];
            
            for(ihp in 1:nh){
              if((Policyh[it, iy, ia, im, ih, id]+1) == ihp){
                dist_h[ihp] = dist_h[ihp] + Puncond[it, iy, ia, im, ih, id];
              }
            }
            for(ir in 1:nr){
              if((Policyr[it, iy, ia, im, ih, id]+1) == ir){
                dist_r[ir] = dist_r[ir] + Puncond[it, iy, ia, im, ih, id];
              }
            }
          }
        }
      }
    }
  }
}

for(iy in 1:ny){
  cross_m[iy] = cross_m[iy]/people_y[iy]
}


# ----------------#
#      values     #
# ----------------#

value_base= 0.0;
value_def = 0.0
for(it in 1:T){
  for(iy in 1:ny){
    for(id in 1:nd){
      for(ih in 1:nh){
        for(im in 1:nm){
          for(ia in 1:na){
            
            value_base = value_base + Puncond_base[it, iy, ia, im, ih, id]*Value_base[it, iy, ia, im, ih, id]
            value_def = value_def + Puncond_def[it, iy, ia, im, ih, id]*Value_def[it, iy, ia, im, ih, id]
          }
        }
      }
    }
  }
}


# ----------------#
#      Data       #
# ----------------#


data3 = read.dta(CESFold)
data3 = subset(data3, age >=22 & age <=80)

# Non durable consumption
a = loess.ancova(x = data3$age, y = log(data3$c_nondurablesp), group = data3$year, degree = 1, criterion = "aicc",
                 family = "gaussian", method="Speckman",
                 iter = 100, tol = 0.01, plot = FALSE)

mat = data.frame(age = data3$age, year = data3$year, fitted = exp(a$linear.fit[1] + a$linear.fit[length(a$linear.fit)] + a$smooth.fit$fitted))
emp_cons = subset(mat$fitted, mat$year == 2014)

# Rellenar con NA las edades que no coinciden
emp_cons = c(NA, NA, emp_cons)
for(i in 62:70){
  emp_cons[i] = NA
}
C_emp = array(0,T)
for(i in 1:T){
  co = 0
  for(j in (1+(i-1)*yearspp):(i*yearspp)){
    if(is.na(emp_cons[j]) == FALSE){
      C_emp[i] = C_emp[i] + emp_cons[j]
      co = co+1
    }
  }
  C_emp[i] = C_emp[i] / co
}


# Labor
a = loess.ancova(x = data3$age, y = data3$labor_unc, group = data3$year, degree = 1, criterion = "aicc",
                 family = "gaussian", method="Speckman",
                 iter = 100, tol = 0.01, plot = FALSE)

mat = data.frame(age = data3$age, year = data3$year, fitted = a$linear.fit[1] + a$linear.fit[length(a$linear.fit)] + a$smooth.fit$fitted)
emp_labor = subset(mat$fitted, mat$year == 2014)

# Rellenar con NA las edades que no coinciden
emp_labor = c(NA, NA, emp_labor)
for(i in 62:70){
  emp_labor[i] = NA
}
L_emp = array(0,T)
for(i in 1:T){
  co = 0
  for(j in (1+(i-1)*yearspp):(i*yearspp)){
    if(is.na(emp_labor[j]) == FALSE){
      L_emp[i] = L_emp[i] + emp_labor[j]
      co = co+1
    }
  }
  L_emp[i] = L_emp[i] / co
}

# Mortgages lifecycle - CEX - NEW!!!
data333 = data3[data3$year<2008,]
a = loess.ancova(x = data333$age, y = log(data333$c_orgmrtx_unc), group = data333$year, degree = 1, criterion = "aicc",
                 family = "gaussian", method="Speckman",
                 iter = 100, tol = 0.01, plot = FALSE)

mat = data.frame(age = data333$age, year = data333$year, fitted = exp(a$linear.fit[1] + a$linear.fit[length(a$linear.fit)] + a$smooth.fit$fitted))
emp_mortgage = subset(mat$fitted, mat$year == 2007)

# Rellenar con NA las edades que no coinciden
emp_mortgage = c(NA, NA, emp_mortgage)
for(i in 62:70){
  emp_mortgage[i] = NA
}
Mortgage_emp = array(0,T)
for(i in 1:T){
  co = 0
  for(j in (1+(i-1)*yearspp):(i*yearspp)){
    if(is.na(emp_mortgage[j]) == FALSE){
      Mortgage_emp[i] = Mortgage_emp[i] + emp_mortgage[j]
      co = co+1
    }
  }
  Mortgage_emp[i] = Mortgage_emp[i] / co
}



# Housing size
a = loess.ancova(x = data3$age, y = data3$c_rentequivalent, group = data3$year, degree = 1, criterion = "aicc",
                 family = "gaussian", method="Speckman",
                 iter = 100, tol = 0.01, plot = FALSE)

mat = data.frame(age = data3$age, year = data3$year, fitted = a$linear.fit[1] + a$linear.fit[length(a$linear.fit)] + a$smooth.fit$fitted)

housing = subset(mat$fitted, mat$year == 2014)

housing = c(NA, NA, housing)
for(i in 62:70){
  housing[i] = NA
}
REqv_emp = array(0,T)
for(i in 1:T){
  co = 0
  for(j in (1+(i-1)*yearspp):(i*yearspp)){
    if(is.na(housing[j]) == FALSE){
      REqv_emp[i] = REqv_emp[i] + housing[j]
      co = co+1
    }
  }
  REqv_emp[i] = REqv_emp[i] / co
}

# House value - unconstrained ACS
data = read.dta(ACSFold)
a = loess.ancova(x = data$agep, y = data$c_value_unc, group = data$year, degree = 1, criterion = "aicc",
                 family = "gaussian", method="Speckman",
                 iter = 100, tol = 0.01, plot = FALSE)

mat = data.frame(age = data$agep, year = data$year, fitted = a$linear.fit[1] + a$linear.fit[length(a$linear.fit)] + a$smooth.fit$fitted)

owner_value = subset(mat$fitted, mat$year == 2014)
owner_value = owner_value[2:length(owner_value)]

H_emp_ACS = array(0,T)
for(i in 1:T){
  co = 0
  for(j in (1+(i-1)*yearspp):(i*yearspp)){
    if(is.na(owner_value[j]) == FALSE){
      H_emp_ACS[i] = H_emp_ACS[i] + owner_value[j]
      co = co+1
    }
  }
  H_emp_ACS[i] = H_emp_ACS[i] / co
}

# Rent - unconstrained ACS
a = loess.ancova(x = data$agep, y = data$c_rent_unc, group = data$year, degree = 1, criterion = "aicc",
                 family = "gaussian", method="Speckman",
                 iter = 100, tol = 0.01, plot = FALSE)

mat = data.frame(age = data$agep, year = data$year, fitted = a$linear.fit[1] + a$linear.fit[length(a$linear.fit)] + a$smooth.fit$fitted)

rent_unc = subset(mat$fitted, mat$year == 2014)
rent_unc = rent_unc[2:length(rent_unc)]

R_emp = array(0,T)
for(i in 1:T){
  co = 0
  for(j in (1+(i-1)*yearspp):(i*yearspp)){
    if(is.na(rent_unc[j]) == FALSE){
      R_emp[i] = R_emp[i] + rent_unc[j]
      co = co+1
    }
  }
  R_emp[i] = R_emp[i] / co
}

# House value - unconstrained SCF
data2 = read.dta(SCFFold)
data2 = subset(data2, age > 23)

a = loess.ancova(x = data2$age, y = data2$value_unc, group = data2$year, degree = 1, criterion = "aicc",
                 family = "gaussian", method="Speckman",
                 iter = 100, tol = 0.01, plot = FALSE)

mat = data.frame(age = data2$age, year = data2$year, fitted = a$linear.fit[1] + a$linear.fit[length(a$linear.fit)] + a$smooth.fit$fitted)
mat = mat[order(mat$year,mat$age),]

owner_value_SCF = subset(mat$fitted, mat$year == 2013)
owner_value_SCF = c(NA, NA, NA, owner_value_SCF)

H_emp_SCF = array(0,T)
for(i in 1:T){
  co = 0
  for(j in (1+(i-1)*yearspp):(i*yearspp)){
    if(is.na(owner_value_SCF[j]) == FALSE){
      H_emp_SCF[i] = H_emp_SCF[i] + owner_value_SCF[j]
      co = co+1
    }
  }
  H_emp_SCF[i] = H_emp_SCF[i] / co
}



# Mortgage value - first mortgage, SCF
data2 = read.dta(SCFFold)
data2 = subset(data2, (age > 23 & age < 81))

a = loess.ancova(x = data2$age, y = data2$mort_val1, group = data2$year, degree = 1, criterion = "aicc",
                 family = "gaussian", method="Speckman",
                 iter = 100, tol = 0.01, plot = FALSE)

mat = data.frame(age = data2$age, year = data2$year, fitted = a$linear.fit[1] + a$linear.fit[length(a$linear.fit)] + a$smooth.fit$fitted)
mat = mat[order(mat$year,mat$age),]

mortgage_SCF = subset(mat$fitted, mat$year == 2013)
mortgage_SCF = c(NA, NA, NA, mortgage_SCF)

M_SCF = array(0,T)
for(i in 1:T){
  co = 0
  for(j in (1+(i-1)*yearspp):(i*yearspp)){
    if(is.na(mortgage_SCF[j]) == FALSE){
      M_SCF[i] = M_SCF[i] + mortgage_SCF[j]
      co = co+1
    }
  }
  M_SCF[i] = M_SCF[i] / co
}


# % of renters
a = loess.ancova(x = data2$age, y = data2$renter, group = data2$year, degree = 1, criterion = "aicc",
                 family = "gaussian", method="Speckman",
                 iter = 100, tol = 0.01, plot = FALSE)

mat = data.frame(age = data2$age, year = data2$year, fitted = a$linear.fit[1] + a$linear.fit[length(a$linear.fit)] + a$smooth.fit$fitted)
mat = mat[order(mat$year,mat$age),]

emp_renters = subset(mat$fitted, mat$year == 2013)

emp_renters = c(NA, NA, NA, emp_renters)

Rp_emp = array(0,T)
for(i in 1:T){
  co = 0
  for(j in (1+(i-1)*yearspp):(i*yearspp)){
    if(is.na(emp_renters[j]) == FALSE){
      Rp_emp[i] = Rp_emp[i] + emp_renters[j]
      co = co+1
    }
  }
  Rp_emp[i] = Rp_emp[i] / co
}



# Home equity distribution - SCF 2007 (as in Chatterjee and Eyigungor)

data_homeequity = read.dta(SCFequity)
dist_equity = array(0,100)
dist_pdf = array(0,130)

initial = 1

for(j in 1:100){
  for(i in initial:dim(data_homeequity)[1]){
    if(is.na(data_homeequity$homeeq[i])==FALSE & data_homeequity$homeeq[i] <= j/100){
      dist_equity[j] = data_homeequity$cdf[i]
    } else{
      initial = i+1
      break
    }
  }
}

# data_homeequity$bin = NA
# for(j in -30:100){
#   for(i in 1:dim(data_homeequity)[1]){
#     if(i == 10 && j == -17){
#       print(data_homeequity$bin[i])
#     }
#     if(is.na(data_homeequity$homeeq[i])==FALSE & data_homeequity$homeeq[i] > (j-1)/100 & data_homeequity$homeeq[i] <= j/100){
#       data_homeequity$bin[i] = j
#       dist_pdf[j+31] = dist_pdf[j+31] + 1
#     }
#   }
# }

# Mortgage debt on life cycle - SCF 2007 (as in Chatterjee and Eyigungor)

data_debt = read.dta(SCFmortdebt)
data_debt = subset(data_debt, (age > 23 & age < 81))

# debt_cycle = array(0,T)
# debt_cycle = data_debt$debt_perc
debt_cycle = array(0,11)

eq_cycle = array(0,11)
# eq_cycle = data_debt$homeeq

for(j in 3:57){
  eq_cycle[ceiling((data_debt$age[j]-25)/5)] = eq_cycle[ceiling((data_debt$age[j]-25)/5)] + data_debt$homeeq[j]/5
  debt_cycle[ceiling((data_debt$age[j]-25)/5)] = debt_cycle[ceiling((data_debt$age[j]-25)/5)] + data_debt$debt_perc[j]/5
}

# ------------------#
#      Graphs       #
# ------------------#

size_factor = 1
size_line = 1.5
size_line_doc = 1
theme_plot = theme(panel.background = element_rect(fill='white', colour='black'),
                   panel.grid.major = element_line(colour = "lightgray"),
                   plot.title = element_text(hjust = 0.5))

age_vector = age = seq(from = 20, to = 79.9, by = yearspp)

aggregates = data.frame(age = age_vector, 
                        C = C, C_emp=C_emp*C[1]/C_emp[1], 
                        L = L, L_emp=L_emp/90, 
                        Income = Income,
                        RU = NetRent, RUemp = (R_emp*NetRent[1]/R_emp[1]),
                        R = q*R, Reqv = REqv_emp*C[1]/C_emp[1],
                        VU = (OwnerOccupied), VUemp = (H_emp_ACS*OwnerOccupied[10]/H_emp_ACS[10]),
                        NR = NR, Rp_emp = Rp_emp, 
                        A = A)

data_mortgages = read.dta(CESlarge)
data_mortgages$orgmrtx_unc[is.na(data_mortgages$orgmrtx_unc)==TRUE] = 0
data_mortgages$orgmrtx_unc = mmax*data_mortgages$orgmrtx_unc/max(data_mortgages$orgmrtx_unc[is.na(data_mortgages$orgmrtx_unc)==FALSE])

pMort = ggplot() + 
  geom_histogram(data = data_mortgages, aes(x=orgmrtx_unc, y=..count../sum(..count..)), binwidth = 0.2) + theme_plot + 
  geom_line(data = data.frame(dist=dist_m, mgrid=seq(mmin,mmax,length.out =nm)), aes(x=mgrid, y = dist)) + 
  ylab("Frequency") + xlab("Mortgage value ($)") + ggtitle("Mortgage Value ($)") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pC = ggplot(data = aggregates, aes(x=age)) + 
  geom_line(aes(y = R, linetype="Model"), colour="blue", size = size_line) + 
  geom_line(aes(y = Reqv, linetype="Data"), colour="blue", size=size_line)  + 
  geom_line(aes(y = C, linetype="Model"), colour="black", size = size_line) + 
  geom_line(aes(y = C_emp, linetype="Data"), colour="black", size=size_line) +
  ggtitle("Consumption (black) + Rental (blue)") + 
  scale_linetype_manual("", breaks = c("Model", "Data"), values = c("Model"="solid", "Data"="dashed"))+
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("$") +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(pC) 

pCratios = ggplot(data = aggregates, aes(x=age)) + 
  geom_line(aes(y = C/R, linetype="Model"), colour="blue", size = size_line) + 
  geom_line(aes(y = C_emp/Reqv, linetype="Data"), colour="blue", size=size_line)  + 
  ggtitle("Consumption (black) + Rental (blue)") + 
  scale_linetype_manual("", breaks = c("Model", "Data"), values = c("Model"="solid", "Data"="dashed"))+
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))


pR = ggplot(data = aggregates, aes(x=age)) + 
  geom_line(aes(y = NR), size=size_line) + 
  geom_line(aes(y = Rp_emp), colour="black", linetype="dashed", size=size_line) + 
  ggtitle("Net Renters by Age (%)") + ylab("%") + xlab("Age") + 
  scale_y_continuous(limit = c(0, 1.1), expand = c(0, 0))  + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pL = ggplot(data = aggregates, aes(x=age)) + 
  geom_line(aes(y = L), size=size_line) + 
  geom_line(aes(y = L_emp), colour="black", linetype="dashed", size=size_line) + 
  ggtitle("Labor") + ylab("") + xlab("Age") + 
  scale_y_continuous(limit = c(0, 0.5), expand = c(0, 0))  + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pIncome = ggplot(data = aggregates, aes(x=age)) + 
  geom_line(aes(y = Income), size=size_line) + 
  ggtitle("Income") + ylab("$") + xlab("Age") + 
  scale_y_continuous(limit = c(0, 1.2), expand = c(0, 0))  + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pA = ggplot(data = data.frame(A=A, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = A), size=size_line) + 
  ggtitle("Assets")   + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot

# pM = ggplot(data = data.frame(M=M, Mp=Mp, age=age_vector), aes(x=age)) + 
#   geom_line(aes(y = M), size=size_line) + 
#   geom_line(aes(y = M_SCF*M[4]/M_SCF[4]), colour = "black", size=size_line, linetype="dashed") + 
#   scale_y_continuous(expand = c(0, 0)) + 
#   scale_x_continuous(expand = c(0, 0)) + 
#   ggtitle("Mortgage")  + theme_plot

pM = ggplot(data = data.frame(M=M, Mp=Mortgage_emp*M[4]/Mortgage_emp[4], age=age_vector), aes(x=age)) + 
  geom_line(aes(y = M), size=size_line) + 
  geom_line(aes(y = Mp), colour="black", linetype="dashed", size=size_line) + 
  scale_y_continuous(expand = c(0, 0), limit = c(0, 1.2*max(Mp))) + 
  scale_x_continuous(expand = c(0, 0)) + 
  ggtitle("Mortgage")  + theme_plot

pH = ggplot(data = aggregates, aes(x=age)) + 
  geom_line(aes(y = VU), size=size_line) + 
  geom_line(aes(y = VUemp), colour="black", linetype="dashed", size=size_line) + 
  ggtitle("Owner Occupied Housing") + 
  scale_y_continuous(limit = c(0, 1.5*max((H_emp_ACS*OwnerOccupied[5]/H_emp_ACS[5]))), expand = c(0, 0))  + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot


conditionalDefault = array(0, dim=c(T,ny))
for(it in 1:T){
  for(iy in 1:ny){
    ind = (it-1)*ny + iy;
    conditionalDefault[it, iy] = condDD[ind];
  }
}

conditionalDefaultDelta = array(0, dim=c(T,nd))
for(it in 1:T){
  for(id in 1:nd){
    ind = (it-1)*nd + id;
    conditionalDefaultDelta[it, id] = condDDdelta[ind];
  }
}

cDefTxAxY = array(0, dim=c(T,na,ny))
for(it in 1:T){
  for(ia in 1:na){
    for(iy in 1:ny){
      ind = (it-1)*na*ny + (ia-1)*ny + iy;
      cDefTxAxY[it, ia, iy] = DcondTxAxY[ind];
    }
  }
}

cDefAxY = array(0, dim=c(na,ny))
for(ia in 1:na){
  for(iy in 1:ny){
    ind = (ia-1)*ny + iy;
    cDefAxY[ia, iy] = DcondAxY[ind];
  }
}

cPipolAxY = array(0, dim=c(na,ny))
for(ia in 1:na){
  for(iy in 1:ny){
    ind = (ia-1)*ny + iy;
    cPipolAxY[ia, iy] = PeoplecondAxY[ind];
  }
}

cLevAxY = array(0, dim=c(na,ny))
for(ia in 1:na){
  for(iy in 1:ny){
    ind = (ia-1)*ny + iy;
    cLevAxY[ia, iy] = LcondAxY[ind];
  }
}

cHevAxY = array(0, dim=c(na,ny))
for(ia in 1:na){
  for(iy in 1:ny){
    ind = (ia-1)*ny + iy;
    cHevAxY[ia, iy] = HcondAxY[ind];
  }
}

cMevAxY = array(0, dim=c(na,ny))
for(ia in 1:na){
  for(iy in 1:ny){
    ind = (ia-1)*ny + iy;
    cMevAxY[ia, iy] = McondAxY[ind];
  }
}

cDefHxM = array(0, dim=c(nh,nm))
for(ih in 1:nh){
  for(im in 1:nm){
    ind = (ih-1)*nm + im;
    cDefHxM[ih, im] = DcondHxM[ind];
  }
}


pD = ggplot(data = data.frame(D=D, age=age_vector), aes(x=age)) +
  geom_line(aes(y = D), size=size_line) +
  ggtitle("Default (%)") +
  scale_y_continuous(limit = c(0, 0.4), expand = c(0, 0))  +  ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot
# 
# pDcond = ggplot(data = data.frame(D1=cDefTxAxY[,3,1], 
#                                   D2=cDefTxAxY[,3,2], 
#                                   D3=cDefTxAxY[,3,3], 
#                                   D4=cDefTxAxY[,3,4], 
#                                   D5=cDefTxAxY[,3,5], 
#                                   age=age_vector), aes(x=age)) + 
#   geom_line(aes(y = D1), size=size_line, colour="red") + 
#   geom_line(aes(y = D3), size=size_line, colour="yellow") + 
#   geom_line(aes(y = D5), size=size_line, colour="blue") + 
#   ggtitle("Default (%)") + 
#   scale_y_continuous(limit = c(0, 0.4), expand = c(0, 0))  +  ylab("%") +
#   scale_x_continuous(expand = c(0, 0)) + theme_plot
# 
# # Default conditional on income for every a
# ggplot(data = data.frame(D1=cDefAxY[,1], 
#                          D2=cDefAxY[,2], 
#                          D3=cDefAxY[,3], 
#                          D4=cDefAxY[,4], 
#                          D5=cDefAxY[,5], 
#                          assets=grida(na,amin,amax)), aes(x=assets)) + 
#   geom_line(aes(y = D1), size=size_line, colour="red") + 
#   geom_line(aes(y = D2), size=size_line, colour="orange") + 
#   geom_line(aes(y = D3), size=size_line, colour="yellow") + 
#   geom_line(aes(y = D4), size=size_line, colour="green") + 
#   geom_line(aes(y = D5), size=size_line, colour="blue") + 
#   ggtitle("Default (%)") + 
#   scale_y_continuous(limit = c(0, 0.4), expand = c(0, 0))  +  ylab("%") +
#   scale_x_continuous(expand = c(0, 0)) + theme_plot


# Default conditional on a for every y
ygrid = exp(Tgrid(ny, ssigma_y, llambda_y, m_y))
# ggplot(data = data.frame(D1=cDefAxY[1,], 
#                          D2=cDefAxY[7,], 
#                          D3=cDefAxY[18,], 
#                          D4=cDefAxY[26,], 
#                          D5=cDefAxY[35,], 
#                          income=ygrid), aes(x=income)) + 
#   geom_line(aes(y = D1), size=size_line, colour="red") + 
#   geom_line(aes(y = D2), size=size_line, colour="orange") + 
#   geom_line(aes(y = D3), size=size_line, colour="yellow") + 
#   geom_line(aes(y = D4), size=size_line, colour="green") + 
#   geom_line(aes(y = D5), size=size_line, colour="blue") + 
#   ggtitle("Default (%)") + 
#   scale_y_continuous(limit = c(0, 0.4), expand = c(0, 0))  +  ylab("%") +
#   scale_x_continuous(expand = c(0, 0)) + theme_plot


pRR = ggplot(data = data.frame(RR=RR, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = RR), size=size_line) + ggtitle("Refinance (%)")   + 
  scale_y_continuous(limit = c(0, 0.5), expand = c(0, 0)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pPrepay = ggplot(data = data.frame(Prepay=Prepay, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = Prepay), size=size_line) + ggtitle("Prepay (%)")   + 
  scale_y_continuous(limit = c(0, 0.5), expand = c(0, 0)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pPP = ggplot(data = data.frame(PP=PP, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = PP), size=size_line) + ggtitle("Pricing")   + 
  scale_y_continuous(expand = c(0, 0)) + ylab("$") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pimpliedR = ggplot(data = data.frame(r=1/(PP), age=age_vector), aes(x=age)) + 
  geom_line(aes(y = r), size=size_line) + ggtitle("Implied Interest Rate")   + 
  scale_y_continuous(limit = c(0.0, 5), expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pimpliedR = ggplot(data = data.frame(PP=1/PP-1, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = PP), size=size_line) + ggtitle("Pricing")   + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pdista = ggplot(data = data.frame(dist=dist_a, agrid=seq(amin,amax,length.out =na)), aes(x=agrid)) + 
  geom_line(aes(y = dist), size=size_line) + ggtitle("Dist: A") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pdistm = ggplot(data = data.frame(dist=dist_m, mgrid=seq(mmin,mmax,length.out =nm)), aes(x=mgrid)) + 
  geom_line(aes(y = dist), size=size_line) + ggtitle("Dist: M") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pdisth = ggplot(data = data.frame(dist=dist_h, hgrid=seq(hmin,hmax,length.out =nh)), aes(x=hgrid)) + 
  geom_line(aes(y = dist), size=size_line) + ggtitle("Dist: H (black), R (blue)") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot +
  geom_line(data = data.frame(dist=dist_r, rgrid=seq(rmin,rmax,length.out =nr)), aes(x=rgrid, y=dist), size=size_line, colour = "blue")

pdistr = ggplot(data = data.frame(dist=dist_r, rgrid=seq(rmin,rmax,length.out =nr)), aes(x=rgrid)) + 
  geom_line(aes(y = dist), size=size_line) + ggtitle("Dist: R") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pequity = ggplot(data = data.frame(dist=dist_equity, grid=seq(0,1,length.out =100)), aes(x=grid)) + 
  geom_line(aes(y = dist), size=size_line, linetype = "dashed") + ggtitle("Home Equity CDF") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity") + 
  geom_line(data = data.frame(dist=home_equity, grid=seq(-0.3,1,length.out =130)), aes(x=grid, y=dist), size=size_line, linetype = "solid") + 
  geom_line(data = data.frame(dist=default_dist, grid=seq(-0.3,1,length.out =130)), aes(x=grid, y=dist), size=size_line, linetype = "solid")





# pdebt_cycle = ggplot(data = data.frame(dist=debt_cycle, grid=seq(24,80,length.out =81-24)), aes(x=grid)) + 
#   geom_line(aes(y = dist), size=size_line, linetype = "dashed") + ggtitle("Mortgage debt / House value") + 
#   scale_x_continuous(expand = c(0, 0)) + theme_plot +
#   geom_line(data = data.frame(dist=mort_debt, grid=seq(24,80,length.out =T)), aes(x=grid, y=dist), size=size_line, linetype = "solid")

pdebt_cycle = ggplot(data = data.frame(dist=debt_cycle, grid=seq(24,80,length.out =11)), aes(x=grid)) + 
  geom_line(aes(y = dist), size=size_line, linetype = "dashed") + ggtitle("Mortgage debt / House value") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot +
  geom_line(data = data.frame(dist=mort_debt, grid=seq(24,80,length.out =T)), aes(x=grid, y=dist), size=size_line, linetype = "solid")

# peq_cycle = ggplot(data = data.frame(dist=eq_cycle, grid=seq(24,80,length.out =81-24)), aes(x=grid)) + 
#   geom_line(aes(y = dist), size=size_line, linetype = "dashed") + ggtitle("Home equity") + 
#   scale_x_continuous(expand = c(0, 0)) + theme_plot +
#   geom_line(data = data.frame(dist=equity_cycle, grid=seq(24,80,length.out =T)), aes(x=grid, y=dist), size=size_line, linetype = "solid")

peq_cycle = ggplot(data = data.frame(dist=eq_cycle, grid=seq(24,80,length.out =11)), aes(x=grid)) + 
  geom_line(aes(y = dist), size=size_line, linetype = "dashed") + ggtitle("% Home Equity by Age") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Age") + 
  geom_line(data = data.frame(dist=equity_cycle, grid=seq(24,80,length.out =T)), aes(x=grid, y=dist), size=size_line, linetype = "solid")

grid.arrange(arrangeGrob(peq_cycle, ncol = 1, nrow = 1),
             arrangeGrob(pdebt_cycle, ncol = 1, nrow = 1),
             ncol = 1, nrow = 2, heights=c(1,1))

# Table 2
size_factor = .7
tt  = ttheme_default(colhead = list(fg_params=list(parse=TRUE, cex = size_factor)), 
                     rowhead = list(fg_params=list(cex = size_factor)),
                     core = list(fg_params=list(cex = size_factor)))

aggreg_stats = data.frame(age=c("% owners w. mort.", "Mortgages / House value", "Default rate"), 
                          data = c("0.66", "0.42", "1.52%"),
                          mod = c(round(owners_w_mortgage,2), round(mortgage_debt/housing_value,2), paste0(round(default_rate*100,2),"%")))

names(aggreg_stats) = c("Variable", "Data", "Model")

tbl2 = tableGrob(aggreg_stats, rows=NULL, theme=tt)


gra = grid.arrange(arrangeGrob(pC + theme(legend.position="none"), pR, pdista, pdebt_cycle, ncol = 4, nrow = 1),
             arrangeGrob(pH, pM, pMort, peq_cycle, ncol = 4, nrow = 1),
             arrangeGrob(pA, pD, pdisth, pL, ncol = 4, nrow = 1),
             arrangeGrob(pRR, pPrepay, pequity, pPP, ncol = 4, nrow = 1), 
             ncol = 1, nrow = 4, heights=c(1,1,1,1))

print(gra)


XXX=apply(Default, c(2,6), mean)

dgrid = array(0,nd)
size = nd;
dstep = (dmax - dmin) /(size - 1);
for(i in 1:nd){
  dgrid[i] = dmin + (i-1)*dstep;
}
ygrid = exp(Tgrid(ny, ssigma_y, llambda_y, m_y))

# mesh = expand.grid(seq(1,ny,by=1), seq(1,nd,by=1))
# 
# level = data.frame(def = matrix(XXX[,], ny*nd, 1), ygrid = mesh$Var1, dgrid = mesh$Var2)
# 
# ggplot(level, aes(x=dgrid, y=ygrid)) + geom_raster(aes(fill = def), interpolate = TRUE) + 
#   scale_fill_continuous(values=c("gray41", "gray81")) +
#   labs(x = "Net Wealth ($)", y = "log(shock)") + ggtitle("Renting vs Owning") +
#   geom_contour(aes(x=xgrid, y=egrid, z = probability, colour = ..level..), size=2, bins = 8) + 
#   scale_colour_continuous(limits = c(0,0.1), guide=FALSE)


# ------------------------------------------------#
#        Graphs to understand default in eq       #
# ------------------------------------------------#

agrid = grida(na,amin,amax)
defaultAxY = data.frame(default = expand.grid(cDefAxY)$Var1, 
                        assets = expand.grid(agrid,ygrid)$Var1, 
                        income = expand.grid(agrid,ygrid)$Var2)
ggplot(defaultAxY, aes(x=assets, y=log(income))) + geom_raster(aes(fill = default), interpolate = TRUE) 


defaultAxY = data.frame(default = expand.grid(cDefTxAxY[6,,])$Var1, 
                        assets = expand.grid(agrid,ygrid)$Var1, 
                        income = expand.grid(agrid,ygrid)$Var2)
ggplot(defaultAxY, aes(x=assets, y=log(income))) + geom_raster(aes(fill = default), interpolate = TRUE) 



# People over state space
pipolAxY = data.frame(people = expand.grid(cPipolAxY)$Var1, 
                        assets = expand.grid(agrid,ygrid)$Var1, 
                        income = expand.grid(agrid,ygrid)$Var2)
ggplot(pipolAxY, aes(x=assets, y=log(income))) + geom_raster(aes(fill = people), interpolate = TRUE) 


# State space: leverage over Asset x Income
leverAxY = data.frame(leverage = expand.grid(cLevAxY)$Var1, 
                      assets = expand.grid(agrid,ygrid)$Var1, 
                      income = expand.grid(agrid,ygrid)$Var2)
ggplot(leverAxY, aes(x=assets, y=log(income))) + geom_raster(aes(fill = leverage), interpolate = TRUE) 


# State space: Housing over Asset x Income
houseAxY = data.frame(housing = expand.grid(cHevAxY)$Var1, 
                      assets = expand.grid(agrid,ygrid)$Var1, 
                      income = expand.grid(agrid,ygrid)$Var2)
ggplot(houseAxY, aes(x=assets, y=log(income))) + geom_raster(aes(fill = housing), interpolate = TRUE) 


mortAxY = data.frame(mort = expand.grid(cMevAxY)$Var1, 
                      assets = expand.grid(agrid,ygrid)$Var1, 
                      income = expand.grid(agrid,ygrid)$Var2)
ggplot(mortAxY, aes(x=assets, y=log(income))) + geom_raster(aes(fill = mort), interpolate = TRUE) 



hgrid = gridh(nh,hmin,hmax)
mgrid = gridm(nm,mmin,mmax)
defaultHxM = data.frame(default = expand.grid(cDefHxM)$Var1, 
                        housing = expand.grid(hgrid,mgrid)$Var1, 
                        mortgag = expand.grid(hgrid,mgrid)$Var2)
ggplot(defaultHxM, aes(x=housing, y=mortgag)) + geom_raster(aes(fill = default), interpolate = TRUE) 


