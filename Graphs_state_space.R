
Ph = 1;
ddeltabar = 0;
repay_coeff = c(7.084407,7.041049,6.994780,6.942594,6.883242,6.816373,6.741416,6.657291,6.562933,6.457651,6.340766,6.211675,6.070287,5.916850,5.750890,5.571630,5.377996,5.168729,4.942374,4.696728,4.426358,4.127028,3.796141,3.430416,3.025226,2.573416,2.065566,1.485583,0.808801,0.000000)
lgrid = seq(0,0.4,by=0.2)
eq_dist = c()

for(it in 1:T){
  for(im in 1:nm){
    for(ih in 2:nh){
      for(id in 1:nd){
        equity = (Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih] - mgrid[im]*(1+repay_coeff[it]))/(Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih])
        ex = 0
        if(length(eq_dist)>0){
          for(k in 1:length(eq_dist)){
            if(equity == eq_dist[k]){
              ex = 1
            }
          }
        }
        if(ex == 0){
          eq_dist = c(eq_dist, equity)
        }
      }
    }
  }
}
eq_dist = sort(eq_dist)

PTI_dist = c(0)
incshock = array(0, dim=c(T))
for(it in 1:T){
  if(it >= 1 & it < 5){
    incshock[it] = 0.128;
  } else if(it >= 5 && it < 10){
    incshock[it] = 0.111;
  } else if(it >= 10 && it < 15){
    incshock[it] = 0.088;
  } else if(it >= 15 && it < 20){
    incshock[it] = 0.096;
  } else if(it >= 20 && it < 25){
    incshock[it] = 0.044;
  } else{
    incshock[it] = 0.0;
  }
}

for(it in 1:T){
  for(iy in 1:ny){
    for(id in 1:nd){
      for(ih in 1:nh){
        for(im in 1:nm){
          for(ia in 1:na){
            if(Policyl[it, iy, ia, im, ih, id]>0){
              PTI = mgrid[im]/(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl[it, iy, ia, im, ih, id]+1] + agrid[ia]*r + q*hgrid[ih])
              
              ex = 0
              if(length(PTI_dist)>0){
                for(k in 1:length(PTI_dist)){
                  if(PTI == PTI_dist[k]){
                    ex = 1
                  }
                }
              }
              if(ex == 0){
                PTI_dist = c(PTI_dist, PTI)
              }
            }
          }
        }
      }
    }
  }
}
PTI_dist = sort(PTI_dist)

sspace = array(0, dim = c(length(eq_dist), na))
rspace = array(0, dim = c(length(eq_dist), na))
pspace = array(0, dim = c(length(eq_dist), na))

saspace = array(0, dim = c(length(eq_dist), T))
raspace = array(0, dim = c(length(eq_dist), T))
paspace = array(0, dim = c(length(eq_dist), T))

sspacecr = array(0, dim = c(length(eq_dist), length(PTI_dist)))
rspacecr = array(0, dim = c(length(eq_dist), length(PTI_dist)))
pspacecr = array(0, dim = c(length(eq_dist), length(PTI_dist)))

pipol = 0.0
full_equity = 0.0

for(it in 1:T){
  for(iy in 1:ny){
    for(ia in 1:na){
      for(im in 1:nm){
        for(ih in 2:nh){
          for(id in 1:nd){
            equity = (Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih] - mgrid[im]*(1+repay_coeff[it]))/(Ph*(1-dgrid[id] - ddeltabar)*hgrid[ih])
            i = 1
            while(equity != eq_dist[i]){
              i=i+1
            }
            pspace[i, ia] = pspace[i, ia] + Puncond[it, iy, ia, im, ih, id]
            sspace[i, ia] = sspace[i, ia] + Puncond[it, iy, ia, im, ih, id]*Default[it, iy, ia, im, ih, id]
            rspace[i, ia] = rspace[i, ia] + Puncond[it, iy, ia, im, ih, id]*Renew[it, iy, ia, im, ih, id]*(Policyh[it, iy, ia, im, ih, id]>1)

            paspace[i, it] = paspace[i, it] + Puncond[it, iy, ia, im, ih, id]
            saspace[i, it] = saspace[i, it] + Puncond[it, iy, ia, im, ih, id]*Default[it, iy, ia, im, ih, id]
            raspace[i, it] = raspace[i, it] + Puncond[it, iy, ia, im, ih, id]*Renew[it, iy, ia, im, ih, id]*(Policyh[it, iy, ia, im, ih, id]>1)
            
            PTI = mgrid[im]/(ygrid[iy]*(1-incshock[it])*eprocess[it]*lgrid[Policyl[it, iy, ia, im, ih, id]+1] + agrid[ia]*r + q*hgrid[ih])
            if(Policyl[it, iy, ia, im, ih, id] == 0){
              k = 1
            } else{
              k = 1
              while(PTI != PTI_dist[k]){
                k=k+1
              }
            }

            pspacecr[i, k] = pspacecr[i, k] + Puncond[it, iy, ia, im, ih, id]
            sspacecr[i, k] = sspacecr[i, k] + Puncond[it, iy, ia, im, ih, id]*Default[it, iy, ia, im, ih, id]
            rspacecr[i, k] = rspacecr[i, k] + Puncond[it, iy, ia, im, ih, id]*Renew[it, iy, ia, im, ih, id]*(Policyh[it, iy, ia, im, ih, id]>1)
            
            pipol = pipol + Puncond[it, iy, ia, im, ih, id]
            
            if(equity == 1){
              full_equity = full_equity + Puncond[it, iy, ia, im, ih, id]
            }
          }
        }
      }
    }
  }
}

full_equity = full_equity/pipol

ss = array(0, dim=c(13, na))
rr = array(0, dim=c(13, na))
pp = array(0, dim=c(13, na))

sas = array(0, dim=c(13, T))
rar = array(0, dim=c(13, T))
pap = array(0, dim=c(13, T))

sscr = array(0, dim=c(13, 10))
rrcr = array(0, dim=c(13, 10))
ppcr = array(0, dim=c(13, 10))

for(ia in 1:na){
  for(i in 1:length(eq_dist)){
    for(k in 1:13){
      if(eq_dist[i] > (k-4)*0.1 && eq_dist[i] <= (k-3)*0.1){
        ss[k, ia] = ss[k, ia] + sspace[i, ia]
        rr[k, ia] = rr[k, ia] + rspace[i, ia]
        pp[k, ia] = pp[k, ia] + pspace[i, ia]
      }
    }
  }
}

distr_eq = array(0,dim=c(13))
distr_PTI = array(0,dim=c(10))

for(iP in 1:length(PTI_dist)){
  for(kP in 1:10){
    for(ie in 1:length(eq_dist)){
      for(ke in 1:13){
        if(eq_dist[ie] > (ke-4)*0.1 && eq_dist[ie] <= (ke-3)*0.1 && PTI_dist[iP] > kP*0.1 && PTI_dist[iP] <= (kP+1)*0.1){
          sscr[ke, kP] = sscr[ke, kP] + sspacecr[ie, iP]
          rrcr[ke, kP] = rrcr[ke, kP] + rspacecr[ie, iP]
          ppcr[ke, kP] = ppcr[ke, kP] + pspacecr[ie, iP]
          
          distr_eq[ke] = distr_eq[ke] + pspacecr[ie, iP]
          distr_PTI[ke] = distr_PTI[ke] + pspacecr[ie, iP]
        }
      }
    }
  }
}

persp(x=seq(-0.3,0.9,by=0.1), y=seq(0.1,1,by=0.1), z=ppcr, theta = 90, phi = 30, col = "lightblue")


persp(x=seq(-0.3,0.9,by=0.1), y=agrid, z=pp, theta = 310, phi = 30, col = "lightblue")

for(it in 1:T){
  for(i in 1:length(eq_dist)){
    for(k in 1:13){
      if(eq_dist[i] > (k-4)*0.1 && eq_dist[i] <= (k-3)*0.1){
        sas[k, it] = sas[k, it] + saspace[i, it]
        rar[k, it] = rar[k, it] + raspace[i, it]
        pap[k, it] = pap[k, it] + paspace[i, it]
      }
    }
  }
}

pap = pap/pipol
rar = rar/pip
pp = pp/pipol

# for(it in 1:T){
#   for(i in 1:13){
#     sas[i, it] = sas[i, it]/pap[i, it]
#     rar[i, it] = rar[i, it]/pap[i, it]
#   }
# }

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/owners_distribution_asset.pdf",width=6,height=6)
persp(x=seq(-30,90,by=10), y=agrid, z=pp, theta = 310, phi = 15, col = "indianred1", shade = 0.5, zlim=c(0,0.105), ticktype = "detailed", xlab="Equity (%)", ylab="Assets ($)", zlab="%")
dev.off()


pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/owners_distribution.pdf",width=6,height=6)
persp(x=seq(-30,90,by=10), y=seq(22,20+2*T,by=2), z=pap, theta = 310, phi = 15, col = "indianred1", shade = 0.5, zlim=c(0,0.028), ticktype = "detailed", xlab="Equity (%)", ylab="Age", zlab="%")
dev.off()


pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/refinancing_state.pdf",width=7,height=7)
persp(x=seq(-30,90,by=10), y=seq(22,20+2*T,by=2), z=rar, theta = 310, phi = 15, col = "indianred1", shade = 0.5, zlim=c(0,0.135), ticktype = "detailed", xlab="Equity (%)", ylab="Age", zlab="%")
dev.off()

persp(x=seq(-0.3,0.9,by=0.1), y=seq(1,T,by=1), z=sas, theta = 310, phi = 10, col = "lightblue", zlim=c(0,0.5))

# for(ia in 1:na){
#   for(i in 1:13){
#     ss[i, ia] = ss[i, ia]/pp[i, ia]
#     rr[i, ia] = rr[i, ia]/pp[i, ia]
#   }
# }

ss = array(0, dim=c(13, na))

