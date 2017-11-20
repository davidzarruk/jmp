
# Initial states:

sims = 10000

Cpath = array(0, dim=c(sims, T+1-it))
Apath = array(0, dim=c(sims, T+1-it))
Mpath = array(0, dim=c(sims, T+1-it))
Hpath = array(0, dim=c(sims, T+1-it))
Dpath = array(0, dim=c(sims, T+1-it))
Rpath = array(0, dim=c(sims, T+1-it))

Ppath = array(0, dim=c(sims, T-it))
Lpath = array(0, dim=c(sims, T-it))
Ypath = array(0, dim=c(sims, T-it))
Cpath = array(0, dim=c(sims, T-it))

idp = 0
iyp = 0

P = Rtauchen(ny, ssigma_y, llambda_y, m_y)
ygrid = exp(Tgrid(ny, ssigma_y, llambda_y, m_y))
agrid = grida(na, amin, amax)
dgrid = gridd(nd, dmin, dmax)
mgrid = gridm(nm, mmin, mmax)
hgrid = gridh(nh, hmin, hmax)
lgrid = seq(0, 0.4, length.out=3)

for(k in 1:sims){
  
  it = 6
  ia = 15
  im = 2
  ih = 3
  id = 1
  iy = 3
  
  idp = 0
  iyp = 0
  imp = 0
  
  Apath[k, 1] = agrid[ia]
  Mpath[k, 1] = mgrid[im]
  Hpath[k, 1] = hgrid[ih]
  Dpath[k, 1] = Default[it, iy, ia, im, ih, id]
  Rpath[k, 1] = Renew[it, iy, ia, im, ih, id]
  
  for(i in it:(T-1)){
    Apath[k, i-it+2] = agrid[Policya[i, iy, ia, im, ih, id]+1]
    Mpath[k, i-it+2] = mgrid[Policym[i, iy, ia, im, ih, id]+1]
    Hpath[k, i-it+2] = hgrid[Policyh[i, iy, ia, im, ih, id]+1]
    Dpath[k, i-it+2] = Default[i, iy, ia, im, ih, id]
    Rpath[k, i-it+2] = Renew[i, iy, ia, im, ih, id]
    
    Ppath[k, i-it+1] = Pricing[i, iy, Policya[i, iy, ia, im, ih, id]+1, Policym[i, iy, ia, im, ih, id]+1,
                               Policyh[i, iy, ia, im, ih, id]+1, id]
    Ypath[k, i-it+1] = iy
    Cpath[k, i-it+1] = Policyc[i, iy, ia, im, ih, id]
    Lpath[k, i-it+1] = Policyl[i, iy, ia, im, ih, id]
    
    dshock = runif(1,0,1)
    yshock = runif(1,0,1)
    mshock = runif(1,0,1)
    
    for(j in 1:nd){
      if(dshock >= ((j-1)/nd) && dshock < (j/nd)){
        idp = j
      }
    }
  
    if(mshock <= rrho){
      imp = Policym[i, iy, ia, im, ih, id]+1
    } else {
      imp = 1
    }
  
    Pcum = P[iy, 1]
    j = 1
    iyp = 1
    while(yshock > Pcum && j <= (ny-1)){
      j = j+1
      Pcum = Pcum + P[iy, j]
    }
    iyp = j
    
    ia = Policya[i, iy, ia, im, ih, id]+1
    ih = Policyh[i, iy, ia, im, ih, id]+1
    id = idp
    iy = iyp
    im = imp
  }
}

Cmean = apply(Cpath, c(2), mean)
Dmean = apply(Dpath, c(2), mean)
Hmean = apply(Hpath, c(2), mean)
Mmean = apply(Mpath, c(2), mean)
Amean = apply(Apath, c(2), mean)
Ymean = apply(Ypath, c(2), mean)
Lmean = apply(Lpath, c(2), mean)
Pmean = apply(Ppath, c(2), mean)

sH = ggplot(data = data.frame(H=Hmean, age=seq(it, T, length.out = (T+1-it))), aes(x=age)) + 
  geom_line(aes(y = H), size=size_line) + ggtitle("H")   + 
  scale_y_continuous(expand = c(0, 0)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

sM = ggplot(data = data.frame(M=Mmean, age=seq(it, T, length.out = (T+1-it))), aes(x=age)) + 
  geom_line(aes(y = M), size=size_line) + ggtitle("M")   + 
  scale_y_continuous(expand = c(0, 0)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

sA = ggplot(data = data.frame(A=Amean, age=seq(it, T, length.out = (T+1-it))), aes(x=age)) + 
  geom_line(aes(y = A), size=size_line) + ggtitle("A")   + 
  scale_y_continuous(expand = c(0, 0)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

sD = ggplot(data = data.frame(D=Dmean, age=seq(it, T, length.out = (T+1-it))), aes(x=age)) + 
  geom_line(aes(y = D), size=size_line) + ggtitle("D")   + 
  scale_y_continuous(expand = c(0, 0)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

sY = ggplot(data = data.frame(Y=Ymean, age=seq(it, T-1, length.out = (T-it))), aes(x=age)) + 
  geom_line(aes(y = Y), size=size_line) + ggtitle("Y")   + 
  scale_y_continuous(expand = c(0, 0)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

sC = ggplot(data = data.frame(C=Cmean, age=seq(it, T-1, length.out = (T-it))), aes(x=age)) + 
  geom_line(aes(y = C), size=size_line) + ggtitle("C")   + 
  scale_y_continuous(expand = c(0, 0)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

sL = ggplot(data = data.frame(L=Lmean, age=seq(it, T-1, length.out = (T-it))), aes(x=age)) + 
  geom_line(aes(y = L), size=size_line) + ggtitle("L")   + 
  scale_y_continuous(expand = c(0, 0)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

sP = ggplot(data = data.frame(P=Pmean, age=seq(it, T-1, length.out = (T-it))), aes(x=age)) + 
  geom_line(aes(y = P), size=size_line) + ggtitle("P")   + 
  scale_y_continuous(expand = c(0, 0)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

paths = grid.arrange(arrangeGrob(sH, sM, sC, ncol = 3, nrow = 1),
                     arrangeGrob(sD, sA, ncol = 3, nrow = 1), 
                     ncol = 1, nrow = 2)

print(paths)
