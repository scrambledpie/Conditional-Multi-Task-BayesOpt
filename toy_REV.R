##################################################################################################
##################################################################################################
rm(list=ls())
library(MASS)
library(FastGP)
library(Rcpp)
outernew = function (X, Y,Z=1) {
  dX <- length(X)
  dY <- length(Y)
  robj<- rep(Y,rep.int(dX,dY)) - rep(X,dY)
  dim(robj) <- c(dX, dY)
  robj
}
tile     = function(m,n)matrix(rep(t(m),n),ncol=ncol(m),byrow = T)
inflate  = function(m,n)matrix(rep(m,each = n),ncol=ncol(m))
nearestneigthbor = function(XA){
  SQ = lapply(1:ncol(XA), function(d)outernew(XA[,d],XA[,d])^2)
  
  Edist = Reduce("+",SQ)
  
  diag(Edist) = 10000000000
  
  sqrt(range(apply(Edist,1,min)))
}

# set dimensions and hyperparameters
verbose  = TRUE
dimx     = 1
dima     = 1
dims     = dimx + dima
lX       = c(10,10,16,22,28,33)[dims]
TrueH    = c( rep(lX,dims), 100)


#######################################################################
#######################################################################
# Synthetic test functions made with SE kernel

SEkernel = function(XA1,XA2,Hp = TrueH){
  
  if(is.null(dim(XA1)))XA1 = matrix(XA1,1)
  if(is.null(dim(XA2)))XA2 = matrix(XA2,1)
  dims = ncol(XA1)
  
  ilXA = -0.5/Hp[1:dims]^2
  
  SQ = lapply(1:dims,function(i){   outernew(XA1[,i],XA2[,i])^2 *ilXA[i]  })
  
  (Hp[dims+1])*exp( Reduce(function(a,b)a+b,SQ) )
  
}
LHCdim = function(N0,dimss){
  sapply(1:dimss,function(d){  (sample(1:N0) - runif(N0))  })*(100/N0)
}
GenFun=function(seed=1,dims = 2){
  set.seed(seed)
  NN0 = 250*dims
  XA0 = LHCdim(NN0,dims)
  TrueK  = SEkernel(XA0,XA0) + diag(0.01,NN0)
  TrueiK = rcppeigen_invert_matrix(TrueK)
  y = mvrnorm(1,rep(0,NN0),TrueK)
  iKy=TrueiK%*%y
  result=function(X){
    if(is.null(dim(X)))X=matrix(X,1)
    Ks = SEkernel(X,XA0)
    Ks%*%iKy
  }
  return(result)
}

UniDens       =  {list(
  ddens=function(x)rep(1/100^dimx,length(x)/dimx),
  rdens = function(n)100*matrix(runif(n*dimx),n),
  rlhc = function(n)LHCdim(n,dimx)
)}
TriDens       =  {list(
  ddens   = function(x){ x = matrix(x,ncol=dimx,byrow = T);  x[,1]/(0.5*100^(dimx+1)) },
  rdens   = function(n){ X = matrix(runif(dimx*n),n,dimx); X[,1] = sqrt(X[,1]); 100*X },
  rlhc    = function(n){
    X     = LHCdim(n,dimx);
    X[,1] = 100*sqrt(X[,1]*0.01);
    X
  },
  name = "TriDens"
)}
TestDistros = list(UniDens, TriDens)


##################################################################################################
##################################################################################################
# Gaussian Process Regression Model and parameters, all global variables because I like my code to be clunky and hard to reuse
# we assume hyperparameters are known in advance in this code

H        = TrueH
GPmean0  = 0

invK     = matrix(0,1,1)
invKy    = rep(0,1)

UpdateGP = function(fixedpars=1,q=1,r=1){
  KK = SEkernel(Data[,1:dims],Data[,1:dims]) + diag(V2,nrow(Data))
  invK<<-rcppeigen_invert_matrix(KK)
  invKy<<-invK%*%(Data[,dims+1]-GPmean0)
  XAn <<- Data[,1:dims]
}
MU       = function(xai){
  if(is.null(dim(xai)))xai=matrix(xai,1)
  
  Ks  = SEkernel(xai,XAn)
  out = Ks%*% invKy
  GPmean0 + out[1:length(out)]
}
COV      = function(xa1,xa2,r=1){
  
  if(is.null(dim(xa1)))xa1=matrix(xa1,1)
  if(is.null(dim(xa2)))xa2=matrix(xa2,1)
  
  Ks1 = SEkernel(xa1,XAn)
  Ks2 = SEkernel(XAn,xa2)
  Kss = SEkernel(xa1,xa2)
  
  Kss - Ks1 %*% invK %*% Ks2
  
}
COV2     = function(xa1,xa2,r=2){
  
  xa1 = xa1[1:dims]
  
  Kss  = SEkernel(xa1,xa2)
  Neig = Kss > H[dims+1]*exp(-0.5*r^2)
  
  xa2  = xa2[Neig,,drop=FALSE]
  
  Ks1  = SEkernel(xa1, XAn)
  Ks2  = SEkernel(XAn, xa2)
  
  out2 = rep(0,length(Neig))
  out2[Neig] = Kss[Neig] - Ks1 %*% invK %*% Ks2
  out2
}
SD       = function(xa){
  if(!is.null(dim(xa))){
    A=apply(xa,1,function(xai)abs(COV(xai,xai)))
    A=A[1:length(A)]
  }else{
    A=abs(COV(xa,xa))[1]
  }
  A
}

##################################################################################################
##################################################################################################
# Sampling Methods

# utility functions, single task knowledge gradient (KG function from the paper), firstly an R version 
# and an old faster cpp version that does not include the filtering steps at the start of the R version.
KGCB      = function(a,b){
  
  if(all(abs(b)<0.000001))return(0)
  
  big = abs(b)>0.000001
  if(any(!big)>0){
    a   = c( a[big], max(a[!big]) )
    b   = c( b[big], 0)
  }
  
  n=length(a)
  
  if(n==1)return(0)
  
  O=order(b)
  b=b[O]
  a=a[O]
  #sort (a,b) by ascending b
  
  A     = 1
  Alast = 1
  C=-Inf
  
  while(Alast<n){
    s1     = Alast 
    sI     = (s1+1):n
    CI     = -(a[s1]-a[sI])/(b[s1]-b[sI])
    
    bestsI = which.min(CI)
    Alast  = sI[bestsI]
    A      = c(A,Alast)
    C      = c(C,CI[bestsI])
  }
  C=c(C,Inf)
  
  pC=pnorm(C)
  dC=dnorm(C)
  
  sum(  a[A]*(pC[-1]-pC[-length(pC)]) - b[A]*(dC[-1]-dC[-length(pC)])  )-max(a)
  
}
if(T){
  cppFunction('NumericMatrix KGCB_sorted(NumericVector a,NumericVector b) {
              int n = a.size()-1;
              int s = 2;
              double z = 0.0;
              
              NumericVector I(n);
              NumericVector Z(n);
              bool loopagain = true;
              I[1]=1; I[2]=2;
              Z[1]=-10000000000000; Z[2]=(a[1]-a[2])/(b[2]-b[1]);
              
              for(int i=3; i<n+1; i++){
              do{
              z=(a[i]-a[I[s]])/(b[I[s]]-b[i]);
              if(z<Z[s]){
              loopagain = true;
              s=s-1;
              }else{
              loopagain = false;
              s=s+1;
              I[s]=i;
              Z[s]=z;
              }
              }while(loopagain);
              
              }
              
              NumericMatrix out(s,2);
              for(int i =0; i<s; i++){out(i,0)=I[i+1]; out(i,1)=Z[i+1];}
              return out;
}')
}
KGCBcpp   = function(ai,bi){
  if(max(abs(bi))<1e-10)return(0)
  
  # baddies=which(abs(b)<1e-10)
  # if(length(baddies)>0){
  #   b[baddies]=0
  #   maxa=max(a[baddies])
  #   a=c(a[b!=0],maxa)
  #   b=c(b[b!=0],0)
  # }
  
  O=order(bi)
  ai = ai[O]
  bi = bi[O]
  
  # db=which(b[-1]-b[-length(b)]<1e-10)
  # bottoms=a[db]>a[db+1]
  
  
  if(!(length(ai)>0 | length(bi)>0))browser()
  O = KGCB_sorted(c(0,ai),c(0,bi))
  
  Z=c(-Inf,O[-1,2],Inf)
  I=O[,1]
  
  pC=pnorm(Z)
  dC=dnorm(Z)
  
  sum(  ai[I]*(pC[-1]-pC[-length(pC)]) - bi[I]*(dC[-1]-dC[-length(pC)])  ) - max(ai)
}

# Each function returns a new recomended sample point, they define Expected Improvment function over XxA which is passed to EI_optimizer
UNIdesign     = function(N0,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  
  X = Dens1$rlhc(N0)[sample(1:N0),]
  A = sapply(1:dima,function(d){  (sample(1:N0) - runif(N0))*(100/N0)  })
  xa = cbind(X,A)
  y  = TestFun(xa)+rnorm(N0,0,TnoiseSD)
  Data<<-cbind(xa,y)
}
MaxLEVI       = function(Na=NULL){
  if(is.null(Na))Na=nrow(Data)
  
  Am      = LHCdim(Na,dima)
  
  
  LEVI=function(xa){
    X1   = xa[1:dimx]
    A1   = xa[dimx+1:dima]
    X1A  = cbind( matrix(X1,Na+1,dimx,byrow = T) , rbind(A1,Am))
    sigt = as.numeric( COV2(xa,X1A) )
    sigt = sigt / sqrt(abs(sigt[1])+V2)
    MM   = MU(X1A)
    OO   = Dens1$ddens(X1)*KGCB(MM,sigt) 
    
    
    return( OO)
    
  }
  
  topxa = EI_optimizer(LEVI)
  
  return(topxa)
  
}
MaxConvLEVI   = function(Na=NULL){
  if(is.null(Na))Na=nrow(Data)
  
  Am      = LHCdim(Na,dima)
  
  KX = Dens1$rlhc(1000)
  Dens2 = function(x)sum(SEkernel(x,KX))
  
  LEVI=function(xa){
    X1   = xa[1:dimx]
    A1   = xa[dimx+1:dima]
    X1A  = cbind( matrix(X1,Na+1,dimx,byrow = T) , rbind(A1,Am))
    sigt = as.numeric( COV2(xa,X1A) )
    sigt = sigt / sqrt(abs(sigt[1])+V2)
    MM   = MU(X1A)
    OO   = Dens2(X1)*KGCB(MM,sigt) 
    
    
    return( OO)
    
  }
  
  topxa = EI_optimizer(LEVI)
  
  return(topxa)
  
}
MaxSparseREVI = function(Nx=NULL,Na=NULL,r=2){
  
  # Nx defines the number of X values to use for refence points, likewise for Na, therefore there will be NxNa ref points in total
  # r defines the radius around the proposed sample point for which refernce points will be included in REVI calculation
  # delta thus gives the covriance threshold for a given radius for all reference points
  
  delta = H[dims+1]*exp(-0.5*r^2)
  
  if(is.null(Na))Na = nrow(Data)
  if(is.null(Nx))Nx = nrow(Data) 
  
  
  # {Xm}*{Am} defines the reference points over which improvement is measured, precalculate means and parts of the posterior covariance
  Am     = LHCdim(Na,dima)
  Xm     = Dens1$rlhc(Nx)
  XmAm   = cbind( inflate(Xm,Na), tile(Am,Nx)  ) 
  MM     = matrix(MU(XmAm),Na,Nx)
  Ks0    = SEkernel(XAn,XmAm)
  iKKs0  = invK%*%Ks0
  
  
  ################################################################################
  sigt = rep(0,nrow(XmAm))
  REVI=function(xa){
    X1 = xa[1:dimx]
    A1 = xa[dimx + 1:dima]
    
    Ks1  = SEkernel(xa,XAn)
    ysd = 1/sqrt(abs(COV(xa,xa))[1]+V2)
    
    
    Kss    = SEkernel(xa,XmAm)
    Neig   = Kss > delta
    
    sigt[Neig] =  Kss[Neig] - Ks1 %*% iKKs0[,Neig,drop=FALSE]
    
    sigt = matrix( sigt, Na, Nx)*ysd
    
    KGm  = sum(sapply(1:Nx,function(i)KGCB(MM[,i],sigt[,i])))
    
    return(KGm)
    
    
  }
  ########################################################################################
  
  topxa = EI_optimizer(REVI)
  
  return(topxa)
  
}
MaxPEI        = function(Ns=NULL){
  td = max(Data[,dims+1])
  Target=function(x1){
    x1A = cbind(matrix(x1,99,dimx,byrow = T),1:99)
    M1  = MU(x1A)
    aa  = which.max(M1)
    M1  = M1[aa]
    if(M1<td) M1 = optimise(function(a)MU(c(x1,a)),interval = c(aa-1,aa+1),maximum = T,tol=0.001)$objective
    min(td,M1)
  }
  
  PEI=function(xa){
    S    = sqrt(abs(COV(xa,xa)[1]))
    M    = MU(xa)
    
    x1   = xa[1:dimx]
    T1   = Target(x1)
    
    z    = (T1 - M)/S
    Dens1$ddens(x1)*(pnorm(-z)*(M-T1) + dnorm(z)*S)
    
  }
  
  topxa = EI_optimizer(PEI,Ns)
  
  return(topxa)
  
}

EI_optimizer = function(EIfun,Ns=NULL){
  
  if (is.null(Ns)) Ns = (dims+2)*min(nrow(Data)/2,30)
  
  
  sXA = LHCdim(Ns,dims)
  
  topO  = -Inf
  botO  = Inf
  topxa = 0
  t=0
  Obj = function(xa){
    if(any(abs(xa-50)>50)) return(0)
    
    out = EIfun(xa)
    
    if (out > topO){
      topO  <<- out
      topxa <<- xa
    }
    if(out<botO){
      botO <<- out
    }
    out
  }
  
  OO=sapply(1:Ns,function(i)optim(sXA[i,],Obj,method = "Nelder-Mead",control = list(maxit = 50,fnscale = -1)))
  
  topxa
}


##################################################################################################
##################################################################################################
# Benchmark Utilities, measure the cost of the mapping

PredCost=function(x0=Testx0,truth=F){
  
  if(!truth){ # for each task x1 find the peak of the true function
    
    besta0=function(x1){
      x1A = cbind( matrix(x1,99,dimx,byrow=T), 1:99 )
      top = which.max(MU(x1A) )
      top  = optimise(function(a)MU(c(x1,a)),interval = c(top-1,top+1),maximum = T,tol=0.001)$maximum
      TestFun(c(x1,top))
    }
    
  }else{ # for each task x1 find the peak of the GP mean and evaluate the true function in the same place
    
    besta0=function(x1){
      x1A        =  cbind(matrix(x1,99,dimx,byrow=T), 1:99)
      top        =  which.max(TestFun(x1A))
      optimise(function(a)TestFun(c(x1,a)),interval = c(top-1,top+1),maximum = T,tol=0.001)$objective
    }
    
  }
  
  x0 = matrix(x0,ncol = dimx)
  return(mean(apply(x0,1,besta0)))
  
}

plotTestFun = function(TF,showD=F){
  X  = seq(0,100,len=50)
  A  = seq(0,100,len=40)
  XA = cbind(X,rep(A,each=length(X)))

  Y  = matrix(TF(XA),length(X),length(A))

  tops = apply(Y,1,which.max)
  contour(X,A,Y, main="GP mean", xlab="task/state", ylab="tool/action")

  lines(X,A[tops],col="blue",lwd=2)

  if(showD)points(Data[,1],Data[,2],pch=19)

}


##################################################################################################
##################################################################################################
# Initialize benchmark, choose algortihm, task distro, sampling budget


method     = 3                     # sampling methods, ffrom uniform to REVI and Profile-EI


cat("Starting, X-dimensions: ",dimx,"\n")
Dens1      = TestDistros[[1]]      # task space distro, 1: uniform, 2: triangular
Seed       = 100              
TestFun    = GenFun(Seed,dims)     # synthetic test function generator
set.seed(Seed)
N0       = dims*10                 # initial and final sample sizes, test set of task values and true optimal mapping
budget   = N0 + 60

TnoiseSD   = 0.1                   # True noise SD, and variance, set GP hyper parameter to the true value
TV2        = TnoiseSD^2
V2         = TV2

Testx0   = Dens1$rlhc((dims*100)+100); if(dimx==1) Testx0 = sort(Testx0)
TrueP    = PredCost(x0 =Testx0,truth=T)

UNIdesign(N0)                      # initialise the GP with uniform data, update the model
UpdateGP()

Cost     = data.frame(N=N0,OC=TrueP-PredCost()) # performance recording


if(verbose)cat("Uni...",N0,"   ", tail(Cost,1)[[2]],"\n")



################################################################################################################################
################################################################################################################################
# Run benchmark

while(nrow(Data)< budget){
  if(verbose)par(mfrow=c(2,1))
  N0 = nrow(Data)
  
  # select new point to sample and add it to Data  
  
  if(method==1){
    if(verbose)cat("LHC...",nrow(Data)+1,"   ")
    UNIdesign(nrow(Data)+1,Seed+1)
    
  }else if(method==2){
    if(verbose)cat("ConvLEVI...",N0+1,"   ")
    xanew=MaxConvLEVI()
    Data = rbind(Data, c(xanew,TestFun(xanew)+rnorm(1,0,TnoiseSD)) )
    
  }else if(method==3){
    if(verbose)cat("SparseREVI...",N0+1,"   ")
    Na     = nrow(Data) 
    Nx     = (3+dimx)*ceiling(sqrt(nrow(Data) ))
    xanew  = MaxSparseREVI(Nx,Na,r = 2.15)
    Data   = rbind(Data, c(xanew,TestFun(xanew)+rnorm(1,0,TnoiseSD)) )
    
  }else if(method==4){
    if(verbose)cat("ProfileEI...",N0+1,"   ")
    xanew=MaxPEI()
    Data = rbind(Data, c(xanew,TestFun(xanew)+rnorm(1,0,TnoiseSD)) )
    
  }else if(method==5){
    if(verbose)cat("LEVI n...",N0+1,"   ")
    xanew = MaxLEVI()
    Data  = rbind(Data, c(xanew,TestFun(xanew)+rnorm(1,0,TnoiseSD)) )
  }
  
  # update model and measure performance 
  UpdateGP()
  Cost     =  rbind(Cost,c(nrow(Data),TrueP-PredCost()))
  
  if(verbose){
    cat(tail(Cost,1)[[2]],"\n")
    
    plot(Cost,log="y", main="Regret/Opportunity Cost")
    lines(Cost)
    
    plotTestFun(MU)
    points(Data[,1],Data[,2],pch=19)
    xy = Data[nrow(Data),]
    points(xy[1],xy[2],pch=19,col="red",cex=2)
    points(Data[nrow(Data),1],Data[nrow(Data),2],pch=19)
    
    Sys.sleep(0.5)
  }
}



##################################################################################################
##################################################################################################
