rm(list=ls(all=TRUE)) 
library(mnormt)
library(Rlab)
library(nlme)
library(CompQuadForm)
library(glmnet)
library(expm)
library(huge)


library(Formula)
library(gtable)
library(Rcpp)
library(colorspace)
library(munsell)
library(plyr)
library(scales)
library(ggplot2)
library(acepack)
library(latticeExtra)
library(RColorBrewer)
library(gridExtra)
library(Hmisc)
library(snow)
library(snowfall)
library(fdrtool)
library(sfsmisc,lib.loc="/Users/Yin/Downloads/")
library(reshape)
library(rags2ridges)
library(PDSCE)

n1=15;
n2=15
rho1=0.5;
p=200;
q=50
rep=100;
alpha=0.01;
M=2;
ll=(floor(p/M))*(floor(p/M)-1)/2;

I = array(0:0, dim=c(p,p))
for (i in 1:p){
  I [i,i]=1;
}
mu = rep(0,q)



e1 = rep(1,n1*q)
e2 = rep(1,n2*q)
ff = array(0:0, dim=c(rep,20))
fp = array(0:0, dim=c(rep,20))
fdp = rep(0,rep)
fpo = rep(0,rep)

ffp = array(0:0, dim=c(rep,20))
fpp = array(0:0, dim=c(rep,20))
fdpp = rep(0,rep)
fpop = rep(0,rep)

Affp = array(0:0, dim=c(rep,20))
Afpp = array(0:0, dim=c(rep,20))
Afdpp = rep(0,rep)
Afpop = rep(0,rep)

Bffp = array(0:0, dim=c(rep,20))
Bfpp = array(0:0, dim=c(rep,20))
Bfdpp = rep(0,rep)
Bfpop = rep(0,rep)

Wffp = array(0:0, dim=c(rep,20))
Wfpp = array(0:0, dim=c(rep,20))
Wfdpp = rep(0,rep)
Wfpop = rep(0,rep)

ffpO = array(0:0, dim=c(rep,20))
fppO = array(0:0, dim=c(rep,20))
fdppO = rep(0,rep)
fpopO = rep(0,rep)


################### Sigma for time
PHI = array(0:0, dim=c(q,q))
for (ik in 1:q){
	for (jk in 1:q){
		PHI[ik,jk] = 0.4^{abs(ik-jk)}
	}
}

II = array(0:0, dim=c(q,q))
for (i in 1:q){
  II [i,i]=1;
}


PHI2 = array(0:0, dim=c(q,q))
for (ik in 1:q){
	for (jk in 1:q){
		PHI2[ik,jk] = 0.5^{abs(ik-jk)}
	}
}

PHI3 = II
for (ik in 1:(q-5)){
	for (jk in (ik+1):(ik+3)){
		PHI3[ik,jk] = 1/{abs(ik-jk)+1}
		PHI3[jk,ik] = PHI3[ik,jk]
	}
}

eP = eigen (PHI3)
de = abs (min (eP$values))+0.05
PHI3 = (PHI3+de*II)/(1+de)

PHI4 = II
for (ik in 1:(q-6)){
	for (jk in (ik+1):(ik+4)){
		PHI4[ik,jk] = 1/{abs(ik-jk)+1}
		PHI4[jk,ik] = PHI4[ik,jk]
	}
}

eP = eigen (PHI4)
de = abs (min (eP$values))+0.05
PHI4 = (PHI4+de*II)/(1+de)


PHI=PHI3
PHI2=PHI4

PHIINV = sqrtm(solve(PHI))
PHI2INV = sqrtm(solve(PHI2))




Omega1 = huge.generator (d=p,graph = "band",g=3, verbose = FALSE)$omega
Omega1 = Omega1*(abs(Omega1)>0.001)
Sigma1 = solve(Omega1)

Omega2 = huge.generator(d=p, verbose = FALSE)$omega
Omega2 = Omega2*(abs(Omega2)>0.001)
Sigma2 = solve(Omega2)

Omega3 = huge.generator(d=p,prob=0.01, verbose = FALSE)$omega
Omega3 = Omega3*(abs(Omega3)>0.001)
Sigma3 = solve(Omega3)

Omega4 = huge.generator(d=p,graph = "hub",g=20, verbose = FALSE)$omega
Omega4 = Omega4*(abs(Omega4)>0.001)
Sigma4 = solve(Omega4)

Omega5= createS(n = 1, p = p,topology="small-world",precision=TRUE,banded.n=5)
Omega5 = Omega5*(abs(Omega5)>0.001)
eO5=eigen (Omega5)
de = abs (min (eO5$values))+0.05;
Omega5 = (Omega5 + de*I)
Sigma5 = solve(Omega5)




for (l in 1:rep){
	
	
Sigma=Sigma1;
Omega=Omega1
DD = diag(diag(Omega))
O = Omega
for (i in 1:(p)){
	for (j in i:p){
		O[j,i]=0
	}
}
nonzero = which(O!=0)
nn = length(nonzero)*0.5
aa = sample(which(O!=0),nn)
OOmega = Omega
OOmega[aa]=rep(0,nn)
for (i in 1:(p-1)){
	for (j in (i+1):p){
		OOmega[j,i]=OOmega[i,j]
	}

}

eO=eigen (OOmega)
de = abs (min (eO$values))+0.05;
OOmega = (OOmega + de*I)
Omega = (Omega + de*I)
Sigma = solve(Omega)
SSigma = solve(OOmega)


RO = I;
for (i in 1:(p-1)){
	for (j in (i+1):p){
		RO[i,j]= Omega[i,j]/sqrt(Omega[i,i]*Omega[j,j])
		RO[j,i]=RO[i,j]
	}
}
RRO = I;
for (i in 1:(p-1)){
	for (j in (i+1):p){
		RRO[i,j]= OOmega[i,j]/sqrt(OOmega[i,i]*OOmega[j,j])
		RRO[j,i]=RRO[i,j]
	}
}

SQS1 = sqrtm(Sigma)
SQS2 = sqrtm(SSigma)
OOO=RO-RRO;
OOO = OOO*(abs(OOO)>0.001)
OOO2=OOO*OOO;



bx = array(0:0, dim=c(p-1,p))
bxO = array(0:0, dim=c(p-1,p))
by = array(0:0, dim=c(p-1,p))
c1 = array(0:0, dim=c(n1*q,p))
T1 = array(0:0, dim=c(p,p))
nbx = array(0:0, dim=c(p-1,p))
nbxO = array(0:0, dim=c(p-1,p))
nby = array(0:0, dim=c(p-1,p))
c2 = array(0:0, dim=c(n2*q,p))
T2 = array(0:0, dim=c(p,p))
T = array(0:0, dim=c(p,p))
diff = array(0:0, dim=c(p,p))
S1 = array(0:0, dim=c((floor(p/M)),(floor(p/M))))
S = array(0:0, dim=c((floor(p/M)),(floor(p/M))))
OO = array(0:0, dim=c((floor(p/M)),(floor(p/M))))
sss=rep(0,20)
sssp=rep(0,20)
ssspO=rep(0,20)
Asssp=rep(0,20)
Bsssp=rep(0,20)
Wsssp=rep(0,20)
X = array(0:0, dim=c(n1,p,q))
CX = array(0:0, dim=c(n1,p,q))
XX = array(0:0, dim=c(n1,p,q))
XXO = array(0:0, dim=c(n1,p,q))
Y = array(0:0, dim=c(n2,p,q))
CY = array(0:0, dim=c(n1,p,q))
YY = array(0:0, dim=c(n2,p,q))
YYO = array(0:0, dim=c(n2,p,q))


#    x = rmnorm(n = n1, mean = mu, Sigma)

##########################################################################
##########################################################################
#data-driven
##########################################################################
##########################################################################

##for x
#separate kronecker product
for (nk in 1:n1){
	Z=rmnorm(n=p, mean=mu, PHI)
	X[nk,,]=SQS1%*%Z
}
PHI0=array(0:0, dim=c(q,q))
for (nk in 1:n1){
	PHI0=PHI0+t(X[nk,,]) %*% X[nk,,]
}
PHI0=PHI0/n1/p

######################
#adaptive thresholding
######################
# generate samples to do cv
for (nk in 1:n1){
	CZ=rmnorm(n=p, mean=mu, PHI)
	CX[nk,,]=SQS1%*%CZ
}
CPHI0=array(0:0, dim=c(q,q))
for (nk in 1:n1){
	CPHI0=CPHI0+t(CX[nk,,]) %*% CX[nk,,]
}
CPHI0X=CPHI0/n1/p

# get the adaptive threshold est
xxx=X[,1,]
for (i in 2:p){
	xxx=rbind(xxx,X[,i,])
}

PHI0INV = sqrtm(solve(PHI0))

for (nk in 1:n1){
	XX[nk,,]=X[nk,,] %*% PHI0INV
}
x = XX[,,1]
for(np in 2:q){
	x=rbind(x,XX[,,np])
}

##for y
#separate kronecker product
for (nk in 1:n2){
	Z=rmnorm(n=p, mean=mu, PHI2)
	Y[nk,,]=SQS2%*%Z
}
PHI0=array(0:0, dim=c(q,q))
for (nk in 1:n2){
	PHI0=PHI0+t(Y[nk,,]) %*% Y[nk,,]
}
PHI0=PHI0/n2/p

######################
#adaptive thresholding
######################
# generate samples to do cv
for (nk in 1:n2){
	CZ=rmnorm(n=p, mean=mu, PHI2)
	CY[nk,,]=SQS2%*%CZ
}
CPHI0=array(0:0, dim=c(q,q))
for (nk in 1:n2){
	CPHI0=CPHI0+t(CY[nk,,]) %*% CY[nk,,]
}
CPHI0Y=CPHI0/n2/p

# get the adaptive threshold est
yyy=Y[,1,]
for (i in 2:p){
	yyy=rbind(yyy,Y[,i,])
}

PHI0INV = sqrtm(solve(PHI0))

for (nk in 1:n2){
	YY[nk,,]=Y[nk,,] %*% PHI0INV
}
y = YY[,,1]
for(np in 2:q){
	y=rbind(y,YY[,,np])
}



#inverse regression

Shatmimix=cov(x);
sigmax=diag(Shatmimix);
Shatmimiy=cov(y);
sigmay=diag(Shatmimiy);
for (ii in 3:20){
      for (i in 1:p){
        newx=t(x);
        newx=newx[-i,];
         newy=t(y);
        newy=newy[-i,];
        
         glm2=glmnet(t(newx),x[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmax[i]*log(p)/(n1*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         bx[,i] = coef[-1]
         
         glm2=glmnet(t(newy),y[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmay[i]*log(p)/(n2*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         nbx[,i] = coef[-1]
       
            
        #calculate \xi
        c1[,i] = x[,i]-mean(x[,i])*e1- t(newx-(colMeans(t(newx)))%*%t(e1))%*%bx[,i];
        c2[,i] = y[,i]-mean(y[,i])*e2- t(newy-(colMeans(t(newy)))%*%t(e2))%*%nbx[,i];

    }
      R1=t(c1)%*%c1/n1/q;
      R2=t(c2)%*%c2/n2/q;
      s1=colMeans(c1*c1);
      s2=colMeans(c2*c2);
      for (i in 1:(p-1)){
          for (j in (i+1):p){
              T1[i,j]=-(R1[i,j]+s1[i]*bx[i,j]+s1[j]*bx[j-1,i]);
              T2[i,j]=-(R2[i,j]+s2[i]*nbx[i,j]+s2[j]*nbx[j-1,i]);
              diff[i,j]=(T1[i,j]/sqrt(R1[i,i]*R1[j,j])-T2[i,j]/sqrt(R2[i,i]*R2[j,j]))^2/(1/n1/q*(1+bx[i,j]^2*R1[i,i]/R1[j,j])+1/n2/q*(1+nbx[i,j]^2*R2[i,i]/R2[j,j]));

              
          }
      }


      z=seq(0:((floor(sqrt(4*log(p)))+2)*100))/100
      le=length(z);
      rr=rep(0,le)
      for (k in 1:le){
          rr[k]=p*(p-1)/2*(2-2*pnorm(z[k]))/max(sum(sqrt(diff)>=z[k]),1);
      }
      if (rr[le]<=alpha){
          t=min(which((rr-alpha)<=0))}
    #  else
    #      {t=le};
      

      
      ffp[l,ii]=sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t]))/max(sum(sum(sqrt(diff)>=z[t])),1);
      sssp[ii]=0;
      for (kk in 1:10){
          sssp[ii]=((sum(sum(sqrt(diff)>=qnorm(1-kk/floor(10/(1-pnorm(sqrt(log(p)))))))))/(kk/floor(10/(1-pnorm(sqrt(log(p)))))*(p*(p-1)))-1)^2+sssp[ii];
      }
      fpp[l,ii]=(max(sum(sum(sqrt(diff)>=z[t])),1)-sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t])))/(sum(sum(OOO!=0))/2);
      }
      

##########################################################################
##########################################################################
#using banded estimate by Bickel and Levina (2008) using package "PDSCE" 
##########################################################################
##########################################################################

qq=10
lambda = seq(1:qq)
dif = rep(0,qq)
for (j in 1:qq){
    dif[j] = sum(sum(abs(band.chol(xxx, j, centered = TRUE, method = "safe")-CPHI0X)))
    }
thr = max(lambda[which(dif==min(dif))])
BPHI0 = band.chol(xxx, thr, centered = TRUE, method = "safe")

BPHI0INV = sqrtm(solve(BPHI0))
XX = array(0:0, dim=c(n1,p,q))
for (nk in 1:n1){
	XX[nk,,]=X[nk,,] %*% BPHI0INV
}
x = XX[,,1]
for(np in 2:q){
	x=rbind(x,XX[,,np])
}

lambda = seq(1:qq)
dif = rep(0,qq)
for (j in 1:qq){
    dif[j] = sum(sum(abs(band.chol(yyy, j, centered = TRUE, method = "safe")-CPHI0Y)))
    }
thr = max(lambda[which(dif==min(dif))])
BPHI0 = band.chol(yyy, thr, centered = TRUE, method = "safe")

BPHI0INV = sqrtm(solve(BPHI0))
YY = array(0:0, dim=c(n2,p,q))
for (nk in 1:n2){
	YY[nk,,]=Y[nk,,] %*% BPHI0INV
}
y = YY[,,1]
for(np in 2:q){
	y=rbind(y,YY[,,np])
}



#inverse regression

Shatmimix=cov(x);
sigmax=diag(Shatmimix);
Shatmimiy=cov(y);
sigmay=diag(Shatmimiy);
for (ii in 3:20){
      for (i in 1:p){
        newx=t(x);
        newx=newx[-i,];
         newy=t(y);
        newy=newy[-i,];
        
         glm2=glmnet(t(newx),x[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmax[i]*log(p)/(n1*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         bx[,i] = coef[-1]
         
         glm2=glmnet(t(newy),y[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmay[i]*log(p)/(n2*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         nbx[,i] = coef[-1]
       
            
        #calculate \xi
        c1[,i] = x[,i]-mean(x[,i])*e1- t(newx-(colMeans(t(newx)))%*%t(e1))%*%bx[,i];
        c2[,i] = y[,i]-mean(y[,i])*e2- t(newy-(colMeans(t(newy)))%*%t(e2))%*%nbx[,i];

    }
      R1=t(c1)%*%c1/n1/q;
      R2=t(c2)%*%c2/n2/q;
      s1=colMeans(c1*c1);
      s2=colMeans(c2*c2);
      for (i in 1:(p-1)){
          for (j in (i+1):p){
              T1[i,j]=-(R1[i,j]+s1[i]*bx[i,j]+s1[j]*bx[j-1,i]);
              T2[i,j]=-(R2[i,j]+s2[i]*nbx[i,j]+s2[j]*nbx[j-1,i]);
              diff[i,j]=(T1[i,j]/sqrt(R1[i,i]*R1[j,j])-T2[i,j]/sqrt(R2[i,i]*R2[j,j]))^2/(1/n1/q*(1+bx[i,j]^2*R1[i,i]/R1[j,j])+1/n2/q*(1+nbx[i,j]^2*R2[i,i]/R2[j,j]));

              
          }
      }


      z=seq(0:((floor(sqrt(4*log(p)))+2)*100))/100
      le=length(z);
      rr=rep(0,le)
      for (k in 1:le){
          rr[k]=p*(p-1)/2*(2-2*pnorm(z[k]))/max(sum(sqrt(diff)>=z[k]),1);
      }
      if (rr[le]<=alpha){
          t=min(which((rr-alpha)<=0))} else
      {t=le};      

      
      Bffp[l,ii]=sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t]))/max(sum(sum(sqrt(diff)>=z[t])),1);
      Bsssp[ii]=0;
      for (kk in 1:10){
          Bsssp[ii]=((sum(sum(sqrt(diff)>=qnorm(1-kk/floor(10/(1-pnorm(sqrt(log(p)))))))))/(kk/floor(10/(1-pnorm(sqrt(log(p)))))*(p*(p-1)))-1)^2+Bsssp[ii];
      }
      Bfpp[l,ii]=(max(sum(sum(sqrt(diff)>=z[t])),1)-sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t])))/(sum(sum(OOO!=0))/2);
      }

##########################################################################
##########################################################################
#oracle
##########################################################################
##########################################################################

for (nk in 1:n1){
	XXO[nk,,]=X[nk,,] %*% PHIINV
}
xO = XXO[,,1]
for(np in 2:q){
	xO=rbind(xO,XXO[,,np])
}

for (nk in 1:n2){
	YYO[nk,,]=Y[nk,,] %*% PHI2INV
}
yO = YYO[,,1]
for(np in 2:q){
	yO=rbind(yO,YYO[,,np])
}


      ShatmimixO=cov(xO);
  sigmaxO=diag(ShatmimixO);
        ShatmimiyO=cov(yO);
  sigmayO=diag(ShatmimiyO);

for (ii in 3:20){
      for (i in 1:p){
        newxO=t(xO);
        newxO=newxO[-i,];
        newyO=t(yO);
        newyO=newyO[-i,]
                
         glm2=glmnet(t(newxO),xO[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmaxO[i]*log(p)/(n1*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         bxO[,i] = coef[-1]
         glm2=glmnet(t(newyO),yO[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmayO[i]*log(p)/(n2*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         nbxO[,i] = coef[-1]
                
            
        #calculate \xi
        c1[,i] = xO[,i]-mean(xO[,i])*e1- t(newxO-(colMeans(t(newxO)))%*%t(e1))%*%bxO[,i];
        c2[,i] = yO[,i]-mean(yO[,i])*e2- t(newyO-(colMeans(t(newyO)))%*%t(e2))%*%nbxO[,i];
    }
      R1=t(c1)%*%c1/n1/q;
      R2=t(c2)%*%c2/n2/q;
      s1=colMeans(c1*c1);
      s2=colMeans(c2*c2);
     for (i in 1:(p-1)){
          for (j in (i+1):p){
              T1[i,j]=-(R1[i,j]+s1[i]*bxO[i,j]+s1[j]*bxO[j-1,i]);
              T2[i,j]=-(R2[i,j]+s2[i]*nbxO[i,j]+s2[j]*nbxO[j-1,i]);              
              diff[i,j]=(T1[i,j]/sqrt(R1[i,i]*R1[j,j])-T2[i,j]/sqrt(R2[i,i]*R2[j,j]))^2/(1/n1/q*(1+bxO[i,j]^2*R1[i,i]/R1[j,j])+1/n2/q*(1+nbxO[i,j]^2*R2[i,i]/R2[j,j]));
              
              
          }
      }


      z=seq(0:((floor(sqrt(4*log(p)))+2)*100))/100
      le=length(z);
      rr=rep(0,le)
      for (k in 1:le){
          rr[k]=p*(p-1)/2*(2-2*pnorm(z[k]))/max(sum(sqrt(diff)>=z[k]),1);
      }
      if (rr[le]<=alpha){
          t=min(which((rr-alpha)<=0))} else
      {t=le};      

      
      ffpO[l,ii]=sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t]))/max(sum(sum(sqrt(diff)>=z[t])),1);
      ssspO[ii]=0;
      for (kk in 1:10){
          ssspO[ii]=((sum(sum(sqrt(diff)>=qnorm(1-kk/floor(10/(1-pnorm(sqrt(log(p)))))))))/(kk/floor(10/(1-pnorm(sqrt(log(p)))))*(p*(p-1)))-1)^2+ssspO[ii];
      }
      fppO[l,ii]=(max(sum(sum(sqrt(diff)>=z[t])),1)-sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t])))/(sum(sum(OOO!=0))/2);
      }




##########################################################################
##########################################################################
#compare with whitening using banded estimate by Bickel and Levina (2008) using package "PDSCE" + Xia et al.
##########################################################################
##########################################################################

for (i in 1:p){

#####################

CPHI0=t(CX[,i,]) %*% (CX[,i,])
#for (nk in 1:n1){
#	CPHI0=CPHI0+(CX[nk,i,]) %*% t(CX[nk,i,])
#}
CPHI0=CPHI0/n1

qq=q
lambda = seq(1:qq)
dif = rep(0,qq)
for (j in 1:qq){
    dif[j] = sum(sum(abs(band.chol(X[,i,], j, centered = TRUE, method = "fast")-CPHI0)))
    }
thr = max(lambda[which(dif==min(dif))])
APHI0 = band.chol(X[,i,], thr, centered = TRUE, method = "safe")
de = abs (min (eigen(APHI0)$values))+0.05
APHI0 = (APHI0+de*II)/(1+de)

#####################
XX[,i,]=X[,i,] %*% sqrtm(solve(APHI0))
    #for (nk in 1:n1){
	#     XX[nk,i,]=X[nk,i,] %*% sqrtm(solve(APHI0))
    #}  
}

x = XX[,,1]
for(np in 2:q){
	x=rbind(x,XX[,,np])
}

for (i in 1:p){

#####################

CPHI0=t(CY[,i,]) %*% (CY[,i,])
#for (nk in 1:n1){
#	CPHI0=CPHI0+(CX[nk,i,]) %*% t(CX[nk,i,])
#}
CPHI0=CPHI0/n2

qq=q
lambda = seq(1:qq)
dif = rep(0,qq)
for (j in 1:qq){
    dif[j] = sum(sum(abs(band.chol(Y[,i,], j, centered = TRUE, method = "fast")-CPHI0)))
    }
thr = max(lambda[which(dif==min(dif))])
APHI0 = band.chol(Y[,i,], thr, centered = TRUE, method = "safe")
de = abs (min (eigen(APHI0)$values))+0.05
APHI0 = (APHI0+de*II)/(1+de)

#####################
YY[,i,]=Y[,i,] %*% sqrtm(solve(APHI0))
    #for (nk in 1:n1){
	#     XX[nk,i,]=X[nk,i,] %*% sqrtm(solve(APHI0))
    #}  
}

y = YY[,,1]
for(np in 2:q){
	y=rbind(y,YY[,,np])
}


#inverse regression

Shatmimix=cov(x);
sigmax=diag(Shatmimix);
Shatmimiy=cov(y);
sigmay=diag(Shatmimiy);
for (ii in 3:20){
      for (i in 1:p){
        newx=t(x);
        newx=newx[-i,];
         newy=t(y);
        newy=newy[-i,];
        
         glm2=glmnet(t(newx),x[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmax[i]*log(p)/(n1*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         bx[,i] = coef[-1]
         
         glm2=glmnet(t(newy),y[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmay[i]*log(p)/(n2*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         nbx[,i] = coef[-1]
       
            
        #calculate \xi
        c1[,i] = x[,i]-mean(x[,i])*e1- t(newx-(colMeans(t(newx)))%*%t(e1))%*%bx[,i];
        c2[,i] = y[,i]-mean(y[,i])*e2- t(newy-(colMeans(t(newy)))%*%t(e2))%*%nbx[,i];

    }
      R1=t(c1)%*%c1/n1/q;
      R2=t(c2)%*%c2/n2/q;
      s1=colMeans(c1*c1);
      s2=colMeans(c2*c2);
      for (i in 1:(p-1)){
          for (j in (i+1):p){
              T1[i,j]=-(R1[i,j]+s1[i]*bx[i,j]+s1[j]*bx[j-1,i]);
              T2[i,j]=-(R2[i,j]+s2[i]*nbx[i,j]+s2[j]*nbx[j-1,i]);
              diff[i,j]=(T1[i,j]/sqrt(R1[i,i]*R1[j,j])-T2[i,j]/sqrt(R2[i,i]*R2[j,j]))^2/(1/n1/q*(1+bx[i,j]^2*R1[i,i]/R1[j,j])+1/n2/q*(1+nbx[i,j]^2*R2[i,i]/R2[j,j]));

              
          }
      }


      z=seq(0:((floor(sqrt(4*log(p)))+2)*100))/100
      le=length(z);
      rr=rep(0,le)
      for (k in 1:le){
          rr[k]=p*(p-1)/2*(2-2*pnorm(z[k]))/max(sum(sqrt(diff)>=z[k]),1);
      }
      if (rr[le]<=alpha){
          t=min(which((rr-alpha)<=0))}  else
      {t=le};
      

      
      Wffp[l,ii]=sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t]))/max(sum(sum(sqrt(diff)>=z[t])),1);
      Wsssp[ii]=0;
      for (kk in 1:10){
          Wsssp[ii]=((sum(sum(sqrt(diff)>=qnorm(1-kk/floor(10/(1-pnorm(sqrt(log(p)))))))))/(kk/floor(10/(1-pnorm(sqrt(log(p)))))*(p*(p-1)))-1)^2+Wsssp[ii];
      }
      Wfpp[l,ii]=(max(sum(sum(sqrt(diff)>=z[t])),1)-sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t])))/(sum(sum(OOO!=0))/2);
      }
      



  iihatp=min(which(sssp==min(sssp[3:20])));
  Aiihatp=min(which(Asssp==min(Asssp[3:20])));
  Biihatp=min(which(Bsssp==min(Bsssp[3:20])));
  Wiihatp=min(which(Wsssp==min(Wsssp[3:20])));
  iihatpO=min(which(ssspO==min(ssspO[3:20])));
  iihat=min(which(sss==min(sss[3:20])));
  fdpp[l]=ffp[l,iihatp]
  fpop[l]=fpp[l,iihatp]
  Afdpp[l]=Affp[l,Aiihatp]
  Afpop[l]=Afpp[l,Aiihatp]
  Bfdpp[l]=Bffp[l,Biihatp]
  Bfpop[l]=Bfpp[l,Biihatp]
  Wfdpp[l]=Wffp[l,Wiihatp]
  Wfpop[l]=Wfpp[l,Wiihatp]
  fdppO[l]=ffpO[l,iihatpO]
  fpopO[l]=fppO[l,iihatpO]
  fdp[l]=ff[l,iihat]
  fpo[l]=fp[l,iihat]
  
  print(l)
  print(fdpp[l])
  print(fpop[l])
  print(Bfdpp[l])
  print(Bfpop[l])
  print(fdppO[l])
  print(fpopO[l])
  print(Wfdpp[l])
  print(Wfpop[l])

}


sum(fdpp)/rep
sum(fpop)/rep
sum(Bfdpp)/rep
sum(Bfpop)/rep
sum(fdppO)/rep
sum(fpopO)/rep
sum(Wfdpp)/rep
sum(Wfpop)/rep

c(sd(fdpp),sd(Bfdpp),sd(fdppO),sd(Wfdpp))

