rm(list=ls(all=TRUE)) 
library(mnormt)
library(Rlab)
library(nlme)
library(CompQuadForm)
library(glmnet)
library(expm)
#library(lars)

n1=50;
rho1=0.5;
p=200;
q=30
rep=100;
alpha=0.01;
M=2;
ll=(floor(p/M))*(floor(p/M)-1)/2;

I = array(0:0, dim=c(p,p))
for (i in 1:p){
  I [i,i]=1;
}
mu = rep(0,p)



e1 = rep(1,n1*q)
ff = array(0:0, dim=c(rep,20))
fp = array(0:0, dim=c(rep,20))
fdp = rep(0,rep)
fpo = rep(0,rep)

ffp = array(0:0, dim=c(rep,20))
fpp = array(0:0, dim=c(rep,20))
fdpp = rep(0,rep)
fpop = rep(0,rep)

ffpO = array(0:0, dim=c(rep,20))
fppO = array(0:0, dim=c(rep,20))
fdppO = rep(0,rep)
fpopO = rep(0,rep)

for (l in 1:rep){

# creat Sigma1

Omega1 = I;
for (i in 1:(p-1)){
    Omega1[i,(i+1)]=0.6;
    Omega1[i+1,i]=0.6;
}
for (i in 1:(p-2)){
    Omega1[i,i+2]=0.3;
    Omega1[i+2,i]=0.3;
}

#Omega1 = I;
#for (i in 1:(p/2-1)){
#	for (j in (i+1):(p/2)){
#		Omega1[i,j]=rbinom(1,1,2/p)*0.6
#		Omega1[j,i]=Omega1[i,j]
#	}
#}

#eO1 = eigen(Omega1)
#de = abs(min(eO1$values))+0.05;
#Omega1 = (Omega1 + de*I)/(1+de);

#creat Omega2
Omega2 = I
for (k in 1:(p/10)){
    for (j in (10*(k-1)+2):(10*(k-1)+10)){
        Omega2[10*(k-1)+1,j]=0.5;
        Omega2[j,10*(k-1)+1]=Omega2[10*(k-1)+1,j];
    }
}
eO2 = eigen(Omega2)
de = abs(min(eO2$values))+0.05;
Omega2 = (Omega2 + de*I)/(1+de);


    

D = array(0:0, dim=c(p,p))
for (i in 1:p){
  D[i,i] = runif(1,1,3)
  D[i,i] = 1
#  D[i,i] = D[i,i]^(1/2)
  }


Omega1 = D%*%Omega1%*%D;
Sigma1=solve(Omega1);


Omega2 = D%*%Omega2%*%D;
Sigma2 = solve(Omega2);


vec = sample(1:p-1,size=p-1,replace=FALSE)
rd = vec[1:floor(p/M)];
rd = sort(rd);
rd =c(0,rd);

#Sigma3
Omega3=I;
for (i in 1:(p-1)){
	for (j in (i+1):p){
		Omega3[i,j] = rbinom(1,1,2/p)*0.8
		Omega3[j,i]=Omega3[i,j]
	}
}
eO3=eigen (Omega3)
de = abs (min (eO3$values))+0.05;
Omega3 = (Omega3 + de*I)/(1+de);
Omega3=D%*%Omega3%*%D;

Sigma3=solve(Omega3);


KK=2;

Sigma5=I
for (k in 1:(p/KK)){
	for (i in ((k-1)*KK+1):(k*KK-1)){
		for (j in (i+1):(k*KK)){
			Sigma5[i,j]=0.8
			Sigma5[j,i]=Sigma5[i,j]
		}
	}
}
eO5=eigen (Sigma5)
de = abs (min (eO5$values))+0.05;
Sigma5 = (Sigma5 + de*I)/(1+de)
Omega5=solve(Sigma5)
Omega5=D%*%Omega5%*%D;
Sigma5=solve(Omega5);


Sigma=Sigma1;
SQS = sqrtm(Sigma)
OOO=Omega1-I;
OOO2=OOO*OOO;

################### Sigma for time
PHI = array(0:0, dim=c(q,q))
for (ik in 1:q){
	for (jk in 1:q){
		PHI[ik,jk] = 0.4^{abs(ik-jk)}
	}
}

bx = array(0:0, dim=c(p-1,p))
bxO = array(0:0, dim=c(p-1,p))
by = array(0:0, dim=c(p-1,p))
c1 = array(0:0, dim=c(n1*q,p))
T1 = array(0:0, dim=c(p,p))
T = array(0:0, dim=c(p,p))
diff = array(0:0, dim=c(p,p))
S1 = array(0:0, dim=c((floor(p/M)),(floor(p/M))))
S = array(0:0, dim=c((floor(p/M)),(floor(p/M))))
OO = array(0:0, dim=c((floor(p/M)),(floor(p/M))))
sss=rep(0,20)
sssp=rep(0,20)
ssspO=rep(0,20)
X = array(0:0, dim=c(n1,p,q))
XX = array(0:0, dim=c(n1,p,q))
XXO = array(0:0, dim=c(n1,p,q))

#    x = rmnorm(n = n1, mean = mu, Sigma)

##########################################################################
##########################################################################
#data-driven
##########################################################################
##########################################################################

#separate kronecker product
for (nk in 1:n1){
	Z=rmnorm(n=p, mean=0, PHI)
	X[nk,,]=SQS%*%Z
}
PHI0=array(0:0, dim=c(q,q))
for (nk in 1:n1){
	PHI0=PHI0+t(X[nk,,]) %*% X[nk,,]
}
PHI0=PHI0/n1/p
for (nk in 1:n1){
	XX[nk,,]=X[nk,,] %*% sqrtm(solve(PHI0))
}
x = XX[,,1]
for(np in 2:q){
	x=rbind(x,XX[,,np])
}


#inverse regression

Shatmimix=cov(x);
sigmax=diag(Shatmimix);
for (ii in 3:20){
      for (i in 1:p){
        newx=t(x);
        newx=newx[-i,];
        
         glm2=glmnet(t(newx),x[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmax[i]*log(p)/(n1*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         bx[,i] = coef[-1]
       
            
        #calculate \xi
        c1[,i] = x[,i]-mean(x[,i])*e1- t(newx-(colMeans(t(newx)))%*%t(e1))%*%bx[,i];
    }
      R1=t(c1)%*%c1/n1/q;
      s1=colMeans(c1*c1);
      for (i in 1:(p-1)){
          for (j in (i+1):p){
              T1[i,j]=-(R1[i,j]+s1[i]*bx[i,j]+s1[j]*bx[j-1,i]);
              diff[i,j]=(T1[i,j]/(R1[i,i]*R1[j,j]))^2/(1/(R1[i,i]*R1[j,j])/n1/q*(1+bx[i,j]^2*R1[i,i]/R1[j,j]));
              
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
#oracle
##########################################################################
##########################################################################

for (nk in 1:n1){
	XXO[nk,,]=X[nk,,] %*% sqrtm(solve(PHI))
}
xO = XXO[,,1]
for(np in 2:q){
	xO=rbind(xO,XXO[,,np])
}


      ShatmimixO=cov(xO);
  sigmaxO=diag(ShatmimixO);
for (ii in 3:20){
      for (i in 1:p){
        newxO=t(xO);
        newxO=newxO[-i,];
        
         glm2=glmnet(t(newxO),xO[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmaxO[i]*log(p)/(n1*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         bxO[,i] = coef[-1]
       
            
        #calculate \xi
        c1[,i] = xO[,i]-mean(xO[,i])*e1- t(newxO-(colMeans(t(newxO)))%*%t(e1))%*%bxO[,i];
    }
      R1=t(c1)%*%c1/n1/q;
      s1=colMeans(c1*c1);
      for (i in 1:(p-1)){
          for (j in (i+1):p){
              T1[i,j]=-(R1[i,j]+s1[i]*bxO[i,j]+s1[j]*bxO[j-1,i]);
              diff[i,j]=(T1[i,j]/(R1[i,i]*R1[j,j]))^2/(1/(R1[i,i]*R1[j,j])/n1/q*(1+bxO[i,j]^2*R1[i,i]/R1[j,j]));
              
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
      

      
      ffpO[l,ii]=sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t]))/max(sum(sum(sqrt(diff)>=z[t])),1);
      ssspO[ii]=0;
      for (kk in 1:10){
          ssspO[ii]=((sum(sum(sqrt(diff)>=qnorm(1-kk/floor(10/(1-pnorm(sqrt(log(p)))))))))/(kk/floor(10/(1-pnorm(sqrt(log(p)))))*(p*(p-1)))-1)^2+ssspO[ii];
      }
      fppO[l,ii]=(max(sum(sum(sqrt(diff)>=z[t])),1)-sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t])))/(sum(sum(OOO!=0))/2);
      }



##########################################################################
##########################################################################
#compare with the case ignoring the kronecker product structure
##########################################################################
##########################################################################

y = X[,,1]
for(np in 2:q){
	y=rbind(y,X[,,np])
}

      Shatmimiy=cov(y);
  sigmay=diag(Shatmimiy);
for (ii in 3:20){
      for (i in 1:p){
        newy=t(y);
        newy=newy[-i,];
        
         glm2=glmnet(t(newy),y[,i],family="gaussian",lambda=2*ii/40*sqrt(sigmay[i]*log(p)/(n1*q)));
         coef = as.matrix(t(coef(glm2,mode="lambda")))
         by[,i] = coef[-1]
       
            
        #calculate \xi
        c1[,i] = y[,i]-mean(y[,i])*e1- t(newy-(colMeans(t(newy)))%*%t(e1))%*%by[,i];
    }
      R1=t(c1)%*%c1/n1/q;
      s1=colMeans(c1*c1);
      for (i in 1:(p-1)){
          for (j in (i+1):p){
              T1[i,j]=-(R1[i,j]+s1[i]*by[i,j]+s1[j]*by[j-1,i]);
              diff[i,j]=(T1[i,j]/(R1[i,i]*R1[j,j]))^2/(1/(R1[i,i]*R1[j,j])/n1/q*(1+by[i,j]^2*R1[i,i]/R1[j,j]));
              
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
      

      
      ff[l,ii]=sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t]))/max(sum(sum(sqrt(diff)>=z[t])),1);
      sss[ii]=0;
      for (kk in 1:10){
          sss[ii]=((sum(sum(sqrt(diff)>=qnorm(1-kk/floor(10/(1-pnorm(sqrt(log(p)))))))))/(kk/floor(10/(1-pnorm(sqrt(log(p)))))*(p*(p-1)))-1)^2+sss[ii];
      }
      fp[l,ii]=(max(sum(sum(sqrt(diff)>=z[t])),1)-sum(sum(sqrt(diff[which(OOO==0,arr.ind = T)])>=z[t])))/(sum(sum(OOO!=0))/2);

 
  }
  
  iihatp=min(which(sssp==min(sssp[3:20])));
  iihatpO=min(which(ssspO==min(ssspO[3:20])));
  iihat=min(which(sss==min(sss[3:20])));
  fdpp[l]=ffp[l,iihatp]
  fpop[l]=fpp[l,iihatp]
  fdppO[l]=ffpO[l,iihatpO]
  fpopO[l]=fppO[l,iihatpO]
  fdp[l]=ff[l,iihat]
  fpo[l]=fp[l,iihat]
  
  print(l)
  print(fdpp[l])
  print(fpop[l])
  print(fdppO[l])
  print(fpopO[l])
  print(fdp[l])
  print(fpo[l])

}


sum(fdpp)/rep
sum(fpop)/rep
sum(fdppO)/rep
sum(fpopO)/rep
sum(fdp)/rep
sum(fpo)/rep

