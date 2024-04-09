#----------------------------------------------------------------------------------------#
# R Code for groupwise dimension reduction                                               #
# Author: Lexin Li and Bing Li                                                           # 
# Reference: Li, L., Li, B., and Zhu, L.X. (2010). Groupwise Dimension Reduction.        #
#            Journal of the American Statistical Association. In press.                  #
# Last updated: 2010/02/09                                                               #
#----------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------
# supporting functions 
#-----------------------------------------------------------------------------------------

# normalize a vector
norm<-function(v)  
{ 
   sumv2<-sum(v^2)
   if(sumv2 == 0) sumv2<-1
   v/sqrt(sumv2)
}

# Gram-Schmidt orthonormalization
orthnorm<-function(X)
{
   X<-as.matrix(X)
   n<-nrow(X)
   p<-ncol(X)

   W<-NULL
   if(p > 1) {
      W<-cbind(W, X[,1])
      for(k in 2:p) {
         gw<-rep(0, n)
         for(i in 1:(k-1)) {
            gki<-as.vector((t(W[,i]) %*% X[,k])/(t(W[,i]) %*% W[,i]))
            gw<-gw + gki * W[,i]
         }
         W<-cbind(W, X[,k] - gw)
      }
   } else {
      W<-cbind(W, X[,1])
   }

   W<-apply(W, 2, norm)
   W
}

# angle between two spaces
angles<-function(B1, B2)
{
   if(!is.matrix(B1)) B1<-as.matrix(B1)
   if(!is.matrix(B2)) B2<-as.matrix(B2)

   if(ncol(B1) >= ncol(B2)) {
      B<-B1; B.hat<-B2
   } else {
      B<-B2; B.hat<-B1
   }

   P1<-B %*% solve(t(B) %*% B) %*% t(B)
   if(ncol(B.hat) == 1) {
      nume<-as.vector(t(B.hat) %*% P1 %*% B.hat)
      deno<-as.vector(t(B.hat) %*% B.hat)
      ratio<-nume / deno
   } else {
      BtB<-t(B.hat) %*% B.hat
      ei<-eigen(BtB)
      BtB2<-ei$vectors %*% diag(1/sqrt(ei$values)) %*% t(ei$vectors)
      M<-BtB2 %*% t(B.hat) %*% P1 %*% B.hat %*% BtB2
      ratio<-abs(eigen(M)$values[nrow(M)])
   }
   ans<-acos(sqrt(ratio))/pi * 180
   if(ans > 90) ans<-180 - ans
   return(ans)
}

# vector correlation and trace correlation between two spaces
eval.space<-function(A, B, orthnm=TRUE) 
{
   if(!is.matrix(A)) A<-as.matrix(A)
   if(!is.matrix(B)) B<-as.matrix(B)
   if(orthnm) { 
      A<-orthnorm(A)
      B<-orthnorm(B) 
   }

   mat<-t(B) %*% A %*% t(A) %*% B
   d<-eigen(mat)$values
   d<-(d+abs(d))/2
   q<-sqrt(prod(d))
   r<-sqrt(mean(d))
   ans<-list(q=q, r=r)
   return(ans)
}



#-----------------------------------------------------------------------------------------
# mave and groupwise mave 
#-----------------------------------------------------------------------------------------

VAL.LARGE<-1e+9
VAL.SMALL<-1e-9

std<-function(v) { (v - mean(v)) / sd(v) }
snorm<-function(v) { v / sum(v) }

wls<-function(x,y,w) { ginv(t(x*w)%*%(x))%*%(t(x*w)%*%y) }

# input:  U = an n x kd matrix, kd is the kernel dimension
#         ktype = type of kernel function to use
# output: an n x n weight matrix
comp.wts<-function(U, ktype="gaussian", h=NULL)
{
   # parameters
   n<-nrow(U)
   kd<-ncol(U)
   if(is.null(h)) { h<-(4/(kd+2))^(1/(kd+4)) * n^(-1/(kd+4)) }

   # standardize U
   Us<-apply(U, 2, std)
   
   # compute weight matrix
   W<-matrix(NA, n, n)
   for(i in 1:n) {
      for(j in i:n) {
         wij<-exp(-0.5*(t(Us[i,] - Us[j,])%*%(Us[i,] - Us[j,])/(h^2))) / (sqrt(2*pi))
         W[i,j]<-W[j,i]<-wij / (h^kd)
      }
   }
   W<-apply(W, 2, snorm)

   # return
   ans<-list(W=W, h=h)
   return(ans)
}


# mave with known d
mave<-function(X, y, d=1, ktype="gaussian", bdwd=NULL, max.iter=100, eps.conv=1e-6, nslices=10)
{
   # parameters
   n<-nrow(X)
   p<-ncol(X)

   # initialization
   beta.init<-as.matrix(dr(y~X, method="sir", nslices=nslices)$evectors[,1:d])
   beta<-beta.init
   G<-delta.G<-VAL.LARGE
   Gs<-NULL
   iter<-0

   # loop
   while((delta.G >= eps.conv) & (iter <= max.iter)) {
      # compute weights
      U<-X %*% beta
      out.wts<-comp.wts(U, ktype, h=bdwd) 
      W<-out.wts$W; bdwd<-out.wts$h

      # solve (a, b) given beta
      a<-b<-NULL
      for(j in 1:n) {
         y.star<-y
         x.star<-cbind(1, t(t(X) - X[j,]) %*% beta)
         cf<-as.vector(wls(x.star, y.star, W[,j]))        
         a<-c(a, cf[1])
         b<-rbind(b, cf[-1])
      } 

      # solve beta given (a, b)
      y.star2<-x.star2<-NULL
      for(j in 1:n) {
         y.star2<-c(y.star2, (y - a[j]))
         x.star2<-rbind(x.star2, t(b[j,]) %x% t(t(X) - X[j,]))
      }
      bt.new<-as.vector(wls(x.star2, y.star2, as.vector(W)))
      beta.new<-matrix(bt.new, nrow=p, byrow=FALSE)

      # compute objective
      rsd<-(y.star2 - x.star2 %*% bt.new)
      G.new<-sum(as.vector(W) * rsd^2) 
      delta.G<-G - G.new

      # update
      if(delta.G > 0) {
         beta<-beta.new
         G<-G.new
      }
      Gs<-c(Gs, G)
      iter<-iter + 1
   }
   
   # normalize
   b.est<-apply(beta, 2, norm)

   # return
   ans<-list(b.est=b.est, b.init=beta.init, Gs=Gs, iter=iter, rss=G, h=bdwd)
   return(ans)
}

# mave BIC criterion to estimate d 
mave.comp.d<-function(X, y, d.max=p-1, ktype="gaussian", bdwd=NULL, max.iter=100, eps.conv=1e-6, nslices=10)
{
   # parameters
   n<-nrow(X)
   p<-ncol(X)
   d.max<-min(d.max, p-1)

   # BIC criterion
   outs<-list()
   crit<-log(sum(y^2) / n)
   for(d in 1:d.max) {
      out.d<-mave(X, y, d=d, ktype=ktype, bdwd=bdwd, max.iter=max.iter, eps.conv=eps.conv, nslices=nslices)
      crit.d<-log(out.d$rss / n) + log(n) * d / (n * (out.d$h)^d)
      crit<-c(crit, crit.d)
      outs[[d]]<-out.d
   }
   d<-order(crit)[1] - 1

   # return
   ans<-list(d=d, crit=crit, out.d=outs[[d]])
   return(ans)  
}


# groupwise mave with known d's
gmave<-function(X, y, grp.id, ds=rep(1,ngrp), ktype="gaussian", bdwd=NULL, max.iter=100, eps.conv=1e-6, nslices=10)
{
   # parameters
   n<-length(y)
   ngrp<-length(grp.id)

   # initialization
   beta.init<-list()
   for(g in 1:ngrp) {
      beta.init[[g]]<-unname(as.matrix(dr(y~X[,grp.id[[g]]], method="sir", nslices=nslices)$evectors[,1:ds[g]]))
   }
   beta<-beta.init
   G<-delta.G<-VAL.LARGE
   Gs<-NULL
   iter<-0

   # loop
   while((delta.G >= eps.conv) & (iter <= max.iter)) {
      # compute weights
      U<-NULL
      for(g in 1:ngrp) { 
         U<-cbind(U, X[,grp.id[[g]]] %*% beta[[g]])
      }
      out.wts<-comp.wts(U, ktype, h=bdwd) 
      W<-out.wts$W; bdwd<-out.wts$h

      # solve (a, b) given beta
      a<-b<-NULL
      for(j in 1:n) {
         y.star<-y
         x.star<-rep(1, n)
         for(g in 1:ngrp) {
            x.star<-cbind(x.star, t(t(X[,grp.id[[g]]]) - X[j,grp.id[[g]]]) %*% beta[[g]])
         }
         cf<-as.vector(wls(x.star, y.star, W[,j]))        
         a<-c(a, cf[1])
         b<-rbind(b, cf[-1])
      } 
      bg<-list(); bst<-1
      for(g in 1:ngrp) {
         ben<-bst+ds[g]-1
         bg[[g]]<-as.matrix(b[, bst:ben])
         bst<-ben+1
      }

      # solve beta given (a, b)
      beta.new<-list()
      for(g in 1:ngrp) {
         b.g<-bg[[g]]
         y.star2<-x.star2<-NULL
         for(j in 1:n) {
            bbx<-rep(0, n) 
            for(gn in 1:ngrp) {
               if(gn != g) {
                  b.gn<-bg[[gn]]
                  bbx<-bbx+t(t(X[,grp.id[[gn]]]) - X[j,grp.id[[gn]]]) %*% beta[[gn]] %*% b.gn[j,]
               }
            }
            y.star2<-c(y.star2, (y - a[j] - bbx))
            x.star2<-rbind(x.star2, t(b.g[j,]) %x% t(t(X[,grp.id[[g]]]) - X[j,grp.id[[g]]]))
         }
         bt.new<-as.vector(wls(x.star2, y.star2, as.vector(W)))
         beta.new[[g]]<-matrix(bt.new, nrow=length(grp.id[[g]]), byrow=FALSE)

         rsd<-(y.star2 - x.star2 %*% bt.new)
      }
      G.new<-sum(as.vector(W) * rsd^2) 

      # compute objective
      delta.G<-G - G.new

      # update
      if(delta.G > 0) {
         beta<-beta.new
         G<-G.new
      }
      Gs<-c(Gs, G)
      iter<-iter + 1
   }
   
   # normalize
   b.est<-list()
   for(g in 1:ngrp) {
      b.est[[g]]<-apply(beta[[g]], 2, norm)
   }

   # return
   ans<-list(b.est=b.est, b.init=unname(beta.init), Gs=Gs, iter=iter, rss=G, h=bdwd)
   return(ans)
}

# groupwise mave BIC criterion to estimate d
gmave.comp.d<-function(X, y, grp.id, d.max=3, ktype="gaussian", bdwd=NULL, max.iter=100, eps.conv=1e-6, nslices=10)
{
   # parameters
   n<-nrow(X)
   ngrp<-length(grp.id)

   # generate d series
   dsa<-NULL
   for(i in 1:ngrp) {
      dsa<-cbind(dsa, rep(seq(0, d.max), each=(d.max+1)^(i-1), times=(d.max+1)^(ngrp-i)))
   }
   dsa<-dsa[,ngrp:1]

   # BIC criterion
   outs<-list()
   ys2<-sum((y - mean(y))^2)
   crit<-log(ys2/n)
   for(k in 1:(nrow(dsa)-1)) {
      dsw<-dsa[k+1,]
 
      pos.d<-seq(1, ngrp)[dsw != 0]
      grp.id.w<-grp.id[pos.d]
      ds.w<-dsw[pos.d]

      if(length(ds.w)==1) { 
         out.d<-mave(X[,grp.id.w[[1]]], y, d=ds.w, ktype=ktype, bdwd=bdwd, max.iter=max.iter, eps.conv=eps.conv, nslices=nslices) 
      } else {   
         out.d<-gmave(X, y, grp.id.w, ds=ds.w, ktype=ktype, bdwd=bdwd, max.iter=max.iter, eps.conv=eps.conv, nslices=nslices)
      }
      
      crit.d<-log(out.d$rss / n) + log(n) * sum(ds.w) / (n * (out.d$h)^sum(ds.w))
      crit<-c(crit, crit.d)

      outs[[k]]<-out.d
   }

   # d estimate
   d<-dsa[order(crit)[1],]

   # return
   ans<-list(d=d, crit=crit, dsa=dsa, outs=outs)
   return(ans)  
}



#-----------------------------------------------------------------------------------------
# simulations
#-----------------------------------------------------------------------------------------

# simulation model 1
gen.data1<-function(n, ps)
{
   ngrp<-length(ps)
   p<-sum(ps)

   b1<-c(1, -1, 0, rep(0, ps[1]-3))
   b2<-c(1, 1, 1,  rep(0, ps[2]-3))
   b3<-c(1, 0, 0,  rep(0, ps[3]-3))
   b4<-c(0, 1, 0,  rep(0, ps[4]-3))
   b.true<-list(b1, b2, b3, b4)

   X<-xb<-NULL
   for(g in 1:ngrp) { 
      X.g<-matrix(rnorm(n * ps[g]), nrow=n)
      xb.g<-as.vector(X.g %*% b.true[[g]])
      X<-cbind(X, X.g)
      xb<-cbind(xb, xb.g)
   }

   y<-xb[,1] + 0.1*(xb[,2]+2)^2 + exp(0.5*xb[,3]) + sin(0.2*pi*xb[,4]) + 0.5*rnorm(n) 

   ans<-list(X=X, y=y, b.true=b.true, xb=xb)
   return(ans)
}


# evaluation
eval.d<-function(b.true, b.est) 
{
   ngrp<-min(length(b.true), length(b.est))
   ang<-vcr<-NULL
   for(g in 1:ngrp) {
      ang<-c(ang, angles(b.true[[g]], b.est[[g]]))
      vcr<-c(vcr, eval.space(b.true[[g]], b.est[[g]])$q)
   }
   ans<-rbind(ang, vcr)
   return(ans)
}



#-----------------------------------------------------------------------------------------
# sample code
#-----------------------------------------------------------------------------------------

#library(dr)
#source("group.R")

# generate data
#set.seed(100)
#ps<-c(10, 10, 5, 5)
#data<-gen.data1(n=100, ps=ps)
#X<-data$X; y<-data$y; b.true<-data$b.true

# groupwise mave estimation
#grp.id<-list(seq(1,ps[1]), seq(ps[1]+1,sum(ps[1:2])), seq(sum(ps[1:2])+1,sum(ps[1:3])), seq(sum(ps[1:3])+1,sum(ps[1:4])))
#ds=rep(1, length(grp.id))
#outg<-gmave(X, y, grp.id, ds, ktype="gaussian", bdwd=NULL, max.iter=100, eps.conv=1e-6, nslices=5)

# evaluation of gRMAVE and aSIR
#eval.d(b.true, outg$b.est)
#eval.d(b.true, outg$b.init)

# a closer look at the groupwise mave output
#outg$iter
#outg$Gs; plot(outg$Gs)




























































