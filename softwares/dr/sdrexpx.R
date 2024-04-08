#---------------------------------------------------------------------------------------------# 
# R functions for dimnesion reduction with exponential predictors                             #
# Author: Lexin Li and R. Dennis Cook                                                         # 
# Reference: Cook, R.D., and Li, L. Dimension Reduction in Regressions with Exponential       #
#            Family Predictors. Journal of Computational and Graphical Statistics,            #
#            18, 774-791.                                                                     #
#---------------------------------------------------------------------------------------------#



#----------------------------------------------------------------------------------------------
# general supporting functions and libraries
#----------------------------------------------------------------------------------------------

library(Matrix)
   
VAL.LARGE<-1e+9
VAL.SMALL<-1e-9


# center a vector
center<-function(v)  v - mean(v)

# standardize a vector
stand<-function(v)  
{ 
   de<-sd(v) 
   if(de == 0) de<-1
   (v - mean(v)) / de
}

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

# power of a matrix
mat.power<-function(A, a)
{
   ei<-eigen(A)
   d<-ei$values
   d<-(d+abs(d))/2
   d2<-d^a
   d2[d == 0]<-0
   ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
   return(ans)
}

# covariance matrix 
cov.x<-function(X)
{
   Xc<-apply(X, 2, center)
   t(Xc) %*% Xc / nrow(Xc)
}

# covariance vector
cov.xy<-function(X, y)
{
   Xc<-apply(X, 2, center)
   yc<-center(y)
   t(Xc) %*% matrix(yc, ncol=1) /nrow(Xc) 
}

# angle (in degree) between two spaces
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



#----------------------------------------------------------------------------------------------
# supporting functions for exponential family predictors
#----------------------------------------------------------------------------------------------

# matrix exponential
mat.exp<-function(A)
{
   out<-expm(Matrix(A))
   matrix(attributes(out)$x, nrow=nrow(A), ncol=ncol(A)) 
}

# negative log likelihood for all binary predictors
obj.func.bin<-function(X, Gm, mu, nu) 
{ 
   loglik.j<-NULL
   for(j in 1:ncol(X)) {
      eta.j<-mu[j] + nu %*% Gm[j,]
      p.j<-binomial()$linkinv(eta.j)
      y.j<-X[,j]
      loglik.j<-c(loglik.j, sum(y.j * log(p.j) + (1 - y.j) * log(1 - p.j))) 
   }
   negloglik<-sum(loglik.j) * -1
 
   return(negloglik)  
}

# negative log likelihood for binary (the first q predictors) and 
# gaussian (the last p-q predictors) predictors
obj.func.mix<-function(X, Gm, mu, nu, q) 
{ 
   logliks<-NULL
   for(j in 1:ncol(X)) {
      eta.j<-mu[j] + nu %*% Gm[j,]
      y.j<-X[,j]

      if(j <= q) {
         mu.j<-binomial()$linkinv(eta.j)
         loglik.j<-sum(y.j * log(mu.j) + (1 - y.j) * log(1 - mu.j))
      }
      if(j > q) {              
         mu.j<-gaussian()$linkinv(eta.j)
         loglik.j<--0.5*sum((y.j - mu.j)^2 + log(2*pi)) 
      }

      logliks<-c(logliks, loglik.j) 
   }
   negloglik<-sum(logliks) * -1
   
   return(negloglik)  
}

# generate the n by r matrix Fy
gen.Fy<-function(y, term.y=c("poly", "inter", "slice", "spline"), degree.y=1) 
{
   # generate basis
   term.y<-match.arg(term.y)
   if(term.y == "poly") {
      Fy<-NULL
      for(i in 1:degree.y) 
         Fy<-cbind(Fy, y^i)
   }
   if(term.y == "inter") {
      Fy<-cbind(y, y[,1]*y[,2])
   }
   if(term.y == "slice") {
      library(dr)
      sy<-dr.slices(y, degree.y)
      Fy<-matrix(0, nrow=length(y), ncol=sy$nslices)
      for(i in 1:sy$nslices) 
         Fy[sy$slice.indicator==i, i]<-1
   }
   if(term.y == "spline") {
     order<-4
     inknots<-degree.y
     knots<-quantile(y,probs = seq(0, 1,length=(inknots+2)))
     Fy<-spline.des(c(rep(knots[1],order-1),knots,rep(knots[inknots+2],order-1)), y, ord=order)$design
   }

   # center Fy
   Fy.c<-apply(Fy, 2, center)

   # drop one level if y is categorical
   if((term.y == "slice") & (length(unique(y)) <= degree.y[1])) {
      Fy.c<-as.matrix(Fy.c[, -ncol(Fy.c)])
   } 

   # return
   return(Fy.c)
}



#----------------------------------------------------------------------------------------------
# basis estimation
#----------------------------------------------------------------------------------------------

# all binary predictors
bs.bin<-function(X, y, d=1, term.y="poly", degree.y=3, max.iter.full=50, eps.conv=1e-3, max.iter.gn=50, delta.init=1, epsilon=1e-3, Gm.init=NULL)
{
   # parameters 
   n<-nrow(X)
   p<-ncol(X)
   est.full<-list()

   # generate Fy
   Fy<-gen.Fy(y, term.y=term.y, degree.y=degree.y) 
   r<-ncol(Fy)
   
   # initialization
   if(is.null(Gm.init)) {
      Gm<-matrix(norm(rep(1, p)), ncol=1)
      d.w<-d   
      while(d.w > 1) {
         Gm<-orthnorm(cbind(Gm, rnorm(p*1)))
         d.w<-d.w-1
      }
   } else {
      Gm<-Gm.init
   }
   Qt<-orthnorm(cbind(Gm, diag(1,p)[,(d+1):p]))

   # iterative estimation
   iter.full<-0
   objvs.full<-NULL
   est.full<-list()
   crit.conv<-VAL.LARGE
   objv.full<--VAL.LARGE
   while((iter.full < max.iter.full) & (crit.conv > eps.conv)){
      # update mu and beta using GLM
      x.star<-NULL
      for(j in 1:p) {
         lj<-rep(0, p); lj[j]<-1
         Gmj<-matrix(Gm[j,], nrow=1)
         x.star.j<-cbind(matrix(rep(lj, n), nrow=n, byrow=TRUE), Fy %x% Gmj)
         x.star<-rbind(x.star, x.star.j)
      }
      y.star<-as.vector(X)
      fit.glm<-glm(y.star ~ x.star-1, family=binomial)
      mu<-unname(coef(fit.glm)[1:p])
      beta<-matrix(unname(coef(fit.glm)[-(1:p)]), nrow=d, ncol=r)
      nu<-Fy %*% t(beta)      

      # update Gamma with fixed max.iter.gn steps
      objv<-obj.func.bin(X, Gm, mu, nu)*-1
      delta<-delta.init
      objvs<-objv; deltas<-delta
      for(iter in 1:max.iter.gn) {
         # alpha
         alpha<-matrix(0, p, p)
         for(i in 1:d) {
            for(j in (d+1):p) {
               Gm.tilde<-Gm
               Gm.tilde[,i]<-cos(epsilon) * Qt[,i] - sin(epsilon) * Qt[,j]
               objv.tilde<-obj.func.bin(X, Gm.tilde, mu, nu)*-1 
               alpha[i,j]<- (objv.tilde - objv) / epsilon
            }
         }

         # A
         A<-matrix(0, p, p)
         for(i in 1:d) {
            for(j in (d+1):p) {
               Eij<-matrix(0, p, p)
               Eij[i, j]<-1
               Eij[j, i]<--1
               A<-A + alpha[i, j] * Eij
            }
         }

         # matrix exponential
         e2<-mat.exp(-1 * delta * A)
     
         # update Gm and Qt
         Qt.new<-Qt %*% t(e2)
         Gm.new<-as.matrix(Qt.new[, 1:d])

         # update objective value
         objv.new<-obj.func.bin(X, Gm.new, mu, nu)*-1
      
         # re-assign
         if(objv.new > objv) {
            Gm<-Gm.new
            Qt<-Qt.new
            objv<-objv.new
         } else {
            delta<-delta * -0.5
         }

         # re-assign delta
         if(abs(delta) < 1e-6) delta<-delta.init

         # record
         objvs<-c(objvs, objv); deltas<-c(deltas, delta)
      }

      # convergence criterion
      objv.full.new<-obj.func.bin(X, Gm, mu, nu)*-1 
      crit.conv<-objv.full.new - objv.full

      # update
      iter.full<-iter.full + 1
      objv.full<-objv.full.new

      # record
      objvs.full<-c(objvs.full, objv.full)
      est.full[[iter.full]]<-list(Gm=Gm, mu=mu, nu=nu)
   }

   # final estimate
   if(crit.conv < 0) {
      Gm<-est.full[[iter.full - 1]]$Gm
      mu<-est.full[[iter.full - 1]]$mu
      nu<-est.full[[iter.full - 1]]$nu
   } else {
      Gm<-est.full[[iter.full]]$Gm
      mu<-est.full[[iter.full]]$mu
      nu<-est.full[[iter.full]]$nu
   }

   # -2log-likelihood
   neglik2<-2 * obj.func.bin(X, Gm, mu, nu)

   # return
   ans<-list(Gamma.est=Gm, mu=mu, nu=nu, neglik2=neglik2, objvs=objvs.full, iter=iter.full)
   return(ans)
}

# binary (the first q predictors) and gaussian (the last p-q predictors) predictors
bs.mix<-function(X, y, d=1, q, term.y="poly", degree.y=3, max.iter.full=50, eps.conv=1e-3, max.iter.gn=50, delta.init=1, epsilon=1e-3, Gm.init=NULL)
{
   # parameters 
   n<-nrow(X)
   p<-ncol(X)
   est.full<-list()

   # generate Fy
   Fy<-gen.Fy(y, term.y=term.y, degree.y=degree.y) 
   r<-ncol(Fy)

   # initialization
   if(is.null(Gm.init)) {
      Gm<-matrix(norm(rep(1, p)), ncol=1)
      d.w<-d   
      while(d.w > 1) {
         Gm<-orthnorm(cbind(Gm, rnorm(p*1)))
         d.w<-d.w-1
      }   
   } else {
      Gm<-Gm.init
   }
   Qt<-orthnorm(cbind(Gm, diag(1,p)[,(d+1):p]))

   # iterative estimation
   iter.full<-0
   objvs.full<-NULL
   est.full<-list()
   crit.conv<-VAL.LARGE
   objv.full<--VAL.LARGE
   while((iter.full < max.iter.full) & (crit.conv > eps.conv)){
      # update mu and beta using GLM
      x.star<-NULL
      for(j in 1:p) {
         lj<-rep(0, p); lj[j]<-1
         Gmj<-matrix(Gm[j,], nrow=1)
         x.star.j<-cbind(matrix(rep(lj, n), nrow=n, byrow=TRUE), Fy %x% Gmj)
         x.star<-rbind(x.star, x.star.j)
      }
      y.star<-as.vector(X)
      fit.glm<-glm(y.star[-seq(1,(q*n))] ~ x.star[-seq(1,(q*n)),]-1, family=gaussian)
      beta<-matrix(unname(coef(fit.glm)[-(1:p)]), nrow=d, ncol=r)
      nu<-Fy %*% t(beta)      
      mu<-unname(coef(fit.glm)[1:p])
      for(j in 1:q) {
         fit.glm<-glm(y.star[((j-1)*n+1):(j*n)] ~ x.star[((j-1)*n+1):(j*n),]-1, family=binomial)
         mu[j]<-unname(coef(fit.glm)[j])
      }

      # update Gamma with fixed max.iter.gn steps
      objv<-obj.func.mix(X, Gm, mu, nu, q=q)*-1
      delta<-delta.init
      objvs<-objv; deltas<-delta
      for(iter in 1:max.iter.gn) {
         # alpha
         alpha<-matrix(0, p, p)
         for(i in 1:d) {
            for(j in (d+1):p) {
               Gm.tilde<-Gm
               Gm.tilde[,i]<-cos(epsilon) * Qt[,i] - sin(epsilon) * Qt[,j]
               objv.tilde<-obj.func.mix(X, Gm.tilde, mu, nu, q=q)*-1 
               alpha[i,j]<- (objv.tilde - objv) / epsilon
            }
         }

         # A
         A<-matrix(0, p, p)
         for(i in 1:d) {
            for(j in (d+1):p) {
               Eij<-matrix(0, p, p)
               Eij[i, j]<-1
               Eij[j, i]<--1
               A<-A + alpha[i, j] * Eij
            }
         }

         # matrix exponential
         e2<-mat.exp(-1 * delta * A)
     
         # update Gm and Qt
         Qt.new<-Qt %*% t(e2)
         Gm.new<-as.matrix(Qt.new[, 1:d])

         # update objective value
         objv.new<-obj.func.mix(X, Gm.new, mu, nu, q=q)*-1
      
         # re-assign
         if(objv.new > objv) {
            Gm<-Gm.new
            Qt<-Qt.new
            objv<-objv.new
         } else {
            delta<-delta * -0.5
         }

         # re-assign delta
         if(abs(delta) < 1e-6) delta<-delta.init

         # record
         objvs<-c(objvs, objv); deltas<-c(deltas, delta)
      }

      # convergence criterion
      objv.full.new<-obj.func.mix(X, Gm, mu, nu, q=q)*-1 
      crit.conv<-objv.full.new - objv.full

      # update
      iter.full<-iter.full + 1
      objv.full<-objv.full.new

      # record
      objvs.full<-c(objvs.full, objv.full)
      est.full[[iter.full]]<-list(Gm=Gm, mu=mu, nu=nu)
   }

   # final estimate
   if(crit.conv < 0) {
      Gm<-est.full[[iter.full - 1]]$Gm
      mu<-est.full[[iter.full - 1]]$mu
      nu<-est.full[[iter.full - 1]]$nu
   } else {
      Gm<-est.full[[iter.full]]$Gm
      mu<-est.full[[iter.full]]$mu
      nu<-est.full[[iter.full]]$nu
   }

   # -2log-likelihood
   neglik2<-2 * obj.func.mix(X, Gm, mu, nu, q=q)

   # return
   ans<-list(Gamma.est=Gm, mu=mu, nu=nu, neglik2=neglik2, objvs=objvs.full, iter=iter.full)
   return(ans)
}



#----------------------------------------------------------------------------------------------
# prediction
#----------------------------------------------------------------------------------------------

# conditional mean predicton for all binary predictors
y.pred.bin<-function(x.new, y, Gm, mu, nu)
{
   n<-length(y)
   p<-length(x.new)

   w.ij<-NULL
   for(j in 1:p) {
      eta.j<-mu[j] + nu %*% Gm[j,]
      a.j<-1 - exp(eta.j) / (1 + exp(eta.j))
      w.ij<-cbind(w.ij, a.j * exp(x.new[j] * eta.j))
   }
   w.j<-apply(w.ij, 1, prod)
   w<-w.j / sum(w.j)
   y.hat<-sum(y * w)

   ans<-y.hat
   return(ans)
}

# conditional mean predicton for binary (the first q predictors) and 
# gaussian (the last p-q predictors) predictors
y.pred.mix<-function(x.new, y, Gm, mu, nu, q)
{
   n<-length(y)
   p<-length(x.new)

   w.ij<-NULL
   for(j in 1:p) {
      eta.j<-mu[j] + nu %*% Gm[j,]
 
      if(j <= q) {
         a.j<-1 - exp(eta.j) / (1 + exp(eta.j))
      }
      if(j > q) {
         a.j<-exp(-0.5*eta.j^2) / sqrt(2*pi)
      }

      w.ij<-cbind(w.ij, a.j * exp(x.new[j] * eta.j))
   }
   w.j<-apply(w.ij, 1, prod)
   w<-w.j / sum(w.j)
   y.hat<-sum(y * w)

   ans<-y.hat
   return(ans)
}



#----------------------------------------------------------------------------------------------
# dimension estimation
#----------------------------------------------------------------------------------------------

dm.bin<-function(X, y, d.max=3, term.y="poly", degree.y=3, max.iter.full=50, eps.conv=1e-3, Gm.init.s=NULL) 
{
   # parameters 
   n<-nrow(X)
   p<-ncol(X)
   Fy<-gen.Fy(y, term.y=term.y, degree.y=degree.y) 
   r<-ncol(Fy)
   outs<-list()

   # d=0
   mu<-NULL
   for(j in 1:p) { 
      mu<-c(mu, unname(coef(glm(X[,j] ~ 1, family=binomial))))
   }
   Gm<-matrix(norm(rep(1, p)), ncol=1)
   nu<-matrix(0, nrow=n, ncol=1)
   neglik2<-2 * obj.func.bin(X, Gm, mu, nu) 
   p.e<-p
   crit.lik<-neglik2
   crit.aic<-neglik2 + 2 * p.e
   crit.bic<-neglik2 + log(n) * p.e

   # d = 1 to min(p, r) 
   d.max<-min(d.max, p, r)
   for(d in 1:d.max) {
      if(is.null(Gm.init.s)) {
         if(d == 1) {
            Gm.init<-matrix(norm(rep(1, p)), ncol=1)
         } else {
            Gm.init<-orthnorm(cbind(outs[[d - 1]]$Gamma.est, norm(rep(1, p))))
         }
      } else {
         Gm.init<-Gm.init.s[[d]]
      }

      out<-bs.bin(X=X, y=y, d=d, term.y=term.y, degree.y=degree.y, max.iter.full=max.iter.full, eps.conv=eps.conv, Gm.init=Gm.init)

      neglik2<-out$neglik2
      p.e<-p + d * (p - d) + d * r

      crit.lik<-c(crit.lik, neglik2)
      crit.aic<-c(crit.aic, neglik2 + 2 * p.e)
      crit.bic<-c(crit.bic, neglik2 + log(n) * p.e)

      outs[[d]]<-out
   }

   # full model
   if(d < min(p, r)) {
      mu<-Gmb<-NULL 
      for(j in 1:p) {
         x.star<-Fy; y.star<-X[,j]
         fit.glm<-glm(y.star ~ x.star, family=binomial)
         mu<-c(mu, unname(coef(fit.glm)[1]))
         Gmb<-rbind(Gmb, unname(coef(fit.glm)[-1]))
      }
      loglik.j<-NULL
      for(j in 1:ncol(X)) {
         eta.j<-mu[j] + Fy %*% Gmb[j,]
         p.j<-binomial()$linkinv(eta.j)
         y.j<-X[,j]
         loglik.j<-c(loglik.j, sum(y.j * log(p.j) + (1 - y.j) * log(1 - p.j))) 
      }
      neglik2.full<-2 * sum(loglik.j) * -1
   } else {
      neglik2.full<-neglik2
   }

   # estimate of d based on ic
   d.aic<-order(crit.aic)[1] - 1
   d.bic<-order(crit.bic)[1] - 1

   # hypothesis test
   st<-df<-pv<-rep(0, d.max)
   for(m in 0:(d.max-1)) {
      st[m+1]<-crit.lik[m+1] - neglik2.full
      df[m+1]<-(p - m) * (r - m)
      pv[m+1]<-1 - pchisq(st[m+1], df[m+1])
   }
   asy.test<-data.frame(cbind(st, df, pv))
   rr<-paste(0:(d.max-1),"D vs >= ",1:d.max,"D",sep="")
   dimnames(asy.test)<-list(rr,c("STAT", "DF", "PVAL")) 

   # return
   ans<-list(d.aic=d.aic, d.bic=d.bic, crit.lik=crit.lik, crit.aic=crit.aic, crit.bic=crit.bic, asy.test=asy.test) 
   return(ans)   
}



#----------------------------------------------------------------------------------------------
# predictor selection
#----------------------------------------------------------------------------------------------

# constrained basis estimation
bs.bin.H<-function(X, y, d=1, H0, Gm.star.init, term.y="poly", degree.y=3, max.iter.full=50, eps.conv=1e-3, max.iter.gn=50, delta.init=1, epsilon=1e-3)
{
   # parameters 
   n<-nrow(X)
   p<-ncol(X)
   p0<-ncol(H0)

   # generate Fy
   Fy<-gen.Fy(y, term.y=term.y, degree.y=degree.y) 
   r<-ncol(Fy)

   # initialization
   Gm.star<-Gm.star.init
   Qt<-orthnorm(cbind(Gm.star, diag(1,p0)[,(d+1):p0]))

   # iterative estimation
   iter.full<-0
   objvs.full<-NULL
   est.full<-list()
   crit.conv<-VAL.LARGE
   objv.full<--VAL.LARGE
   while((iter.full < max.iter.full) & (crit.conv > eps.conv)){
      # update mu and beta using GLM
      x.star<-NULL
      for(j in 1:p) {
         lj<-rep(0, p); lj[j]<-1
         Gmj<-matrix((H0 %*% Gm.star)[j,], nrow=1)
         x.star.j<-cbind(matrix(rep(lj, n), nrow=n, byrow=TRUE), Fy %x% Gmj)
         x.star<-rbind(x.star, x.star.j)
      }      
      y.star<-as.vector(X)
      fit.glm<-glm(y.star ~ x.star-1, family=binomial)
      mu<-unname(coef(fit.glm)[1:p])
      beta<-matrix(unname(coef(fit.glm)[-(1:p)]), nrow=d, ncol=r)
      nu<-Fy %*% t(beta)      

      # update Gamma with fixed max.iter.gn steps
      objv<-obj.func.bin(X, H0 %*% Gm.star, mu, nu)*-1
      delta<-delta.init
      objvs<-objv; deltas<-delta
      for(iter in 1:max.iter.gn) {
         # alpha
         alpha<-matrix(0, p0, p0)
         for(i in 1:d) {
            for(j in (d+1):p0) {
               Gm.star.tilde<-Gm.star
               Gm.star.tilde[,i]<-cos(epsilon) * Qt[,i] - sin(epsilon) * Qt[,j]
               objv.tilde<-obj.func.bin(X, H0 %*% Gm.star.tilde, mu, nu)*-1 
               alpha[i,j]<- (objv.tilde - objv) / epsilon
            }
         }

         # A
         A<-matrix(0, p0, p0)
         for(i in 1:d) {
            for(j in (d+1):p0) {
               Eij<-matrix(0, p0, p0)
               Eij[i, j]<-1
               Eij[j, i]<--1
               A<-A + alpha[i, j] * Eij
            }
         }

         # matrix exponential
         e2<-mat.exp(-1 * delta * A)
     
         # update Gm.star and Qt
         Qt.new<-Qt %*% t(e2)
         Gm.star.new<-as.matrix(Qt.new[, 1:d])

         # update objective value
         objv.new<-obj.func.bin(X, H0 %*% Gm.star.new, mu, nu)*-1
      
         # re-assign
         if(objv.new > objv) {
            Gm.star<-Gm.star.new
            Qt<-Qt.new
            objv<-objv.new
         } else {
            delta<-delta * -0.5
         }

         # re-assign delta
         if(abs(delta) < 1e-6) delta<-delta.init

         # record
         objvs<-c(objvs, objv); deltas<-c(deltas, delta)
      }

      # convergence criterion
      objv.full.new<-obj.func.bin(X, H0 %*% Gm.star, mu, nu)*-1 
      crit.conv<-objv.full.new - objv.full

      # update
      iter.full<-iter.full + 1
      objv.full<-objv.full.new

      # record
      objvs.full<-c(objvs.full, objv.full)
      est.full[[iter.full]]<-list(Gm=(H0 %*% Gm.star), mu=mu, nu=nu)
   }

   # final estimate
   if(crit.conv < 0) {
      Gm<-est.full[[iter.full - 1]]$Gm
      mu<-est.full[[iter.full - 1]]$mu
      nu<-est.full[[iter.full - 1]]$nu
   } else {
      Gm<-est.full[[iter.full]]$Gm
      mu<-est.full[[iter.full]]$mu
      nu<-est.full[[iter.full]]$nu
   }

   # -2log-likelihood
   neglik2<-2 * obj.func.bin(X, Gm, mu, nu)

   # return
   ans<-list(Gamma.est=Gm, neglik2=neglik2, objvs.full=objvs.full, iter.full=iter.full)
   return(ans)
}

vs.bin<-function(X, y, d=1, term.y="poly", degree.y=3, max.iter.full=50, eps.conv=1e-3, Gm.init=NULL, ord=NULL)
{
   # parameters 
   n<-nrow(X)
   p<-ncol(X)
   Fy<-gen.Fy(y, term.y=term.y, degree.y=degree.y) 
   r<-ncol(Fy)
   outs<-list()

   # predictor pre-ordering
   if(is.null(ord)) {
      ord<-seq(p, 1, by=-1)
   }
   X1<-X[, ord]
   if(!is.null(Gm.init)) { 
      Gm.init<-as.matrix(Gm.init[ord, ])
   } 

   # full model 
   out0<-bs.bin(X=X1, y=y, d=d, term.y=term.y, degree.y=degree.y, max.iter.full=max.iter.full, eps.conv=eps.conv, Gm.init=Gm.init)
   neglik2<-out0$neglik2
   p.e<-p + d * (p - d) + d * r
   crit.lik<-neglik2
   crit.aic<-neglik2 + 2 * p.e
   crit.bic<-neglik2 + log(n) * p.e
   outs[[1]]<-out0

   # backward selection
   for(pos in 1:(p-d-1)) {
      H0<-as.matrix(diag(1, p)[,-seq(1,pos)])
      Gm.star.init<-t(H0) %*% out0$Gamma.est
      Gm.star.init<-apply(Gm.star.init, 2, norm)
      out.H<-bs.bin.H(X=X1, y=y, d=d, H0=H0, Gm.star.init=Gm.star.init, term.y=term.y, degree.y=degree.y, max.iter.full=max.iter.full, eps.conv=eps.conv)

      neglik2<-out.H$neglik2
      p.e<-p + d * (p - pos - d) + d * r

      crit.lik<-c(crit.lik, neglik2)
      crit.aic<-c(crit.aic, neglik2 + 2 * p.e)
      crit.bic<-c(crit.bic, neglik2 + log(n) * p.e)

      outs[[pos+1]]<-out.H
   }

   # null model
   Gm.star<-diag(1, d, d)
   H0<-as.matrix(diag(1, p)[,-seq(1,p-d)]) 
   x.star<-NULL
   for(j in 1:p) {
      lj<-rep(0, p); lj[j]<-1
      for(i in 1:n) {
         x.star<-rbind(x.star, c(lj, Fy[i,] %x% (H0 %*% Gm.star)[j,]))
      }
   }
   y.star<-as.vector(X1)
   fit.glm<-glm(y.star ~ x.star-1, family=binomial)
   mu<-unname(coef(fit.glm)[1:p])
   beta<-matrix(unname(coef(fit.glm)[-(1:p)]), nrow=d, ncol=r)
   nu<-Fy %*% t(beta)      
   neglik2<-2 * obj.func.bin(X1, (H0 %*% Gm.star), mu, nu) 
   p.e<-p + d * r
   crit.lik<-c(crit.lik, neglik2)
   crit.aic<-c(crit.aic, neglik2 + 2 * p.e)
   crit.bic<-c(crit.bic, neglik2 + log(n) * p.e)
   outs[[p-d+1]]<-list(Gamma.est=(H0 %*% Gm.star))

   # selection
   pos.aic<-order(crit.aic)[1] - 1
   if(pos.aic == 0) { id.aic<-ord } else { id.aic<-ord[-seq(1, pos.aic)] }
   pos.bic<-order(crit.bic)[1] - 1
   if(pos.bic == 0) { id.bic<-ord } else { id.bic<-ord[-seq(1, pos.bic)] }

   # estimate after selection
   Gm.est.aic<-(outs[[pos.aic + 1]]$Gamma.est)
   Gm.est.bic<-(outs[[pos.bic + 1]]$Gamma.est)
   Gm.est<-(outs[[1]]$Gamma.est)

   # back to original order
   Gm.est<-Gm.est[order(ord),]
   Gm.est.aic<-Gm.est.aic[order(ord),]
   Gm.est.bic<-Gm.est.bic[order(ord),]
   
   # return
   ans<-list(Gamma.est=Gm.est, Gamma.est.aic=Gm.est.aic, Gamma.est.bic=Gm.est.bic, 
             id.aic=id.aic, id.bic=id.bic, ord=ord, 
             crit.lik=crit.lik, crit.aic=crit.aic, crit.bic=crit.bic) 
   return(ans)   
}



#----------------------------------------------------------------------------------------------
# functions for simulations
#----------------------------------------------------------------------------------------------

# generate a simulated data with all binary predictors and d=1
gen.data.b1<-function(n, p, p1, sigma.y) 
{
   y<-rnorm(n) * sigma.y
   y<-y - mean(y)    
   b1<-c(rep(1, p1), rep(0, p-p1))
   b1<-norm(b1)
   Gamma<-matrix(b1, ncol=1)   
   X<-NULL
   for(i in 1:n) {
      X.i<-NULL
      for(j in 1:p) { 
         eta.ij<-as.vector(Gamma[j, ] * y[i])
         p.ij<-1 / (1 + exp(-eta.ij))
         X.i<-c(X.i, rbinom(1, 1, p.ij))
      }
      X<-rbind(X, X.i)
   }

   ans<-list(X=unname(X), y=y, b.true=Gamma)
   return(ans)
}

# generate a simulated data with binary and gaussian predictors and d=2
gen.data.m2<-function(n, p, p1, q=5, sigma.y) 
{
   y<-rnorm(n) * sigma.y
   fy<-cbind(0.5*y, 0.1*y^2) 
   fy<-apply(fy, 2, center) 
   b1<-c(rep(1, p1), rep(0, p-p1)); b2<-c(rep(0, p-p1), rep(1, p1)) 
   Gamma<-cbind(b1, b2)
   Gamma<-apply(Gamma, 2, norm)
   X<-NULL
   for(i in 1:n) {
      X.i<-NULL
      for(j in 1:q) { 
         eta.ij<-as.vector(t(Gamma[j, ]) %*% fy[i,])
         p.ij<-1 / (1 + exp(-eta.ij))
         X.i<-c(X.i, rbinom(1, size=1, prob=p.ij))
      }
      for(j in (q+1):p) {
         eta.ij<-as.vector(t(Gamma[j, ]) %*% fy[i,])
         sigma<-1
         X.i<-c(X.i, rnorm(1, mean=eta.ij, sd=sigma))
      }
      X<-rbind(X, X.i)
   }

   ans<-list(X=unname(X), y=y, b.true=Gamma, fy=fy)
   return(ans)
}



#----------------------------------------------------------------------------------------------
# code for demonstration
#----------------------------------------------------------------------------------------------

run.ex1<-function()
{
   set.seed(100)
   n<-200; p<-20; p1<-10; sigma.y<-5
   data<-gen.data.b1(n=n, p=p, p1=p1, sigma.y=sigma.y) 
   X<-data$X; y<-data$y; Gamma.true<-data$b.true

   # basis estimation
   out.pfc<-bs.bin(X, y, d=1)
   angles(out.pfc$Gamma.est, Gamma.true)              
   
   # prediction
   y.hat<-apply(X, 1, y.pred.bin, y=y, Gm=out.pfc$Gamma.est, mu=out.pfc$mu, nu=out.pfc$nu)
   mean((y-y.hat)^2) / (sigma.y^2)
      
   # dimension estimation
   # this procedure will take a while
   out.dm<-dm.bin(X, y, d.max=3)
   
   # predictor selection   
   out.vs<-vs.bin(X, y, d=1, ord=NULL)                 
}

run.ex2<-function()
{
   set.seed(100)
   n<-200; p<-20; p1<-10; q<-5; sigma.y<-5
   data<-gen.data.m2(n=n, p=p, p1=p1, q=q, sigma.y=sigma.y) 
   X<-data$X; y<-data$y; Gamma.true<-data$b.true

   # basis estimation
   out.pfc<-bs.mix(X, y, d=2, q=5)
   angles(out.pfc$Gamma.est, Gamma.true)              
   
   # prediction
   y.hat<-apply(X, 1, y.pred.mix, y=y, Gm=out.pfc$Gamma.est, mu=out.pfc$mu, nu=out.pfc$nu, q=5)
   mean((y-y.hat)^2) / (sigma.y^2)               
}



















