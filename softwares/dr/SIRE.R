#---------------------------------------------------------------------------------------------# 
# R functions for shrinkage inverse regression estimation                                     #
# Author: Lexin Li                                                                            # 
# Reference: Bondell, H.D., and Li, L. (2009). Shrinkage Inverse Regression Estimation        #
#            for Model Free Variable Selection. Journal of the Royal Statistical Society,     #
#            Series B. 71, 287-299.                                                           #
#---------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------
# supporting functions
#----------------------------------------------------------------------------------------------

# center a vector
center<-function(v)  v - mean(v)

# normalize a vector
norm<-function(v)  
{ 
   sumv2<-sum(v^2)
   if(sumv2 == 0) sumv2<-1
   v/sqrt(sumv2)
}

# covariance matrix 
cov.x<-function(X)
{
   Xc<-apply(X, 2, center)
   t(Xc) %*% Xc / nrow(Xc)
}

# square-root of a matrix
mat.sqrt<-function(A)
{
   ei<-eigen(A)
   d<-ei$values
   d<-(d+abs(d))/2
   d2<-sqrt(d)
   ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
   return(ans)
}

# square-root-inverse of a matrix
mat.sqrt.inv<-function(A)
{
   ei<-eigen(A)
   d<-ei$values
   d<-(d+abs(d))/2
   d2<-1 / sqrt(d)
   d2[d == 0]<-0
   ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
   return(ans)
}

# garrote with a given tuning parameter of tau
garrote.tau<-function(U, V, tau, PRECISION=1e-10)
{ 
   p<-ncol(V)
   Dmat<-t(V) %*% V
   dvec<-as.vector(t(V) %*% U)
   Amat<-cbind(diag(1, p), rep(-1, p))
   bvec<-c(rep(0, p), -tau)
   out<-solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)

   omega.hat<-out$solution
   omega.hat[omega.hat < PRECISION]<-0

   ans<-list(omega.hat=omega.hat, obj.val=out$value)
   return(ans)
}

# garrote with a sequence of tuning parameters of tau 
garrote<-function(U, V, tau.range)
{
   n.e<-nrow(V)
   p.e<-ncol(V)
   h<-as.integer(n.e / p.e)
   sigma2<-sum(resid(lm(U ~ V - 1))^2) / (n.e - p.e)

   ests<-NULL
   crit<-objv<-dfes<-NULL
   for(k in 1:length(tau.range)) {
      tau<-tau.range[k]
      out.tau<-garrote.tau(U, V, tau)

      omega.hat<-out.tau$omega.hat
      ests<-cbind(ests, omega.hat)

      objv.k<-sum((U - V%*%omega.hat)^2)
      df.e<-2 * sum(omega.hat > 0) + sum(omega.hat * (h - 2))
      crit.k<-objv.k / sigma2 + log(n.e) * df.e
      crit<-c(crit, crit.k)
      objv<-c(objv, objv.k)
      dfes<-c(dfes, df.e)
   }
   omega.hat<-ests[, order(crit)[1]]
   colnames(ests)<-as.character(tau.range)
   names(crit)<-names(objv)<-names(dfes)<-as.character(tau.range)

   ans<-list(omega.hat=omega.hat, tau=tau.range[order(crit)[1]], tau.range=tau.range, crit=crit, objv=objv, dfes=dfes, sigma2=sigma2, ests=ests)
   return(ans)
}

# alternating least squares optimization
also<-function(zeta, V, B.init, max.iter=1000, eps=1e-6)
{
   # parameters
   p<-nrow(B.init)
   d<-ncol(B.init)
   h1<-as.integer(ncol(V) / nrow(zeta))

   # initialization
   B<-B.init
   vec.zeta<-vec.mat(zeta)
   vec.C<-solve((diag(1,h1) %x% t(B)) %*% V %*% (diag(1,h1) %x% B)) %*% (diag(1,h1) %x% t(B)) %*% V %*% vec.zeta
   C<-matrix(vec.C, nrow=d, ncol=h1)
   comp<-as.vector(zeta) - as.vector(B %*% C)
   F<-t(comp) %*% V %*% comp
   Fvals<-F
   delta.F<-F
   iter<-0

   # loop
   while((delta.F >= eps) & (iter <= max.iter)) {
      B.old<-B; C.old<-C; F.old<-F

      if(d == 1) {
         alpha.k<-vec.mat(zeta)
         ck<-matrix(C, ncol=1)
         bk<-ginv((t(ck) %x% diag(1,p)) %*% V %*% (ck %x% diag(1,p))) %*% (t(ck) %x% diag(1,p)) %*% V %*% alpha.k   # Moore-Penrose generalized inverse
         bk<-norm(bk)
         B.new<-matrix(bk, ncol=1)
         B<-B.new

         vec.C<-solve((diag(1,h1) %x% t(B)) %*% V %*% (diag(1,h1) %x% B)) %*% (diag(1,h1) %x% t(B)) %*% V %*% vec.zeta
         C<-matrix(vec.C, nrow=d, ncol=h1)
      } else {
         for(k in 1:d) {
            B.k<-B[, -k]; if(d == 2) B.k<-matrix(B.k, ncol=1)
            C.k<-C[-k, ]; if(d == 2) C.k<-matrix(C.k, nrow=1)
            alpha.k<-vec.mat(zeta - B.k %*% C.k)
            Q.B.k<-Q.mat(B.k)
            ck<-matrix(C[k, ], ncol=1)
            bk<-Q.B.k %*% ginv(Q.B.k %*% (t(ck) %x% diag(1,p)) %*% V %*% (ck %x% diag(1,p)) %*% Q.B.k) %*% Q.B.k %*% (t(ck) %x% diag(1,p)) %*% V %*% alpha.k
            bk<-norm(bk)
            B.new<-B; B.new[, k]<-bk
            B<-B.new
   
            vec.C<-solve((diag(1,h1) %x% t(B)) %*% V %*% (diag(1,h1) %x% B)) %*% (diag(1,h1) %x% t(B)) %*% V %*% vec.zeta
            C<-matrix(vec.C, nrow=d, ncol=h1)
         }
      }

      comp<-as.vector(zeta) - as.vector(B %*% C)
      F.new<-t(comp) %*% V %*% comp
      delta.F<-F - F.new
      F<-F.new
      Fvals<-c(Fvals, F)

      iter<-iter + 1
   }

   # disgard the last iteration if F goes up
   if(delta.F < 0) {
      B<-B.old; C<-C.old; F<-F.old
   }

   ans<-list(B.est=B, C.est=C, F=F, iter=iter, Fvals=Fvals)
   return(ans)
}



#------------------------------------------------------------------------------------------------------
# shrinkage cire
#------------------------------------------------------------------------------------------------------

# cire components
delta.cire<-function(X, y, sy)
{
   h<-sy$nslices
   YJ<-NULL
   for(s in 1:h) {
      y.s<-y
      y.s[sy$slice.indicator != s]<-0
      YJ<-cbind(YJ, y.s) 
   }
   YJ.c<-apply(YJ, 2, center)
   beta.hat<-beta.cire(X, y, sy)
   X.c<-apply(X, 2, center)
   delta<-YJ.c - X.c %*% beta.hat
   return(delta)
}
beta.cire<-function(X, y, sy)
{
   h<-sy$nslices
   Sigma.x.inv<-solve(cov.x(X))

   beta<-NULL
   for(s in 1:h) {
      y.s<-y
      y.s[sy$slice.indicator != s]<-0
      beta.s<-Sigma.x.inv %*% cov.xy(X, y.s)
      beta<-cbind(beta, beta.s)
   }
   return(beta)
}
Gamma.cire<-function(X, y, sy)
{
   n<-nrow(X)
   p<-ncol(X)
   h<-sy$nslices
   Sigma.x.inv<-solve(cov.x(X))
   delta.hat<-delta.cire(X, y, sy)
   X.c<-apply(X, 2, center)

   Gamma.n<-matrix(0, p*h, p*h)
   for(i in 1:n) {
      Gamma.n<-Gamma.n + (delta.hat[i,] %*% t(delta.hat[i,])) %x% (Sigma.x.inv %*% X.c[i,] %*% t(X.c[i,]) %*% Sigma.x.inv)      
   }
   Gamma.n<-Gamma.n / n
   return(Gamma.n)
}


# cire
cire<-function(X, y, nslices=5, d=1) 
{
   # parameters
   n<-nrow(X)
   p<-ncol(X)

   # slice y
   sy<-dr.slices(y, nslices)
   h<-sy$nslices

   # components of cire 
   beta.hat<-beta.cire(X, y, sy)
   Gamma.hat<-Gamma.cire(X, y, sy)

   # initial estimate of B
   B.init<-sir(X, y, nslices=nslices, d=d)

   # alternating least squares optimization
   V<-solve(Gamma.hat)
   out<-also(beta.hat, V, B.init)

   # return
   ans<-list(eta.hat=out$B.est, gamma.hat=out$C.est, theta.hat=beta.hat, Gamma.hat=Gamma.hat, obj.val=out$F, nslices=h, iter=out$iter, Fvals=out$Fvals)
   return(ans)
}


# shrinkage cire 
cire.garot<-function(X, y, nslices=5, d=1, tau.range) 
{
   # cire estimation
   out<-cire(X, y, nslices=nslices, d=d) 
   h<-out$nslices

   # form data
   Gamma.hat.inv2<-mat.sqrt.inv(out$Gamma.hat)
   U<-Gamma.hat.inv2 %*% vec.mat(out$theta.hat)
   V<-NULL
   for(s in 1:h) {
      V<-rbind(V, diag(as.vector(out$eta.hat %*% out$gamma.hat[,s])))
   }
   V<-Gamma.hat.inv2 %*% V

   # non-negative garrote
   outs<-garrote(U=U, V=V, tau.range=tau.range)

   # basis estimate
   b.est<-apply(diag(outs$omega.hat) %*% out$eta.hat, 2, norm)

   # variable list
   vlist<-seq(1, ncol(X))[outs$omega.hat != 0]

   # return
   ans<-list(b.est=b.est, eta.hat=out$eta.hat, omega.hat=outs$omega.hat, vlist=vlist, out.ns=out, out.s=outs)
   return(ans)
}


