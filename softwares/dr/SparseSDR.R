#---------------------------------------------------------------------------------------------# 
# R functions for sparse sufficient dimension reduction                                       #
# Author: Lexin Li                                                                            # 
# Reference: Li, L. (2007). Sparse sufficient dimension reduction. Biometrika. 94, 603-613.   #                                                          #
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

# Gram-Schmidt orthonormalization
orthnormal<-function(X)
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



#----------------------------------------------------------------------------------------------
# functions for sparse principal components wrt I inner product 
# (copied from Hui Zou's elasticnet library)
#----------------------------------------------------------------------------------------------

updateRR<-
function(xnew, R = NULL, xold, lambda,eps = .Machine$double.eps)
{
	xtx <- (sum(xnew^2)+lambda)/(1+lambda)
	norm.xnew <- sqrt(xtx)
	if(is.null(R)) {
		R <- matrix(norm.xnew, 1, 1)
		attr(R, "rank") <- 1
		return(R)
	}
	Xtx <- drop(t(xnew) %*% xold)/(1+lambda)
	r <- backsolvet(R, Xtx)
	rpp <- norm.xnew^2 - sum(r^2)
	rank <- attr(R, "rank")	
	if(rpp <= eps)
		rpp <- eps
	else {
		rpp <- sqrt(rpp)
		rank <- rank + 1
	}
	R <- cbind(rbind(R, 0), c(r, rpp))
	attr(R, "rank") <- rank
	R
}


solvebeta<-function(x, y, paras, max.steps, sparse=c("penalty","varnum"), eps = .Machine$double.eps)
{
       sparse <- match.arg(sparse)
       if (missing(sparse)) sparse <- "penalty"
	nm <- dim(x)
	n <- nm[1]
	m <- nm[2]
	im <- seq(m)
	one <- rep(1, n)
	vn <- dimnames(x)[[2]]
        lambda<-paras[1]
        if(lambda>0){
	   maxvars <- m
        }
        if (lambda==0) {
           maxvars <- min(m,n-1)
           if (m==n){
             maxvars<-m
           }
        }
	d1 <- sqrt(lambda)
	d2 <- 1/sqrt(1 + lambda)
	Cvec <- drop(t(y) %*% x) * d2
        ssy <- sum(y^2)
	residuals <- c(y, rep(0, m))
        if(missing(max.steps)) {max.steps <- 50 * maxvars}
        penalty<-max(abs(Cvec))
        if (sparse=="penalty" && penalty*2/d2<=paras[2])
          { 
           beta<-rep(0,m)
          }
        else {
               beta <- rep(0, m)
               first.in <- integer(m)
               active <- NULL
               ignores <- NULL
               actions <- as.list(seq(max.steps))
               drops <- FALSE
               Sign <- NULL
               R <- NULL
	       k <- 0
	       while((k < max.steps) & (length(active) < maxvars - length(ignores))) 
		{
		 action <- NULL
		 k <- k + 1
        	 inactive <- if(k == 1) im else im[ - c(active, ignores)]
		 C <- Cvec[inactive]
		 Cmax <- max(abs(C))
		 if(!any(drops)) {
			new <- abs(C) == Cmax 
			C <- C[!new]
			new <- inactive[new]       
			for(inew in new) {
                        R <- updateRR(x[, inew], R, x[, active], lambda) 
                               if(attr(R, "rank") == length(active)) {
					nR <- seq(length(active))
					R <- R[nR, nR, drop = FALSE]
					attr(R, "rank") <- length(active)
					ignores <- c(ignores, inew)
					action <- c(action,  - inew)
				}
				else {
					if(first.in[inew] == 0)
						first.in[inew] <- k
					active <- c(active, inew)
					Sign <- c(Sign, sign(Cvec[inew]))
					action <- c(action, inew)
				}
			}
		}
		else action <-  - dropid
                Gi1 <- backsolve(R, backsolvet(R, Sign))
            	A <- 1/sqrt(sum(Gi1 * Sign))
            	w <- A * Gi1
                u1<-drop(x[,active,drop=FALSE]%*%w*d2)  
                u2<-rep(0,m)
                u2[active]<-d1*d2*w
                u<-c(u1,u2)
                if(lambda>0){
	            maxvars <- m-length(ignores)
                           }
                if (lambda==0){
                    maxvars <- min(m-length(ignores),n-1)
                           }
       		if(length(active) == maxvars - length(ignores)) {
			gamhat <- Cmax/A
		}
		else {
			a <- (drop(u1 %*% x[,  - c(active, ignores)]) + d1 *
				u2[ - c(active, ignores)]) * d2
      			gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
				gamhat <- min(gam[gam > eps], Cmax/A)                        
			Cdrop <- c(C - gamhat * a,  - C + gamhat * a) - (Cmax -
				gamhat * A)
		}
		dropid <- NULL
		b1 <- beta[active]
		z1 <-  - b1/w
		zmin <- min(z1[z1 > eps], gamhat)
		if(zmin < gamhat) {
			gamhat <- zmin
			drops <- z1 == zmin
		}
               else drops <- FALSE
               beta2<-beta          
               beta[active] <- beta[active] + gamhat * w
               residuals <- residuals - (gamhat*u)
               Cvec <- (drop(t(residuals[1:n]) %*% x) + d1 * residuals[ - (1:n)]) * d2
               penalty <- c(penalty,penalty[k]-abs(gamhat*A))
               if(sparse=="penalty" && rev(penalty)[1]*2/d2<=paras[2]){
                   s1<-rev(penalty)[1]*2/d2
                   s2<-rev(penalty)[2]*2/d2
                   beta<-(s2-paras[2])/(s2-s1)*beta+(paras[2]-s1)/(s2-s1)*beta2
                   beta<-beta*d2
                   break
               }
		if(any(drops)) {
			dropid <- seq(drops)[drops]
			for(id in rev(dropid)) {
				R <- downdateR(R, id)
      			}
			dropid <- active[drops]
                        beta[dropid] <- 0
			active <- active[!drops]
			Sign <- Sign[!drops]
                      }
               if(sparse=="varnum" && length(active)>=paras[2]){
                 break
               }
	       if(!is.null(vn))
			names(action) <- vn[abs(action)]
                	actions[[k]] <- action 
	}
             }
     
        return (beta)
}



#----------------------------------------------------------------------------------------------
# functions for sparse sdr wrt G inner product
#----------------------------------------------------------------------------------------------

comp.sdr<-function(X, y, method, d, nslices)
{
   out.m<-switch(method,
                 pc     = comp.sdr.pc(X, y, d), 
                 sir    = comp.sdr.sir(X, y, d, nslices),
                 save   = comp.sdr.save(X, y, d, nslices),
                 phdres = comp.sdr.phdres(X, y, d)
                 )

   ans<-list(m=out.m$m, G=out.m$G, beta.sdr=out.m$beta.sdr)
   return(ans)
}



comp.sdr.pc<-function(X, y, d)
{
   # m matrix
   Xc<-apply(X, 2, center)
   m<-mat.sqrt(t(Xc) %*% Xc)

   # G matrix 
   G<-diag(1, ncol(X))

   # pc
   v<-eigen(cov.x(X))$vectors[, 1:d]
   if(d == 1) v<-matrix(v, ncol=1)

   # return
   ans<-list(m=m, G=G, beta.sdr=v)
   return(ans)
}



comp.sdr.sir<-function(X, y, d, nslices)
{
   # parameters
   n<-nrow(X)
   p<-ncol(X)
   Sigma.x<-cov.x(X)
   Sigma.x2<-mat.sqrt(Sigma.x)
   Sigma.x.inv2<-mat.sqrt.inv(Sigma.x)

   # standardize X
   Z<-apply(X, 2, center) %*% Sigma.x.inv2

   # slice y
   sy<-dr.slices(y, nslices)
   nslices<-sy$nslices

   # compute sdr kernel matrix
   M.sir.z<-matrix(0, nrow=p, ncol=p)
   for(s in 1:nslices) {
      Z.s<-Z[sy$slice.indicator == s, ]
      if(sy$slice.sizes[s] == 1) Z.s<-matrix(Z.s, nrow=1)
      Z.sm<-as.vector(apply(Z.s, 2, mean))
      M.sir.z<-M.sir.z + (sy$slice.sizes[s]/n) * Z.sm %*% t(Z.sm)
   }
   M.sir<-Sigma.x2 %*% M.sir.z %*% Sigma.x2
   #m<-mat.sqrt(n*M.sir)
   m<-mat.sqrt(M.sir)

   # compute sdr estimate w/o shrinkage
   v<-eigen(M.sir.z)$vectors[,1:d]
   if(d == 1) v<-matrix(v, ncol=1)
   beta.sdr<-Sigma.x.inv2 %*% v
   beta.sdr<-apply(beta.sdr, 2, norm)

   # return
   ans<-list(m=m, G=Sigma.x, beta.sdr=beta.sdr)
   return(ans)
}



comp.sdr.save<-function(X, y, d, nslices)
{
   # parameters
   n<-nrow(X)
   p<-ncol(X)
   Sigma.x<-cov.x(X)
   Sigma.x2<-mat.sqrt(Sigma.x)
   Sigma.x.inv2<-mat.sqrt.inv(Sigma.x)

   # standardize X
   Z<-apply(X, 2, center) %*% Sigma.x.inv2

   # slice y
   sy<-dr.slices(y, nslices)
   nslices<-sy$nslices

   # compute sdr kernel matrix
   M.save.z<-matrix(0, nrow=p, ncol=p)
   for(s in 1:nslices) {
      Z.s<-Z[sy$slice.indicator == s, ]
      if(sy$slice.sizes[s] == 1) Z.s<-matrix(Z.s, nrow=1)
      iVz<-diag(1, p) - cov.x(Z.s) 
      M.save.z<-M.save.z + (sy$slice.sizes[s]/n) * iVz %*% iVz 
   }
   M.save<-Sigma.x2 %*% M.save.z %*% Sigma.x2
   m<-mat.sqrt(M.save)

   # compute sdr estimate w/o shrinkage
   v<-eigen(M.save.z)$vectors[,1:d]
   if(d == 1) v<-matrix(v, ncol=1)
   beta.sdr<-Sigma.x.inv2 %*% v
   beta.sdr<-apply(beta.sdr, 2, norm)

   # return
   ans<-list(m=m, G=Sigma.x, beta.sdr=beta.sdr)
   return(ans)
}



comp.sdr.phdres<-function(X, y, d)
{
   # parameters
   n<-nrow(X)
   p<-ncol(X)
   Sigma.x<-cov.x(X)
   Sigma.x2<-mat.sqrt(Sigma.x)
   Sigma.x.inv2<-mat.sqrt.inv(Sigma.x)

   # standardize X
   Z<-apply(X, 2, center) %*% Sigma.x.inv2

   # residual
   e<-resid(lm(y~Z))

   # compute sdr kernel matrix
   M.phd.z<-matrix(0, nrow=p, ncol=p)
   for(s in 1:n) {
      M.phd.z<-M.phd.z + e[s] * Z[s,] %*% t(Z[s,])
   }
   M.phd.z<-M.phd.z / n
   M.phd<-Sigma.x2 %*% (M.phd.z %*% t(M.phd.z)) %*% Sigma.x2
   m<-mat.sqrt(M.phd)

   # compute sdr estimate w/o shrinkage
   v<-eigen(M.phd.z %*% t(M.phd.z))$vectors[,1:d]
   if(d == 1) v<-matrix(v, ncol=1)
   beta.sdr<-Sigma.x.inv2 %*% v
   beta.sdr<-apply(beta.sdr, 2, norm)

   # return
   ans<-list(m=m, G=Sigma.x, beta.sdr=beta.sdr)
   return(ans)
}



# this function serves as a wrapper to call Hui Zou's solvebeta function
solve.beta<-function(x, y, G2, lambda1, lambda2)
{  
   # transform data to an L1 problem
   x.star<-rbind(x, sqrt(lambda2) * G2) 
   y.star<-c(y, rep(0, ncol(x)))

   # call solvebeta() with lambda2=0
   beta.est<-solvebeta(x.star, y.star, para=c(0, lambda1), sparse="penalty")

   # return
   return(beta.est)
}



ssdr.lambda<-function(X, y, method=c("pc", "sir", "save", "phdres"), d=1, nslices=5, lambda1=NULL, lambda2=1e-6, max.iter=200, eps.conv=1e-3)
{
   # parameters
   n<-nrow(X)
   p<-ncol(X)
   method<-match.arg(method)
 
   # compute sdr components
   out.m<-comp.sdr(X, y, method, d, nslices)
   m<-out.m$m
   M<-t(m) %*% m   
   G<-out.m$G
   G2<-mat.sqrt(G)
   G2.inv<-mat.sqrt.inv(G2)
 
   # initial estimate of alpha and beta
   alpha<-out.m$beta.sdr      
   beta<-alpha      
   for(i in 1:d) {
      ym<-m %*% alpha[,i]
      beta[,i]<-solve.beta(m, ym, G2, lambda1[i], lambda2)
   }

   # iteration
   iter<-0
   beta.n<-apply(beta, 2, norm)
   diff.conv<-1
   while((iter < max.iter) & (diff.conv[iter+1] > eps.conv)){
      z<-svd(G2.inv %*% M %*% beta)
      alpha<-G2.inv %*% (z$u) %*% t(z$v)
      for(i in 1:d) {
         ym<-m %*% alpha[,i]
         beta[,i]<-solve.beta(m, ym, G, lambda1[i], lambda2)
      }
  
      beta.n.new<-apply(beta, 2, norm)
      diff.conv<-c(diff.conv, max(abs(beta.n.new - beta.n)))
      beta.n<-beta.n.new

      iter<-iter + 1
   }
  
   # compute objective value
   comp1<-m %*% solve(G) - m %*% beta %*% t(beta)
   rss1<-sum(diag(comp1 %*% G %*% t(comp1)))

   z<-svd(G2.inv %*% M %*% beta)
   alpha<-G2.inv %*% (z$u) %*% t(z$v)
   comp2<-m %*% solve(G) - m %*% beta %*% t(alpha)
   rss2<-sum(diag(comp2 %*% G %*% t(comp2)))

   p.e<-sum(as.vector(beta) != 0)

   # normalize beta
   beta<-apply(beta, 2, norm)

   # return
   ans<-list(beta=beta, beta0=out.m$beta.sdr, alpha=alpha, diff.conv=diff.conv, iter=iter, p.e=p.e, rss1=rss1, rss2=rss2)
   return(ans)
}



ssdr.wrap<-function(X, y, method=c("pc", "sir", "save", "phdres"), d=1, nslices=5, s1.range=NULL, s2.range=NULL, max.iter=200, eps.conv=1e-3)
{
   # parameters 
   n<-nrow(X)

   # compute criteria for all s1 and s2
   crit.all<-list()
   for(j in 1:length(s2.range)) {
      s2<-s2.range[j]

      crit<-NULL
      for(k in 1:length(s1.range)) {
         s1<-s1.range[k]

         out1<-ssdr.lambda(X, y, method=method, d=d, nslices=nslices, lambda1=rep(s1, d), lambda2=s2)

         aic1<-n*out1$rss1 + 2 * out1$p.e
         bic1<-n*out1$rss1 + log(n) * out1$p.e
         out.crit<-c(aic1, bic1)
         crit<-cbind(crit, out.crit) 
      }
  
      crit.all[[j]]<-crit
   }

   # locate optimal s 
   s1.min<-NULL
   ct.min<-NULL
   for(j in 1:length(s2.range)) {
      s1.min.j<-ct.min.j<-NULL
      for(l in 1:length(out.crit)) {
         s1.min.j<-c(s1.min.j, s1.range[order(crit.all[[j]][l,])[1]])
         ct.min.j<-c(ct.min.j, min(crit.all[[j]][l,], na.rm=T))
      }
      s1.min<-rbind(s1.min, s1.min.j)
      ct.min<-rbind(ct.min, ct.min.j)
   }
   rownames(s1.min)<-as.character(s2.range); colnames(s1.min)<-c("aic1", "bic1")
   rownames(ct.min)<-as.character(s2.range); colnames(ct.min)<-c("aic1", "bic1")

   # beta estimate with given optimal s
   beta.est.all<-NULL
   s12.est.all<-NULL
   for(l in 1:length(out.crit)) {
      pos<-order(ct.min[, l])[1]
      s1<-s1.min[pos, l]
      s2<-s2.range[pos]
      out1<-ssdr.lambda(X, y, method=method, d=d, nslices=nslices, lambda1=rep(s1, d), lambda2=s2)
      beta.est.all<-cbind(beta.est.all, out1$beta)
      s12.est.all <-cbind(s12.est.all,  c(s1, s2))
   }
   rownames(s12.est.all)<-c("s1", "s2"); colnames(s12.est.all)<-c("aic1", "bic1")

   # return
   ans<-list(beta.est.all=beta.est.all, beta.est0=out1$beta0, s12.est.all=s12.est.all, crit.all=crit.all, s1.min=s1.min, ct.min=ct.min)
   return(ans)
}













