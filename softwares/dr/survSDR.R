#-----------------------------------------------------------------------------------------------------------------#
# Code for: 
#    "Dimension Reduction Methods for Microarrays with Application to Censored Survival Data"
#    by Lexin Li, Hongzhe Li
#
# Author: 
#    Lexin Li (lexli@ucdavis.edu)
#
# Last Modified:
#    Nov 26, 2003
# 
# Comment:
#    library(mva); library(dr); library(survival)
#    tested on R 1.7.0
#-----------------------------------------------------------------------------------------------------------------#





#----------------------------------------------------------------------
# Code of SIR for survival data, double slicing
#----------------------------------------------------------------------

dr.fit.M.kir <-function(object,z,y,w=NULL,nslices=NULL,
                        slice.info=NULL, mincl=2, del.cluster=FALSE,...) {
# get slice information
    h <- if (!is.null(nslices)) nslices else max(8, NCOL(z)+3)
    library(mva)
    slices<- if(is.null(slice.info)) dr.clust.surv(y, h, mincl=mincl, del.cluster=del.cluster) else slice.info
# initialize slice means matrix
    zmeans <- matrix(0,slices$nslices,NCOL(z))
    slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z))
# make sure weights add to n
    wts <- if(is.null(w)) rep(1,NROW(z)) else NROW(z) * w /sum(w)
# compute weighted means within slice 
    wmean <- function (x, wts) { sum(x * wts) / sum (wts) }
    for (j in 1:slices$nslices){
 ######### print(j)
      sel <- slices$slice.indicator==j
      zmeans[j,]<- apply(z[sel,],2,wmean,wts[sel])
      slice.weight[j]<-sum(wts[sel])}
# get M matrix for sir
    M <- t(zmeans) %*% apply(zmeans,2,"*",slice.weight)/ sum(slice.weight)
    return (list (M=M,slice.info=slices))
}



dr.test.kir <- function(...) {dr.test.sir(...)}



dr.clust.surv<-function(y, nslices, mincl=2, del.cluster=FALSE)
{
  u<-sort(unique(y[,2]))
  h<-nslices
  ind<-sizes<-0
  for (j in 1:length(u)) {
      pos<-which(y[,2]==u[j])
      s2<-dr.slice2(y[pos, 1], h)
      ind[pos]<-(j - 1)*h + s2$slice.indicator 
      sizes[((j - 1)*h+1) : (j * h)]<-s2$slice.sizes 
  } 
  list(slice.indicator=ind, nslices=length(u)*h, slice.sizes=sizes)
}



dr.slice.1d <- function(y,h) {
  z<-unique(y)
  if (length(z) >= h) dr.slice2(y,h) else dr.slice1(y,length(z),sort(z))
}


dr.slice1 <- function(y,h,u){
  z <- sizes <- 0
  for (j in 1:length(u)) {
      temp <- which(y==u[j])
      z[temp] <- j
      sizes[j] <- length(temp) }
  list(slice.indicator=z, nslices=length(u), slice.sizes=sizes)
}


dr.slice2<-function(y,h)
{
  or <- order(y)
  n <- length(y)
  m<-floor(n/h)
  r<-n-m*h
  start<-sp<-ans<-0
  j<-1
  while((start+m)<n) {
      if (r==0)
        start<-start
      else
        {start<-start+1
         r<-r-1
        }
       while (y[or][start+m]==y[or][start+m+1])
          start<-start+1
       sp[j]<-start+m
       start<-sp[j]
       j<-j+1
  }
  sp[j]<-n
  ans[or[1:sp[1]]] <- 1
  for (k in 2:j){ans[ or[(sp[k-1]+1):sp[k] ] ] <- k}
  list(slice.indicator=ans, nslices=j, slice.sizes=c(sp[1],diff(sp)))
}





#----------------------------------------------------------------------
# Misc functions
#----------------------------------------------------------------------

# selection based on univariate cox model
select.gene.uc<-function(data.x, data.y, alpha=0.05)
{
   stime<-data.y[,1]
   status<-data.y[,2]
   gindex<-NULL
   for(i in 1:ncol(data.x)) {
      x<-data.x[,i]
      fit<-coxph(Surv(stime, status)~x, method="breslow")
      pval<-1 - pchisq((fit$coefficients/sqrt(fit$var))^2, 1)
      if(pval < alpha) 
         gindex<-c(gindex, i)
   }

   ans<-list(index=gindex, num=length(gindex))
   return(ans)
}



# accumlative percentage of variance for principal components
apv.pc<-function(data.x, gindex=seq(1, ncol(data.x)), pc.center=F, pc.scale=F, plot=F)
{
   # data
   X0<-data.x[, gindex]
   
   # pc
   pc.out<-prcomp(X0, center=pc.center, scale=pc.scale)

   # accumulative percentage of variance
   pc.var<-pc.out$sdev^2
   total.var<-sum(pc.var)
   perc.var<-NULL
   for(i in 1:length(pc.var))
       perc.var<-c(perc.var, sum(pc.var[1:i])/total.var)
   
   if(plot)
      plot(seq(1, length(pc.var)), perc.var)

   ans<-perc.var
   return(ans)
}





#----------------------------------------------------------------------
# Code of filling missing values in gene data
#----------------------------------------------------------------------

fill.missing.1<-function(data, num.nn=8) 
{
   n<-nrow(data)
   p<-ncol(data)

   pos.no.mis.gene<-seq(1:p)[!is.na(apply(data, 2, sum))]    
   data.no.mis.gene<-data[, pos.no.mis.gene]

   for(i in 1:p) {
      if(! i %in% pos.no.mis.gene) {
         u<-data[, i]
         pos.na<-seq(1, n)[is.na(u)]

         dist.func<-function(v, u) {
            diff<-(u - v)^2
            temp<-is.na(diff)
            count.na<-sum(temp)
            diff[temp]<-0
            sum(diff)/(n - count.na)
         }
         dist<-as.vector(apply(data.no.mis.gene, 2, dist.func, u=u), mode="numeric")

         o<-order(dist)
         pos.i<-o[1:num.nn]

         u[pos.na]<-apply(data.no.mis.gene[, pos.i], 1, mean)[pos.na]
         data[, i]<-u
      }
   }
   return(data)
}



fill.missing.2<-function(data, num.nn=8) 
{
   n<-nrow(data)
   p<-ncol(data)

   for(i in 1:p) {
      u<-data[, i]
      pos.na<-seq(1, n)[is.na(u)]

      dist.func<-function(v, u) {
         diff<-(u - v)^2
         temp<-is.na(diff)
         count.na<-sum(temp)
         diff[temp]<-0
         sum(diff)/(n - count.na)
      }
      dist<-as.vector(apply(data, 2, dist.func, u=u), mode="numeric")

      o<-order(dist)
      pos.no.mis.gene<-seq(1:p)[!is.na(apply(data[pos.na, o], 2, sum))]    
      pos.i<-o[pos.no.mis.gene][1:num.nn]

      u[pos.na]<-apply(data[, pos.i], 1, mean)[pos.na]
      data[, i]<-u
   }
   
   return(data)
}





#----------------------------------------------------------------------
# Code of evaluation of survival fit
#----------------------------------------------------------------------

score.plot<-function(scores, data.y, plot=F)
{
   stime<-data.y[,1]
   status<-data.y[,2]

   o<-order(scores)
   scores.o<-scores[o]

   signs<-rep(1, length(status))
   signs[status == 0]<- -1
   stime.sign<-stime * signs
   stime.sign.o<-stime.sign[o]

   if(plot) {
      plot(scores.o, stime.sign.o, xlab="Scores", ylab="Time to Death or Censoring", type="h")
      lines(range(scores)*1.2, c(0, 0)) 
      points(0, 0, pch=19)
   }
   
   # number of discrepancies
   ans<-sum(scores * signs < 0)
   return(ans)
}



# concordant score
score.cc<-function(scores, data.y)
{
   stime<-data.y[,1]
   status<-data.y[,2]
   # retrieve uncersored data / 'dead'
   scores.1<-scores[status == 1]
   stime.1<-stime[status==1]

   m<-length(stime.1)
   cc<-0
   for(i in 1:(m-1)) {
      for(j in (i+1):m) {
         if( xor((scores.1[i] >= scores.1[j]), (stime.1[i] <= stime.1[j])) )
            cc<-cc+1
      }
   }
   tot<-m * (m-1) / 2
   ans<-cc/tot
   return(ans)
}



surv.plot<-function(scores, data.y, sep="zero", num=4, legend=T, legend.pos=c(13,1), xmax=max(stime))
{
   stime<-data.y[,1]
   status<-data.y[,2]

   if(sep == "zero") {
      group<-rep(1, length(status))
      group[scores > 0]<- 2
   }

   if(sep == "slice") {      
      group<-dr.slice.1d(scores, num)$slice.indicator
   }

   # (log rank) test for survival curve differences
   test.out<-survdiff(Surv(stime, status) ~ group)

   if(legend) {
      out<-plot(survfit(Surv(stime, status) ~ strata(group)), xlab="Time to death", ylab="Death-free survival", xmax=xmax, 
           lty=c(1, 2), legend.bty="o", legend.pos=legend.pos, legend.text=c("low-risk patients", "high-risk patients"))
   } else {
      out<-plot(survfit(Surv(stime, status) ~ strata(group)), xlab="Time to death", ylab="Death-free survival", xmax=xmax)
   }

   return(list(test=test.out, x=out$x, y=out$y))
}





#----------------------------------------------------------------------
# Code of computing survival scores
#----------------------------------------------------------------------

score.q<-function(data.x.tr, data.y.tr, data.x.te, data.y.te, q, pc.center=F, pc.scale=F, nslices=5)
{
   # training data
   y<-data.y.tr
   X0<-data.x.tr

   # pc  (q selected components)
   pc.out<-prcomp(X0, center=pc.center, scale=pc.scale)
   X<-pc.out$x[, seq(1,q)]   
   rotation<-pc.out$rotation
   colm<-as.vector(apply(X0, 2, mean))
   cols<-as.vector(apply(X0, 2, sd))

   # sdr
   obj.kir<-dr(y ~ X, method="kir", nslices=nslices)

   sdr1<-as.vector(X %*% obj.kir$evectors[,1])
   sdr12<-sdr1^2
   d0<-data.frame(cbind(y, sdr1, sdr12))
   colnames(d0)<-c("stime", "status", "sdr1", "sdr12")

   # survival
   fit<-coxph(Surv(stime, status)~sdr1 + sdr12, data=d0, method="breslow")
   scores<-as.vector(predict(fit, type="lp"))

   # scores for testing data
   for(i in 1:nrow(data.y.te)) {
      # produce 'signatures' for testing data
      X0.test<-data.x.te[i, ]
      X0.test.s<-as.vector(X0.test)
      if(pc.center)  X0.test.s<-X0.test.s - colm
      if(pc.scale)   X0.test.s<-X0.test.s / cols     
      X.test<-as.vector((X0.test.s %*% rotation)[1:q])
      sdr1.test<-as.vector(X.test %*% obj.kir$evectors[,1])
      sdr12.test<-sdr1.test^2
      d.test<-data.frame(matrix(c(sdr1.test, sdr12.test),nrow=1))
      colnames(d.test)<-c("sdr1", "sdr12")
   
      # score for the i-th observation
      score.i<-as.vector(predict(fit, newdata=d.test, type="lp"))
      scores<-c(scores, score.i)
   }
   return(scores)
}



score.pc.q<-function(data.x.tr, data.y.tr, data.x.te, data.y.te, q, pc.center=F, pc.scale=F)
{
   # training data
   y<-data.y.tr
   X0<-data.x.tr

   # pc  (q selected components)
   pc.out<-prcomp(X0, center=pc.center, scale=pc.scale)
   X<-pc.out$x[, seq(1,q)]   
   rotation<-pc.out$rotation
   colm<-as.vector(apply(X0, 2, mean))
   cols<-as.vector(apply(X0, 2, sd))
   xnames<-"PC1"
   for(j in 2:q) xnames<-c(xnames, paste("PC", j, sep=""))

   # survival
   d0<-data.frame(cbind(y, X))
   colnames(d0)<-c("stime", "status", xnames)
   xform<-paste("Surv(stime, status) ~", xnames[1])
   for(j in 2:q) xform<-paste(xform, "+", xnames[j])
   xform<-as.formula(xform)
   fit<-coxph(xform, data=d0, method="breslow")
   scores<-as.vector(predict(fit, type="lp"))

   # scores for testing data
   for(i in 1:nrow(data.y.te)) {
      # produce 'signatures' for testing data
      X0.test<-data.x.te[i, ]
      X0.test.s<-as.vector(X0.test)
      if(pc.center)  X0.test.s<-X0.test.s - colm
      if(pc.scale)   X0.test.s<-X0.test.s / cols     
      X.test<-as.vector((X0.test.s %*% rotation)[1:q])
      d.test<-data.frame(matrix(c(X.test),nrow=1))
      colnames(d.test)<-xnames

      # score for the i-th observation
      score.i<-as.vector(predict(fit, newdata=d.test, type="lp"))
      scores<-c(scores, score.i)
   }
   return(scores)
}
     










