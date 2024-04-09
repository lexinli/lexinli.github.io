#######################################################################
#### Code for "Mixed-Effect Time-Varying Network Model and         ####
#### Application in Brain Connectivity Analysis"                   ####
#######################################################################

library(lme4)
library(Iso)
library(splines)
library(genlasso)

#######################################################################
#### Truncate a vector x to keep s components of largest magnitude ####
#######################################################################

mytruncate<-function(x,s){
	xtruncate = rep(0, length(x))
	if(length(which(x!=0)) ==0 ){
		xtruncate = x 
	}else{
		xtruncate[ which(rank(-abs(x))<= s)] = x[which(rank(-abs(x))<= s)]
		if(length(which(xtruncate!=0)) < s){
			## the s+1-th rank is tied, pick one randomly.
			tie_index = which(s<rank(-abs(x)) & rank(-abs(x))<= s+1)
			if(length(tie_index) >= 1){
				xtruncate[tie_index[1]] = x[tie_index[1]]
			}
		}
	}
	xtruncate
}

###############################################################
#### Fuse a vector beta to keep at most f different values ####
###############################################################

myfuse<-function(beta, f){
	d = length(beta)
	D = getD1dSparse(d)
	fusetmp = mytruncate( D %*% beta, f-1)
	fuseindex = rep(0, d)
	fuseindex[1] = 1
	for(i in 1:(d-1)){

		if(fusetmp[i] == 0){
			fuseindex[i+1] = fuseindex[i]
		}else{
			fuseindex[i+1] = fuseindex[i] + 1
		}
	}
	xfuse = beta
	for(i in 1: length(unique(fuseindex))){
		xfuse[which(fuseindex == i)] = mean(beta[which(fuseindex == i)])
	}
	xfuse
}


##############################################
######### somg auxiliary functions ###########
##############################################
logit<-function(x){
	log(x/(1-x))
}

logistic<-function(x){
	exp(x)/(1+exp(x))
}

bprob<-function(seg,y,t){
	n<-length(seg)
	unlist(lapply(t,function(x) {y[max(c(1,c(1:n)[seg<x]))]}))
}


###################################################################
######### functions for simulating data ###########################
## n: the number of nodes                                        ##    
## K: the number of communities                                  ##
## iter(intra): between(within) community connecting probability ##
## sigma: random effect sd                                       ##
###################################################################

generateMatrix<-function(n,K,inter,intra,sigma){
    npc<-n/K
	C<-diag(K)%x%matrix(1,npc,1)
	r<-rnorm(1,0,sd=sigma)
    p1<-logistic(logit(inter)+r)
    r1<-logistic(logit(intra)+r)
    P<-matrix(p1,K,K);diag(P)<-r1
    A<-1*(matrix(runif(n*n),n,n)<C%*%P%*%t(C))
	A[lower.tri(A)]<-t(A)[lower.tri(A)]; diag(A)<-0
	A
}



################################################################
#################### main function  ############################
# network.list: a list of N binary networks,each with n nodes ##
# age: ages of the N subjects                                 ##
# t: binning of the age range                                 ##
# subnet: known community structure                           ##
################################################################

REdnet<-function(network.list,age,t,subnet){
	#step 1
	coef<-matrix(0,0,K*(K+1)/2) 
	coef_se<-matrix(0,0,K*(K+1)/2)
	sighat<-c()
	for(s in 1:(length(t)-1)){
		glmmatrix<-matrix(0,0,4)
		for (k in c(1:length(age))[age>=t[s] & age<t[s+1]]){
			ind<-0
			for(i in 1:K){
				for (j in i:K){
					ind<-ind+1
				    sub<-network.list[[k]][subnet==i,subnet==j]
			        n1<-dim(sub)[1]
			        n2<-dim(sub)[2]
			        if(i==j){
			            glmmatrix<-rbind(glmmatrix,c(k,sum(sub)/2,n1*(n1-1)/2,ind))
			        }
			        if(i!=j){
			            glmmatrix<-rbind(glmmatrix,c(k,sum(sub),n1*n1,ind))
			        }
			    }
			}
		}
		colnames(glmmatrix)<-c("subject","success","trial","block")
		glmmatrix<-data.frame(glmmatrix)
		glmmatrix$block<-as.factor(glmmatrix$block)
		glmmatrix$subject<-as.factor(glmmatrix$subject)
		glm1<-glmer(cbind(success,trial-success)~ block+ (1|subject)-1, data=glmmatrix, family = binomial,control=glmerControl(optCtrl=list(maxfun=50000)))
		coef<-rbind(coef,getME(glm1,"beta"))
		coef_se<-rbind(coef_se,sqrt(diag(vcov(glm1))))
		sighat<-c(sighat,getME(glm1,"theta"))
	} 
	# coef is unconstrained estimator 
	# sighat is random effect sd estimate
    
    #step 2
	coef_shape<-c()
	for(p in 1:(K*(K+1)/2)){
		z <- ufit(coef[,p],x=t[-1],type="b")
		coef_shape<-cbind(coef_shape,z$y)
	}
	# coef_shape is shape contrained estimator

    #step 3
	allsub<-list()
	for(s in 1:(length(t)-1)){
		glmmatrix<-matrix(0,0,4)
		for (k in c(1:length(age))[age>=t[s] & age<t[s+1]]){
			ind<-0
			for(i in 1:K){
				for (j in i:K){
					ind<-ind+1
					sub<-network.list[[k]][subnet==i,subnet==j]
					n1<-dim(sub)[1]
					n2<-dim(sub)[2]
					if(i==j){ 
						diag(sub)<-0
				        glmmatrix<-rbind(glmmatrix,c(k,sum(sub)/2,n1*(n1-1)/2,ind))
				    }
				    if(i!=j){ 
				        glmmatrix<-rbind(glmmatrix,c(k,sum(sub),n1*n1,ind))
				    }
				}
		    }
		}
		allsub[[s]]<-glmmatrix
	}
    # lam_vec is the tuning parameter for fusion, selected using BIC
	lam_vec<-c()
	for(ind in 1:(K*(K+1)/2)){
		BIC_single<-function(coef,sighat,ind){
		r5 <- GHrule(50, asMatrix=FALSE)
		logsum<-0
		for(s in 1:(length(t)-1)){
		    glmmatrix<-allsub[[s]]
		    for(k in c(1:length(age))[age>=t[s] & age<t[s+1]]){
		    	submatrix<-glmmatrix[glmmatrix[,1]==k,]
		    	noname<-(submatrix[ind,2])*log(logistic(r5$z*sighat[s]+coef[s,ind]))+(submatrix[ind,3]-submatrix[ind,2])*log(1-logistic(r5$z*sighat[s]+coef[s,ind]))
		    	logsum<-logsum+log(exp(noname-max(noname))%*%r5$w)+max(noname)
		    }
		}
		-2*logsum+log(length(age))*sum(!duplicated(coef[,ind]))
	    }
	    BICvec_single<-c()
	    for(f in 1:25){
	    	BICvec_single<-c(BICvec_single,BIC_single(apply(coef_shape,2,FUN=myfuse,f=f),sighat,ind))
	    }
	    lam_vec[ind]<-c(1:25)[which.min(BICvec_single)]
	}
    coef_shape_fuse<-c()
	for(ind in 1:(K*(K+1)/2)){
		coef_shape_fuse<-cbind(coef_shape_fuse,myfuse(coef_shape[,ind],lam_vec[ind]))
	}
	# coef_shape_fuse is shape and fuse contrained estimator

    ## the columns in the coefficent output have length equal to the bins selected (i.e., length of t)
	return(list(coef,coef_shape,coef_shape_fuse,lam_vec))

}


###########################
######### example 1 #######
###########################


## setting the parameters
age<-seq(0,1,length=201)[-201]  #subject age (N=200)
n=100; K=5; r1=0.1; r2=0.2 #network information
subnet<-as.vector(outer(rep(1,n/K),(1:K)))  #known community membership
seg<-c(0,0.2,0.6,1) # seg: vector with end points of the step function intervals 
y<-c(0.05,0.1,0.2)  # y: values of the step function  


## network simulation
network.list<-list()
for(s in 1:200){
	time<-seq(0,1,length=201)[s]
	network.list[[s]]<-generateMatrix(n=n,K=K,inter=bprob(seg,y,time)+r2,intra=bprob(seg,y,time)+r1+r2,sigma=0.1)
}

## estimation with 20 equal sized bins
t<-c(age[c((0:19)*10+1)],1)    
output<-REdnet(network.list=network.list,age=age,t=t,subnet=subnet)



###########################
######### example 2 #######
###########################
intra_prob<-function(t){
	0.2/(1+exp(-25*(t-0.4)))+0.1
}

inter_prob<-function(t){
	tt<-1-t
	0.2/(1+exp(-25*(tt-0.4)))+0.1
}


sigma<-function(t){
	seg<-c(0,0.5,1)
	n<-length(seg)
	y<-c(0.2,0.1)
	unlist(lapply(t,function(x) {y[max(c(1,c(1:n)[seg<x]))]}))
}
   
age<-seq(0,1,length=201)[-201]       #subject age (N=200)
n=100; K=5                           #network information
subnet<-as.vector(outer(rep(1,n/K),(1:K)))  #known community membership


## network simulation
network.list<-list()
for(s in 1:200){
	time<-age[s] 
	network.list[[s]]<-generateMatrix(n=n,K=K,inter=inter_prob(time),intra=intra_prob(time),sigma=sigma(time))
}

## estimation with 20 equal sized bins
t<-c(age[c((0:19)*10+1)],1)          
output<-REdnet(network.list=network.list,age=age,t=t,subnet=subnet)


