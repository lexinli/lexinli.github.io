## Implementation of Dynamic Tensor Clustering
## Author: Wei Sun



library(mnormt)
library(MASS)
library(rTensor)
#library(caret) # in case there is an error in "library(caret)": install.packages("~/Dropbox/2016-SunWei/code/pbkrtest_0.4-5.tar.gz")
library(genlasso) #generalized lasso
library(cluster)

EPS = 1e-4
norm = base::norm

#####################
## Utility Functions
#####################
fnorm <- function(T){
	### Tensor Frobenius Norm \| T \|_F ###
	return(sqrt(sum(T^2)))
}

mysd = function(x){
	if(norm(x,"2") == 0){
		return(x)
	}else{
		return(x/sqrt(t(x)%*%x))
	}	
}


mytruncate = function(x,s){
	## Truncate a vector x to keep at most s components of largest magnitude

	if(s > length(x)){
		print("Wanring: s is larger than the size of the vector; use the size of the vector to substitute s!")
		s = length(x)
	}
	if(s <= 0 ){
		stop("s should be a positive integer!")
	}
	set.seed(1)
	xtruncate = rep(0, length(x))
	Rank_vec = rank(-abs(x))
	xtruncate[ which(Rank_vec<= s)] = x[which(Rank_vec<= s)]
	nonzero = length(which(xtruncate!=0))
	if(nonzero > s){
		## means some values in xtruncate are tied, remove some randomly.
		tie_value = -sort(-abs(xtruncate))[s]
		remove_index = sample(which(abs(xtruncate) ==tie_value), nonzero - s)
		xtruncate[remove_index] = 0
		nonzero = length(which(xtruncate!=0))
	}
	t = 0
	while(nonzero < s && t < s-nonzero ){
		## if the s+1-th rank is tied, pick additional ones randomly.
		t = t + 1
		if(max(abs(xtruncate)) == 0){
			tie_max_value = max(abs(x))
		}else{
			tie_max_value = max(abs(x)[which(abs(x) < min(abs(xtruncate[which(xtruncate!=0)])))])
		}
		tie_index = which(abs(x) == tie_max_value)
		if(length(tie_index) > s-nonzero){
			## only need to randomly pick s-nonzero elements.
			use_index = sample(tie_index, s-nonzero)
			xtruncate[use_index] = x[use_index]
			nonzero = length(which(xtruncate!=0))
		}else{
			use_index = sample(tie_index, s-nonzero)
			xtruncate[use_index] = x[use_index]
			nonzero = length(which(xtruncate!=0))
		}
	}
	xtruncate
}

myerrorlist = function(beta.old, beta.new){
	## compute the F-norm distance of two lists of recovered rank1-components, if there have different Rank, use the smaller rank to compute error!

	nlist = length(beta.old)
	for(i in 1:nlist){
		beta.old[[i]] = as.matrix(beta.old[[i]])
		beta.new[[i]] = as.matrix(beta.new[[i]])
	}

	Rank_old = ncol(beta.old[[1]])
	Rank_new = ncol(beta.new[[1]])

	if(Rank_old != Rank_new){
		Rank = min(Rank_old, Rank_new)
		error = 0
		for(i in 1:nlist){
			for(j in 1:Rank){
				error = error + min(fnorm(beta.old[[i]][,j] - beta.new[[i]][,j])^2, fnorm(beta.old[[i]][,j] + beta.new[[i]][,j])^2)
			}
		}
		out = sqrt(error/(nlist*Rank))
	}else{
		Rank = ncol(beta.old[[1]])
		error = 0
		for(i in 1:nlist){
			for(j in 1:Rank){
				error = error + min(fnorm(beta.old[[i]][,j] - beta.new[[i]][,j])^2, fnorm(beta.old[[i]][,j] + beta.new[[i]][,j])^2)
			}
		}
		out = sqrt(error/(nlist*Rank))
	}
	out

}


myerror_betaw = function(Beta.est, Beta.true, w.est, w.true){
	#compute the final estimation error of components and weight w.
	
	A.est = Beta.est[[1]]
	B.est = Beta.est[[2]]
	C.est = Beta.est[[3]]
	A = truebeta[[1]]
	B = truebeta[[2]]
	C = truebeta[[3]]
	K 	  = dim(A)[2]
	K.est = dim(A.est)[2]	
	if(K.est > K){
		# more latent decomposition. truncate it to first K columns.
		A.est = A.est[,1:K]
		B.est = B.est[,1:K] 
		C.est = C.est[,1:K]
		w.est = w.est[1:K]
		warning("K.est = ", K.est, " is larger than true K = ", K, ". Use first K latent estimations!")
	}else if(K.est < K){
		# less latent decomposition. add columns with zeros to fill out K columns.
		virtual.A = c(1, rep(0, nrow(A.est) - 1))
		virtual.B = c(1, rep(0, nrow(B.est) - 1))
		virtual.C = c(1, rep(0, nrow(C.est) - 1))
		A.est = cbind(A.est, matrix(rep(virtual.A, K-K.est),nrow(A.est), K-K.est))
		B.est = cbind(B.est, matrix(rep(virtual.B, K-K.est),nrow(B.est), K-K.est)) 
		C.est = cbind(C.est, matrix(rep(virtual.C, K-K.est),nrow(C.est), K-K.est))
		w.est = c(w.est,rep(0,K-K.est))
		warning("K.est = ", K.est, " is smaller than true K = ", K, ". Fill up zeros to get K latent estimations!")
	}

	## estimation error 
	ERROR = matrix(NA,K,K)
	j.best = rep(0,K)

	for(i in 1:K){
		for(j in 1:K){	
			ERROR[i,j] = min(norm(A.est[,i] - A[,j],type="2"),norm(A.est[,i] + A[,j],type="2"))
		}
		j.best[i] = which.min(ERROR[i,])	
	}

	error.abc = rep(0, K)
	error.w = rep(0,K)
	for(ii in 1:K){

		error.a = min(norm(A.est[,ii] - A[,j.best[ii]],type="2"), norm(A.est[,ii] + A[,j.best[ii]],type="2"))
		error.b = min(norm(B.est[,ii] - B[,j.best[ii]],type="2"), norm(B.est[,ii] + B[,j.best[ii]],type="2"))
		error.c = min(norm(C.est[,ii] - C[,j.best[ii]],type="2"), norm(C.est[,ii] + C[,j.best[ii]],type="2"))
		error.abc[ii] = mean(error.a, error.b, error.c)
		error.w[ii] = abs((w.est[ii] - w.true[j.best[ii]])/w.true[j.best[ii]])
	
	}

	P = list()
	P$error.beta = mean(error.abc)
	P$error.w = mean(error.w)
	P
}


myerror_tensor = function(Beta.est, Beta.true, w.est, w.true){
	#compute the tensor recovery error

	num_m = length(Beta.true)

	if(num_m == 3){
		d1 = dim(Beta.true[[1]])[1]
		d2 = dim(Beta.true[[2]])[1]
		d3 = dim(Beta.true[[3]])[1]
		R.true = length(w.true)
		R.hat = length(w.est)
		T.true = array(0, c(d1,d2,d3))
		T.hat = array(0, c(d1,d2,d3))
		for(rrr in 1:R.true){
			T.true = T.true + w.true[rrr] * outer(outer(as.numeric(Beta.true[[1]][,rrr]), as.numeric(Beta.true[[2]][,rrr])), as.numeric(Beta.true[[3]][,rrr]))
		}
		for(rrr in 1:R.hat){
			T.hat = T.hat + w.est[rrr] * outer(outer(as.numeric(Beta.est[[1]][,rrr]), as.numeric(Beta.est[[2]][,rrr])), as.numeric(Beta.est[[3]][,rrr]))
		}
		error = sqrt(sum((T.hat - T.true)^2))/ sqrt(sum((T.true)^2))

	}else if(num_m == 4){

		d1 = dim(Beta.true[[1]])[1]
		d2 = dim(Beta.true[[2]])[1]
		d3 = dim(Beta.true[[3]])[1]
		d4 = dim(Beta.true[[4]])[1]
		R.true = length(w.true)
		R.hat = length(w.est)
		T.true = array(0, c(d1,d2,d3,d4))
		T.hat = array(0, c(d1,d2,d3,d4))
		for(rrr in 1:R.true){
			T.true = T.true + w.true[rrr] * outer(outer(outer(as.numeric(Beta.true[[1]][,rrr]), as.numeric(Beta.true[[2]][,rrr])), as.numeric(Beta.true[[3]][,rrr])), as.numeric(Beta.true[[4]][,rrr]))
		}
		for(rrr in 1:R.hat){
			T.hat = T.hat + w.est[rrr] * outer(outer(outer(as.numeric(Beta.est[[1]][,rrr]), as.numeric(Beta.est[[2]][,rrr])), as.numeric(Beta.est[[3]][,rrr])), as.numeric(Beta.est[[4]][,rrr]))
		}
		error = sqrt(sum((T.hat - T.true)^2))/ sqrt(sum((T.true)^2))
	}
	error
}




Tabc = function(T, a, b, c){
	#compute tensor vector product T(a,b,c) for tensor T and vectors, a,b,c with different dimensions.
	
	d1 = dim(T)[1]	
	d2 = dim(T)[2]
	d3 = dim(T)[3]
	tmp = 0
	for(i in 1:d1){
		for(j in 1:d2){	
			for(k in 1:d3){		
				tmp = tmp + a[i]*b[j]*c[k]*T[i,j,k]
			}
		}
	}
	tmp
}

Tabcd = function(T, a, b, c, d){
	#compute tensor vector product T(a,b,c,d) for tensor T and vectors, a,b,c,d with different dimensions.
	
	d1 = dim(T)[1]	
	d2 = dim(T)[2]
	d3 = dim(T)[3]
	d4 = dim(T)[4]
	tmp = 0
	for(i in 1:d1){
		for(j in 1:d2){	
			for(k in 1:d3){	
				for(l in 1:d4){
					tmp = tmp + a[i]*b[j]*c[k]*d[l]*T[i,j,k,l]
				}	
			}
		}
	}
	tmp
}


clus.dist = function(clus1, clus2) {
	## Input two vectors of membership, output the disagreement between them. if one of them is truth, output is the clustering error.
	if (length(clus1)!=length(clus2)) return("cluster sizes don't match!")
	n=length(clus1)
	s=0
	for ( i in 2:n ) {
		for ( j in 1:(i-1)) {
			s=s+as.numeric(clus1[i]==clus1[j] & clus2[i]!=clus2[j])
			s=s+as.numeric(clus1[i]!=clus1[j] & clus2[i]==clus2[j])
		}
	}
	return(2*s/(n*(n-1)))
}


TuneK_gap = function(gap_result, K.max){
	# gap_result: an object returns from "clusGap" function
	# K.max: the maximal K considered
	# reference: https://web.stanford.edu/~hastie/Papers/gap.pdf
	gap_se = gap_result$Tab[,3] - gap_result$Tab[,4]
	satisfied_index = NULL
	if(K.max == 1){
		output = K.max
	}else{

		for(kk in 2:(K.max-1)){
			if(gap_result$Tab[kk,3] >= gap_se[kk+1]){
				satisfied_index = append(satisfied_index, kk)
			}
		}
		if(is.null(satisfied_index) == TRUE){
			output = K.max
		}else{
			output = min(satisfied_index)
		}
	}
	output
}


Count_change = function(x){
	# x: a vector
	# output: the number of changes in vector x.

	n = length(x)
	count = 0
	for(i in 1:(n-1)){
		if(x[i] != x[i+1]){count = count + 1}
	}
	count
}


count.rows <-function(x){
	if(ncol(x) == 1){
		value = unique(x)
	}else{

		order.x <- do.call(order,as.data.frame(x))
		equal.to.previous <-
		rowSums(x[tail(order.x,-1),] != x[head(order.x,-1),])==0
		tf.runs <- rle(equal.to.previous)
		counts <- c(1,
		         unlist(mapply( function(x,y) if (y) x+1 else (rep(1,x)),
		                       tf.runs$length, tf.runs$value )))
		counts <- counts[ c(diff(counts) <= 0, TRUE ) ]
		unique.rows <- which( c(TRUE, !equal.to.previous ) )
		value = cbind( counts, x[order.x[ unique.rows ], ,drop=F] )
	}
}


mydf = function(beta){
	## compute the degree of freedom when there are both sparsity and fuse structures.
	beta_nonzero = beta[which(beta!=0)]
	length(unique(beta_nonzero))
}

myBIC = function(W, Beta, T){
	## W: the list of weights (w_1, ..., w_R) outputed from stdtruncate() function
	## Beta: the list of beta's outputed from stdtruncate() function
	## T: the original tensor
	## output: BIC value

	rank = length(W)
	num_m = length(Beta)
	num_dim = rep(0, num_m)
	df.est = 0
	for(j in 1:num_m){
		num_dim[j] = nrow(Beta[[j]])
		df.est = df.est + sum(apply(Beta[[j]], 2, mydf))
	}

	if(num_m == 3){
		# mode-3 tensor
		d1 = num_dim[1]
		d2 = num_dim[2]
		d3 = num_dim[3]
		T.est = array(0, c(d1, d2, d3))
		for(i in 1:rank){
			T.est = T.est + W[i] * outer(outer(Beta[[1]][,i], Beta[[2]][,i]), Beta[[3]][,i])
		}
	    bic = log(fnorm(T - T.est)^2 / (d1 * d2 * d3)) + log(d1 * d2 * d3) / (d1 * d2 * d3) * df.est
	    bic2 = log(fnorm(T - T.est)^2 / (d1 * d2 * d3)) + log(d1 * d2 * d3) / (d1 * d2 * d3) * df.est / num_m
	}else if(num_m == 4){
		# mode-4 tensor
		d1 = num_dim[1]
		d2 = num_dim[2]
		d3 = num_dim[3]
		d4 = num_dim[4]
		T.est = array(0, c(d1, d2, d3, d4))
		for(i in 1:rank){
			T.est = T.est + W[i] * outer(outer(outer(Beta[[1]][,i], Beta[[2]][,i]), Beta[[3]][,i]), Beta[[4]][,i])
		}
	    bic = log(fnorm(T - T.est)^2 / (d1 * d2 * d3 * d4)) + log(d1 * d2 * d3 * d4) / (d1 * d2 * d3 * d4) * df.est
	    bic2 = log(fnorm(T - T.est)^2 / (d1 * d2 * d3 * d4)) + log(d1 * d2 * d3 * d4) / (d1 * d2 * d3 * d4) * df.est/ num_m
	}

	out = list()
	out$bic = bic
	out$bic2 = bic2
	out
}



### sparse (Sun et al., 2016) and smoothed (Tibshirani et al., 2005) tensor decomposition ###
stdtruncate <- function(T, sparse_para, smooth_para, Rank, niter = 20, use_latest_update = TRUE, is_first2modes_symmetric = TRUE) {
	# T: a three-order or forth-order tensor
	# sparse_para = c(s1,...,sm): number of non-zero slements in each columns of components. If sj = dj, it means no sparsity for component j.
	# smooth_para = c(l_1, ..., l_m): keep l_j different parameters for the j-th components. if l_j = d_j, it means no smoothness in the component j.
	# Rank: rank of tensor decomposition.
	# niter: number of iterations 
	# use_latest_update: = TRUE means we use beta_1^{t}, ..., beta_{j-1}^{t} in the update of beta_j^{t}, otherwise use beta_1^{t-1},...beta_m^{t-1} in the update of beta_j^{t}.
	# is_first2modes_symmetric: = TRUE if first 2 modes are symmetric (used in real data), = FALSE otherwise. 

	T.original = T
	DD = dim(T)
	num_m = length(DD) # tensor mode. current code can handle num_m = 3 and num_m = 4.

	Beta = list()
	for(j in 1:num_m){
		Beta[[j]] = matrix(0, DD[j], Rank)
	}
	TER.ERROR = NULL
	w.tmp = NULL
	t.tmp = NULL
	W = rep(1, Rank)

	for(rrr in 1:Rank){

		print(paste("Start Estimation of Rank ", rrr, sep = ""))

		## random initialization
		set.seed(1+rrr)
		beta.new = list()
		beta.old = list()
		for(j in 1:num_m){
			beta.new[[j]] = rnorm(DD[j])
			beta.old[[j]] = rep(0, DD[j])
		}
		t = 0

		## start iteration update
		while(myerrorlist(beta.old, beta.new) >= EPS && t < niter)
		{	
			#debug
			print(paste("Start iteration ", t+1, " of Rank ", rrr, sep = ""))
			TER.ERROR = append(TER.ERROR, myerrorlist(beta.old, beta.new))
			
			t = t + 1
			beta.old = beta.new
			beta.tmp = list()
			for(j in 1:num_m){
				beta.tmp[[j]] = rep(0, DD[j])
			}
			if(num_m == 3){
				## handle 3-rd mode tensor

				if(is_first2modes_symmetric == TRUE){
					## beta_1 = beta_2 ##
					for (j in 1 : DD[2]) {
						for (k in 1 : DD[3]) {
						  beta.tmp[[1]] = beta.tmp[[1]] + beta.old[[2]][j] * beta.old[[3]][k] * T[, j, k]
						}
					}
					if(smooth_para[1]!=DD[1]){
						## apply fuselasso
						D.mat = getD1dSparse(DD[1])
						fuseoutput = fusedlasso(beta.tmp[[1]], D = D.mat)
						beta.tmp[[1]] = as.numeric(coef(fuseoutput, fuseoutput$lambda[smooth_para[1]])$beta)
					}
					if(sparse_para[1]!=DD[1]){
						## apply truncation
						beta.tmp[[1]] = mytruncate(mysd(beta.tmp[[1]]),sparse_para[1])
					}
					beta.new[[1]] = mysd(beta.tmp[[1]])
					beta.new[[2]] = mysd(beta.tmp[[1]])
				}else{
					## beta_1 ##
					for (j in 1 : DD[2]) {
						for (k in 1 : DD[3]) {
						  beta.tmp[[1]] = beta.tmp[[1]] + beta.old[[2]][j] * beta.old[[3]][k] * T[, j, k]
						}
					}
					if(smooth_para[1]!=DD[1]){
						## apply fuselasso
						D.mat = getD1dSparse(DD[1])
						fuseoutput = fusedlasso(beta.tmp[[1]], D = D.mat)
						beta.tmp[[1]] = as.numeric(coef(fuseoutput, fuseoutput$lambda[smooth_para[1]])$beta)
					}
					if(sparse_para[1]!=DD[1]){
						## apply truncation
						beta.tmp[[1]] = mytruncate(mysd(beta.tmp[[1]]),sparse_para[1])
					}
					beta.new[[1]] = mysd(beta.tmp[[1]])
					## beta_2 ##
					for (i in 1 : DD[1]) {
						for (k in 1 : DD[3]) {
							if(use_latest_update == TRUE){
								beta.tmp[[2]] = beta.tmp[[2]] + beta.new[[1]][i] * beta.old[[3]][k] * T[i, , k]
							}else{
								beta.tmp[[2]] = beta.tmp[[2]] + beta.old[[1]][i] * beta.old[[3]][k] * T[i, , k]
							}
						  
						}
					}
					if(smooth_para[2]!=DD[2]){
						## apply fuselasso
						D.mat = getD1dSparse(DD[2])
						fuseoutput = fusedlasso(beta.tmp[[2]], D = D.mat)
						beta.tmp[[2]] = as.numeric(coef(fuseoutput, fuseoutput$lambda[smooth_para[2]])$beta)
					}
					if(sparse_para[2]!=DD[2]){
						## apply truncation
						beta.tmp[[2]] = mytruncate(mysd(beta.tmp[[2]]),sparse_para[2])
					}
					beta.new[[2]] = mysd(beta.tmp[[2]])
				}

				## beta_3 ##
				for (i in 1 : DD[1]) {
					for (j in 1 : DD[2]) {
						if(use_latest_update == TRUE){
							beta.tmp[[3]] = beta.tmp[[3]] +  beta.new[[1]][i] * beta.new[[2]][j] * T[i, j,]
						}else{
							beta.tmp[[3]] = beta.tmp[[3]] +  beta.old[[1]][i] * beta.old[[2]][j] * T[i, j,]							
						}
					  
					}
				}
				if(smooth_para[3]!=DD[3]){
					## apply fuselasso
					D.mat = getD1dSparse(DD[3])
					fuseoutput = fusedlasso(beta.tmp[[3]], D = D.mat)
					beta.tmp[[3]] = as.numeric(coef(fuseoutput, fuseoutput$lambda[smooth_para[3]])$beta)
				}
				if(sparse_para[3]!=DD[3]){
					## apply truncation
					beta.tmp[[3]] = mytruncate(mysd(beta.tmp[[3]]),sparse_para[3])
				}
				beta.new[[3]] = mysd(beta.tmp[[3]])

			}else if(num_m == 4){
				## handle 4-th mode tensor

				if(is_first2modes_symmetric == TRUE){
					## beta_1 = beta_2 ##
					for (j in 1 : DD[2]) {
						for (k in 1 : DD[3]) {
							for (l in 1 : DD[4]) {
						  		beta.tmp[[1]] = beta.tmp[[1]] + beta.old[[2]][j] * beta.old[[3]][k] * beta.old[[4]][l] * T[, j, k, l]
							}
						}
					}
					if(smooth_para[1]!=DD[1]){
						## apply fuselasso
						D.mat = getD1dSparse(DD[1])
						fuseoutput = fusedlasso(beta.tmp[[1]], D = D.mat)
						beta.tmp[[1]] = as.numeric(coef(fuseoutput, fuseoutput$lambda[smooth_para[1]])$beta)
					}
					if(sparse_para[1]!=DD[1]){
						## apply truncation
						beta.tmp[[1]] = mytruncate(mysd(beta.tmp[[1]]),sparse_para[1])
					}
					beta.new[[1]] = mysd(beta.tmp[[1]])
					beta.new[[2]] = mysd(beta.tmp[[1]])
				}else{
					## beta_1 ##
					for (j in 1 : DD[2]) {
						for (k in 1 : DD[3]) {
							for (l in 1 : DD[4]) {
						  		beta.tmp[[1]] = beta.tmp[[1]] + beta.old[[2]][j] * beta.old[[3]][k] * beta.old[[4]][l] * T[, j, k, l]
							}
						}
					}
					if(smooth_para[1]!=DD[1]){
						## apply fuselasso
						D.mat = getD1dSparse(DD[1])
						fuseoutput = fusedlasso(beta.tmp[[1]], D = D.mat)
						beta.tmp[[1]] = as.numeric(coef(fuseoutput, fuseoutput$lambda[smooth_para[1]])$beta)
					}
					if(sparse_para[1]!=DD[1]){
						## apply truncation
						beta.tmp[[1]] = mytruncate(mysd(beta.tmp[[1]]),sparse_para[1])
					}
					beta.new[[1]] = mysd(beta.tmp[[1]])
					## beta_2 ##
					for (i in 1 : DD[1]) {
						for (k in 1 : DD[3]) {
							for (l in 1 : DD[4]) {
								if(use_latest_update == TRUE){
									beta.tmp[[2]] = beta.tmp[[2]] + beta.new[[1]][i] * beta.old[[3]][k] * beta.old[[4]][l] * T[i, , k, l]
								}else{
									beta.tmp[[2]] = beta.tmp[[2]] + beta.old[[1]][i] * beta.old[[3]][k] * beta.old[[4]][l] * T[i, , k, l]
								}
							}					  
						}
					}
					if(smooth_para[2]!=DD[2]){
						## apply fuselasso
						D.mat = getD1dSparse(DD[2])
						fuseoutput = fusedlasso(beta.tmp[[2]], D = D.mat)
						beta.tmp[[2]] = as.numeric(coef(fuseoutput, fuseoutput$lambda[smooth_para[2]])$beta)
					}
					if(sparse_para[2]!=DD[2]){
						## apply truncation
						beta.tmp[[2]] = mytruncate(mysd(beta.tmp[[2]]),sparse_para[2])
					}
					beta.new[[2]] = mysd(beta.tmp[[2]])
				}

				## beta_3 ##
				for (i in 1 : DD[1]) {
					for (j in 1 : DD[2]) {
						for (l in 1 : DD[4]) {
							if(use_latest_update == TRUE){
								beta.tmp[[3]] = beta.tmp[[3]] +  beta.new[[1]][i] * beta.new[[2]][j] * beta.old[[4]][l] * T[i, j, , l]
							}else{
								beta.tmp[[3]] = beta.tmp[[3]] +  beta.old[[1]][i] * beta.old[[2]][j] * beta.old[[4]][l] * T[i, j, , l]							
							}
						}
					}
				}
				if(smooth_para[3]!=DD[3]){
					## apply fuselasso
					D.mat = getD1dSparse(DD[3])
					fuseoutput = fusedlasso(beta.tmp[[3]], D = D.mat)
					beta.tmp[[3]] = as.numeric(coef(fuseoutput, fuseoutput$lambda[smooth_para[3]])$beta)
				}
				if(sparse_para[3]!=DD[3]){
					## apply truncation
					beta.tmp[[3]] = mytruncate(mysd(beta.tmp[[3]]),sparse_para[3])
				}
				beta.new[[3]] = mysd(beta.tmp[[3]])

				## beta_4 ##
				for (i in 1 : DD[1]) {
					for (j in 1 : DD[2]) {
						for (k in 1 : DD[3]) {
							if(use_latest_update == TRUE){
								beta.tmp[[4]] = beta.tmp[[4]] +  beta.new[[1]][i] * beta.new[[2]][j] * beta.new[[3]][k] * T[i, j, k, ]
							}else{
								beta.tmp[[4]] = beta.tmp[[4]] +  beta.old[[1]][i] * beta.old[[2]][j] * beta.old[[3]][k] * T[i, j, k, ]							
							}
						}
					}
				}
				if(smooth_para[4]!=DD[4]){
					## apply fuselasso
					D.mat = getD1dSparse(DD[4])
					fuseoutput = fusedlasso(beta.tmp[[4]], D = D.mat)
					beta.tmp[[4]] = as.numeric(coef(fuseoutput, fuseoutput$lambda[smooth_para[4]])$beta)
				}
				if(sparse_para[4]!=DD[4]){
					## apply truncation
					beta.tmp[[4]] = mytruncate(mysd(beta.tmp[[4]]),sparse_para[4])
				}
				beta.new[[4]] = mysd(beta.tmp[[4]])
			}

			#debug
			w.tmp = append(w.tmp, ifelse(num_m ==3, Tabc(T, as.matrix(beta.new[[1]]), as.matrix(beta.new[[2]]), as.matrix(beta.new[[3]])), Tabcd(T, as.matrix(beta.new[[1]]), as.matrix(beta.new[[2]]), as.matrix(beta.new[[3]]), as.matrix(beta.new[[4]]))))
			print(paste("W is ", w.tmp, sep = ""))			
			print(paste("End iteration ", t, " of Rank ", rrr, sep = ""))

		}
		#debug	
		t.tmp = append(t.tmp, t)
		
		for(j in 1:num_m){
			beta.new[[j]] = as.matrix(beta.new[[j]])
			Beta[[j]][,rrr] = beta.new[[j]]
		}
		W[rrr] = ifelse(num_m ==3, Tabc(T, beta.new[[1]], beta.new[[2]], beta.new[[3]]), Tabcd(T, beta.new[[1]], beta.new[[2]], beta.new[[3]], beta.new[[4]]))

		## update residual tensor
		if(num_m == 3){
			T = T - W[rrr] * outer(outer(as.numeric(beta.new[[1]]), as.numeric(beta.new[[2]])), as.numeric(beta.new[[3]]))
		}else if(num_m == 4){
			T = T - W[rrr] * outer(outer(outer(as.numeric(beta.new[[1]]), as.numeric(beta.new[[2]])), as.numeric(beta.new[[3]])),as.numeric(beta.new[[4]]))
		}
		print(paste("End Estimation of Rank ", rrr, sep = ""))
	}


	OUT = list()
	OUT$Beta = Beta
	OUT$W = W
	return(OUT)
}



#############################################################
###### Tuning via BIC for Structured tensor factorization ###
#############################################################
mytune <- function(T, rank_list = seq(1,5,1), sparse_list = seq(0.1, 0.9, 0.1), fuse_list = seq(0.1, 0.9, 0.1), tune_method = "together"){  
	# T: a tensor
	# rank_list: list of possible ranks for tuning. default rank_list = seq(1,5,1)
	# sparse_list: list of possible sparsity parametters. default sparse_list = seq(0.1, 0.9, 0.1). Do not add sparsity on the last mode!
	# fuse_list: list of possible fuse parametters. default fuse_list = seq(0.1, 0.9, 0.1). Do not add fuse on the last mode!
	# tune_method: {"sequential", "together"}. "sequential" tunes Rank --> Sparse --> Fuse. "together" tunes all three parameters together. default = "sequential"

	DD = dim(T)
	num_m = length(DD) # tensor mode. current code can handle num_m = 3 and num_m = 4.
	N = DD[num_m]
	len_rank = length(rank_list)
	len_tune_sparse = length(sparse_list)
	len_tune_fuse = length(fuse_list)
	TIME = array(0, c(len_rank, len_tune_sparse, len_tune_fuse))

	if(tune_method == "together"){
		## tune_method = "together"
		bic_value_mat = rep(NA, 4)
		for(irank in 1:len_rank){
			for(isparse in 1:len_tune_sparse){
				for(ifuse in 1:len_tune_fuse){
					tt2 = proc.time()
					Rank = rank_list[irank]
					sparse_para = c(floor(sparse_list[isparse] * DD[1:(num_m-1)]), N)
					smooth_para = c(floor(fuse_list[ifuse] * DD[1:(num_m-1)]), N)
					out_stdtruncate = stdtruncate(T, sparse_para, smooth_para, Rank, niter = 10)
					BIC_OUT = myBIC(out_stdtruncate$W, out_stdtruncate$Beta, T)
					bic_value_mat = rbind(bic_value_mat, c(BIC_OUT$bic, Rank, sparse_list[isparse], fuse_list[ifuse]))
					TIME[irank, isparse, ifuse] = ((proc.time() - tt2)[3])
				}
			}
		}
		time.rank = sum(apply(TIME, 1, mean))
		time.sparse = sum(apply(TIME, c(1,2), mean)) - time.rank
		time.fuse = sum(apply(TIME, c(1,3), mean)) - time.rank
		bic_value_mat = bic_value_mat[-1,]
		colnames(bic_value_mat) = c("bic", "Rank", "Sparse", "Fuse")
		min_BIC_index = which.min(bic_value_mat[,1])
		Rank.opt = bic_value_mat[min_BIC_index, 2]
		sparse.opt = c(floor(bic_value_mat[min_BIC_index, 3] * DD[1:(num_m-1)]), N) 
		fuse.opt = c(floor(bic_value_mat[min_BIC_index, 4] * DD[1:(num_m-1)]), N)

	}else{

		## tune_method = "sequential"
		tt1 = proc.time()
		bic_value_mat_rank = rep(NA, 2)
		for(irank in 1:len_rank){
			Rank = rank_list[irank]
			out_stdtruncate = stdtruncate(T, sparse_para = DD, smooth_para = DD, Rank, niter = 10)
			BIC_OUT = myBIC(out_stdtruncate$W, out_stdtruncate$Beta, T)
			bic_value_mat_rank = rbind(bic_value_mat_rank, c(BIC_OUT$bic, Rank))
		}
		bic_value_mat_rank = bic_value_mat_rank[-1,]
		min_BIC_index = which.min(bic_value_mat_rank[,1])
		Rank.opt = bic_value_mat_rank[min_BIC_index, 2]
		time.rank = (proc.time() - tt1)[3]

		tt1 = proc.time()		
		bic_value_mat_fuse = rep(NA, 2)
		for(ifuse in 1:len_tune_fuse){
			smooth_para = c(floor(fuse_list[ifuse] * DD[1:(num_m-1)]), N)
			out_stdtruncate = stdtruncate(T, sparse_para = DD, smooth_para, Rank.opt, niter = 10)
			BIC_OUT = myBIC(out_stdtruncate$W, out_stdtruncate$Beta, T)
			bic_value_mat_fuse = rbind(bic_value_mat_fuse, c(BIC_OUT$bic, fuse_list[ifuse]))
		}
		bic_value_mat_fuse = bic_value_mat_fuse[-1,]
		min_BIC_index = which.min(bic_value_mat_fuse[,1])
		fuse.opt = c(floor(bic_value_mat_fuse[min_BIC_index, 2] * DD[1:(num_m-1)]), N)
		time.fuse = (proc.time() - tt1)[3]
		
		tt1 = proc.time()
		bic_value_mat_sparse = rep(NA, 2)
		for(isparse in 1:len_tune_sparse){
			sparse_para = c(floor(sparse_list[isparse] * DD[1:(num_m-1)]), N)
			out_stdtruncate = stdtruncate(T, sparse_para, smooth_para = fuse.opt, Rank.opt, niter = 10)
			BIC_OUT = myBIC(out_stdtruncate$W, out_stdtruncate$Beta, T)
			bic_value_mat_sparse = rbind(bic_value_mat_sparse, c(BIC_OUT$bic, sparse_list[isparse]))
		}
		bic_value_mat_sparse = bic_value_mat_sparse[-1,]
		min_BIC_index = which.min(bic_value_mat_sparse[,1])
		sparse.opt = c(floor(bic_value_mat_sparse[min_BIC_index, 2] * DD[1:(num_m-1)]), N) 
		time.sparse = (proc.time() - tt1)[3]

		bic_value_mat = rbind(bic_value_mat_rank, bic_value_mat_fuse, bic_value_mat_sparse)
	}

	P <- list()
	P$Rank.opt = Rank.opt
	P$sparse.opt = sparse.opt
    P$fuse.opt = fuse.opt
    P$time.rank = time.rank
    P$time.fuse = time.fuse
    P$time.sparse = time.sparse
    P$bic_value_mat = bic_value_mat
    return(P)
}




###############################################################
###### Our main function of Dynamic Tensor Clustering (DTC) ###
###############################################################
mydtc <- function(out_stdtruncate, K){  
    ## out_stdtruncate: object of sparse tensor decomposition. it is a list consisting of Beta and W.
    ## K: number of clusters in the 3-rd mode for a third-mode tensor; if length(K) == 2, refer to # of clusters in the third and forth mode for a forth order tensor.

	Beta = out_stdtruncate$Beta
	W = out_stdtruncate$W

	if(ncol(Beta[[3]]) == 1 && ncol(Beta[[1]]) > 1 && K[1] == 1){
		## window_size = 1 and we do not need to cluster over time. Need to transfer Beta[[3]] from a column vec to a row vec 
		Beta[[3]] = t(Beta[[3]])
	}

	## Run clustering based on decomposed components
	num_m = length(Beta) # tensor mode. current code can handle num_m = 3 and num_m = 4.
	DD = NULL
	for(ii in 1:num_m){
		DD = append(DD, nrow(Beta[[ii]]))
	}
	Rank = ncol(Beta[[1]])
	membership = list()
	if(num_m == 3){
    	out.kmeans = kmeans(Beta[[3]], K, nstart = 100)
    	membership = out.kmeans$cluster
	}else if(num_m == 4 & length(K) == 1){
		out.kmeans1 = kmeans(Beta[[4]], K, nstart = 100)
		membership = out.kmeans1$cluster
	}else if(num_m == 4 & length(K) == 2){
		out.kmeans1 = kmeans(Beta[[3]], K[1], nstart = 100)
		out.kmeans2 = kmeans(Beta[[4]], K[2], nstart = 100)
		membership$mode3 = out.kmeans1$cluster
		membership$mode4 = out.kmeans2$cluster
	}

	## Computer new cluster centers via components
	centers = list()
	if(num_m == 3){
    	for(k in 1:K){
    		tmp = matrix(0, DD[1], DD[2])
    		for(rr in 1:Rank){
    			tmp = tmp + outer(Beta[[1]][,rr], Beta[[2]][,rr]) * W[rr] * mean(Beta[[3]][which(membership == k),rr])
    		}
    		centers[[k]] = tmp
    	}
	}else if(num_m == 4 & length(K) == 1){
		for(k in 1:K){
    		tmp = array(0, c(DD[1], DD[2], DD[3]))
    		for(rr in 1:Rank){
    			tmp = tmp + outer(outer(Beta[[1]][,rr], Beta[[2]][,rr]), Beta[[3]][,rr])* W[rr] * mean(Beta[[4]][which(membership == k),rr])
    		}
    		centers[[k]] = tmp
		}
	}else if(num_m == 4 & length(K) == 2){
		for(k1 in 1:K[1]){
			for(k2 in 1:K[2]){
				index = (k1-1)*K[2] + k2
	    		tmp = matrix(0, DD[1], DD[2])
	    		for(rr in 1:Rank){
	    			tmp = tmp + outer(Beta[[1]][,rr], Beta[[2]][,rr]) * W[rr] * mean(Beta[[3]][which(membership$mode3 == k1),rr]) * mean(Beta[[4]][which(membership$mode4 == k2),rr])
	    		}
	    		centers[[index]] = tmp
			}
		}
	}

	P <- list()
	P$W = W
	P$Beta = Beta
    P$centers = centers
    P$membership = membership
    return(P)
}






##################################################################################
###### Alternative Method 1: Kmeans clustering on the original tensor space ######
##################################################################################
tenkmeans <- function(data, K, n.iter = 20, EPS = 1e-2, seed = 1){  
    ## k-means clustering for tensor-variate data
    ## data: contains all X1, ..., Xn each of dimension m1 * m2 * m3
    ## K: number of clusters
    ## n.iter: maximal number of iterations allowed.
    ## EPS: tolerance of distance between old and new centers for iteration termination.

    nmode = length(dim(data)) - 1

    set.seed(seed)
    ## initialization via random tensors ##   
    if(nmode == 2){
    ## matrix samples
        m1 = dim(data)[1]
        m2 = dim(data)[2]
        n  = dim(data)[3]
        center_ind = sample(1:n, K)
        cluster_center = data[,,center_ind]
        cluster_center.old = array(0, c(m1,m2,K))
    }else if(nmode == 3){
    ## tensor samples
        m1 = dim(data)[1]
        m2 = dim(data)[2]
        m3 = dim(data)[3]
        n  = dim(data)[4]
        center_ind = sample(1:n, K)
        cluster_center = data[,,,center_ind]
        cluster_center.old = array(0, c(m1,m2,m3,K))
    }
      
    t = 0
    #debug_distance = fnorm(cluster_center - cluster_center.old)
    #debug_distance_normalized = fnorm(cluster_center - cluster_center.old)/fnorm(cluster_center)
    while(fnorm(cluster_center - cluster_center.old)/fnorm(cluster_center) >= EPS & t < n.iter)
    {  
        t = t + 1
        cluster_center.old = cluster_center
       
        ## update cluster assignment
        DIST = matrix(0, n, K)
        L.mat = matrix(0, n, K)
        for(i in 1:n){
            for(k in 1:K){
                if(nmode == 2){
                    DIST[i,k] = fnorm(data[,,i] - cluster_center[,,k])
                }else if(nmode == 3){
                    DIST[i,k] = fnorm(data[,,,i] - cluster_center[,,,k])
                }
            }
            L.mat[i,which.min(DIST[i,])] = 1
        }

        ## update cluster centers
        for(k in 1:K){

            if(nmode == 2){
                cluster_center[,,k] = apply(data[,,which(L.mat[,k] == 1)], c(1,2), mean)
            }else if(nmode == 3){
                cluster_center[,,,k] = apply(data[,,,which(L.mat[,k] == 1)], c(1,2,3), mean)
            }
        }

        #debug_distance = append(debug_distance, fnorm(cluster_center - cluster_center.old))
        #debug_distance_normalized = append(debug_distance_normalized, fnorm(cluster_center - cluster_center.old)/fnorm(cluster_center))

    }   
	
	P <- list()
    P$cluster_center = cluster_center
    P$membership <- apply(L.mat, 1, function(x) which(x==1))
    return(P)
}




####################################################################################
###### Alternative Method 2: vectorize tensor first, and then kmeans clustering  ###
####################################################################################
veckmeans <- function(data, K){  
    ## vectorize tensor, then k-means clustering
    ## data: contains all X1, ..., Xn each of dimension m1 * m2 * m3
    ## K: number of clusters

    nmode = length(dim(data)) - 1

    if(nmode == 2){
    ## matrix samples
        Comp.data = t(apply(DATA,3,as.vector))

    }else if(nmode == 3){
    ## tensor samples
    	Comp.data = t(apply(DATA,4,as.vector))
    }

    out.kmeans = kmeans(Comp.data, K, nstart = 100)

	P <- list()
    P$cluster_center = out.kmeans$centers
    P$membership <- out.kmeans$cluster
    return(P)
}




