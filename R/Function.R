
#
# Get Parameter for the Liu et. al 
#

Get_Satterthwaite<-function(muQ, varQ){

	a1<-varQ/muQ /2
	a2<-muQ/a1

	re<-list(df=a2, a=a1)
	return(re)
}

Get_Liu_Params<-function(c1){
  ## Helper function for getting the parameters for the null approximation
  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    a = 1/s1
    d = 0
    l = 1/s1^2
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

Get_Liu_Params_Mod<-function(c1){
  ## Helper function for getting the parameters for the null approximation
  
  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}


Get_Liu_Params_Mod_Lambda<-function(lambda, df1=NULL){
  ## Helper function for getting the parameters for the null approximation

	if(is.null(df1)){
		df1=rep(1, length(lambda))
	}
  	c1<-rep(0,4)
  	for(i in 1:4){
  		# changed... multiply df1
		c1[i]<-sum(lambda^i * df1)
  	}

  	muQ<-c1[1]
  	sigmaQ<-sqrt(2 *c1[2])
  	s1 = c1[3] / c1[2]^(3/2)
  	s2 = c1[4] / c1[2]^2

  	beta1<-sqrt(8)*s1
  	beta2<-12*s2
  	type1<-0

  	#print(c(s1^2,s2))
  	if(s1^2 > s2){
    	a = 1/(s1 - sqrt(s1^2 - s2))
    	d = s1 *a^3 - a^2
    	l = a^2 - 2*d
  	} else {
    	type1<-1
    	l = 1/s2
    	a = sqrt(l)
    	d = 0
  	}
  	muX <-l+d
  	sigmaX<-sqrt(2) *a

  	re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  	return(re)
}

Get_Liu_PVal<-function(Q, W, Q.resampling = NULL){
	##added by Zhangchen, for sparse matrix, 12.17.2018
    	Q=as.matrix(Q)
	W=as.matrix(W)
	if (length(Q.resampling)>0){Q.resampling=as.matrix(Q.resampling)}

	Q.all<-c(Q,Q.resampling)

	A1<-W/2
	A2<-A1 %*% A1

	c1<-rep(0,4)
	c1[1]<-sum(diag(A1))
	c1[2]<-sum(diag(A2))
	c1[3]<-sum(A1*t(A2))
	c1[4]<-sum(A2*t(A2))
	param<-Get_Liu_Params(c1)

	Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
	Q.Norm1<-Q.Norm * param$sigmaX + param$muX
	p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE)

	p.value.resampling = NULL
	if(length(Q.resampling) > 0){
		p.value.resampling<-p.value[-1]
	}

	re<-list(p.value = p.value[1], param=param, p.value.resampling = p.value.resampling )  
	return(re)
}

Get_Liu_PVal.MOD<-function(Q, W, Q.resampling = NULL){
	##added by Zhangchen, for sparse matrix, 12.17.2018
    	Q=as.matrix(Q)
	W=as.matrix(W)
	if (length(Q.resampling)>0){Q.resampling=as.matrix(Q.resampling)}
	
	Q.all<-c(Q,Q.resampling)

	A1<-W/2
	A2<-A1 %*% A1

	c1<-rep(0,4)
	c1[1]<-sum(diag(A1))
	c1[2]<-sum(diag(A2))
	c1[3]<-sum(A1*t(A2))
	c1[4]<-sum(A2*t(A2))
	param<-Get_Liu_Params_Mod(c1)

	Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
	Q.Norm1<-Q.Norm * param$sigmaX + param$muX
	p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE)

	p.value.resampling = NULL
	if(length(Q.resampling) > 0){
		p.value.resampling<-p.value[-1]
	}

	re<-list(p.value = p.value[1], param=param, p.value.resampling = p.value.resampling ) 

	return(re)
}

Get_Liu_PVal.MOD.Lambda<-function(Q.all, lambda, df1=NULL, log.p=FALSE){

	param<-Get_Liu_Params_Mod_Lambda(lambda, df1=df1)

	Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
	Q.Norm1<-Q.Norm * param$sigmaX + param$muX
	p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)

	return(p.value)

}


Get_Liu_PVal.MOD.Lambda.Zero<-function(Q, muQ, muX, sigmaQ, sigmaX, l, d){


	Q.Norm<-(Q - muQ)/sigmaQ
	Q.Norm1<-Q.Norm * sigmaX + muX
	
	temp<-c(0.05,10^-10, 10^-20,10^-30,10^-40,10^-50, 10^-60, 10^-70, 10^-80, 10^-90, 10^-100)
	#qchisq(temp, df=1000000000,lower.tail=FALSE)	
	out<-qchisq(temp,df = l,ncp=d, lower.tail=FALSE)
	#cat(c(Q.Norm1,l,d, out))
	#cat("\n")
	IDX<-max(which(out < Q.Norm1))
	
	pval.msg<-sprintf("Pvalue < %e", temp[IDX])
	return(pval.msg)

}


Get_Davies_PVal<-function(Q, W, Q.resampling = NULL, isFast=FALSE){
    
    ##added by Zhangchen, for sparse matrix, 12.17.2018
    Q=as.matrix(Q)
	W=as.matrix(W)
	if (length(Q.resampling)>0){Q.resampling=as.matrix(Q.resampling)}
	
	K<-W/2
	
	Q.all<-c(Q,Q.resampling)

	re<-Get_PValue(K,Q.all, isFast=isFast)
	param<-list()
	param$liu_pval<-re$p.val.liu[1]
	param$Is_Converged<-re$is_converge[1]


	p.value.resampling = NULL
	if(length(Q.resampling) > 0){
		p.value.resampling<-re$p.value[-1]
		param$liu_pval.resampling<-re$p.val.liu[-1]
		param$Is_Converged.resampling<-re$is_converge[-1]

	}
	

	re<-list(p.value = re$p.value[1], param=param,p.value.resampling = p.value.resampling
	, pval.zero.msg=re$pval.zero.msg )  
	return(re)
}



Get_Lambda_Approx<-function(K, maxK=100){

	#maxK=200
	lambda<-NULL
	p1<-ncol(K)

	lambda_sum = sum(diag(K))
	lambda2_sum = sum(K*K)
		
	lambda<-rep(0, p1)
	k1<-floor(p1/10)
	if(k1 > maxK){
		k1=maxK
	}
	
	out.s<-RSpectra::eigs_sym(K, k=k1, which = "LM")
	lambda1<-out.s$values
		
	IDX1<-which(lambda1 >= 0)
	IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)
	if(length(IDX2) == 0){
		stop("No Eigenvalue is bigger than 0!!")
	} 
		
	if(length(lambda) <=k1){
		lambda1<-lambda1[IDX2]
	}  
		
	lambda_sum = lambda_sum - sum(lambda1)
	lambda2_sum = lambda2_sum - sum(lambda1^2)
		
	
	df1 = lambda_sum^2/ lambda2_sum
	df1 = round(df1)
	
	# make df1 as integer
	
	lambda_last = lambda_sum/df1
	
	lambda = c(lambda1, lambda_last)
	
	df1 = c(rep(1, length(lambda1)), df1)
	out_lambda = list(lambda=lambda, df1=df1)
	
	
	return(out_lambda)

}


Get_Lambda<-function(K, isFast=FALSE, maxK=100){

	lambda<-NULL

	out.s<-eigen(K,symmetric=TRUE, only.values = TRUE)

	lambda1<-out.s$values
	IDX1<-which(lambda1 >= 0)

	# eigenvalue bigger than sum(eigenvalues)/1000
	IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)
	
	if(length(IDX2) == 0){
		stop("No Eigenvalue is bigger than 0!!")
	}
	lambda<-lambda1[IDX2]

	return(lambda)

}




Get_Lambda_U_From_Z<-function(Z1){
	
	if(dim(Z1)[2]==1){
		Is.OnlyOne = TRUE
		lambda<-sum(Z1^2)
		U<-Z1/sqrt(lambda)		
		return(list( lambda = lambda, U = cbind(U)))
	}

	#########################################
	#try1<-try(svd(Z1, LINPACK = TRUE),silent = TRUE)
	#
	#if(class(try1) == "try-error"){
	#	# try LAPACK
	#	try1<-try(svd(Z1, LINPACK = FALSE),silent = TRUE)
	#}

	try1<-try(svd(Z1),silent = TRUE)
	
	if(Is_TryError(try1)){
		stop("SVD error!");
	} else {
		out.svd = try1
	}
	
	lambda.org<-out.svd$d^2
	IDX<-which(lambda.org > mean(lambda.org)/100000)
	if(length(IDX) <= 1){
		Is.OnlyOne = TRUE
	}

	if(length(IDX) == 0){
		return(list(lambda=NULL, U=NULL))
	}
	return(list( lambda = lambda.org[IDX], U = cbind(out.svd$u[,IDX])))
}


Get_PValue<-function(K,Q, isFast=FALSE){
	
	
	df1=NULL
	p1<-ncol(K)
	if(!isFast || p1 < 2000){
		lambda<-Get_Lambda(K)
		re<-Get_PValue.Lambda(lambda,Q)
		
	} else {

		out_lambda=Get_Lambda_Approx(K)
		lambda = out_lambda$lambda
		df1=out_lambda$df1
		re<-Get_PValue.Lambda(lambda,Q, df1=df1)
	}
	
	return(re)
}


Get_PValue.Lambda<-function(lambda,Q, df1=NULL){
	
	#print(lambda)
	n1<-length(Q)

	p.val<-rep(0,n1)
	p.val.liu<-rep(0,n1)
	is_converge<-rep(0,n1)
	p.val.liu<-Get_Liu_PVal.MOD.Lambda(Q, lambda, df1)
	
	for(i in 1:n1){
		
		if(is.null(df1)){
			out<-SKAT_davies(Q[i],lambda, acc=10^(-6))
		} else {
			out<-SKAT_davies(Q[i],lambda, h=df1, acc=10^(-6))
		}
		p.val[i]<-out$Qq
		#p.val.liu[i]<-SKAT_liu(Q[i],lambda)

		is_converge[i]<-1
		
		# check convergence
		if(length(lambda) == 1){
			p.val[i]<-p.val.liu[i]
		} else if(out$ifault != 0){
			is_converge[i]<-0
		}
	
		# check p-value
		if(p.val[i] > 1 || p.val[i] <= 0 ){
			is_converge[i]<-0
			p.val[i]<-p.val.liu[i]
		}
	}
	
	p.val.msg = NULL
	p.val.log=NULL
	#cat(p.val[1])
	if(p.val[1] == 0){

		param<-Get_Liu_Params_Mod_Lambda(lambda, df1)
		p.val.msg<-Get_Liu_PVal.MOD.Lambda.Zero(Q[1], param$muQ, param$muX, param$sigmaQ, param$sigmaX, param$l, param$d)
		p.val.log<-Get_Liu_PVal.MOD.Lambda(Q[1], lambda, log.p=TRUE)[1]

	}

	return(list(p.value=p.val, p.val.liu=p.val.liu, is_converge=is_converge, p.val.log=p.val.log, pval.zero.msg=p.val.msg))

}



SKAT_Get_MAF<-function(Z,id_include=NULL, Is.chrX=FALSE, SexVar=NULL){
	
	nSNP = ncol(Z)
	
	if(Is.chrX==FALSE){
		if(is.null(id_include)){
			MAF<-colMeans(Z, na.rm = TRUE)/2
		} else {
			MAF<-colMeans(as.matrix(Z[id_include,]),na.rm=TRUE)/2
		}
	} else {
		if(is.null(SexVar)){
			stop("Error SexVar!")
		}
		
		id.male<-which(SexVar==1)
		id.female<-which(SexVar==2)
		
		if(!is.null(id_include)){
			id.male<-intersect(id.male, id_include)
			id.female<-intersect(id.female, id_include)
		}

		n.male.v= 0
		n.female.v=0
		
		if(length(id.male)>0){
		  Z.nomiss<-matrix(1, length(id.male), ncol(Z))
		  id.na<-which(is.na(cbind(Z[id.male,])))
		  Z.nomiss[id.na]<-0
		  n.male.v<-colSums(Z.nomiss)
		}
		if(length(id.female)>0){
		  Z.nomiss<-matrix(1, length(id.female), ncol(Z))
		  id.na<-which(is.na(cbind(Z[id.female,])))
		  Z.nomiss[id.na]<-0
		  
		  n.female.v<-colSums(Z.nomiss)
		}
		id.all<-union(id.male, id.female)
		MAF<-colSums(cbind(Z[id.all,]), na.rm = TRUE)/(2*n.female.v+n.male.v)
				
	}	
	
	
	return(MAF)
}



# Simple Imputation
# Z : an n x p genotype matrix with n samples and p SNPs
# Missing : a missing genotype value. Default is 9

Impute<-function(Z, impute.method){
	
	p<-dim(Z)[2]

	if(impute.method =="random"){
		for(i in 1:p){
			IDX<-which(is.na(Z[,i]))
			if(length(IDX) > 0){
				maf1<-mean(Z[-IDX,i])/2
				Z[IDX,i]<-rbinom(length(IDX),2,maf1)
			}
		}
	} else if(impute.method =="fixed"){
		for(i in 1:p){
			IDX<-which(is.na(Z[,i]))
			if(length(IDX) > 0){
				maf1<-mean(Z[-IDX,i])/2
				Z[IDX,i]<-2 * maf1
			}
		}
	} else if(impute.method =="bestguess") {
		
		for(i in 1:p){
			IDX<-which(is.na(Z[,i]))
			if(length(IDX) > 0){
				maf1<-mean(Z[-IDX,i])/2
				Z[IDX,i]<-round(2 * maf1)
			}
		}
	
	} else {
		stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
	}

	return(Z)
}


##################################################
# Get polymorphic SNP
SKAT_Get_Polymorphic_SNP<-function(Z){

	temp<-apply(Z,2,var)
	ID<-which(temp == 0)
	return(ID)
}

###########################################
#
# Functions related to weights


# Get Beta Weights
# Z : an n x p genotype matrix with n samples and p SNPs

Beta.Weights<-function(MAF,weights.beta){

	n<-length(MAF)
	weights<-rep(0,n)	
	IDX_0<-which(MAF == 0)
	if(length(IDX_0) == n){
		stop("No polymorphic SNPs")
	} else if( length(IDX_0) == 0){
		weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
	} else {
		weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
	}

	
	#print(length(IDX_0))
	#print(weights[-IDX_0])
	return(weights)
	
}



Get_MAF<-function(Z){

	is.missing<-which(Z == 9)
	Z[is.missing]<-NA

	maf<-colMeans(Z,na.rm = TRUE)/2
	return(maf)

}


Get_Logistic_Weights_MAF<-function(MAF,par1= 0.07, par2=150){

	n<-length(MAF)
	weights<-rep(0,n)	
	IDX<-which(MAF > 0)
	if(length(IDX) == 0){
		stop("No polymorphic SNPs")
	} else {

		x1<-(MAF[IDX] - par1) * par2
		weights[IDX]<-exp(-x1)/(1+exp(-x1))
	} 

	return(weights)
	
}



Get_Logistic_Weights<-function(Z,par1=0.07, par2=150){

	MAF<-Get_MAF(Z)
	re<-Get_Logistic_Weights_MAF(MAF,par1,par2)
	return(re)
}




Get_Matrix_Square.1<-function(A){
	
	out<-eigen(A,symmetric=TRUE)
	ID1<-which(out$values > 0)
	if(length(ID1)== 0){
		stop("Error to obtain matrix square!")
	}
	out1<-t(out$vectors[,ID1]) * sqrt(out$values[ID1])
	return(out1)
}

Get_Resampling_Bin<-function(ncase, prob, n.Resampling){

	n<-length(prob)
	err<-0
	temp<-.C("SL_Binary_Boot", as.integer(n), as.integer(n.Resampling), 
	as.integer(ncase), as.double(prob), integer(n), integer(n), integer( n * n.Resampling), as.integer(err) )
	
	if(temp[[8]] != 1){
		return(NULL)
	}
	
	out<-matrix(temp[[7]],byrow=FALSE, ncol=n.Resampling)
	return(out)

}

Check_Class<-function(obj, class_type){
  re<-TRUE
  
  #change to use inherits
  #if(!any(class(obj) %in% class_type)){
  if(!inherits(obj, class_type)){
    re<-FALSE
  }
  return(re)
}

Is_TryError<-function(obj){
  
  re<-Check_Class(obj, "try-error")
  return(re)
}

# obj<-Z; class_type=c("matrix", "temp1"); Check_Class(obj, class_type)




#
# x is either y or SKAT_NULL_Model 
#
SKAT.SSD.OneSet_SetIndex_OLD = function(SSD.INFO, SetIndex, obj, ..., obj.SNPWeight=NULL){
  
  id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
  #id1 = SetIndex
  if(length(id1) == 0){
    MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
    stop(MSG)
  }	
  SetID<-SSD.INFO$SetInfo$SetID[id1]
  
  is_ID = FALSE
  if(!is.null(obj.SNPWeight)){
    is_ID = TRUE
  }
  try1<-try(Get_Genotypes_SSD(SSD.INFO, SetIndex, is_ID=is_ID),silent = TRUE)
  if(!Is_TryError(try1)){
    Z<-try1
    Is.Error<-FALSE	
  } else {
    err.msg<-geterrmessage()
    msg<-sprintf("Error to get genotypes of %s: %s",SetID, err.msg)
    stop(msg)
  }
  
  
  if(is.null(obj.SNPWeight)){
    
    re<-SKAT(Z, obj, ...)
  } else {
    
    SNP_ID<-colnames(Z)
    p<-ncol(Z)
    weights<-rep(0, p)
    for(i in 1:p){
      val1<-SNP_ID[i]			
      val2<-obj.SNPWeight$hashset[[val1]]
      
      if(is.null(val2)){
        msg<-sprintf("SNP %s is not found in obj.SNPWeight!", val1)
        stop(msg)
      }
      
      weights[i]<-val2
    }
    re<-SKAT(Z, obj, weights=weights, ...)
  }
  
  return(re)
}

