

#
#	Modified by Seunggeun Lee - Ver 0.1
#
KMTest.logistic.Linear.VarMatching = function(res, Z, X1, kernel, weights = NULL, pi_1, method, res.out,n.Resampling, r.corr, mu, res.moments = NULL){

	n<-length(pi_1)
  	D  = diag(pi_1)   

	# Weighted Linear Kernel 
	if (kernel == "linear.weighted") {
	    Z = t(t(Z) * (weights))
	}

	# r.corr
	if(r.corr == 1){
		Z<-cbind(rowSums(Z))
	} else if(r.corr > 0){

		p.m<-dim(Z)[2]	
		R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Z<- Z %*% t(L) 
	}

  	Q.Temp = t(res)%*%Z
  	Q = Q.Temp %*% t(Q.Temp)/2

  	Q.res = NULL
  	if(n.Resampling > 0){
  		Q.Temp.res = t(res.out)%*%Z
  		Q.res = rowSums(rbind(Q.Temp.res^2))/2
  	}
  	Z1 = (Z * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z * pi_1))

	Q.sim = NULL
	if(!is.null(res.moments)){

		Q.Temp.res1 = t(cbind(res.moments))%*%Z
  		Q.sim = rowSums(rbind(Q.Temp.res1^2))/2

	} 

	Q.all<-c(Q,Q.res)
	p_all<-mu

	type = "Other"
	if(method =="ECP"){
		type = "OnlySim"
	}
	re<-SKAT_PValue_Logistic_VarMatching(Q.all, Z1 /sqrt(2), p_all, Q.sim, type)

	
	# re$p.value is p-values of aSKAT
	
	p.value.resampling = NULL
	p.value.noadj.resampling = NULL

	
	p.value= re$p.value[1]
	if(length(Q.all) > 1){
		p.value.resampling<-re$p.value[-1]
	}

	p.value.noadj= re$p.value.noadj[1]
	if(length(Q.all) > 1){
	
		p.value.noadj.resampling<-re$p.value.noadj[-1]
	}

	re<-list(p.value = p.value, p.value.resampling = p.value.resampling
, p.value.noadj = p.value.noadj, p.value.noadj.resampling = p.value.noadj.resampling
, Test.Type = method, Q = Q,  Q.resampling = Q.res, param=NULL)  
  
  	return(re)
}




SKAT_PValue_Logistic_VarMatching<-function(Q, Z1, p_all, Q.sim, type="Other"){


	param<-SKAT_Logistic_VarMatching_GetParam1(Z1, p_all, Q.sim, type)

	# aSKAT pvalue
	Q.Norm<-(Q - param$muQ)/sqrt(param$varQ)
	Q.Norm1<-Q.Norm * sqrt(2*param$df) + param$df
	p.value<- 1-pchisq(Q.Norm1,  df = param$df, ncp=0)


	# SKAT pvalue
	param.noadj<-param$param.noadj
	Q.Norm<-(Q - param.noadj$muQ)/param.noadj$sigmaQ
	Q.Norm1<-Q.Norm * param.noadj$sigmaX + param.noadj$muX
	p.value.noadj<- 1-pchisq(Q.Norm1,  df = param.noadj$l,ncp=0)

	out<-list(p.value=p.value, p.value.noadj=p.value.noadj, param=param)	

	return(out)

}

SKAT_GET_skewness <-  function(x) {
	m3 <- mean((x-mean(x))^3)
	skew <- m3/(sd(x)^3)
	return(skew)
}


SKAT_GET_kurtosis <- function(x) {  

	if(sd(x) == 0){
		return(-100)
	}
	m4 <- mean((x-mean(x))^4) 
	kurt <- m4/(sd(x)^4)-3  
	kurt
}


SKAT_Get_DF_Sim<-function(Q.sim){


	s2.sim<-SKAT_GET_kurtosis(Q.sim)
	df.sim<-12/s2.sim

	if(s2.sim <= 0){
		
		df.sim=100000
	} else if(df.sim < 0.01 ){
		s1.sim<-SKAT_GET_skewness(Q.sim)
		df.sim<-8/s1.sim^2
	}

	return(df.sim)
}

SKAT_Logistic_VarMatching_GetParam1<-function(Z1, p_all, Q.sim, type="Other"){


	if(type != "OnlySim"){

		try1<-try(Get_Lambda_U_From_Z(Z1),silent = TRUE)
		if(class(try1) == "try-error"){
			type="OnlySim"
		} else {
			out.svd = try1
			lambda<-out.svd$lambda
    			U<-out.svd$U
			param<-SKAT_Logistic_VarMatching_GetParam(lambda, U, p_all, Q.sim)
		}

	}

	if(type == "OnlySim"){
		param<-SKAT_Logistic_VarMatching_GetParam1_OnlySim(Z1, p_all, Q.sim)
		
	} 

	
	return(param)

}


SKAT_Logistic_VarMatching_GetParam1_OnlySim<-function(Z1, p_all, Q.sim){


	out.svd = Get_Lambda_U_From_Z(Z1)
	lambda<-out.svd$lambda
	
	muQ = sum(lambda)
	varQ.sim<-var(Q.sim)

	#print(c(varQ, varQ.sim))
	df.sim<-SKAT_Get_DF_Sim(Q.sim)

	# No adjustment
	c1<-rep(0,4)	
	for(i in 1:4){
		c1[i]<-sum(lambda^i)
	}	
	param<-Get_Liu_Params_Mod(c1)

	return(list(muQ = muQ, varQ = varQ.sim, df=df.sim, lambda.new=NULL, param.noadj = param))

}


#
#	lambda : eigenvalues
#	U : eigenvectors
#	p_all : 
#	If Q.resample == NULL, it does not estimate kurtosis and df.sim < 0
SKAT_Logistic_VarMatching_GetParam<-function(lambda, U, p_all, Q.sim){

	# Var match
	re<-SKAT_Get_Cov_Param(lambda, p_all, U)

	# New Lambda 
	lambda.new<- re$lambda.new
	
	# new parameters
	muQ<-re$muQ

	# new var
	varQ<-re$varQ

	# df
	s2 = sum(lambda.new^4) / sum(lambda.new^2)^2
	df<-1/s2

	if(!is.null(Q.sim)){
		df<-SKAT_Get_DF_Sim(Q.sim)
	}

	
	# No adjustment
	c1<-rep(0,4)	
	for(i in 1:4){
		c1[i]<-sum(lambda^i)
	}	
	param<-Get_Liu_Params_Mod(c1)

	return(list(muQ = muQ, varQ = varQ, df=df, lambda.new=lambda.new, param.noadj = param))

}


SKAT_Get_Var_Elements<-function(m4,p_all,u1,u2){

	temp1<-u1^2 * u2^2


	a1<-sum(m4 * temp1)
	a2<-sum(u1^2) * sum(u2^2) - sum(temp1)
	a3<-sum(u1*u2)^2 - sum(temp1)

	a3<-a3*2
	
	
	a1+a2+a3 
}


SKAT_Get_Cov_Param<-function(lambda,p_all,U){

	#p_all<-obj$mu
	#U<-out$U
	#lambda<-out$lambda

	p.m<-length(lambda)
	m4<-p_all*(1-p_all)*(3*p_all^2-3*p_all +1) / (p_all*(1-p_all))^2
		
	zeta<-rep(0,p.m)
	var_i<-rep(0,p.m)
	varQ<-0	

	for(i in 1:p.m){
		temp.M1<-sum(U[,i]^2)^2 - sum(U[,i]^4) 
		zeta[i]<-sum(m4 * U[,i]^4) + 3* temp.M1 # because ( \sum .)^4, not ^2
		var_i[i]<-zeta[i] - 1	
	}

	if(p.m == 1){
		Cov_Mat<-matrix(zeta* lambda^2, ncol=1,nrow=1)
	} else if(p.m > 1){

		Cov_Mat<-diag(zeta* lambda^2)
		for(i in 1:(p.m-1)){
			for(j in (i+1):p.m){
				Cov_Mat[i,j]<-SKAT_Get_Var_Elements(m4,p_all,U[,i],U[,j])
				Cov_Mat[i,j]<-Cov_Mat[i,j]* lambda[i]* lambda[j]
			}
		}
	}
	
	Cov_Mat<-Cov_Mat + t(Cov_Mat)
	diag(Cov_Mat)<-diag(Cov_Mat)/2

	varQ<-sum(Cov_Mat) - sum(lambda)^2	
	muQ=sum(lambda)
	lambda.new<-lambda * sqrt(var_i)/sqrt(2)
	return(list(zeta=zeta, var_i=var_i, varQ = varQ, muQ=muQ, lambda.new=lambda.new))
	
}


