#
#	Method.Bin = "Hybrid", "ER", "Adaptive", "QuantileAdj, "MomentAdj", "NoAdj"
#
#		Hybrid:
#	
#		Adaptive: 10^5 initial, max to 10^7
#		QuantileAdj: Only effective when MAF > 20, if not warning message
#		MomentAdj: any ....
#		NoAdj: NoAdj


SKATBinary<-function(Z, obj, kernel = "linear.weighted", method="SKAT", method.bin="Hybrid", 
weights.beta=c(1,25), weights = NULL, r.corr=0, impute.method = "fixed", 
is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1, N.Resampling=2 *10^6, 
seednum=100, epsilon=10^-6){


	Is.Seed=FALSE
	if(!is.null(seednum)){
		Is.Seed=TRUE
		seed.store=.Random.seed
		set.seed(seednum)
	}
	obj.res = SKATExactBin_CheckObj(obj)
	
	
	out.method<-SKAT_Check_Method(method,r.corr)
	method= out.method$method
	r.corr= out.method$r.corr
	Is.Hybrid=FALSE
	ExactMax=20000
	
	y = round(obj.res$mu + obj.res$res)
	case_ratio = abs(0.5 - sum(y)/ length(y))
	
	
	
	out.Z = SKATExactBin_Check(Z=Z, obj=obj.res, kernel = kernel, weights.beta=weights.beta, weights = weights, impute.method = impute.method, r.corr=r.corr, is_dosage = is_dosage, estimate_MAF=estimate_MAF,
		missing_cutoff=missing_cutoff, SetID = NULL, Is.Single=FALSE, Is.MakeZ1=FALSE)
	
	if(out.Z$return==1){
		out.Z$MAC=0
		out.Z$m=0
		out.Z$MAP=1
		out.Z$method.bin="NA"
		
		if(Is.Seed){
			# Restore seed
			.Random.seed = seed.store
		}
		
		return(out.Z)
	}
	
	MAC = out.Z$MAC
	m = out.Z$m
		
	# use MAC for select test if MAC < m. If not, use m. 
	# MAC > m when hard call genotyping was used. 
	MAC1 = max(MAC, m)
	
	if(method.bin == "Hybrid"){
		
		if(MAC1 <=20){
			method.bin="ER"	
			N.Resampling=2 *10^6 # get p-value from exact dist.
				
		} else if(MAC1 <=40){
			method.bin="ER.A"
		} else if(m <= 500){ # considering computation time, we currently use QA and MA up to m <=500
			method.bin="QA"
			if(case_ratio < 0.1){
				method.bin="MA"
			}
		}  else {
			method.bin="UA"
		}
		Is.Hybrid =TRUE
	}	
	
	if(method.bin == "UA"){
		
		re = SKAT(Z=Z, obj=obj.res, kernel = kernel, method=method, weights.beta=weights.beta, weights=weights, 
       	impute.method=impute.method, r.corr=r.corr, is_check_genotype=is_check_genotype,
       	is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF)	
		
		re$MAP = -1

	
	} else if(method.bin == "MA" || method.bin == "QA"){
		
		re = SKATExactBin(Z=Z, obj=obj.res, kernel = kernel, method.bin=method.bin, weights.beta=weights.beta, weights=weights,  
		impute.method=impute.method, r.corr=r.corr, is_check_genotype=is_check_genotype,
		is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF,
		N.Resampling=N.Resampling, epsilon=epsilon)
		

	
	} else if(method.bin == "ER" ){
	
		re = SKATExactBin(Z=Z, obj=obj.res, kernel = kernel, method.bin=method.bin, weights.beta=weights.beta, weights=weights,  
		impute.method=impute.method, r.corr=r.corr, is_check_genotype=is_check_genotype,
		is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF,
		N.Resampling=N.Resampling, epsilon=epsilon)

		if(!re$is.accurate && Is.Hybrid){
			p.value.exact = re$p.value
			method.bin = "QA"
			if(case_ratio < 0.1){
				method.bin="MA"
			}
			
			re = SKATExactBin(Z=Z, obj=obj.res, kernel = kernel, method.bin=method.bin, weights.beta=weights.beta, weights=weights,  
			impute.method=impute.method, r.corr=r.corr, is_check_genotype=is_check_genotype,
			is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF,
			N.Resampling=N.Resampling, epsilon=epsilon)
			re$p.value.exact = p.value.exact
		}	
	
	} else if(method.bin == "ER.A"){
	
		
		re = SKATExactBin.Adaptive(Z=Z, obj=obj.res, kernel = kernel, weights.beta=weights.beta, weights=weights, 
		impute.method=impute.method, r.corr=r.corr, is_check_genotype=is_check_genotype,
		is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF,
		N.Iter=50000, N.Resampling=N.Resampling, epsilon=epsilon)
		
		if(!re$is.accurate && Is.Hybrid){
			p.value.exact = re$p.value
			method.bin = "QA"
			if(case_ratio < 0.1){
				method.bin="MA"
			}
			
			re = SKATExactBin(Z=Z, obj=obj.res, kernel = kernel, method.bin=method.bin, weights.beta=weights.beta, weights=weights,  
			impute.method=impute.method, r.corr=r.corr, is_check_genotype=is_check_genotype,
			is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF,
			N.Resampling=N.Resampling, epsilon=epsilon)
			re$p.value.exact = p.value.exact
		}	

	} else if(method.bin == "Firth"){
	
		
		re = SKATExactBin.Firth(Z=Z, obj=obj.res, kernel = kernel, weights.beta=weights.beta, weights=weights, 
		impute.method=impute.method, r.corr=r.corr, is_check_genotype=is_check_genotype,
		is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF)

	} else {
		if(Is.Seed){
			# Restore seed
			.Random.seed = seed.store
		}
		
		msg<-sprintf("Error! %s is not correct method.bin type !", method.bin)
		stop(msg)
		
	}
	
	if(Is.Seed){
		# Restore seed
		.Random.seed = seed.store
	}
		
	re$MAC = MAC
	re$m = m
		
	re$method.bin = method.bin
	return(re)
	 
}



SKATBinary_Single<-function(Z, obj, method.bin="Hybrid", impute.method = "bestguess", 
is_check_genotype=TRUE, is_dosage = FALSE, estimate_MAF=1, N.Resampling=2*10^6, seednum=100, epsilon=10^-6){


	Is.Seed=FALSE
	if(!is.null(seednum)){
		Is.Seed=TRUE
		seed.store=.Random.seed
		set.seed(seednum)
	}

	obj.res = SKATExactBin_CheckObj(obj)
	out.Z=SKATBinary.Single.CheckZ(Z=Z,id_include=obj.res$id_include, impute.method=impute.method, is_check_genotype=is_check_genotype
	, is_dosage=is_dosage, estimate_MAF=estimate_MAF)

	# MAC, m and MAC1
	MAC = out.Z$MAC
	m = out.Z$m
	MAC1 = max(MAC, m)
		
	
	if(out.Z$return==1){
		
		if(Is.Seed){
			# Restore seed
			.Random.seed = seed.store
		}
		return(out.Z)
	}

	
	Is.Hybrid=FALSE
	Z = out.Z$Z.test 	
	if(method.bin=="Hybrid"){

		if(MAC1 <=20){
			
			method.bin="ER"
			N.Resampling=2 *10^6 # get p-value from exact dist.
			
		} else {
			method.bin="Firth"
		}
		Is.Hybrid =TRUE
	}	
	
	if(method.bin == "ER" ){
	
		re = SKATExactBin(Z=cbind(Z), obj=obj.res, kernel = "linear", method.bin=method.bin, weights.beta=c(1,1), weights=NULL,  
		impute.method=impute.method, r.corr=1, is_check_genotype=is_check_genotype,
		is_dosage = is_dosage, missing_cutoff=1, estimate_MAF=estimate_MAF,
		N.Resampling=N.Resampling, Is.Single=TRUE, epsilon=epsilon)

	
	} else if(method.bin == "Firth"){
	
		
		re = SKATExactBin.Firth(Z=cbind(Z), obj=obj.res, kernel = "linear", weights.beta=c(1,1), weights=NULL, 
		impute.method=impute.method, r.corr=1, is_check_genotype=is_check_genotype,
		is_dosage = is_dosage, missing_cutoff=1, estimate_MAF=estimate_MAF, Is.Single=TRUE)

	} else {
	
		if(Is.Seed){
			# Restore seed
			.Random.seed = seed.store
		}
	
		msg<-sprintf("Error! %s is not correct method.bin type !", method.bin)
		stop(msg)
		
	}
	
	if(Is.Seed){
		# Restore seed
		.Random.seed = seed.store
	}
	
	re$method.bin = method.bin
	return(re)
	 
}


###################################
# SSD

SKATBinary.SSD.OneSet = function(SSD.INFO, SetID, obj, ..., obj.SNPWeight=NULL){
	
	id1<-which(SSD.INFO$SetInfo$SetID == SetID)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set id [%s] from SSD!", SetID)
		stop(MSG)
	}	
	SetIndex<-SSD.INFO$SetInfo$SetIndex[id1]

	re = SKATBinary.SSD.OneSet_SetIndex(SSD.INFO, SetIndex, obj, ...)
	
	return(re)
}

#
# x is either y or SKAT_NULL_Model 
#
SKATBinary.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, obj, ..., obj.SNPWeight=NULL){

	re1 = SKAT.SSD.GetSNP_Weight(SSD.INFO, SetIndex, obj.SNPWeight=obj.SNPWeight)

	if(!re1$Is.weights){
		re<-SKATBinary(re1$Z, obj, ...)
	} else {
	
		re<-SKATBinary(re1$Z, obj, weights=re1$weights, ...)
	}
	
	return(re)

}


#
# Only SKAT_Null_Model obj can be used
#
SKATBinary.SSD.All = function(SSD.INFO, obj, ..., obj.SNPWeight=NULL){
	
	N.Set<-SSD.INFO$nSets
	OUT.Pvalue<-rep(NA,N.Set)
	OUT.Marker<-rep(NA,N.Set)
	OUT.Marker.Test<-rep(NA,N.Set)
	OUT.Error<-rep(-1,N.Set)
	OUT.Pvalue.Resampling<-NULL

	OUT.MAC<-rep(NA,N.Set)
	OUT.MAP<-rep(NA,N.Set)
	OUT.m<-rep(NA,N.Set)
	OUT.Method.bin<-rep("",N.Set)

	Is.Resampling = FALSE
	n.Resampling = 0
	
	obj.res = SKATExactBin_CheckObj(obj)
	
	if(obj.res$n.Resampling > 0){
		Is.Resampling = TRUE
		n.Resampling = obj$n.Resampling
		OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
	}
	

	for(i in 1:N.Set){
		Is.Error<-TRUE
		try1<-try(SKATBinary.SSD.OneSet_SetIndex(SSD.INFO, i, obj.res, ..., obj.SNPWeight=obj.SNPWeight) ,silent = TRUE)
		if(class(try1) == "try-error"){
			
				err.msg<-geterrmessage()
				msg<-sprintf("Error to run SKATBinary for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
				warning(msg,call.=FALSE)
			
		} else {
			re =try1
			OUT.Pvalue[i]<-re$p.value
			OUT.Marker[i]<-re$param$n.marker
			OUT.Marker.Test[i]<-re$param$n.marker.test
			OUT.MAC[i]<-re$MAC
			OUT.Method.bin[i]<-re$method.bin
			OUT.MAP[i]<-re$MAP
			OUT.m[i]<-re$m
			if(Is.Resampling){
				OUT.Pvalue.Resampling[i,]<-re$p.value.resampling
			}
		}
	}

	
	out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue
	, N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test
	, MAC=OUT.MAC, m=OUT.m
	, Method.bin = OUT.Method.bin
	, MAP = OUT.MAP )
	re<-list(results=out.tbl,P.value.Resampling=OUT.Pvalue.Resampling)
	class(re)<-"SKATBinary_SSD_ALL"

	return(re)	
}

