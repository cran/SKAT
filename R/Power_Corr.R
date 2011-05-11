
#
#	W is beta^2, so make it sqrt(w)
Get_Power_Corr_WR_rho<-function(W,Pi,r.corr){
	
	n<-length(W)
	R<-matrix(rep(1,n*n),ncol=n) * r.corr
	diag(R)<-rep(1,n)

	temp<-sqrt(W) * Pi
	
	W_rho<- t(t(R * temp) * temp)
	diag(W_rho)<-W * Pi

	return(W_rho)

}

Get_Power_Corr.GetPower<-function(K,Mu,alpha.ALL,N.Sample){

	n.a<-length(alpha.ALL)
	OUT.Power<-rep(0,n.a)

	c1<-rep(0,4)
	c2<-rep(0,4)
	c_a<-rep(0,4)

	K1<-K
	K2<-K %*% K
	K3<-K2 %*% K
	
	c1[1] = trace.SKAT(K) * N.Sample
	c1[2] = sum(K *t(K)) * N.Sample^2
	c1[3] = sum(K2 * t(K) ) * N.Sample^3
	c1[4] = sum(K2 * t(K2)) * N.Sample^4

	c2[1] = trace.SKAT(  Mu) * N.Sample^2
	c2[2] = 2 * sum(  K * t(Mu)) * N.Sample^3
	c2[3] = 3 * sum(  K2 * t(Mu)) * N.Sample^4
	c2[4] = 4 * sum(  K3 * t(Mu)) * N.Sample^5

	for(i in 1:4){
		c_a[i]<-c1[i] + c2[i]
	}
		

	for(k in 1:n.a){
		alpha<-alpha.ALL[k]

		out<-Get_Critical_Value(c1, alpha)
		param<-Get_Liu_Params(c_a)

		t = (param$sigmaX/param$sigmaQ)*(out$q.crit-param$muQ) + param$muX
		power<-1-pchisq(t, df = param$l, ncp = param$d)
			
		OUT.Power[k]<-power
	}
	return(OUT.Power)

}

Get_Power_Logistic_R<-function(Z,eta,Beta0,ratio,alpha.ALL,N.Sample.ALL,Weight.Param=c(1,25),r.corr){

	Dist<-Dist_Case_Control(eta,Beta0,ratio)
	MAF<-colSums(Z * Dist) /2

	A<-Get_A(Z,eta,Dist,ratio)
	B<-Get_B(Z,eta, Dist,ratio)
	
	W<-Get_W_New(MAF,Weight.Param)		
	
	OUT.Power<-matrix(rep(0,length(alpha.ALL)*length(N.Sample.ALL)),ncol=length(alpha.ALL))
	rownames(OUT.Power)<-N.Sample.ALL
	colnames(OUT.Power)<-alpha.ALL

	for(j in 1:length(N.Sample.ALL)){
			
		N.Sample<-N.Sample.ALL[j]
		Pi1<-Get_PI_New(MAF,N.Sample)
		W_R<-Get_Power_Corr_WR_rho(W,Pi1,r.corr)

		K<-A %*% W_R
		Mu1<-B %*% W_R
		
		OUT.Power[j,]<-Get_Power_Corr.GetPower(K,Mu1,alpha.ALL,N.Sample)
		
	}

	return(OUT.Power)
}

Get_Power_Continuous_R<-function(Z,eta,alpha.ALL,N.Sample.ALL,Weight.Param=c(1,25),r.corr){

	A<-Get_A_Q(Z)
	B<-Get_B_Q(Z,eta)

	MAF<-colMeans(Z ) /2	
	W<-Get_W_New(MAF,Weight.Param)	

	#print(MAF)
	#print(Weight.Param)
	OUT.Power<-matrix(rep(0,length(alpha.ALL)*length(N.Sample.ALL)),ncol=length(alpha.ALL))
	rownames(OUT.Power)<-N.Sample.ALL
	colnames(OUT.Power)<-alpha.ALL
	
	for(j in 1:length(N.Sample.ALL)){
		
		N.Sample<-N.Sample.ALL[j]
		Pi1<-Get_PI_New(MAF,N.Sample)
		W_R<-Get_Power_Corr_WR_rho(W,Pi1,r.corr)

		#print(sum(Pi1))
		#print(sum(W_R))

		K<-A %*% W_R
		Mu1<-B %*% W_R
		
		OUT.Power[j,]<-Get_Power_Corr.GetPower(K,Mu1,alpha.ALL,N.Sample)
	}
	return(OUT.Power)

}

Power_Logistic_R<-function(Haplotypes=NULL, SNP.Location=NULL, SubRegion.Length=-1, Prevalence=0.01, Case.Prop=0.5, Causal.Percent=5, Causal.MAF.Cutoff=0.03, alpha =c(0.01,10^(-3),10^(-6)),N.Sample.ALL = 500 * (1:10), Weight.Param=c(1,25),N.Sim=100, OR.Type = "Log", MaxOR=5, Negative.Percent=0,r.corr=0){
	
	if(is.null(Haplotypes)){
	
		MSG_SKAT_Example()
		data(SKAT.haplotypes)

		Haplotypes<-SKAT.haplotypes$Haplotype
		SNP.Location<-SKAT.haplotypes$SNPInfo$CHROM_POS

	}
	if(is.null(SNP.Location)){
		stop("Error : SNP.Location is NULL")
	}

	Marker.MAF.ALL<-colMeans(Haplotypes) 


	# approximated Beta0
	Beta0<-log(Prevalence/(1-Prevalence))
	
	#####################################
	# 	Compute Power
	######################################	
	OUT.ALL<-NULL
	n1<-dim(Haplotypes)[1]

	for(i in 1:N.Sim){
		
		IDX.Marker<-Get_RandomRegion(SNP.Location,SubRegion.Length)	
		if(n1 >= 5000){		
			H1<-sample(1:n1,replace=FALSE)
			H2<-sample(1:n1,replace=FALSE)
		} else {
			H1<-sample(1:10000,replace=TRUE)
			H2<-sample(1:10000,replace=TRUE)
			H1[1:n1]<-sample(1:n1,replace=FALSE)
			H2[1:n1]<-sample(1:n1,replace=FALSE)
		}
	
			
		X1<-Haplotypes[H1,IDX.Marker] + Haplotypes[H2,IDX.Marker]
		Marker.MAF<-Marker.MAF.ALL[IDX.Marker]
		#print(Marker.MAF)

		Causal.Idx<-Get_CausalSNPs(Marker.MAF,  Causal.Percent/100, Causal.MAF.Cutoff)
		Marker.Causal.MAF<-Marker.MAF[Causal.Idx]
		#print(length(Causal.Idx))
		#print(Marker.Causal.MAF)
		Beta = Get_Beta(OR.Type, Marker.Causal.MAF, log(MaxOR),Negative.Percent/100)

		#print(Marker.Causal.MAF)
		#print(Beta)

		Causal.Idx1<-IDX.Marker[Causal.Idx]

		# Seunggeun Change
		#eta<-(Haplotypes[,Causal.Idx1] %*% Beta)[,1] - (t(Marker.Causal.MAF *2)  %*% Beta)[1,1]
		#eta1<-eta[H1] + eta[H2]
		eta1<-(X1[,Causal.Idx] %*% Beta)[,1] - (t(Marker.Causal.MAF *2)  %*% Beta)[1,1]

		#print(Beta)
		#####################################
		#
		#	Power

		OUT<-Get_Power_Logistic_R(X1,eta1,Beta0, Case.Prop, alpha,N.Sample.ALL,Weight.Param,r.corr)
		if(i==1){
			OUT.ALL<-OUT/N.Sim
		} else {
			OUT.ALL<-OUT.ALL + OUT/N.Sim
		}
		if(floor(i/10) * 10 == i){
			msg<-sprintf("%d/%d",i,N.Sim)
			print(msg)
		}
	}
	re<-list(Power = OUT.ALL)
	class(re)<-"SKAT_Power"

	return(re)
}


Power_Continuous_R<-function(Haplotypes=NULL, SNP.Location=NULL, SubRegion.Length=-1, Causal.Percent=5, Causal.MAF.Cutoff=0.03, alpha =c(0.01,10^(-3),10^(-6)),N.Sample.ALL = 500 * (1:10)
,Weight.Param=c(1,25),N.Sim=100,BetaType = "Log", MaxBeta=1.6, Negative.Percent=0,r.corr=0){
	
	if(is.null(Haplotypes)){
	
		MSG_SKAT_Example()
		data(SKAT.haplotypes)

		Haplotypes<-SKAT.haplotypes$Haplotype
		SNP.Location<-SKAT.haplotypes$SNPInfo$CHROM_POS

	}
	if(is.null(SNP.Location)){
		stop("Error : SNP.Location is NULL")
	}


	Marker.MAF.ALL<-colMeans(Haplotypes) 

	#####################################
	# 	Compute Power
	######################################	
	OUT.ALL<-NULL
	n1<-dim(Haplotypes)[1]
	out.r_2<-rep(0,N.Sim)

	for(i in 1:N.Sim){
		
		IDX.Marker<-Get_RandomRegion(SNP.Location,SubRegion.Length)
		
		if(n1 >= 5000){		
			H1<-sample(1:n1,replace=FALSE)
			H2<-sample(1:n1,replace=FALSE)
		} else {
			H1<-sample(1:10000,replace=TRUE)
			H2<-sample(1:10000,replace=TRUE)
			H1[1:n1]<-sample(1:n1,replace=FALSE)
			H2[1:n1]<-sample(1:n1,replace=FALSE)
		}
				
		X1<-Haplotypes[H1,IDX.Marker] + Haplotypes[H2,IDX.Marker]
		Marker.MAF<-Marker.MAF.ALL[IDX.Marker]

		Causal.Idx<-Get_CausalSNPs(Marker.MAF, Causal.Percent/100, Causal.MAF.Cutoff)
		Marker.Causal.MAF<-Marker.MAF[Causal.Idx]
		Beta = Get_Beta(BetaType, Marker.Causal.MAF, MaxBeta,Negative.Percent/100)
		Causal.Idx1<-IDX.Marker[Causal.Idx]

		out.r_2[i]<-sum(Beta^2*2*Marker.Causal.MAF*(1-Marker.Causal.MAF))

		
		# Seunggeun Change
		#eta<-(Haplotypes[,Causal.Idx1] %*% Beta)[,1] - (t(Marker.Causal.MAF *2)  %*% Beta)[1,1]
		#eta1<-eta[H1] + eta[H2]
		eta1<-(X1[,Causal.Idx] %*% Beta)[,1] - (t(Marker.Causal.MAF *2)  %*% Beta)[1,1]

		#print(Causal.Idx)
		#print(Beta)
		#####################################
		#
		#	Power
		
		OUT<-Get_Power_Continuous_R(X1,eta1,alpha,N.Sample.ALL,Weight.Param,r.corr)
		
		if(i==1){
			OUT.ALL<-OUT/N.Sim
		} else {
			OUT.ALL<-OUT.ALL + OUT/N.Sim
		}

		if(floor(i/10) * 10 == i){
			msg<-sprintf("%d/%d",i,N.Sim)
			print(msg)
		}
	}
	r_sq.v<-mean(out.r_2 /(out.r_2 +1))
	re<-list(Power = OUT.ALL, R.sq = r_sq.v)
	class(re)<-"SKAT_Power"

	return(re)
}






