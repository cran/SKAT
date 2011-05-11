
#
#	Read Fam File
#
Read_Plink_FAM<-function(Filename, Is.binary=TRUE, flag1=0){

	Check_File_Exists(Filename)
	Fam.Obj<-read.table(Filename, header=FALSE)
	colnames(Fam.Obj)<-c("FID","IID","PID","MID", "Sex", "Phenotype")

	id.missing<-which(Fam.Obj$Phenotype == -9)
	if(length(id.missing) > 0){
		Fam.Obj$Phenotype[id.missing]<-NA
	}

	if(Is.binary && flag1 == 0){
		
		id_0<-which(Fam.Obj$Phenotype == 0)
		id_1<-which(Fam.Obj$Phenotype == 1)
		id_2<-which(Fam.Obj$Phenotype == 2)

		if(length(id_0)> 0){
			Fam.Obj$Phenotype[id_0] <-NA
		}
		if(length(id_1)> 0){
			Fam.Obj$Phenotype[id_0] <-0
		}
		if(length(id_2)> 0){
			Fam.Obj$Phenotype[id_0] <-1
		}

	}

	return(Fam.Obj)

}

#
# x is either y or SKAT_NULL_Model 
#
SKAT.SSD.OneSet = function(SSD.INFO, SetID, obj, ...){
	
	id1<-which(SSD.INFO$SetInfo$SetID == SetID)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set id [%s] from SSD!", SetID)
		stop(MSG)
	}	
	Set_Index<-SSD.INFO$SetInfo$SetIndex[id1]

	Z<-Get_Genotypes_SSD(SSD.INFO, Set_Index)
	re<-SKAT(Z, obj, ...)
	
	return(re)
}

#
# x is either y or SKAT_NULL_Model 
#
SKAT.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, obj, ...){
	
	id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
		stop(MSG)
	}	
	SetID<-SSD.INFO$SetInfo$SetID[id1]


	Z<-Get_Genotypes_SSD(SSD.INFO, SetIndex)
	re<-SKAT(Z, obj, ...)
	return(re)
}


#
# Only SKAT_Null_Model obj can be used
#
SKAT.SSD.All = function(SSD.INFO, obj, ...){
	
	N.Set<-SSD.INFO$nSets
	OUT.Pvalue<-rep(NA,N.Set)
	OUT.Marker<-rep(NA,N.Set)
	OUT.Marker.Test<-rep(NA,N.Set)
	OUT.Error<-rep(-1,N.Set)
	OUT.Pvalue.Resampling<-NULL

	Is.Resampling = FALSE
	n.Resampling = 0
	
	if(class(obj) == "SKAT_NULL_Model"){
		if(obj$n.Resampling > 0){
			Is.Resampling = TRUE
			n.Resampling = obj$n.Resampling

			OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
		}
	} else if(class(obj) == "SKAT_NULL_Model_ADJ"){
		if(obj$re1$n.Resampling > 0){
			Is.Resampling = TRUE
			n.Resampling = obj$re1$n.Resampling

			OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
		}
	}


	for(i in 1:N.Set){
		Is.Error<-TRUE
		try1<-try(Get_Genotypes_SSD(SSD.INFO, i),silent = TRUE)
		if(class(try1) != "try-error"){
			Z<-try1
			Is.Error<-FALSE
			
			
		} else {
			err.msg<-geterrmessage()
			msg<-sprintf("Error to get genotypes of %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
			warning(msg,call.=FALSE)
		}
	
		if(!Is.Error){
			Is.Error<-TRUE
			try2<-try(SKAT_1(Z,obj, ...),silent = TRUE)
			
			if(class(try2) != "try-error"){
				re<-try2
				Is.Error<-FALSE
			} else {

				err.msg<-geterrmessage()
				msg<-sprintf("Error to run SKAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
				warning(msg,call.=FALSE)
			}
		}
		
		if(!Is.Error){

			OUT.Pvalue[i]<-re$p.value
			OUT.Marker[i]<-re$param$n.marker
			OUT.Marker.Test[i]<-re$param$n.marker.test
			if(Is.Resampling){
				OUT.Pvalue.Resampling[i,]<-re$p.value.resampling
			}
		}
	}

	
	out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test)
	re<-list(results=out.tbl,P.value.Resampling=OUT.Pvalue.Resampling)
	class(re)<-"SKAT_SSD_ALL"

	return(re)	
}

