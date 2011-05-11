Print_Error_SSD<-function(code){

	if(code == 0 ){
		return(0);
	} else if(code == 1){
		stop("Error Can't open BIM file")
	} else if(code == 2){
		stop("Error Can't open FAM file")
	} else if(code == 3){
		stop("Error Can't open BED file")
	} else if(code == 4){
		stop("Error Can't open SETID file")
	} else if(code == 5){
		stop("Error Can't write SSD file")
	} else if(code == 6){
		stop("Error Can't read SSD file")
	} else if(code == 7){
		stop("Error Can't write INFO file")
	} else if(code == 8){
		stop("Error Can't read INFO file")
	} else if(code == 9){
		stop("Error Can't write INFO file")
	} else if(code == 13){
		stop("Error Wrong SNP or Individual sizes")
	} else if(code == 14){
		stop("Error SetID not found")
	} else {
		MSG<-sprintf("Error [%d]\n",code)
		stop(MSG)
	}
	
	return(1)
}


Check_File_Exists<-function(FileName){
	
	if(!file.exists(FileName)){
		Msg<-sprintf("File %s does not exist\n",FileName)
		stop(Msg)
	}

}


Check_ID_Length<-function(FileName){
	

	SSD.Info<-try(read.table(FileName, header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
	if(class(SSD.Info)=="try-error"){
		stop("Error in SetID file!") 
	}

	n1<-length(which(nchar(SSD.Info[,1]) > 25))
	n2<-length(which(nchar(SSD.Info[,2]) > 25))

	if(n1 > 0){
		stop("Some SetIDs have more than 25 characters!") 
	}
	if(n2 > 0){
		stop("Some SNP_IDs have more than 25 characters!") 
	}	

	nSets<-length(unique(SSD.Info[,1]))
	return(list(nSets=nSets))
}





Print_File_Info<-function(INFO){
	
	MSG<-sprintf("%d Samples, %d Sets, %d Total SNPs\n",INFO$nSample, INFO$nSets, INFO$nSNPs)
	cat(MSG)

}

Read_File_Info<-function(File.Info){

	Check_File_Exists(File.Info)
	
	info1<-read.table(File.Info, header=FALSE, nrows= 6, sep='\t')
	info2<-read.table(File.Info, header=FALSE, skip=7, sep='\t', stringsAsFactors=FALSE,comment.char = "")

	INFO<-list()
	INFO$WindowSize<-info1[1,1]
	INFO$OverlapSize<-info1[2,1]
	INFO$nSNPs<-info1[3,1]	
	INFO$nSample<-info1[4,1]
	INFO$nDecodeSize<-info1[5,1]
	INFO$nSets<-info1[6,1]
	INFO$SetInfo<-data.frame(SetIndex= info2[,1], SetID = as.character(info2[,3])
	, SetSize = info2[,4], Offset = info2[,2],stringsAsFactors=FALSE)

	Print_File_Info(INFO)
	
	return(INFO)


}


Read_File_Info_Head<-function(File.Info, Is.Print=TRUE){

	Check_File_Exists(File.Info)
	
	info1<-read.table(File.Info, header=FALSE, nrows= 6, sep='\t')


	INFO<-list()
	INFO$WindowSize<-info1[1,1]
	INFO$OverlapSize<-info1[2,1]
	INFO$nSNPs<-info1[3,1]	
	INFO$nSample<-info1[4,1]
	INFO$nDecodeSize<-info1[5,1]
	INFO$nSets<-info1[6,1]
	
	if(Is.Print == TRUE){
		Print_File_Info(INFO)
	} 

	return(list(nSets=INFO$nSets))

}

Check_SetID_LOG<-function(File.SSD){

	File.LOG<-sprintf("%s_LOG.txt",File.SSD)
	if(!file.exists(File.LOG)){
		MSG<-sprintf("Error %s file did not generated!\n",File.LOG)
		stop(MSG)
	}
	
	temp<-read.delim(File.LOG, header = FALSE, sep = "\n",comment.char = "")
	N.Miss<-dim(temp)[1] -1

	if(N.Miss > 0){
		MSG1<-sprintf("Warning: %d SNPs in the SetID file were not found in Bim file!\n",N.Miss)
		MSG2<-sprintf("Please check %s file!\n",File.LOG)
		
		cat(MSG1)
		cat(MSG2)
	}

}


Create_Temporaly_SetID<-function(FileName){
	

	FileName.out<-sprintf("%s.temp",FileName)
	SSD.Info<-try(read.table(FileName, header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
	if(class(SSD.Info)=="try-error"){
		stop("Error in SepID file!") 
	}

	ord1<-order(SSD.Info[,1])
	SSD.Info1<-SSD.Info[ord1,]

	write.table(SSD.Info1, file=FileName.out, quote=FALSE, row.names=FALSE, col.names=FALSE)
	
	return(FileName.out)
}



##################################################################
#
#	Generate SSD Files

Generate_SSD_SetID_Work<-function(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info){


	File.Bed<-normalizePath(File.Bed ,mustWork =FALSE)
	File.Bim<-normalizePath(File.Bim ,mustWork =FALSE)
	File.Fam<-normalizePath(File.Fam ,mustWork =FALSE)
	File.SetID<-normalizePath(File.SetID ,mustWork =FALSE)
	File.SSD<-normalizePath(File.SSD ,mustWork =FALSE)
	File.Info<-normalizePath(File.Info ,mustWork =FALSE)


	Check_File_Exists(File.Bed)
	Check_File_Exists(File.Bim)
	Check_File_Exists(File.Fam)
	Check_File_Exists(File.SetID)

	nSet1 = Check_ID_Length(File.SetID)$nSet
	
	err_code<-0


	temp<-.C("R_Generate_MWA_SetID_File"
	, as.character(File.Bed), as.character(File.Bim), as.character(File.Fam)
	, as.character(File.SetID), as.character(File.SSD), as.character(File.Info)
	, as.integer(err_code))
	
	
	error_code<-temp[[7]]
	Print_Error_SSD(error_code)

        # Check SetID_LOG file 
	Kill_SSD_SetID()	
	Check_SetID_LOG(File.SSD)

	nSet2 = Read_File_Info_Head(File.Info, Is.Print=FALSE)$nSet

	re=1;
	if(nSet1 < nSet2){
		re=-1;	# re=-1 not matching numbers
	}

	return(re)	
}


Generate_SSD_SetID<-function(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info){


	re = Generate_SSD_SetID_Work(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)

	if(re == -1){
		MSG1<-sprintf("Warning: SSD file has more SNP sets then SetID file.It happens when SNPs in sets are not contiguous!\n")
		MSG2<-sprintf("SKAT generates a temporary SetID file with contiguous SNP sets. However the order of SNP sets will be different from the order of SNP sets in the original SetID file! \n\n")
		cat(MSG1)
		cat(MSG2)

		File.SetID.temp<-Create_Temporaly_SetID(File.SetID)
		re = Generate_SSD_SetID_Work(File.Bed, File.Bim, File.Fam, File.SetID.temp, File.SSD, File.Info)
		if(re == -1){
			stop("Keep having more SetIDs in the SSD file!")
		}
	}	

	Read_File_Info_Head(File.Info, Is.Print=TRUE)
	print("SSD and Info files are created!")

}

Generate_SSD_MovingWindow<-function(File.Bed, File.Bim, File.Fam, File.SSD, File.Info, WindowSize,Overlap){

	File.Bed<-normalizePath(File.Bed ,mustWork =FALSE)
	File.Bim<-normalizePath(File.Bim ,mustWork =FALSE)
	File.Fam<-normalizePath(File.Fam ,mustWork =FALSE)
	File.SSD<-normalizePath(File.SSD ,mustWork =FALSE)
	File.Info<-normalizePath(File.Info ,mustWork =FALSE)


	Check_File_Exists(File.Bed)
	Check_File_Exists(File.Bim)
	Check_File_Exists(File.Fam)


	err_code<-0

	temp<-.C("R_Generate_MWA_MovingWindow"
	, as.character(File.Bed), as.character(File.Bim), as.character(File.Fam)
	, as.character(File.SSD), as.integer(WindowSize), as.integer(Overlap), as.character(File.Info)
	, as.integer(err_code))

	Kill_SSD_MovingWindow()

	error_code<-temp[[8]]
	Print_Error_SSD(error_code)


	Read_File_Info_Head(File.Info)
	print("SSD and Info files are created!")

}

Kill_SSD_SetID<-function(){
	temp<-.C("R_Kill_MWA_SetID_File")
}

Kill_SSD_MovingWindow<-function(){
	temp<-.C("R_Kill_MWA_MovingWindow")
}

##################################################################
#
#	Open and Close the SSD Files

SSD_FILE_OPEN.isOpen = 0
SSD_FILE_OPEN.FileName = ""

#assign("SSD_FILE_OPEN.isOpen", 0, envir=.GlobalEnv)
#assign("SSD_FILE_OPEN.FileName","", envir=.GlobalEnv)

Close_SSD<-function(){

	if(get("SSD_FILE_OPEN.isOpen", envir=.GlobalEnv) == 1){
		temp<-.C("R_Close_MWA")
		Msg<-sprintf("Close the opened SSD file: %s\n"
		,get("SSD_FILE_OPEN.FileName", envir=.GlobalEnv));
		cat(Msg)
		assign("SSD_FILE_OPEN.isOpen", 0, envir=.GlobalEnv);
	} else{
		Msg<-sprintf("No opened SSD file!\n");
		cat(Msg)		
	}
}

Open_SSD<-function(File.SSD, File.Info){

	err_code<-0
	File.SSD<-normalizePath(File.SSD ,mustWork =FALSE)
	File.Info<-normalizePath(File.Info ,mustWork =FALSE)

	Check_File_Exists(File.SSD)
	Check_File_Exists(File.Info)

	if(get("SSD_FILE_OPEN.isOpen", envir=.GlobalEnv) == 1){
		Close_SSD();
	}

	# Read Info File
	INFO<-Read_File_Info(File.Info)

	# Read SSD File
	temp<-.C("R_Open_MWA", as.character(File.SSD), as.character(File.Info)
	, as.integer(err_code))

	error_code<-temp[[3]]
	Print_Error_SSD(error_code)


	Msg<-sprintf("Open the SSD file\n");
	cat(Msg)

	#SSD_FILE_OPEN.isOpen<<-1
	#SSD_FILE_OPEN.FileName<<-File.SSD

	assign("SSD_FILE_OPEN.isOpen", 1, envir=.GlobalEnv)
	assign("SSD_FILE_OPEN.FileName",File.SSD, envir=.GlobalEnv)


	return(INFO)
	
}

#######################################################



##################################################################
#
#	Get Genotype Matrix

Get_Genotypes_SSD<-function(SSD_INFO, Set_Index){

	Is_MakeFile=0
	if(get("SSD_FILE_OPEN.isOpen", envir=.GlobalEnv) == 0){
		stop("SSD file is not opened. Please open it first!")
	}

	id1<-which(SSD_INFO$SetInfo$SetIndex == Set_Index)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set index [%d] from SSD!\n", Set_Index)
		stop(MSG)
	}	
	Set_Index<-SSD_INFO$SetInfo$SetIndex[id1]

	err_code<-0
	N.SNP<-SSD_INFO$SetInfo$SetSize[id1]
	N.Sample<-SSD_INFO$nSample
	size<-N.SNP * N.Sample

	Z<-rep(9,size)


	temp<-.C("R_Get_Genotypes",as.integer(Set_Index),as.integer(Z),as.integer(size)
	,as.integer(Is_MakeFile), as.integer(err_code))

	error_code<-temp[[5]]
	Print_Error_SSD(error_code)


	Z.out.t<-matrix(temp[[2]],byrow=TRUE, nrow=N.SNP)
	return(t(Z.out.t))

}






