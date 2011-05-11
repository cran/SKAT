 \name{Generate_SSD_SetID}
 \alias{Generate_SSD_SetID}
 \title{Generate SNP set data file (SSD) }
 \description{
	Generate a SNP set data file (SSD) from binary plink formated data files using user specified SNP sets. If you want to use plink formated data files, you must generate the SSD files first. 

 }
 \usage{
	Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID,
	 File.SSD, File.Info)
 }
\arguments{
      \item{File.Bed}{ the name of the binary ped file (BED).}
      \item{File.Bim}{ the name of the binary map file (BIM).}
      \item{File.Fam}{ the name of the FAM file (FAM).}
      \item{File.SetID}{ the name of the Snp set ID file which defines SNP sets. The first column of the file must be Set ID, and the second column must be SNP ID. There should be no header!! }
      \item{File.SSD}{ the name of the SSD file generated. }
      \item{File.Info}{ the name of the SSD info file generated. }
      
}

\details{
 The SetID file is white-space (space or tab) delimitered file with 2 columns:
  SetID and SNP_ID.

 Please keep in mind that there should be no header!!
 The SNP_IDs and SetIDs should be less than 25 characters, otherwise, it will return error message.
      
 The SSD file is a binary formated file with genotype information. 
 The SSD info file is a text file. The first 6 rows have general information of data and SNP sets. The information of each set can be found from the 8th row. 

                                                        
}


\author{Seunggeun Lee, Larisa Miropolsky}

