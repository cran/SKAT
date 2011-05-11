 \name{Read_Plink_FAM}
 \alias{Read_Plink_FAM}
 \title{Read a Plink FAM file}
 \description{
     Read a Plink FAM file
 }
 \usage{
	Read_Plink_FAM(Filename, Is.binary=TRUE, flag1=0)
 }
\arguments{
      \item{Filename}{an input file name of plink FAM file}
      \item{Is.binary}{if TRUE, the phenotype is binary. If phenotype is continuous, it should be FALSE}
      \item{flag1}{0 represents the default coding of unaffected/affected (1/2) (default=0), and 1 represents 0/1 coding. flag1=1 is the same as --1 flag. Please see the plink manual. }      
}
\value{
	A data frame of Family ID (FID), Individual ID (IID), Paternal ID (PID), Maternal ID(MID), Sex, and Phenotype. 
  	
}



\author{Seunggeun Lee}




