 \name{Resampling_FWER}
 \alias{Resampling_FWER}
 \alias{Resampling_FWER_1}
 \title{Return significant SNP sets after controlling family wise error rate (FWER)}
 \description{

      This function returns significant SNP sets after controlling for family wise error rate (FWER) 
      using resampled residuals. To use it, you need to obtain resampling residuals using SKAT_Null_Model, 
      and then conduct the SKAT repeatedly for all genes/SNP sets or use SKAT.SSD.All function.
 }
 \usage{

	Resampling_FWER(obj,FWER=0.05)

	Resampling_FWER_1(P.value, P.value.Resampling, FWER=0.05)

 }
\arguments{
      \item{obj}{an object returned from SKAT.SSD.All function.}
      \item{P.value}{a vector of the SKAT p-value. If you test 100 genes, this vector should have 100 p-values.}
      \item{P.value.Resampling}{a matrix of p-values of the resampled residuals. Each row represents each gene/snp set, and each column represents resampling set. For example, if you have 100 genes, and conducted resampling 1000 times ( ex.n.Resampling=1000 in SKAT_Null_Model), then it should be a 100 x 1000 matrix.}  
      \item{FWER}{a numeric value of FWER rate to control (default=0.05)}
}
\value{
	\item{results}{If you use the returned object from SKAT.SSD.all function, 
	it is a sub-table of significant snp sets of the result table in the obj. 
	If you use P.value and P.value.Resampling, it is a vector of significant p-values. 
	If there is no significant snp set, it has NULL value. }
	\item{n}{a numeric value of the number of significant snp sets.}
  	\item{ID}{a vector of indexes of significant snp sets.}
}

\author{Seunggeun Lee}

