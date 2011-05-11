 \name{Get_Resampling_Pvalue}
 \alias{Get_Resampling_Pvalue}
 \alias{Get_Resampling_Pvalue_1}
 \title{Compute the resampling p-value}
 \description{
     To compute a resampling p-value using the resampled residuals. 
     To use it, you need to obtain resampling residuals using SKAT_Null_Model, and then run SKAT.
 }
 \usage{
	Get_Resampling_Pvalue(obj)

	Get_Resampling_Pvalue_1(p.value, p.value.resampling)


 }
\arguments{
      \item{obj}{a SKAT outcome object.}
      \item{p.value}{a numeric value of the SKAT p-value.}
      \item{p.value.resampling}{a vector of p-values of the resampled residuals.}  
}
\value{
	\item{p.value}{the resampling p-value. It is computed as (n1 +1)/(n+1), where n is the number of resampling, 
	and n1 is the number of resampled residual p-values smaller than the original sample p-value. }
  	\item{is_smaller}{a logical value which indicates whether the resampling p-value should be smaller. 
  	If n1=0, then it has TRUE, otherwise it has FALSE. }
  	
}
\details{
	See SKAT_Null_Model
}

\author{Seunggeun Lee}

