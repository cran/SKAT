 \name{SKATBinary}
 \alias{SKATBinary}
 \alias{SKATBinary.SSD.OneSet}
 \alias{SKATBinary.SSD.OneSet_SetIndex}
 \title{SNP set test for binary traits with asymptotic and efficient resampling methods}
 \description{
     This function computes p-values of Burden test, SKAT, and SKAT-O for binary traits using asymptotic and efficient resampling methods.   
 }
 \usage{
 
	SKATBinary(Z, obj, kernel = "linear.weighted", method="SKAT"
	, method.bin="Hybrid", weights.beta=c(1,25), weights = NULL
	, r.corr=0, impute.method = "fixed", is_check_genotype=TRUE
	, is_dosage = FALSE, missing_cutoff=0.15
	, estimate_MAF=1, N.Resampling=2 *10^6, seednum=100, epsilon=10^-6)
	
	SKATBinary.SSD.OneSet(SSD.INFO, SetID, obj, \dots ,obj.SNPWeight=NULL)

	SKATBinary.SSD.OneSet_SetIndex(SSD.INFO, SetIndex, obj, \dots ,obj.SNPWeight=NULL)


 }
\arguments{
      \item{Z}{a numeric genotype matrix with each row as a different individual and each column as a separate gene/snp. 
      Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing, 
      where A is a major allele and a is a minor allele. }
      \item{obj}{an output object of the SKAT_Null_Model function. }
      \item{kernel}{a type of kernel (default= "linear.weighted"). See details in SKAT function. }
      \item{method}{a type of gene based test (default= "SKAT"). The possible choices are
      "SKAT", "Burden" and "SKATO", which represents SKAT, Burden and SKAT-O tests, respectively. 
      This parameter differs from the "method" parameter in SKAT function. }
      \item{method.bin}{a type of method to compute a p-value (default="Hybrid"). Possible choices are
      "Hybrid", "ER", "ER.A", "QA", "MA" and "UA". See details}
       \item{weights.beta}{a numeric vector of parameters of beta weights. 
      It is only used for weighted kernels. 
      If you want to use your own  weights, please specify the ``weights'' parameter.}
      \item{weights}{a numeric vector of weights for the weighted kernels. See details in SKAT function.}
      \item{impute.method}{a method to impute missing genotypes (default= "fixed"). "bestguess" imputes missing genotypes as most likely 
      values (0,1,2), "random" imputes missing genotypes by generating binomial(2,p) random variables (p is the MAF), 
      and "fixed" imputes missing genotypes by assigning the mean genotype value (2p).}
      \item{r.corr}{the \eqn{\rho} parameter of new class of kernels with compound symmetric correlation structure for genotype effects (default= 0). 
      If it is a vector, SKAT will conduct the optimal test. See details.}
      \item{is_check_genotype}{a logical value indicating whether to check the validity of the genotype matrix Z (default= TRUE). If you use non-SNP type data and want to run kernel machine test, please set it FALSE, otherwise you will get an error message. If you use SNP data or imputed data, please set it TRUE. If it is FALSE, and you use weighted kernels, the weights should be given through ``weights'' parameter.}
      \item{is_dosage}{a logical value indicating whether the matrix Z is a dosage matrix. If it is TRUE, SKAT will ignore ``is_check_genotype''. }
      \item{missing_cutoff}{a cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.}
      \item{estimate_MAF}{a numeric value indicating how to estimate MAFs for the weight calculation and 
      the missing genotype imputation. If estimate_MAF=1 (default), SKAT uses all samples to estimate MAFs. 
      If estimate_MAF=2, only samples with non-missing phenotypes and covariates are used to estimate MAFs. }
      \item{N.Resampling}{a number of resampling to be conducted to get p-values (default=2 *10^6).}
      \item{seednum}{a seed number for random number generation (default=100). If NULL, no seed number will be assigned.}
      \item{epsilon}{a precision level  (default=10^-6).}
      
      \item{SSD.INFO}{an SSD_INFO object returned from Open_SSD. }
      \item{SetID}{a character value of Set ID. You can find a set ID of each set from SetInfo object of SSD.INFO}
      \item{SetIndex}{a numeric value of Set index. You can find a set index of each set from SetInfo object of SSD.INFO  }
      \item{\dots}{further arguments to be passed to ``SKATBinary'' }
      \item{obj.SNPWeight}{an output object of Read_SNP_WeightFile (default=NULL). 
      If NULL, the beta weight with the ``weights.beta'' parameter is used.  }
	  
}
\value{
	\item{p.value}{the p-value. It will be the mid p-value if ER or ER.A is used to compute the p-value.}
	\item{p.value.standard}{(ER and ER.A only)the standard p-value.}
	\item{p.value.resampling}{the p-value from resampled outcome. You can obtain it when n.Resampling in SKAT_Null_Model was > 0. See the SKAT_Null_Model. }
	\item{p.value.standard.resampling}{(ER and ER.A only)the standard p-value from resampled outcome.}
	\item{m}{the number of individuals with minor alleles.}
	\item{MAP}{the minimum possible p-values. It is available when the method.bin="ER" and m is sufficiently small.}
	\item{MAC}{the total minor allele count (MAC).}
	\item{n.total}{(ER only) the number of resampling to be generated to get the p-value. 
	It can be smaller than N.Resampling when the total number of configurations of case-controls among individuals with minor alleles are smaller than
	N.Resampling.}
  	\item{is.accurate}{logical value for the accuracy of the p-value. If it is false, more resampling is needed to accurately estimate the p-value. }
	\item{param$n.marker}{a number of SNPs in the genotype matrix}  
	\item{param$n.marker.test}{a number of SNPs used for the test. It can be different from param$n.marker when 
	some markers are monomorphic or have higher missing rates than the missing_cutoff. } 
	\item{method.bin}{a type of method to be used to compute the p-value.}
}
\details{

This function implements six methods (method.bin) to compute p-values: 1) Efficient resampling (ER);
2) Quantile adjusted moment matching (QA); 3) Moment matching adjustment (MA);
4) No adjustment (UA);  5) Adaptive ER (ER.A); and 6) Hybrid. 
"Hybrid" selects a method based on the total minor allele count (MAC), the number of individuals with minor 
alleles (m), and the degree of case-control imbalance. When method.bin="ER" or "ER.A", SKATBinary compute mid-p-values and minimum achievable 
mid p-values. 

}


\author{Seunggeun Lee}

\examples{


data(SKATBinary.example)
attach(SKATBinary.example)


obj<-SKAT_Null_Model(y ~ x1 + x2, out_type="D")

# run SKAT (default method) with Hybrid
out = SKATBinary(Z, obj)

# p-value
out$p.value

# MAP
out$MAP

# method used to compute p-value (method.bin)
out$method.bin


#
#	Run burden and SKAT-O with Hybrid

SKATBinary(Z, obj, method="Burden")$p.value
SKATBinary(Z, obj, method="SKATO")$p.value

#
#	Run with SKAT-QA, -MA and -UA

SKATBinary(Z, obj, method.bin="QA")$p.value

SKATBinary(Z, obj, method.bin="MA")$p.value

SKATBinary(Z, obj, method.bin="UA")$p.value

# UA from SKAT function
SKAT(Z, obj)$p.value


#
#	Run with Adaptive ER

out =SKATBinary(Z, obj, method.bin="ER.A")

out$p.value

# the number of total resampling is smaller than 2*10^6 (default value)
out$n.total 

}


