 \name{SKAT_ChrX}
 \alias{SKAT_ChrX}
 \title{SNP-set (Sequence) Kernel Association Test for  X chromosome variables}
 \description{
     Test for association between a set of SNPS/genes in the X chromosome and continuous or dichotomous outcomes using the kernel machine.      
 }
 \usage{

SKAT_ChrX(Z, obj, is_X.inact =TRUE
, kernel = "linear.weighted", method="davies", weights.beta=c(1,25)
, weights = NULL, impute.method = "fixed", r.corr=0, is_check_genotype=TRUE
, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1, SetID=NULL)


 }
\arguments{
      \item{Z}{a numeric genotype matrix with each row as a different individual and each column as a separate gene/snp. 
      Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing, where A is a major allele and a is a minor allele. Missing genotypes will be imputed by the simple Hardy-Weinberg equilibrium (HWE) based imputation. }
      \item{obj}{an output object of the SKAT_Null_Model_ChrX function. }
      \item{is_X.inact}{an indicator variable for the X-inactivation coding (default=TRUE). 
      Male genotypes are coded as g=(0,2) when it is TRUE, and  g=(0,1) when it is false.}
      \item{kernel}{a type of kernel (default= "linear.weighted"). }
      \item{method}{a method to compute the p-value (default= "davies"). 
      "davies" represents an exact method that  computes the p-value by inverting the characteristic function of the mixture chisq, 
      "liu" represents an approximation method that matches the first 3 moments, 
      "liu.mod" represents modified "liu" method that matches kurtosis instead of skewness 
      to improve tail probability approximation, "optimal.adj" represents a SKAT-O based on an unified approach, 
      and "optimal" is an old version of the implementation of SKAT-O.}
      \item{weights.beta}{a numeric vector of parameters of beta weights. 
      It is only used for weighted kernels. 
      If you want to use your own  weights, please specify the ``weights'' parameter.}
      \item{weights}{a numeric vector of weights for the weighted kernels. 
      It is \eqn{\sqrt{w}} in the SKAT paper. 
      So if you want to use the Madsen and Browning (2009) weight, you should set each element of weights as \eqn{1/ \sqrt{p(1-p)}}, 
      not \eqn{1/ p(1-p)}. When it is NULL, the beta weight with the ``weights.beta'' parameter is used. }
      \item{impute.method}{a method to impute missing genotypes (default= "fixed"). "fixed" imputes missing genotypes by assigning the mean genotype value (2p), and "bestguess" uses best guess genotype values. }
      \item{r.corr}{the \eqn{\rho} parameter of new class of kernels with compound symmetric correlation structure for genotype effects (default= 0). If you give a vector value, SKAT will conduct the optimal test. }
      \item{is_check_genotype}{a logical value indicating whether to check the validity of the genotype matrix Z (default= TRUE). If you use non-SNP type data and want to run kernel machine test, please set it FALSE, otherwise you will get an error message.
      If you use SNP data or imputed data, please set it TRUE. If it is FALSE, and you use weighted kernels, the weights should be given through ``weights'' parameter.}
      \item{is_dosage}{a logical value indicating whether the matrix Z is a dosage matrix. If it is TRUE, SKAT will ignore ``is_check_genotype''. }
      \item{missing_cutoff}{a cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.}
      \item{estimate_MAF}{a numeric value indicating how to estimate MAFs for the weight calculation and 
      the missing genotype imputation. If estimate_MAF=1 (default), SKAT uses all samples to estimate MAFs. 
      If estimate_MAF=2, only samples with non-missing phenotypes and covariates are used to estimate MAFs. }
      \item{SetID}{Internal use only. }
     
      
}
\value{
	\item{p.value}{the p-value of SKAT. }
	\item{p.value.resampling}{the p-value from resampled outcome. You can get it when you use obj from SKAT_Null_Model function with resampling. See the SKAT_Null_Model. }
	\item{p.value.noadj}{the p-value of SKAT without the small sample adjustment. It only appears when small sample adjustment is applied.}
	\item{p.value.noadj.resampling}{the p-value from resampled outcome without the small sample adjustment. It only appears when small sample adjustment is applied. }
  	\item{pval.zero.msg}{(only when p.value=0) text message that shows how small the p.value is. ex. "Pvalue < 1.000000e-60" when p.value is smaller than \eqn{10^{-60}} } 
  	\item{Q}{the test statistic of SKAT. It has NA when method="optimal.adj" or "optimal".}
	\item{param}{estimated parameters of each method.}   
	\item{param$Is_Converged}{ (only with method="davies") an indicator of the convergence. 1 indicates the method is converged, and 0 indicates the method is not converged. When 0 (not converged), "liu" method is used to compute p-value. }  
	\item{param$n.marker}{a number of SNPs in the genotype matrix}  
	\item{param$n.marker.test}{a number of SNPs used for the test. It can be different from param$n.marker when 
	some markers are monomorphic or have higher missing rates than the missing_cutoff. } 
	
}
\details{


For details of parameters, please see a page for the SKAT function. 
                                   

}


\author{Clement Ma and Seunggeun Lee}


\examples{



data(SKAT.example.ChrX)
attach(SKAT.example.ChrX)

#############################################################
#	Compute the P-value of SKAT 

# binary trait
obj.x<-SKAT_Null_Model_ChrX(y ~ x1 +x2 + Gender, SexVar="Gender", out_type="D")
SKAT_ChrX(Z, obj.x, kernel = "linear.weighted", r.corr=0)

# Burden
SKAT_ChrX(Z, obj.x, kernel = "linear.weighted", r.corr=1)

# SKAT-O
SKAT_ChrX(Z, obj.x, kernel = "linear.weighted", method="optimal.adj")




}


