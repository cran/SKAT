 \name{SKAT_Null_Model}
 \alias{SKAT_Null_Model}
 \title{Get parameters and residuals from the H0 model}
 \description{
     Compute model parameters and residuals for SKAT. You also can obtain resampled residuals that can be used to compute resampling p-value or to control family-wise error rate.
 }
 \usage{

SKAT_Null_Model(formula, data=NULL, out_type="C", n.Resampling=0
, type.Resampling="bootstrap", Adjustment=TRUE)

 }
\arguments{
      \item{formula}{an object of class ``formula'': a symbolic description of the NULL model to be fitted.}
      \item{data}{an optional data frame containing the variables in the model (default=NULL).  If it is NULL, the variables are taken from 'environment(formula)'}
      \item{out_type}{an indicator of the outcome type. "C" for the continuous outcome and "D" for the dichotomous outcome.}
      \item{n.Resampling}{a numeric value of the number of resampling (default=0). If you don't want resampling, please set n.Resampling=0. }
      \item{type.Resampling}{ resampling methods (default="bootstrap"). see details.}
      \item{Adjustment}{If TRUE, a small sample adjustment will be applied when the sample size < 2000 and the trait is binary (default=TRUE). See details}      
}
\value{
	This function returns an object that has model parameters and residuals of the NULL model of no association between genetic variables and outcome phenotypes. After obtaining it, please use SKAT function to conduct the association test.

}
\details{

There are 2 different methods to get resampled residuals.
"bootstrap" conducts the parametric bootstrap to resample residuals under the NULL model with considering covariates. 
If there is no covariate, "bootstrap" is equivalent to the permutation method.
"perturbation" perturbs the residuals by multiplying mean zero and variance one normal random variable. The default method is "bootstrap".

We no longer provide "perturbation" method! 

When the trait is binary, the SKAT can produce conservative results when the sample size is small. To address this, we recently developed a small sample adjustment method, which adjust asymptotic null distribution by estimating small sample moments. See also SKAT_Null_Model_MomentAdjust.

                                                                      
}


\author{Seunggeun Lee}

\examples{


data(SKAT.example)
attach(SKAT.example)

#############################################################
#	Compute the P-value of SKAT 

# binary trait
obj<-SKAT_Null_Model(y.b ~ X, out_type="D")
SKAT(Z, obj, kernel = "linear.weighted")$p.value


#############################################################
# 	When you have no covariate to adjust.

# binary trait
obj<-SKAT_Null_Model(y.b ~ 1, out_type="D")
SKAT(Z, obj, kernel = "linear.weighted")$p.value



#########################################################
# Small sample adjustment
IDX<-c(1:100,1001:1100)

# With-adjustment
obj<-SKAT_Null_Model(y.b[IDX] ~ X[IDX,],out_type="D")
SKAT(Z[IDX,], obj, kernel = "linear.weighted")$p.value

# Without-adjustment
obj<-SKAT_Null_Model(y.b[IDX] ~ X[IDX,],out_type="D", Adjustment=FALSE)
SKAT(Z[IDX,], obj, kernel = "linear.weighted")$p.value


}


