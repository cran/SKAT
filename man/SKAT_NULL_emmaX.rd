 \name{SKAT_NULL_emmaX}
 \alias{SKAT_NULL_emmaX}
 \title{Get parameters and residuals from the null model with kinship matrix}
 \description{
     Compute model parameters and residuals for SKAT with incorporating kinship structure. 
 }
 \usage{

SKAT_NULL_emmaX (formula, data=NULL, K=NULL, 
Kin.File=NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10) 
  
 }
\arguments{
      \item{formula}{an object of class ``formula'': a symbolic description of the NULL model to be fitted.}
      \item{data}{an optional data frame containing the variables in the model (default=NULL).  If it is NULL, the variables are taken from 'environment(formula)'}
      \item{K}{a kinship matrix. If K=NULL, the function reads Kin.File to get a kinship matrix.}
      \item{Kin.File}{an emmax-kin output file name. If K=NULL, the function reads this file. }
  	  \item{ngrids}{Number of grids to search optimal variance component}
      \item{llim}{Lower bound of log ratio of two variance components}
      \item{ulim}{Upper bound of log ratio of two variance components}
      \item{esp}{Tolerance of numerical precision error}

}
\value{
	This function returns an object that has model parameters and residuals of the NULL model of no association between genetic variables and outcome phenotypes. 
	After obtaining it, use SKAT_emmaX function to carry out the association test.

}
\details{

The Emma package code was used to implement this function.

Resampling is not implemented. 
                                                                      
}
\examples{


data(SKAT.fam.example)
attach(SKAT.fam.example)


obj<-SKAT_NULL_emmaX(y ~ X, K=K)
SKAT(Z, obj)$p.value
	

}



\author{Seunggeun Lee}
