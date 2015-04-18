 \name{Get_EffectiveNumberTest}
 \alias{Get_EffectiveNumberTest}
 \title{Estimate the effective number of tests for Bonferroni correction}
 \description{
     Estimate the effective number of tests using the minimum achievable p-values (MAP).
 }
 \usage{
 
 
	Get_EffectiveNumberTest(MAP, alpha=0.05, Is.MidP=TRUE) 


 }
\arguments{
      \item{MAP}{a vector of the minimum achievable p-values (MAP).}
      \item{alpha}{a significant level.}
      \item{Is.MidP}{a logical value indicating whether the p-values are mid-p-values.}
}
\value{
	 Effective number of test for Bonferroni correction at the level alpha. 
}

\author{Seunggeun Lee}
