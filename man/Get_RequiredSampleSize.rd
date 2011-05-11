 \name{Get_RequiredSampleSize}
 \alias{Get_RequiredSampleSize}
 \title{Get the required sample size to achieve the given power}
 \description{
     Get the sample sizes required to achieve the given power.
 }
 \usage{
	Get_RequiredSampleSize(obj, Power=0.8)
 }
\arguments{
      \item{obj}{an object returned from Power_Continuous or Power_Logistic.}
      \item{Power}{a value of the power to be achieved (default= 0.8).}
      
}
\value{
	A list object of the required sample sizes.
  	
}
\details{
	It computes required sample sizes by simple interpolation. 
                                                            
}


\author{Seunggeun Lee}




