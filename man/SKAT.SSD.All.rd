 \name{SKAT.SSD.All}
 \alias{SKAT.SSD.All}
 \title{SNP-set Kernel Association Test}
 \description{
	Iteratively conduct association tests with phenotypes and SNP sets in 
	SSD file. 
 }
 \usage{

	SKAT.SSD.All(SSD.INFO, obj, \dots)
 }
\arguments{

      \item{SSD.INFO}{an SSD_INFO object returned from Open_SSD.   }
      \item{obj}{an output object of the SKAT_Null_Model function. }
      \item{\dots}{ furthuer arguments to be passed to ``SKAT''. }
}
\value{
	\item{results}{the dataframe that contains SetID, p-values (P.value), the number of markers in the SNP sets (N.Marker.All), and the number of markers to test for an association after excluding non-polymorphic or high missing rates markers (N.Marker.Test).  }
	\item{P.value.Resampling}{the matrix that contains p-values of resampled phenotypes. }
}
\details{
Please see SKAT for details.                     

}


\author{Seunggeun Lee}

