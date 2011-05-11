 \name{Get_Genotypes_SSD}
 \alias{Get_Genotypes_SSD}
 \title{Get Genotype data from SSD file}
 \description{
	Read SSD file and return a genotype matrix.
 }
 \usage{
	Get_Genotypes_SSD(SSD_INFO, Set_Index)
 }
\arguments{
      \item{SSD_INFO}{a SSD_INFO object returned from Open_SSD.   }
      \item{Set_Index}{a numeric value of Set index. You can find a set index of each set from SetInfo object of SSD.INFO. }
      
}
\value{
 	The genotype matrix with n rows and m columns, where n is the number of samples, and m is the number of SNPs.
}
\author{Seunggeun Lee, Larisa Miropolsky}

