% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allele_cluster.R
\name{allele_diff}
\alias{allele_diff}
\title{Alleles nucleotide position difference}
\usage{
allele_diff(
  reference_allele,
  sample_allele,
  position_threshold = 0,
  snps = TRUE
)
}
\arguments{
\item{reference_allele}{The nucleotide sequence of the reference allele, character object.}

\item{sample_allele}{The nucleotide sequence of the sample allele, character object.}

\item{position_threshold}{A position from which to check for differential positions. If zero checks all position. Default to zero.}

\item{snps}{If to return the SNP with the position (e.g., A2G where A is for the reference and G is for the sample.). If false returns just the positions. Default to True}
}
\value{
A \code{character} vector of the differential nucleotide positions of the sample allele.
}
\description{
Compare the sequences of two alleles (reference and sample alleles)
and returns the differential nucleotide positions of the sample allele.
}
\details{
The function utilizes c++ script to optimize the run time for large comparisons.
}
\examples{
{
reference_allele = "AAGG"
sample_allele = "ATGA"

# setting position_threshold = 0 will return all differences
diff <- allele_diff(reference_allele, sample_allele)
# "A2T", "G4A"
print(diff)

# setting position_threshold = 3 will return the differences from position three onward
diff <- allele_diff(reference_allele, sample_allele, position_threshold = 3)
# "G4A"
print(diff)

# setting snps = FALSE will return the differences as indices
diff <- allele_diff(reference_allele, sample_allele, snps = FALSE)
# 2, 4
print(diff)

}
}
