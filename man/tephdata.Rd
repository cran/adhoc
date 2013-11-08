\name{tephdata}
\alias{tephdata}
\docType{data}
\title{DNA sequences of Tephritids to use as an example dataset for package adhoc.}
\description{DNA sequences (326 sequences of 656 base pairs) of the cytochrome c oxidase subunit I (COI) gene of Tephritids (Diptera) (see Virgilio et al. 2012) to use as an example dataset for package adhoc (Sonet et al. 2013).
}
\usage{data(tephdata)}

\format{
326 DNA sequences in binary format stored in a matrix.
}

\references{
Sonet G, Jordaens K, Nagy ZT, Breman FC, De Meyer M, Backeljau T, Virgilio M (2013) Adhoc: an R package to calculate ad hoc distance thresholds for DNA barcoding identification. Zookeys.

Virgilio M, Jordaens K, Breman FC, Backeljau T, De Meyer M. (2012) Identifying insects with incomplete DNA barcode libraries, African Fruit flies (Diptera: Tephritidae) as a test case. PLoS ONE 7(2) e31581. doi: 10.1371/journal.pone.0031581.
}

\examples{
data(tephdata);
labels(tephdata);
}

\keyword{datasets}
