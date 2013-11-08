\name{checkDNAbcd}
\alias{checkDNAbcd}
\title{ Evaluation of a reference library of aligned DNA barcodes 
}
\description{
This function provides an overview of the content of a reference library of aligned DNA barcodes (Sonet et al. 2013). It calculates all pairwise distances and delivers an output that can be used by the function adhocTHR.
}
\usage{
checkDNAbcd(seq, DistModel = "K80")
}

\arguments{
  \item{seq}{ an object of class \code{"DNAbin"}. }
  \item{DistModel}{ "K80" (for Kimura two-parameter) or "raw" (for p-distances) or any other nucleotide substitution model available in the function \code{"dist.dna"} (Paradis et al. 2004). }
}

\value{
checkDNAbcd returns a list of 6 components:
\item{mylabels}{ a data.frame providing both parts of the species names and the complete label of each sequence (as extracted from the first argument). }
\item{listsp}{ a data.frame listing the number of sequences (Nseq) and haplotypes (Nhap) for each species of the reference library. }
\item{DNAlength}{ a numeric vector of the sequence lengths of each DNA sequence. }
\item{dist}{ a matrix of all distances obtained by pairwise comparison. }
\item{spdist}{ a list of all pairwise interspecific distances (inter) and all pairwise intraspecific distances (intra). }
\item{seq}{ an object of class \code{"DNAbin"} with all sequences in the reference library (= first argument). }
}
\references{
Paradis E, Claude J, Strimmer K (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.

Sonet G, Jordaens K, Nagy ZT, Breman FC, De Meyer M, Backeljau T, Virgilio M (2013) Adhoc: an R package to calculate ad hoc distance thresholds for DNA barcoding identification. Zookeys.
}
\author{
Gontran Sonet
}

\seealso{ \code{\link{adhocTHR}} }

\examples{
data(tephdata);
out1<-checkDNAbcd(tephdata);
layout(matrix(c(1:4),1,4));
hist(out1$DNAlength,main="Seq. lengths",xlab="Seq. length (bp)");
hist(out1$spdist$inter,main="Intersp. dist",xlab="Distance",col="#0000ff99");
hist(out1$spdist$intra, main="Intrasp. dist.",xlab="Distance",col="#0000ff22");
hist(out1$spdist$inter,main="Intra- & intersp. dist",xlab="Distance",col="#0000ff99");
hist(out1$spdist$intra, add=TRUE,col="#0000ff22");
}
