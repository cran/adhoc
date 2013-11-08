\name{adhocTHR}
\alias{adhocTHR}
\title{Calculation of an ad hoc distance threshold for DNA barcoding identification}
\description{This function utilises the output of \code{\link{checkDNAbcd}} to quantify the relative identification errors obtained for a particular reference library at set of arbitrary distance thresholds. It then performs linear (or polynomial) regression and calculates the ad hoc distance threshold corresponding to an expected relative identification error. The method is described by Virgilio et al. (2012).}
\usage{adhocTHR(a, NbrTh = 30, ErrProb = 0.05, Ambig = "incorrect", Reg = "linear")}
\arguments{
  \item{a}{ the output of the function checkDNAbcd. }
  \item{NbrTh}{ the number of arbitrary distance thresholds used for the fitting (30 by default). }
  \item{ErrProb}{ a value between 0 and 1 setting the relative identification error probability (5\% by default). }
  \item{Ambig}{ "incorrect", "correct" or "ignore" to treat ambiguous identifications as incorrect, correct, or to ignore them (default is "incorrect"). }
  \item{Reg}{ "linear" or " polynomial" to perform a linear or a polynomial fitting (default is linear). }
}
\details{
NbrTh: By default, the function samples 30 distance thresholds equidistantly distributed between 0.0 and the largest distance between all pairs of query - best match observed in the reference library. Using more arbitrary distance thresholds may improve the fitting but will also require more computing time.

ErrProb: By default, the function sets an estimated error probability of 5\%. In some cases, the selected relative identification error (e.g. 5\%) cannot be reached, even when using the most restrictive distance threshold (viz. when setting distance threshold = 0.00). See below.

Ambig: We recommend to use this option with caution. If "correct", it will treat all red-flagged species involved in the same ambiguous identification as a single species. Setting this option to "ignore" will allow the user to evaluate the consequences of the ambiguous identifications on the relative identification errors and hence, on the estimated ad hoc distance threshold.

Reg: The user has the possibility to estimate the ad hoc distance threshold on the basis of polynomial regression (rather than linear). For this, check that the package polynom (Venables et al. 2013) is installed.
}
\value{\code{"adhocTHR"} returns a list of 7 components:
\item{BM}{ a data.frame describing the best matches obtained for each query: distance (distBM), identification (idBM), evaluation of the identification (IDeval) according to Sonet et al. (2013). }
\item{IDcheck}{ a data.frame describing the identification success at each arbitrary threshold: numbers of true positives (TP), false positives (FP), true negatives (TN) and false negatives (FN), the relative identification error (RE), the overall identification error (OE), the accuracy and the precision (Sonet et al. 2013). }
\item{reg}{ an object of class \code{"lm"} describing the regression performed by \code{"adhocTHR"}. }
\item{ErrProb}{ the value of the relative identification error probability fixed as an argument by the user (default is 5\%). }
\item{THR}{ the value of the ad hoc distance threshold producing the desired probability (ErrProb) of relative identification error. }
\item{redflagged}{ lists all ambiguous matches obtained for each ambiguous identification: the number of species involved in the ambiguous identification (Nb_species), the number of reference sequences with the same species name as the query (Nb_conspecific_seq), the number of reference sequences with a different species name than the query (Nb_heteroallospecific_seq) and the labels of all best matches involved in the ambiguous identification (List_of_best-matches). }
\item{redflaggedSP}{ lists all species names that are involved in the same ambiguous identification. }
}
\references{ 
Sonet G, Jordaens K, Nagy ZT, Breman FC, De Meyer M, Backeljau T, Virgilio M (2013) Adhoc: an R package to calculate ad hoc distance thresholds for DNA barcoding identification. Zookeys.

Venables B, Hornik K, Maechler M (2013) polynom: a collection of functions to implement a class for univariate polynomial manipulations. R package version 1.3-7. http://cran.r-project.org/web/packages/polynom/index.html.

Virgilio M, Jordaens K, Breman FC, Backeljau T, De Meyer M. (2012) Identifying insects with incomplete DNA barcode libraries, African Fruit flies (Diptera: Tephritidae) as a test case. PLoS ONE 7(2) e31581. doi: 10.1371/journal.pone.0031581.
}
\author{
Gontran Sonet
}
\note{ 
In particular cases, (e.g. reference libraries with low taxon coverage) all best matches might result in correct identifications, with RE = 0.0 at all distance thresholds and the regression line being parallel to the x axis. In this situation, adhocTHR will give the following message: "All identifications are correct when using the best match method (no distance threshold considered). An ad hoc threshold for best close match identification cannot be calculated".

In other cases, reaching an estimated RE of 5\% might not be possible, even at the most restrictive distance threshold (distance threshold = 0.00) and regression fitting will intercept the y axis above the RE value. In this case, adhocTHR will give the following message: "The estimated RE cannot be reached using this reference library". The user should then either increase the relative error probability (RE) or try and increase the taxon coverage of the library (see Virgilio et al. 2012).
}
\seealso{
\code{\link{checkDNAbcd}}
}

\examples{
data(tephdata);
out1<-checkDNAbcd(tephdata);
out2<-adhocTHR(out1);
layout(matrix(1,1,1));
par(font.sub=8);
plot(RE~thres,out2$IDcheck,ylim=c(0,max(c(out2$IDcheck$RE,out2$ErrProb))),xlab=NA,ylab=NA);
title(main="Ad hoc threshold",xlab="Distance", ylab="Relative identification error (RE)");
title(sub=paste("For a RE of",round(out2$ErrProb,4), "use a threshold of", round(out2$THR,4)));
regcoef<-out2$reg$coefficient;
curve(regcoef[1] + regcoef[2]*x + regcoef[3]*x^2 + regcoef[4]*x^3, add=TRUE);
segments(-1,out2$ErrProb,out2$THR,out2$ErrProb);segments(out2$THR,-1,out2$THR,out2$ErrProb);
}
