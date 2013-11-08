\name{adhoc-package}
\alias{adhoc-package}
\alias{adhoc}
\docType{package}
\title{
Calculation of ad hoc distance thresholds for DNA barcoding identification.
}
\description{
The ad hoc distance threshold (see Virgilio et al. 2012) can be calculated for each particular reference library of DNA barcodes. Using this distance threshold in the "best close match method" (a method where a species name assignment is rejected when the distance between the query and its best match in a reference dataset is above the threshold) will provide identifications with an estimated relative error probability that is fixed by the user (e.g. 5\%).
}
\details{
\tabular{ll}{
Package: \tab adhoc\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-08\cr
License: \tab GPL-2 | GPL-3\cr
}
The procedure consists of two consecutive steps (Sonet et al. 2013). A first function (checkDNAbcd) evaluates a reference library in order to allow some quality control (sequence labelling, sequence lengths, etc.). The output of the first function is used by a second function (adhocTHR), which calculates an ad hoc distance threshold. For this, the second function takes each sequence of the reference library as a query and finds its best match(es) among all other DNA barcodes in the library. It then quantifies the relative identification errors (see Virgilio et al. 2012) obtained for a set of different arbitrary distance thresholds (30 values by default), performs a linear regression (or optionally a polynomial regression) and calculates the ad hoc distance threshold corresponding to an expected identification error probability (5\% by default).
}

\author{
Gontran Sonet
Maintainer: Gontran Sonet <gosonet@gmail.com>
}

\references{
Sonet G, Jordaens K, Nagy ZT, Breman FC, De Meyer M, Backeljau T, Virgilio M (2013) Adhoc: an R package to calculate ad hoc distance thresholds for DNA barcoding identification. Zookeys.

Virgilio M, Jordaens K, Breman FC, Backeljau T, De Meyer M. (2012) Identifying insects with incomplete DNA barcode libraries, African Fruit flies (Diptera: Tephritidae) as a test case. PLoS ONE 7(2) e31581. doi: 10.1371/journal.pone.0031581.
}

\examples{
data(tephdata);
out1<-checkDNAbcd(tephdata);
out2<-adhocTHR(out1);
layout(matrix(1,1,1));
par(font.sub=8);
plot(RE~thres,out2$IDcheck,ylim=c(0,max(c(out2$IDcheck$RE,out2$ErrProb))),xlab=NA,ylab=NA);
title(main="Ad hoc threshold",xlab="Distance", ylab="Relative identification error (RE)")
title(sub=paste("For a RE of",round(out2$ErrProb,4), "use a threshold of", round(out2$THR,4)));
regcoef<-out2$reg$coefficient;
curve(regcoef[1] + regcoef[2]*x + regcoef[3]*x^2 + regcoef[4]*x^3, add=TRUE);
segments(-1,out2$ErrProb,out2$THR,out2$ErrProb);segments(out2$THR,-1,out2$THR,out2$ErrProb);
}
