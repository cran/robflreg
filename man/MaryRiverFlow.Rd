\name{MaryRiverFlow}
\alias{MaryRiverFlow}
\title{Hourly River Flow Measurements in the Mery River}
\description{Hourly river flow measurements obtained from January 2009 to December 2014 (6 years in total) in the Mery River, Australia.}
\usage{
data(MaryRiverFlow)
}
\author{Ufuk Beyaztas and Han Lin Shang}
\examples{
data(MaryRiverFlow)
# Plot
library(fda.usc)
fflow <- fdata(MaryRiverFlow, argvals = 1:24)
plot(fflow, lty = 1, ylab = "", xlab = "Hour",
main = "", mgp = c(2, 0.5, 0), ylim = range(fflow))
}
