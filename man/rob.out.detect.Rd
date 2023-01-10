\name{rob.out.detect}
\alias{rob.out.detect}
\title{Outlier detection in the functional response}
\description{
This function is used to detect outliers in the functional response based on a fitted function-on-function regression model in the output of \code{\link{rob.ff.reg}}.
}
\usage{
rob.out.detect(object, alpha = 0.01, B = 200, fplot = FALSE)
}
\arguments{
  \item{object}{An output object obtained from \code{\link{rob.ff.reg}}.}
  \item{alpha}{Percentile of the distribution of the functional depth. The default value is 0.01.}
  \item{B}{The number of bootstrap samples. The default value is 200.}
  \item{fplot}{If \code{TRUE}, then the outlying points flagged by the method is plotted along with the values of functional response \eqn{Y(t)}.}
}
\details{
The functional depth-based outlier detection method of Febrero-Bande et al. (2008) together with the h-modal depth proposed by Cuaves et al. (2007) is applied to the estimated residual functions obtained from \code{\link{rob.ff.reg}} to determine the outliers in the response variable. This method makes it possible to determine both magnitude and shape outliers in the response variable Hullait et al., (2021).
}
\value{
A vector containing the indices of outlying observations in the functional response.
}
\references{
M. Febrero-Bande and P. Galeano and W. Gonzalez-Mantelga (2008), "Outlier detection in functional data by depth measures, with application to identify abnormal NOx levels", \emph{Environmetrics}, \bold{19}(4), 331-345.

A. Cuaves and M. Febrero and R Fraiman (2007), "Robust estimation and classification for functional data via projection-based depth notions", \emph{Computational Statistics}, \bold{22}(3), 481-496.

H. Hullait and D. S. Leslie and N. G. Pavlidis and S. King (2021), "Robust function-on-function regression", \emph{Technometrics}, \bold{63}(3), 396-409.
}
\author{
Ufuk Beyaztas and Han Lin Shang
}
\examples{
\donttest{
sim.data <- generate.ff.data(n.pred = 5, n.curve = 200, n.gp = 101, out.p = 0.1)
out.indx <- sim.data$out.indx
Y <- sim.data$Y
X <- sim.data$X
gpY = seq(0, 1, length.out = 101) # grid points of Y
gpX <- rep(list(seq(0, 1, length.out = 101)), 5) # grid points of Xs

model.MM <- rob.ff.reg(Y = Y, X = X, model = "full", emodel = "robust", fmodel = "MM", 
                       gpY = gpY, gpX = gpX)
rob.out.detect(object = model.MM, fplot = TRUE)
sort(out.indx)
}
}
