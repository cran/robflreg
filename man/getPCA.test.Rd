\name{getPCA.test}
\alias{getPCA.test}
\title{Get the functional principal component scores for a given test sample}
\description{
This function is used to compute the functional principal component scores of a test sample based on outputs obtained from \code{\link{getPCA}}.
}
\usage{
getPCA.test(data, bs_basis, PCAcoef, gp, emodel = c("classical", "robust"))
}
\arguments{
  \item{data}{An \eqn{n \times p}-dimensional data matrix for functional data \eqn{X(s)} (test sample), where \eqn{n} denotes the sample size and \eqn{p} denotes the number of grid points for \eqn{X(s)}.}
  \item{bs_basis}{A functional data object for B-spline basis expansion obtained from \code{\link{getPCA}}.}
  \item{PCAcoef}{A functional data object for the eigenfunctions obtained from \code{\link{getPCA}}.}
  \item{gp}{A vector containing the grid points for the functional data for \eqn{X(s)}.}
  \item{emodel}{Method to be used for functional principal component decomposition. Possibilities are "classical"" and "robust".}
}
\details{
See \code{\link{getPCA}} for details.
}
\value{A matrix of principal component scores for the functional data.
}
\author{
Ufuk Beyaztas and Han Lin Shang
}
\examples{
sim.data <- generate.ff.data(n.pred = 5, n.curve = 200, n.gp = 101)
Y <- sim.data$Y
Y.train <- Y[1:100,]
Y.test <- Y[101:200,]
gpY = seq(0, 1, length.out = 101) # grid points
rob.fpca <- getPCA(data = Y.train, nbasis = 20, ncomp = 4, gp = gpY, emodel = "robust")
bs_basis <- rob.fpca$bs_basis
PCAcoef <- rob.fpca$PCAcoef
rob.fpca.test <- getPCA.test(Y.test, bs_basis, PCAcoef, gpY, emodel = "robust")
}