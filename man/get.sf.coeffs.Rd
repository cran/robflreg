\name{get.sf.coeffs}
\alias{get.sf.coeffs}
\title{Get the estimated regression coefficient functions for scalar-on-function regression model}
\description{This function is used to obtain the estimated regression coefficient functions \eqn{\beta_m(s)} and the estimated regression coefficients \eqn{\gamma_r} (if \eqn{X.scl \neq NULL}) for scalar-on-function regression model (see the description in \code{\link{rob.sf.reg}} based on output object obtained from \code{\link{rob.sf.reg}}).
}
\usage{
get.sf.coeffs(object)
}
\arguments{
  \item{object}{The output object of \code{\link{rob.sf.reg}}.}
}
\details{In the estimation of regression coefficient functions, the estimated functional principal components of predictor \eqn{\hat{\Psi}_m(s), 1\le m\le M} variables and the estimated regression parameter function obtained from the regression model of scalar response on the principal component scores of the functional predictor variables \eqn{\hat{B}} are used, i.e., \eqn{\hat{\beta}_m(s) = \hat{\Psi}_m^\top(s) \hat{B}}.}
\value{A list object with the following components:
\item{gp}{A list with length \eqn{M}. The \eqn{m}-th element of gp is a vector containing the grid points of  the \eqn{m}-th functional predictor \eqn{X_m(s)}.}
\item{coefficients}{A list with length \eqn{M}. The \eqn{m}-th element of coefficients is a vector of the estimated values of the coefficient function for the \eqn{m}-th functional predictor \eqn{X_m(s)}.}
\item{scl.coefficients}{A vector consisting of the estimated coefficients of the scalar predictor \eqn{X.scl}.}
}
\author{Ufuk Beyaztas and Han Lin Shang}
\examples{
sim.data <- generate.sf.data(n = 400, n.pred = 5, n.gp = 101)
Y <- sim.data$Y
X <- sim.data$X
gp <- rep(list(seq(0, 1, length.out = 101)), 5) # grid points of Xs
model.fit <- rob.sf.reg(Y, X, emodel = "classical", gp = gp)
coefs <- get.sf.coeffs(model.fit)
}