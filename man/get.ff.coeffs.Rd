\name{get.ff.coeffs}
\alias{get.ff.coeffs}
\title{Get the estimated bivariate regression coefficient functions for function-on-function regression model}
\description{This function is used to obtain the estimated bivariate regression coefficient functions \eqn{\beta_m(s,t)} for function-on-function regression model (see the description in \code{\link{rob.ff.reg}} based on output object obtained from \code{\link{rob.ff.reg}}).
}
\usage{
get.ff.coeffs(object)
}
\arguments{
  \item{object}{The output object of \code{\link{rob.ff.reg}}.}
}
\details{In the estimation of bivariate regression coefficient functions, the estimated functional principal components of 
response \eqn{\hat{\Phi}(t)} and predictor \eqn{\hat{\Psi}_m(s)} variables and the estimated regression parameter function obtained from the regression model between the principal component scores of response and predictor variables \eqn{\hat{B}} 
are used, i.e., \eqn{\hat{\beta}_m(s,t) = \hat{\Psi}_m^\top(s) \hat{B} \hat{\Phi}(t)}.}
\value{A list object with the following components:
\item{vars}{A numeric vector specifying the indices of functional predictors used in the function-on-function regression model \code{\link{rob.ff.reg}}.
}
\item{gpY}{A vector containing the grid points of the functional response \eqn{Y(t)}.}
\item{gpX}{A list with length \eqn{M}. The \eqn{m}-th element of gpX is a vector containing the grid points of 
the \eqn{m}-th functional predictor \eqn{X_m(s)}.}
\item{coefficients}{A list with length \eqn{M}. The \eqn{m}-th element of coefficients is a matrix of the estimated values of the coefficient function for the \eqn{m}-th functional predictor \eqn{X_m(s)}.}
}
\author{
Ufuk Beyaztas and Han Lin Shang
}
\examples{
sim.data <- generate.ff.data(n.pred = 5, n.curve = 200, n.gp = 101)
Y <- sim.data$Y
X <- sim.data$X
gpY = seq(0, 1, length.out = 101) # grid points of Y
gpX <- rep(list(seq(0, 1, length.out = 101)), 5) # grid points of Xs
model.fit <- rob.ff.reg(Y, X, model = "full", emodel = "classical", 
                        gpY = gpY, gpX = gpX)
coefs <- get.ff.coeffs(model.fit)
}