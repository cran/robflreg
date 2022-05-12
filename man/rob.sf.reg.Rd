\name{rob.sf.reg}
\alias{rob.sf.reg}
\title{Robust scalar-on-function regression}
\description{This function is used to perform both classical and robust scalar-on-function regression model \deqn{
Y = \sum_{m=1}^M X_m(s) \beta_m(s) ds + \epsilon,} where \eqn{Y} denotes the scalar response, \eqn{X_m(s)} denotes the \eqn{m}-th functional predictor, \eqn{\beta_m(s)} denotes the \eqn{m}-th regression coefficient function, and \eqn{\epsilon} is the error function, which is assumed to follow standard normal distribution.
}
\usage{
rob.sf.reg(Y, X, emodel = c("classical", "robust"),
fmodel = c("LTS", "MM", "S", "tau"), nbasis = NULL, gp = NULL, ncomp = NULL)
}
\arguments{
  \item{Y}{An \eqn{n \times 1}-dimensional matrix containing the observations of scalar response \eqn{Y}, where \eqn{n} denotes the sample size.}
  \item{X}{A list consisting of \eqn{M} functional predictors \eqn{X_m(s), 1\le m\le M}. Each element of \code{X} is an \eqn{n \times p_m}-dimensional matrix containing the observations of \eqn{m}-th functional predictor \eqn{X_m(s)}, where \eqn{n} is the sample size and \eqn{p_m} denotes the number of grid points for \eqn{X_m(s)}.}
  \item{emodel}{Method to be used for functional principal component decomposition. Possibilities are "classical"" and "robust".}
  \item{fmodel}{Fitting model used to estimate the function-on-function regression model. Possibilities are "LTS", "MM", "S", and "tau".}
  \item{nbasis}{A vector with length \eqn{M}. Its \eqn{m}-th value denotes the number of B-spline basis expansion functions to be used to approximate the functional principal components for the \eqn{m}-th functional predictor \eqn{X_m(s)}. If \code{NULL}, then, 20 number of B-spline basis expansion functions are used for each functional predictor.}
  \item{gp}{A list with length \eqn{M}. The \eqn{m}-th element of \code{gp} is a vector containing the grid points of the \eqn{m}-th functional predictor \eqn{X_m(s)}. If \code{NULL}, then, \eqn{p_m} equally spaced time points in the interval [0, 1] are used for the \eqn{m}-th functional predictor.}
  \item{ncomp}{A vector with length \eqn{M}. Its \eqn{m}-th value denotes the number of functional principal components to be computed for the \eqn{m}-th functional predictor \eqn{X_m(s)}. If \code{NULL}, then, four principal components are used for each functional predictor.}
}
\details{
When performing a scalar-on-function regression model based on the functional principal component analysis, first, the functional predictors \eqn{X_m(s), 1\le m\le M} are decomposed by the functional principal component analysis method: \deqn{X_m(s) = \bar{X}_m(s) + \sum_{l=1}^{K_m} \xi_{ml} \psi_{ml}(s),} where \eqn{\bar{X}_m(s)} is the mean function, \eqn{\psi_{ml}(s)} is the weight function, and \eqn{\xi_{ml} = \int (X_m(s) - \bar{X}_m(s)) \psi_{ml}(s)} is the principal component score for the \eqn{m}-th functional predictor. Assume that the \eqn{m}-th regression coefficient function admits the expansion \deqn{\beta_m(s) = \sum_{l=1}^{K_m} b_{ml} \psi_{ml}(s),} where \eqn{b_{ml} = \int \beta_m(s) \psi_{m}(s)  ds}. Then, the following multiple regression model is obtained for the scalar response: \deqn{\hat{Y} = \bar{Y} + \sum_{m=1}^M \sum_{l=1}^{K_m} b_{ml} \xi_{ml}.}

If \code{emodel = "classical"}, then, the least-squares method is used to estimate the scalar-on-function regression model.

If \code{emodel = "robust"}, then, the robust functional principal component analysis of Bali et al. (2011) along with the method specified in \code{fmodel} is used to estimate the  scalar-on-function regression model.

If \code{fmodel = "LTS"}, then, the least trimmed squares robust regression of Rousseeuw (1984) is used to estimate the scalar-on-function regression model.

If \code{fmodel = "MM"}, then, the MM-type regression estimator described in Yohai (1987) and Koller and Stahel (2011) is used to estimate the scalar-on-function regression model.

If \code{fmodel = "S"}, then, the S estimator is used to estimate the scalar-on-function regression model.

If \code{fmodel = "tau"}, then, the tau estimator proposed by Salibian-Barrera et al. (2008) is used to estimate the scalar-on-function regression model.
}
\value{A list object with the following components:
\item{data}{A list of matrices including the original scalar response and functional predictors.}
\item{fitted.values}{An \eqn{n \times 1}-dimensional matrix containing the fitted values of the scalar response.}
\item{residuals}{An \eqn{n \times 1}-dimensional matrix containing the residuals.}
\item{fpca.results}{A list object containing the functional principal component analysis results of the functional predictors variables.}
\item{model.details}{A list object containing model details, such as number of basis functions, number of principal components, and grid points used for each functional predictor variable.}
}
\references{
J. L. Bali and G. Boente and D. E. Tyler and J. -L.Wang (2011), "Robust functional principal components: A projection-pursuit approach", \emph{The Annals of Statistics}, \bold{39}(6), 2852-2882.

P. J. Rousseeuw (1984), "Least median of squares regression", \emph{Journal of the American Statistical Association}, \bold{79}(388), 871-881.

P. J. Rousseeuw and K. van Driessen (1999) "A fast algorithm for the minimum covariance determinant estimator", \emph{Technometrics}, \bold{41}(3), 212-223.

V. J. Yohai (1987), "High breakdown-point and high efficiency estimates for regression", \emph{The Annals of Statistics}, \bold{15}(2), 642-65.

M. Koller and W. A. Stahel (2011), "Sharpening Wald-type inference in robust regression for small samples", \emph{Computational Statistics & Data Analysis}, \bold{55}(8), 2504-2515.

M. Salibian-Barrera and G. Willems and R. Zamar (2008), "The fast-tau estimator for regression", \emph{Journal of Computational and Graphical Statistics}, \bold{17}(3), 659-682 
}
\author{Ufuk Beyaztas and Han Lin Shang}
\examples{
sim.data <- generate.sf.data(n = 400, n.pred = 5, n.gp = 101)
Y <- sim.data$Y
X <- sim.data$X
gp <- rep(list(seq(0, 1, length.out = 101)), 5) # grid points of Xs
model.tau <- rob.sf.reg(Y, X, emodel = "robust", fmodel = "tau", gp = gp)
}