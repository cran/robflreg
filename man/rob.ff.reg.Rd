\name{rob.ff.reg}
\alias{rob.ff.reg}
\title{Robust function-on-function regression}
\description{
This function is used to perform both classical and robust function-on-function regression model \deqn{
Y(t) = \sum_{m=1}^M X_m(s) \beta_m(s,t) ds + \epsilon(t),} where \eqn{Y(t)} denotes the functional response, \eqn{X_m(s)} denotes the \eqn{m}-th functional predictor, \eqn{\beta_m(s,t)} denotes the \eqn{m}-th bivariate regression coefficient function, and \eqn{\epsilon(t)} is the error function.
}
\usage{
rob.ff.reg(Y, X, model = c("full", "selected"), emodel = c("classical", "robust"), 
 fmodel = c("MCD", "MLTS", "MM", "S", "tau"), nbasisY = NULL, nbasisX = NULL, 
 gpY = NULL, gpX = NULL, ncompY = NULL, ncompX = NULL)
}
\arguments{
  \item{Y}{An \eqn{n \times p}-dimensional matrix containing the observations of functional response \eqn{Y(t)}, where \eqn{n} is the sample size and \eqn{p} denotes the number of grid points for \eqn{Y(t)}.}
  \item{X}{A list consisting of \eqn{M} functional predictors \eqn{X_m(s), 1\le m\le M}. Each element of \code{X} is an \eqn{n \times p_m}-dimensional matrix containing the observations of \eqn{m}-th functional predictor \eqn{X_m(s)}, where \eqn{n} is the sample size and \eqn{p_m} denotes the number of grid points for \eqn{X_m(s)}.}
  \item{model}{Model to be fitted. Possibilities are "full" and "selected".}
  \item{emodel}{Method to be used for functional principal component decomposition. Possibilities are "classical"" and "robust".}
  \item{fmodel}{Fitting model used to estimate the function-on-function regression model. Possibilities are "MCD", "MLTS", "MM", "S", and "tau".}
  \item{nbasisY}{An integer value specifying the number of B-spline basis expansion functions to be used to approximate the functional principal components for the response variable \eqn{Y(t)}. If \code{NULL}, then, 20 number of B-spline basis expansion functions are used.}
  \item{nbasisX}{A vector with length \eqn{M}. Its \eqn{m}-th value denotes the number of B-spline basis expansion functions to be used to approximate the functional principal components for the \eqn{m}-th functional predictor \eqn{X_m(s)}. If \code{NULL}, then, 20 number of B-spline basis expansion functions are used for each functional predictor.}
  \item{gpY}{A vector containing the grid points of the functional response \eqn{Y(t)}. If \code{NULL}, then \eqn{p} 
  equally spaced time points in te interval [0, 1] are used.}
  \item{gpX}{A list with length \eqn{M}. The \eqn{m}-th element of \code{gpX} is a vector containing the grid points of 
the \eqn{m}-th functional predictor \eqn{X_m(s)}. If \code{NULL}, then, \eqn{p_m} equally spaced time points in te interval 
[0, 1] are used for the \eqn{m}-th functional predictor.}
  \item{ncompY}{An integer specifying the number of functional principal components to be computed for the functional response \eqn{Y(t)}. If \eqn{NULL}, then, 4 functional principal components are used.}
  \item{ncompX}{A vector with length \eqn{M}. Its \eqn{m}-th value denotes the number of functional principal components to be computed for the \eqn{m}-th functional predictor \eqn{X_m(s)}. If \code{NULL}, then, 4 principal components are used for each functional predictor.}
}
\details{
When performing a function-on-function regression model based on the functional principal component analysis, first, both the functional response \eqn{Y(t)} and functional predictors \eqn{X_m(s), 1\le m\le M} are decomposed by the functional principal component analysis method: \deqn{Y(t) = \bar{Y}(t) + \sum_{k=1}^K \nu_k \phi_k(t),} \deqn{X_m(s) = \bar{X}_m(s) + \sum_{l=1}^{K_m} \xi_{ml} \psi_{ml}(s),} where \eqn{\bar{Y}(t)} and \eqn{\bar{X}_m(s)} are the mean functions, \eqn{\phi_k(t)} and \eqn{psi_{ml}(s)} are the weight functions, and \eqn{\nu_k = \int (Y(t) - \bar{Y}(t)) \phi_k(t)} and \eqn{\xi_{ml} = \int (X_m(s) - \bar{X}_m(s)) \psi_{ml}(s)} are the principal component scores for the functional response and \eqn{m}-th functional predictor, respectively. Assume that the \eqn{m}-th bivariate regression coefficient function admits the expansion \deqn{\beta_m(s,t) = \sum_{k=1}^K \sum_{l=1}^{K_m} b_{mkl} \phi_k(t) \psi_{ml}(s),} where \eqn{b_{mkl} = \int \int \beta_m(s,t) \phi_k(t) \psi_{ml}(s) dt ds}. Then, the following multiple regression model is obtained for the functional response: \deqn{\hat{Y}(t) = \bar{Y}(s) + \sum_{k=1}^K ( \sum_{m=1}^M \sum_{l=1}^{K_m} b_{mkl} \xi_{ml} ) \phi_k(t).}

If \code{model = "full"}, then, all the functional predictor variables are used in the model.

If \code{model = "selected"}, then, only the significant functional predictor variables determined by the forward variable selection procedure of Beyaztas and Shang (2021) are used in the model.

If \code{emodel = "classical"}, then, the least-squares method is used to estimate the function-on-function regression model.

If \code{emodel = "robust"}, then, the robust functional principal component analysis of Bali et al. (2011) along with the method specified in \code{fmodel} is used to estimate the  function-on-function regression model.

If \code{fmodel = "MCD"}, then, the minimum covariance determinant estimator of Rousseeuw et al. (2004) is used to estimate the function-on-function regression model.

If \code{fmodel = "MLTS"}, then, the multivariate least trimmed squares estimator Agullo et al. (2008) is used to estimate the function-on-function regression model.

If \code{fmodel = "MM"}, then, the MM estimator of Kudraszow and Maronna (2011) is used to estimate the function-on-function regression model.

If \code{fmodel = "S"}, then, the S estimator of Bilodeau and Duchesne (2000) is used to estimate the function-on-function regression model.

If \code{fmodel = "tau"}, then, the tau estimator of Ben et al. (2006) is used to estimate the function-on-function regression model.
}
\value{A list object with the following components:
\item{data}{A list of matrices including the original functional response and functional predictors.}
\item{fitted.values}{An \eqn{n \times p}-dimensional matrix containing the fitted values of the functional response.}
\item{residuals}{An \eqn{n \times p}-dimensional matrix containing the residual functions.}
\item{fpca.results}{A list object containing the functional principal component analysis results of the functional predictor and functional predictors variables.}
\item{model.details}{A list object containing model details, such as number of basis functions, number of principal components, and grid points used for each functional variable.}
}
\references{
J. Agullo and C. Croux and S. V. Aelst (2008), "The multivariate least-trimmed squares estimator", \emph{Journal of Multivariate Analysis}, \bold{99}(3), 311-338.

M. G. Ben and E. Martinez and V. J. Yohai (2006), "Robust estimation for the multivariate linear model based on a \eqn{\tau} scale", \emph{Journal of Multivariate Analysis}, \bold{97}(7), 1600-1622.

U. Beyaztas and H. L. Shang (2021), "A partial least squares approach for function-on-function interaction regression", \emph{Computational Statistics}, \bold{36}(2), 911-939.

J. L. Bali and G. Boente and D. E. Tyler and J. -L.Wang (2011), "Robust functional principal components: A projection-pursuit approach", \emph{The Annals of Statistics}, \bold{39}(6), 2852-2882.

M. Bilodeau and P. Duchesne (2000), "Robust estimation of the SUR model", \emph{The Canadian Journal of Statistics}, \bold{28}(2), 277-288.

N. L. Kudraszow and R. A. Moronna (2011), "Estimates of MM type for the multivariate linear model", \emph{Journal of Multivariate Analysis}, \bold{102}(9), 1280-1292.

P. J. Rousseeuw and K. V. Driessen and S. V. Aelst and J. Agullo (2004), "Robust multivariate regression", \emph{Technometrics}, \bold{46}(3), 293-305.
}
\author{
Ufuk Beyaztas and Han Lin Shang
}
\examples{
sim.data <- generate.ff.data(n.pred = 5, n.curve = 200, n.gp = 101)
Y <- sim.data$Y
X <- sim.data$X
gpY <- seq(0, 1, length.out = 101) # grid points of Y
gpX <- rep(list(seq(0, 1, length.out = 101)), 5) # grid points of Xs
model.MM <- rob.ff.reg(Y = Y, X = X, model = "full", emodel = "robust", 
                       fmodel = "MM", gpY = gpY, gpX = gpX)
}

