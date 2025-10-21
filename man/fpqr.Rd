\name{fpqr}
\alias{fpqr}
\title{Function-on-function partial quantile regression}
\description{This function is used to perform function-on-function linear quantile regression model \deqn{
Q_{\tau}[Y(u)|X(v)] = \int X(v) \beta_{\tau}(u,v) dv} based on the functional partial quantile regression.}
\usage{
fpqr(y, x, h, tau, nby, nbx, gpy, gpx, qc.type = c("dodge","choi","li"),
hs, nbys, nbxs, nfold, CV)
}
\arguments{
  \item{y}{An matrix containing the observations of function response \eqn{Y(u)}.}
  \item{x}{A matrix containing the observations of functional predictor \eqn{X(v)}.}
  \item{h}{A numeric value denoting the number of functional partial quantile regression components to be computed.}
  \item{tau}{Quantile level.}
  \item{nby}{A numeric value denoting the number of B-spline basis expansion functions to be used to approximate the functional response variable.}
  \item{nbx}{A numeric value denoting the number of B-spline basis expansion functions to be used to approximate the functional predictor variable.}
  \item{gpy}{A vector containing the grid points of the functional response \eqn{Y(u)}.}
  \item{gpx}{A vector containing the grid points of the functional predictor \eqn{X(v)}.}
  \item{qc.type}{Method type used to estimate the quantile covariance. Possibilities are "dodge", "choi", and "li".}
  \item{hs}{A vector containing the candidate elements for the \code{h}.}
  \item{nbys}{A vector containing the candidate elements for the \code{nby}.}
  \item{nbxs}{A vector containing the candidate elements for the \code{nbx}.}
  \item{nfold}{An integer denoting the number of folds used in the k-fold cross validation. Default value is 5.}
  \item{CV}{Logical. If TRUE, then nfold cross-validation is used to find optimum values of \code{h}, \code{nby}, and \code{nbx}. If FALSE, then the specified \code{h}, \code{nby}, and \code{nbx} values are used in the model.}
}
\details{If qc.type = "dodge", then, the quantile covariance proposed by Dodge and Whittaker (2009) is used in the functional partial quantile regression decomposition.

If qc.type = "choi", then, the quantile covariance proposed by Choi and Shin (2018) is used in the functional partial quantile regression decomposition.

If qc.type = "li", then, the quantile covariance proposed by Li et al. (2015) is used in the functional partial quantile regression decomposition.
}
\value{
\item{fitted.values}{An matrix containing the fitted values of the functional response.}
\item{b0.hat}{A vector containing the estimate of intersecp parameter.}
\item{b.hat}{A matrix containing the estimate of bivariate regression coefficient coducted between the functional response and functional predictor.}
\item{mdts}{A numeric value containing the estimated intercept parameter.}
\item{pqr.coefs}{A vector containing the estimated regression parameter for the regression problem of scalar response on the partial quantile regression components.}
\item{V}{A matrix whose rows are the eigenvectors}
\item{model.details}{A list object containing model details, such as number of basis functions, number of partial quantile components, and grid points used for the functional variables.}
}
\references{
Y. Dodge and J. Whittaker (2009), "Partial quantile regression", \emph{Metrika}, \bold{70}(1), 35-57.
J. E.. Choi and D. W. Shin (2018), "Quantile correlation coefficient: A new tail dependence measure", \emph{Statistical Papers}, \bold{63}, 1075-1104.
G. Li and Y. Li and C. L. Tsai (2015), "Quantile correlations and quantile autoregressive modeling" \emph{Journal of American Statistical Association: Theory and Methods}, \bold{110}(509), 246-261.
}
\author{Muge Mutis, Ufuk Beyaztas, Filiz Karaman, and Han Lin Shang}
\examples{
\dontrun{
  gpx <- (1:50)/50
  gpy <- (1:60)/60
  data <- fpqr_dgp(n = 100, gpy = gpy, gpx = gpx, err.dist = "normal")
  y <- data$y
  x <- data$x

  fpqr.model.dodge <- fpqr(y=y, x=x, tau = 0.5, gpx = gpx, gpy = gpy,
                           qc.type = "dodge")

  fpqr.model.choi <- fpqr(y=y, x=x, tau = 0.5, gpx = gpx, gpy = gpy,
                          qc.type = "choi")

  fpqr.model.li <- fpqr(y=y, x=x, tau = 0.5, gpx = gpx, gpy = gpy,
                        qc.type = "li")
}
}
