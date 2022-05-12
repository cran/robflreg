\name{getPCA}
\alias{getPCA}
\title{Functional principal component analysis}
\description{
This function is used to perform functional principal component analysis.
}
\usage{
getPCA(data, nbasis, ncomp, gp, emodel = c("classical", "robust"))
}
\arguments{
  \item{data}{An \eqn{n \times p}-dimensional data matrix for functional data \eqn{X(s)}, where \eqn{n} denotes the sample size 
  and \eqn{p} denotes the number of grid points for \eqn{X(s)}.}
  \item{nbasis}{An integer specifying the number of B-spline basis expansion functions used to approximate the functional principal components.}
  \item{ncomp}{An integer specifying the number of functional principal components to be computed.}
  \item{gp}{A vector containing the grid points for the functional data for \eqn{X(s)}.}
  \item{emodel}{Method to be used for functional principal component decomposition. Possibilities are "classical"" and "robust".}
}
\details{The functional principal decomposition of a functional data \eqn{X(s)} is given by \deqn{X(s) = \bar{X}(s) + \sum_{k=1}^K \xi_k \psi_k(s),} where \eqn{\bar{X}(s)} is the mean function, \eqn{\psi_k(s)} is the \eqn{k}-th weight function, and \eqn{\xi_k} is the corresponding principal component score which is given by \deqn{\xi_k = \int (X(s) - \bar{X}(s)) \psi_k(s) ds.} When computing the estimated functional principal components, first, the functional data is expressed by a set of B-spline basis expansion. Then, the functional principal components are equal to the principal components extracted from the matrix \eqn{D \varphi^{1/2}}, where \eqn{D} is the matrix of basis expansion coefficients and \eqn{\varphi} is the inner product matrix of the basis functions, i.e., \eqn{\varphi = \int \varphi(s) \varphi^\top(s) ds}. Finally, the \eqn{k}-th weight function is given by \eqn{\psi_k(s) = \varphi^{-1/2} a_k}, where \eqn{a_k} is the \eqn{k}-th eigenvector of the sample covariance matrix of \eqn{D \varphi^{1/2}}.

If \code{emodel = "classical"}, then, the standard functional principal component decomposition is used as given by 
Ramsay and Dalzell (1991).

If \code{emodel = "robust"}, then, the robust principal component algorithm of Hubert, Rousseeuw and Verboven (2002) is used.
}
\value{A list object with the following components:
\item{PCAcoef}{A functional data object for the eigenfunctions.}
\item{PCAscore}{A matrix of principal component scores.}
\item{meanScore}{A functional data object for the mean function.}
\item{bs_basis}{A functional data object for B-spline basis expansion.}
\item{evalbase}{A matrix of the B-spline basis expansion functions.}
}
\references{
J. O. Ramsay and C. J. Dalzell (1991) "Some tools for functional data analysis (with discussion)", \emph{Journal of the Royal Statistical Society: Series B}, \bold{53}(3), 539-572.

M. Hubert and P. J. Rousseeuw and S. Verboven (2002) "A fast robust method for principal components with applications to chemometrics", \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{60}(1-2), 101-111.

P. Filzmoser and H. Fritz and K Kalcher (2021) pcaPP: Robust PCA by Projection Pursuit, R package version 1.9-74, URL: https://cran.r-project.org/web/packages/pcaPP/index.html.

J. L. Bali and G. Boente and D. E. Tyler and J.-L. Wang (2011) "Robust functional principal components: A projection-pursuit approach", \emph{The Annals of Statistics}, \bold{39}(6), 2852-2882.
}
\author{
Ufuk Beyaztas and Han Lin Shang
}
\examples{
sim.data <- generate.ff.data(n.pred = 5, n.curve = 200, n.gp = 101)
Y <- sim.data$Y
gpY <- seq(0, 1, length.out = 101) # grid points
rob.fpca <- getPCA(data = Y, nbasis = 20, ncomp = 4, gp = gpY, emodel = "robust")
}