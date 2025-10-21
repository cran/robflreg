\name{sff_dgp}
\alias{sff_dgp}
\title{Simulate data from a Spatial Function--on--Function Regression model}
\description{
Generates synthetic functional predictors and responses from the spatial
function-on-function regression (SFoFR) data-generating process described in
Beyaztas, Shang and Sezer (2025). The model embeds spatial autoregression on the functional response, Fourier-type
basis structure for the covariate, and user-controlled Gaussian noise.
}
\usage{
sff_dgp(
  n,
  nphi     = 10,
  gpy      = NULL,
  gpx      = NULL,
  rf       = 0.9,
  sd.error = 0.01,
  tol      = 0.001,
  max_iter = 1000
)
}
\arguments{
  \item{n}{Number of spatial units (curves) to generate.}
  \item{nphi}{Number of sine \emph{and} cosine basis functions used to build each functional predictor.  Total latent scores generated are therefore \code{2 * nphi}.}
  \item{gpy}{Numeric vector of evaluation points for the response domain \eqn{t \in [0,1]}. Defaults to an equally-spaced grid of 101 points.}
  \item{gpx}{Numeric vector of evaluation points for the predictor domain \eqn{s \in [0,1]}. Defaults to an equally-spaced grid of 101 points.}
  \item{rf}{Scalar in \eqn{(0,1)} controlling the strength of spatial autocorrelation through the surface \eqn{\rho(t,u)}. Values closer to 1 yield stronger dependence.}
  \item{sd.error}{Standard deviation of the i.i.d. Gaussian noise \eqn{\varepsilon_i(t)} added to the latent regression part.}
  \item{tol}{Absolute tolerance used in the fixed-point iteration that solves the spatial autoregressive operator equation (stopping rule on the sup-norm of successive iterates).}
  \item{max_iter}{Maximum number of fixed-point iterations. Prevents infinite looping when strong spatial feedback and small \code{tol} interact.}
}
\details{
The generator mimics the penalised SFoFR set-up:

\deqn{
  Y_i(t)
  \;=\;
  \sum_{j=1}^{n} w_{ij}\int_0^1 Y_j(u)\,\rho(t,u)\,du
  \;+\;
  \int_0^1 X_i(s)\,\beta(t,s)\,ds
  \;+\;
  \varepsilon_i(t),
}{%
Y_i(t)=\sum_j w_{ij}\int Y_j(u)\rho(t,u)du + \int X_i(s)\beta(t,s)ds + \varepsilon_i(t)
}

where

\itemize{
  \item \eqn{w_{ij}} are row-normalised inverse-distance weights,
  \item \eqn{X_i(s)} is built from Fourier scores
        \eqn{\xi_{ijk} \sim \mathcal{N}(0,1)} and
        damped basis functions
        \eqn{\phi_k^{\cos}(s)=(k^{-3/2})\sqrt{2}\cos(k\pi s)} and
        \eqn{\phi_k^{\sin}(s)=(k^{-3/2})\sqrt{2}\sin(k\pi s)},
  \item the regression surface is
        \eqn{\beta(t,s)=2+s+t+0.5\sin(2\pi s t)},
  \item the spatial autocorrelation surface is
        \eqn{\rho(t,u)=rf\,(1+ut)/(1+|u-t|)},
  \item \eqn{\varepsilon_i(t)\stackrel{\text{i.i.d.}}{\sim}\mathcal{N}(0,\sigma^2)},
        with \eqn{\sigma=\code{sd.error}}.
}

Given the contraction condition
\eqn{\|\rho\|_\infty < 1/\|W\|_\infty}, the Neumann series defining
\eqn{(\mathbb{I}-\mathcal{T})^{-1}} converges and the solution is obtained by
simple fixed-point iterations until the change is below \code{tol}.  Full
details are in Beyaztas, Shang and Sezer (2025).
}
\value{
A named \code{list} with components
\describe{
  \item{Y}{\code{n x length(gpy)} matrix of observed functional responses on the grid \code{gpy}.}
  \item{Y_true}{Same dimension as \code{Y}; noise-free latent responses before adding \eqn{\varepsilon_i(t)}.}
  \item{X}{\code{n x length(gpx)} matrix of functional predictors.}
  \item{W}{\code{n x n} row-normalised spatial weight matrix based on inverse distances.}
  \item{rho}{\code{length(gpy) x length(gpy)} matrix containing \eqn{\rho(t,u)} evaluated on the response grid.}
  \item{beta}{\code{length(gpx) x length(gpy)} matrix containing \eqn{\beta(t,s)} evaluated on the Cartesian product of the predictor and response grids.}
}
}
\references{
Beyaztas, U., Shang, H. L., and Sezer, G. B. (2025). Penalised Spatial Function-on-Function Regression. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, \bold{in press}.
}
\author{
Ufuk Beyaztas, Han Lin Shang, and Gizel Bakicierler Sezer
}
\examples{
\donttest{
# generate a toy data set
dat <- sff_dgp(n = 250, rf = 0.5)
}
}
