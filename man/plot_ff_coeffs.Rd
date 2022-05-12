\name{plot_ff_coeffs}
\alias{plot_ff_coeffs}
\title{3D plot of bivariate regression coefficient functions of a function-on-function regression model}
\description{This function is used to obtain 3D plots of bivariate regression coefficient functions of a function-on-function regression model based on output object obtained from \code{\link{get.ff.coeffs}}.
}
\usage{
plot_ff_coeffs(object, b, phi, theta, cex.axis)
}
\arguments{
  \item{object}{The output object of \code{\link{get.ff.coeffs}}.}
  \item{b}{An integer value indicating which regression parameter function to be plotted.}
  \item{phi}{A numeric value defining the viewing direction.}
  \item{theta}{A numeric value defining the viewing direction.}
  \item{cex.axis}{A numeric value defining the size of the tick label.}
}
\value{
No return value, called for side effects.
}
\references{
K. Soetaert (2021) plot3D: Plotting Multi-Dimensional Data. R package version 1.4, 
URL: https://cran.r-project.org/web/packages/plot3D/index.html.
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
plot_ff_coeffs(object = coefs, b = 1, phi = 5, theta = 40, cex.axis = 0.75)
}