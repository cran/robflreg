\name{plot_sf_coeffs}
\alias{plot_sf_coeffs}
\title{Plot of regression coefficient functions of a scalar-on-function regression model}
\description{This function is used to obtain the plots of regression coefficient functions of a scalar-on-function regression model based on output object obtained from \code{\link{get.sf.coeffs}}.
}
\usage{
plot_sf_coeffs(object, b)
}
\value{
No return value, called for side effects.
}
\arguments{
  \item{object}{The output object of \code{\link{get.sf.coeffs}}.}
  \item{b}{An integer value indicating which regression parameter function to be plotted.}
}
\author{Ufuk Beyaztas and Han Lin Shang}
\examples{
sim.data <- generate.sf.data(n = 400, n.pred = 5, n.gp = 101)
Y <- sim.data$Y
X <- sim.data$X
gp <- rep(list(seq(0, 1, length.out = 101)), 5) # grid points of Xs
model.fit <- rob.sf.reg(Y, X, emodel = "classical", gp = gp)
coefs <- get.sf.coeffs(model.fit)
plot_sf_coeffs(object = coefs, b = 1)
}

