\name{plot_ff_coeffs}
\alias{plot_ff_coeffs}
\title{Image plot of bivariate regression coefficient functions of a function-on-function regression model}
\description{This function is used to obtain image plots of bivariate regression coefficient functions of a function-on-function regression model based on output object obtained from \code{\link{get.ff.coeffs}}.
}
\usage{
plot_ff_coeffs(object, b)
}
\arguments{
  \item{object}{The output object of \code{\link{get.ff.coeffs}}.}
  \item{b}{An integer value indicating which regression parameter function to be plotted.}
}
\value{
No return value, called for side effects.
}
\references{
D. Nychka and R. Furrer and J. Paige and S. Sain (2021) fields: Tools for spatial data. R package version 14.1,
URL: https://github.com/dnychka/fieldsRPackage.
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
plot_ff_coeffs(object = coefs, b = 1)
}