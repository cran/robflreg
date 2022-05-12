\name{generate.sf.data}
\alias{generate.sf.data}
\title{Generate functional data for the scalar-on-function regression model}
\description{This function is used to simulate data for the scalar-on-function regression model \deqn{
Y = \sum_{m=1}^M X_m(s) \beta_m(s) ds + \epsilon,} where \eqn{Y} denotes the scalar response, \eqn{X_m(s)} denotes the \eqn{m}-th functional predictor, \eqn{\beta_m(s)} denotes the \eqn{m}-th regression coefficient function, and \eqn{\epsilon} is the error process. 
}
\usage{
generate.sf.data(n, n.pred, n.gp)
}
\arguments{
  \item{n}{An integer, specifying the number of observations for each variable to be generated.}
  \item{n.pred}{An integer, denoting the number of functional predictors to be generated.}
  \item{n.gp}{An integer, denoting the number of grid points, i.e., a fine grid on the interval [0, 1].}
}
\details{
In the data generation process, first, the  functional predictors are simulated based on the following process: \deqn{
X_m(s) = \sum_{j=1}^5 \kappa_j v_j(s),} where \eqn{ \kappa_j } is a vector generated from a Normal distribution with mean one and variance \eqn{\sqrt{a} j^{-3/2}}, \eqn{a} is a uniformly generated random number between 1 and 4, and \deqn{v_j(s) = \sin(j \pi s) - \cos(j \pi s).} The  regression coefficient functions are generated from a coefficient space that includes ten different functions such as: \deqn{b \sin(2 \pi s)} and \deqn{b \cos(2 \pi s),} where \eqn{b} is generated from a uniform distribution between 1 and 3. The error process is generated from the standard normal distribution. All the functional predictors are generated equally spaced point in the interval \eqn{[0, 1]}.
}
\value{A list object with the following components:
\item{Y}{An \eqn{n \times 1}-dimensional matrix containing the observations of simulated scalar response variable.}
\item{X}{A list with length n.pred. The elements are the \eqn{n \times n.gp}-dimensional matrices containing the observations of simulated functional predictor variables.}
\item{f.coef}{A list with length n.pred. Each element is a vector and contains the generated regression coefficient function.}
}
\author{Ufuk Beyaztas and Han Lin Shang}
\examples{
library(fda.usc)
library(fda)
sim.data <- generate.sf.data(n = 400, n.pred = 5, n.gp = 101)
Y <- sim.data$Y
X <- sim.data$X
coeffs <- sim.data$f.coef
plot(Y, type = "p", pch = 16, xlab = "Index", ylab = "", main = "Response")
fX1 <- fdata(X[[1]], argvals = seq(0, 1, length.out = 101))
plot(fX1, lty = 1, ylab = "", xlab = "Grid point", 
     main = expression(X[1](s)), mgp = c(2, 0.5, 0))
}