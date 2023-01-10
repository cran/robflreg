\name{generate.ff.data}
\alias{generate.ff.data}
\title{Generate functional data for the function-on-function regression model}
\description{This function provides a unified simulation structure for the function-on-function regression model \deqn{
Y(t) = \sum_{m=1}^M \int X_m(s) \beta_m(s,t) ds + \epsilon(t),} where \eqn{Y(t)} denotes the functional response, \eqn{X_m(s)} denotes the \eqn{m}-th functional predictor, \eqn{\beta_m(s,t)} denotes the \eqn{m}-th bivariate regression coefficient function, and \eqn{\epsilon(t)} is the error function.
}
\usage{
generate.ff.data(n.pred, n.curve, n.gp, out.p = 0)
}
\arguments{
  \item{n.pred}{An integer, denoting the number of functional predictors to be generated.}
  \item{n.curve}{An integer, specifying the number of observations for each functional variable to be generated.}
  \item{n.gp}{An integer, denoting the number of grid points, i.e., a fine grid on the interval [0, 1].}
  \item{out.p}{An integer between 0 and 1, denoting the outlier percentage in the generated data.}
}
\details{
In the data generation process, first, the functional predictors are simulated based on the following process: \deqn{
X_m(s) = \sum_{j=1}^5 \kappa_j v_j(s),} where \eqn{ \kappa_j } is a vector generated from a Normal distribution with mean one and variance \eqn{\sqrt{a} j^{-1/2}}, \eqn{a} is a uniformly generated random number between 1 and 4, and \deqn{v_j(s) = \sin(j \pi s) - \cos(j \pi s).} The bivariate regression coefficient functions are generated from a coefficient space that includes ten different functions such as: \deqn{b \sin(2 \pi s) \sin(\pi t)} and \deqn{b e^{-3 (s - 0.5)^2} e^{-4 (t - 1)^2},} where \eqn{b} is generated from a uniform distribution between 1 and 3. The error function \eqn{\epsilon(t)}, on the other hand, is generated from the Ornstein-Uhlenbeck process: \deqn{\epsilon(t) = l + [\epsilon_0(t) - l] e^{-\theta t} + \sigma \int_0^t e^{-\theta (t-u)} d W_u,} where \eqn{l, \theta > 0, \sigma > 0} are constants, \eqn{\epsilon_0(t)} is the initial value of \eqn{\epsilon(t)} taken from \eqn{W_u}, and
\eqn{W_u} is the Wiener process. If outliers are allowed in the generated data, i.e., \eqn{out.p > 0}, then, the randomly selected \eqn{n.curve \times out.p} of the data are generated in a different way from the aforementioned process. In more detail, if \eqn{out.p > 0}, the bivariate regression coefficient functions (possibly different from the previously generated coefficient functions) generated from the coefficient space with \eqn{b^*} (instead of \eqn{b}), where \eqn{b^*} is generated from a uniform distribution between 1 and 2, are used to generate the outlying observations. In addition, in this case, the following process is used to generate functional predictors: \deqn{
X_m^*(s) = \sum_{j=1}^5 \kappa_j^* v_j^*(s),} where \eqn{ \kappa_j^* } is a vector generated from a Normal distribution with mean one and variance \eqn{\sqrt{a} j^{-3/2}} and \deqn{v_j^*(s) = 2 \sin(j \pi s) - \cos(j \pi s).} All the functions are generated equally spaced point in the interval \eqn{[0, 1]}.
}
\value{A list object with the following components:
\item{Y}{An \eqn{n.curve \times n.gp}-dimensional matrix containing the observations of simulated functional response variable.}
\item{X}{A list with length n.pred. The elements are the \eqn{n.curve \times n.gp}-dimensional matrices containing the observations of simulated functional predictor variables.}
\item{f.coef}{A list with length n.pred. Each element is a matrix and contains the generated bivariate regression coefficient function.}
\item{out.indx}{A vector with length \eqn{n.curve \times out.p} denoting the indices of outlying observations.}
}
\author{Ufuk Beyaztas and Han Lin Shang}
\references{
E. Garcia-Portugues and J. Alvarez-Liebana J and G. Alvarez-Perez G and W. Gonzalez-Manteiga W (2021)  "A goodness-of-fit test for the functional linear model with functional response",  \emph{Scandinavian Journal of Statistics}, \bold{48}(2), 502-528.
}
\examples{
library(fda)
library(fda.usc)
set.seed(2022)
sim.data <- generate.ff.data(n.pred = 5, n.curve = 200, n.gp = 101, out.p = 0.1)
Y <- sim.data$Y
X <- sim.data$X
coeffs <- sim.data$f.coef
out.indx <- sim.data$out.indx
fY <- fdata(Y, argvals = seq(0, 1, length.out = 101))
plot(fY[-out.indx,], lty = 1, ylab = "", xlab = "Grid point", 
     main = "Response", mgp = c(2, 0.5, 0), ylim = range(fY))
lines(fY[out.indx,], lty = 1, col = "black") # Outlying functions
}
