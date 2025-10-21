\name{fpqr_dgp}
\alias{fpqr_dgp}
\title{Generate a dataset for the function-on-function regression model}
\description{This function can be used to generate a dataset for the function-on-function regression model \deqn{
Y(u) = \int X(v) \beta(u,v) dv + \epsilon(u),} where \eqn{Y(u)} denotes the functional response, \eqn{X(v)} denotes the functional predictor, \eqn{\beta(u,v)} denotes the bivariate regression coefficient function, and \eqn{\epsilon(u)} is the error function.
}

\usage{
fpqr_dgp(n, nphi, gpy, gpx, err.dist, out.p)
}
\arguments{
  \item{n}{An integer, specifying the number of observations for each variable to be generated.}
  \item{nphi}{An integer, specifying the number of sub-random variables used in the generation of functional predictor. Default value is 10.}
  \item{gpy}{A vector, denoting the grid points of the functional response variable \eqn{Y(u)}.}
  \item{gpx}{A vector, denoting the grid points of the functional predictor variable \eqn{X(v)}.}
  \item{err.dist}{Distribution of the error term. Possibilities are "normal","t","chisq", "cn".}
  \item{out.p}{An integer between 0 and 1, denoting the outlier percentage in the generated data when the error distribution is contaminated normel, that is, "cn".}
}
\details{
When generating the data for the function-on-function regression model, first, the functional predictor is generated using the following process: \deqn{X_i(v) = \sum_{k=1}^{10} \frac{1}{k^{2}} \{ \zeta _{i1,k}\sqrt{2} \sin ( k \pi v ) \zeta _{i2,k} \sqrt{2} \cos (k \pi v) \}}
where \eqn{\zeta _{i1,k}} and \eqn{\zeta _{i2,k}} for \eqn{k=1,\ldots,10} are independent random variables generated from the standard normal distribution. The true intercept and bivariate regression parameter functions are generated as \eqn{\alpha ( u )=2e^\{- ( u-1 )^{2}\}} and \eqn{\beta ( v,u )= 4 \cos ( 2 \pi u ) \sin ( \pi v )}, respectively, Then, the elements of the functional response are generated as follows: \deqn{
Y_{i} ( u ) = \alpha ( u ) +  \int X_{i} ( v ) \beta( v, u) dv + \varepsilon _{i} ( u ),} where \eqn{\varepsilon _{i} ( u )} are generated from one of the following distributions:

When err.dist is "normal", then the error terms are generated from the standard normal distribution.

When err.dist is "t", then the error terms are generated from the t distribution with five degrees of freedom.

When err.dist is "chisq", then the error terms are generated from the chi-square distribution with one degree of freedom.

When err.dist is "cn", then the error terms are generated from the contaminated normal distribution.
}
\value{
\item{x}{A matrix containing the observations of simulated functional predictor variable.}
\item{y}{A matrix containing the observations of simulated functional response variable.}
\item{x.true}{A matrix containing the observations of simulated functional predictor variable without measurement error.}
\item{y.true}{A matrix containing the observations of simulated functional response variable without measurement error.}
\item{coef.a}{A vector containing the intercept function.}
\item{coef.b}{A matrix containing the bivariate regression coefficient vector.}
}
\references{
C. Xiong, X. Liugen, and C. Jiguo (2021), "Robust penalized M-estimation for function-on-function linear regression"", \emph{Stat}, \bold{10}(1), e390.
}
\author{Muge Mutis, Ufuk Beyaztas, Filiz Karaman, and Han Lin Shang
}
\examples{
gpx <- (1:50)/50
gpy <- (1:60)/60
data <- fpqr_dgp(n = 250, gpy = gpy, gpx = gpx, err.dist = "normal")
y <- data$y
x <- data$x
}
