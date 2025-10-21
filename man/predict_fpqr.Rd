\name{predict_fpqr}
\alias{predict_fpqr}
\title{Prediction for a function-on-function linear quantile regression model based on functional partial quantile regression}
\description{This function is used to make prediction for a new set of functional predictors based upon a fitted function-on-function linear quantile regression model in the output of \code{\link{fpqr}}.}
\usage{
predict_fpqr(object, xnew)
}
\arguments{
  \item{object}{An output object obtained from \code{\link{fpqr}}.}
  \item{xnew}{A matricx consisting of the new observations of functional predictor. The argument \code{xnew} must have the same length and the same structure as the input \code{x} of \code{\link{fpqr}}.}
}
\value{
A matrix of predicted values of the functional response variable for the given set of new functional predictor \code{xnew}.
}
\author{Muge Mutis, Ufuk Beyaztas, Filiz Karaman, and Han Lin Shang}
\examples{
\donttest{
  gpx <- (1:25)/25
  gpy <- (1:30)/30
  data <- fpqr_dgp(n = 10, gpy = gpy, gpx = gpx, err.dist = "normal")
  y <- data$y
  x <- data$x
  data.test <- fpqr_dgp(n = 10, gpy = gpy, gpx = gpx, err.dist = "normal")
  x.test <- data.test$x
  y.test <- data.test$y.true
  fpqr.model.li <- fpqr(y=y, x=x, tau = 0.5, gpx = gpx, gpy = gpy, qc.type = "li")
  predictions <- predict_fpqr(object = fpqr.model.li, xnew = x.test)
}
}
