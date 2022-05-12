\name{predict_sf_regression}
\alias{predict_sf_regression}
\title{Prediction for a scalar-on-function regression model}
\description{This function is used to make prediction for a new set of functional predictors based upon a fitted scalar-on-function regression model in the output of \code{\link{rob.sf.reg}}.
}
\usage{
predict_sf_regression(object, Xnew)
}
\arguments{
  \item{object}{An output object obtained from \code{\link{rob.sf.reg}}.}
  \item{Xnew}{A list of matrices consisting of the new observations of functional predictors. The argument \code{Xnew} must have the same length and the same structure as the input \code{X} of \code{\link{rob.sf.reg}}.}
}
\value{An \eqn{n_{test} \times 1}-dimensional matrix of predicted values of the scalar response variable for the given set of new functional predictors \code{Xnew}. Here, \eqn{n_{test}}, the number of rows of the matrix of predicted values, equals to the number of rows of \code{Xnew}.
}
\author{Ufuk Beyaztas and Han Lin Shang}
\examples{
set.seed(987)
sim.data <- generate.sf.data(n = 400, n.pred = 5, n.gp = 101)
Y <- sim.data$Y
X <- sim.data$X
Y.train <- Y[1:200,]
Y.test <- Y[201:400,]
X.train <- X.test <- list()
for(i in 1:5){
  X.train[[i]] <- X[[i]][1:200,]
  X.test[[i]] <- X[[i]][201:400,]
}
gp <- rep(list(seq(0, 1, length.out = 101)), 5) # grid points of Xs

set.seed(123)
sim.data2 <- generate.sf.data(n = 400, n.pred = 5, n.gp = 101)
out.index <- sample(1:200, 20)
Y.train[out.index] <- sim.data2$Y[out.index] + 10
for(i in 1:5)
  X.train[[i]][out.index,] <- sim.data2$X[[i]][out.index,]
model.tau <- rob.sf.reg(Y.train, X.train, emodel = "robust", fmodel = "tau", gp = gp)
pred.tau <- predict_sf_regression(object = model.tau, Xnew = X.test)
round(mean((Y.test - pred.tau)^2), 4)        # 1.9434 (tau method)
}