\name{roll.win}
\alias{roll.win}
\title{Rolling Window Analysis of a Time Series}
\usage{
roll.win(x, window = 3, indicator = c("mean", "var", "cov"), tau = NULL)
}
\arguments{
\item{x}{a vector, usually a time series.}

\item{window}{the width of the window to calculate the indicator.}

\item{indicator}{character string specifying type of indicator: expected value (\code{"mean"}), variance (\code{"var"}) or covariance (\code{"cov"}).}

\item{tau}{number of lags to calculate the covariance. When not specified using \code{"cov"}, the variance is calculated.}
}
\value{
a vector of the calculated indicators.
}
\description{
Helps to (visually) detect whether a time series is stationary or non-stationary. A time series is a data-generating process with every observation - as a random variable - following a distribution. When expectational value, variance, and covariance (between different points in time) are constant, the time series is indicated as weekly dependent and seen as stationary. This desired property is a requirement to overcome the problem of spurious regression. Since there is no distribution but only one observation for each point in time, adjacent observations will be used as stand-in to calculate the indicators. Therefore, the chosen window should not be too large.
}
\note{
Objects generated by \code{roll.win()} can be plotted using the regular \code{plot()} command.
}
\examples{
## Plot the expected values with a window of width 5
exp.values <- roll.win(1:100, window = 5, indicator = "mean")
plot(exp.values)

## Spurious regression example
set.seed(123)
N <- 10^3
p.values <- rep(NA, N)

for (i in 1:N) {
  x <- 1:100 + rnorm(100) # time series with trend
  y <- 1:100 + rnorm(100) # time series with trend
  p.values[i] <- summary(ols(y ~ x))$coef[2,4]
}
sum(p.values < 0.05)/N    # share of significant results (100\%)

for (i in 1:N) {
  x <- rnorm(100)         # time series without trend
  y <- 1:100 + rnorm(100) # time series with trend
  p.values[i] <- summary(ols(y ~ x))$coef[2,4]
}
sum(p.values < 0.05)/N    # share of significant results (~ 5\%)

}
\concept{non-stationarity}
\concept{rolling window}
\concept{stationarity}
\concept{weak dependence}
