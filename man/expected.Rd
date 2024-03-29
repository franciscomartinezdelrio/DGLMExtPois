% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal-probabilities-and-expectations.R
\name{expected}
\alias{expected}
\alias{hP_expected}
\alias{CMP_expected}
\title{Expected Probabilities and Frequencies for the hyper-Poisson and COM-Poisson
Model}
\usage{
hP_expected(object)

CMP_expected(object)
}
\arguments{
\item{object}{a fitted object of class inheriting from \code{"glm_hP"} or
\code{"glm_CMP"}.}
}
\value{
A list containing the following components:

 \item{\code{frequencies}}{the expected counts for the hP or CMP fit.}
 \item{\code{observed_freq}}{the observed distribution.}
 \item{\code{probabilities}}{the expected distribution for the hP or CMP
 fit.} \item{\code{dif}}{sum of the absolute value of differences between
 \code{frequencies} and \code{observed_freq}.} \item{\code{chi2}}{sum of the
 square of differences between \code{frequencies} and \code{observed_freq}.}
}
\description{
The \code{hP_expected} and \code{CMP_expected} functions calculate the
probability distribution of the count response variable Y for each observation
and obtain the corresponding expected frequencies. It is an informal way of
assessing the fit of the hP or CMP model by comparing the predicted
distribution of counts with the observed distribution.
}
\details{
The average expected probabilities are computed as \deqn{\bar(Pr)(y=k) =
\frac{1}{n} \sum_{i=1}^n \widehat{Pr}(y_i = k | x_i)}

The expected frequencies are obtained by multiplying by n.

Two measures are offered for summarizing the comparison between expected and
observed frequencies: the sum of the absolute value of differences and the sum
of the square of differences (similar to the Pearson statistic of goodness of
fit).
}
\examples{
## Fit a hyper-Poisson model

Bids$size.sq <- Bids$size ^ 2
hP.fit <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
                 whtknght + bidprem + insthold + size + size.sq + regulatn,
                 formula.gamma = numbids ~ 1, data = Bids)

## Compute the expected probabilities and the frequencies

hP_expected(hP.fit)
## Estimate a COM-Poisson model

Bids$size.sq <- Bids$size ^ 2
CMP.fit <- glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest +
                   whtknght + bidprem + insthold + size + size.sq + regulatn,
                   formula.nu = numbids ~ 1, data = Bids)

## Compute the expected probabilities and the frequencies

CMP_expected(CMP.fit)
}
\references{
J. M. Hilbe (2011). Negative Binomial Regression. (2nd ed.). Cambridge
University Press.

M. Scott Long and Jeremy Freese (2014). Regression Models for Categorical
Dependent Variables using STATA. (3rd ed.). Stata Press.
}
