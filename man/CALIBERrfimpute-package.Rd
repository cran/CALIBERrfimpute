\name{CALIBERrfimpute-package}
\alias{CALIBERrfimpute-package}
\alias{CALIBERrfimpute}
\docType{package}

\title{
Imputation in MICE using Random Forest
}

\description{
Multivariate Imputation by Chained Equations (MICE) is commonly used to impute
missing values in analysis datasets using full conditional specifications. However,
it requires that the predictor models are specified correctly, including interactions
and nonlinearities. Random Forest is a regression and classification method which can
accommodate interactions and non-linearities without requiring a particular
statistical model to be specified.
}

\details{
CALIBERrfimpute provides functions that can be used as imputation functions with \code{mice}. Note that the mice package itself provides the \code{mice.impute.rf} function for imputation using Random Forest, as of version 2.20. The CALIBERrfimpute package provides different, independently developed imputation functions using Random Forest in MICE.

This package contains reports of two simulation studies:

\describe{
\item{Simulation study (linear regression):}{ A comparison of Random Forest and parametric MICE in a linear regression example. To open this report, run: \code{utils::browseURL(system.file("doc", "simstudy.pdf", package = "CALIBERrfimpute"))}}
\item{Survival analysis with interactions (large sample):}{ Compares the Random Forest MICE algorithm (\code{mice.impute.rfcont}) with parametric MICE and the algorithms of Doove et al. (\code{mice.impute.cart} and \code{mice.impute.rf}). To open this report, run:  \code{utils::browseURL(system.file("doc", "simstudy_survival_largesample.pdf", package = "CALIBERrfimpute"))}}
\item{Survival analysis with interactions (small sample):}{ A small sample version of the above simulation as a vignette. See \code{vignette("simstudy_survival", package = "CALIBERrfimpute")}}
}

\tabular{ll}{
Package: \tab CALIBERrfimpute\cr
Type: \tab Package\cr
Version: \tab 1.0-8\cr
Date: \tab 2026-02-17\cr
License: \tab GPL-3\cr
}

}
\author{
Anoop Shah

Maintainer: anoop@doctors.org.uk
}

\references{
Shah AD, Bartlett JW, Carpenter J, Nicholas O, Hemingway H.
Comparison of Random Forest and parametric imputation models for imputing missing
data using MICE: a CALIBER study. American Journal of Epidemiology 2014; 179(6): 764--774.
\doi{10.1093/aje/kwt312}

Doove LL, van Buuren S, Dusseldorp E.
Recursive partitioning for missing data imputation in the presence of interaction
effects. Computational Statistics and Data Analysis 2014; 72: 92--104.
\doi{10.1016/j.csda.2013.10.025}
}

\keyword{ package }

\seealso{
\code{\link[mice]{mice}},
\code{\link[randomForest]{randomForest}},
\code{\link{mice.impute.rfcont}},
\code{\link{mice.impute.rfcat}},
\code{\link[mice]{mice.impute.rf}}
}


