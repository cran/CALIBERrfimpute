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

The mice package provides the mice.impute.rf function for imputation using Random Forest, as of version 2.20. The CALIBERrfimpute package provides different, independently developed imputation functions using Random Forest in MICE.

This package contains reports of two simulation studies:

\href{../doc/simstudy.pdf}{Simulation study} is a comparison of Random Forest and parametric MICE in a linear regression example.

\href{../doc/simstudy_survival_large.pdf}{Vignette for survival analysis with interactions - large sample} compares the Random Forest MICE algorithm for continuous variables (\code{\link{mice.impute.rfcont}}) with parametric MICE and the algorithm of Doove et al. in the mice package (\code{\link[mice]{mice.impute.cart}} and \code{\link[mice]{mice.impute.rf}}).

\href{../doc/simstudy_survival.pdf}{Vignette for survival analysis with interactions - small sample} is a small sample example of the above, created during compilation of the package.
}

\details{
\tabular{ll}{
Package: \tab CALIBERrfimpute\cr
Type: \tab Package\cr
Version: \tab 1.0-7\cr
Date: \tab 2022-11-17\cr
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
doi:10.1093/aje/kwt312
\url{https://academic.oup.com/aje/article/179/6/764/107562}

Doove LL, van Buuren S, Dusseldorp E.
Recursive partitioning for missing data imputation in the presence of interaction
effects. Computational Statistics and Data Analysis 2014; 72: 92--104.
doi:10.1016/j.csda.2013.10.025
\url{https://www.sciencedirect.com/science/article/abs/pii/S0167947313003939}
}

\keyword{ package }

\seealso{
\code{\link[mice]{mice}},
\code{\link[randomForest]{randomForest}},
\code{\link{mice.impute.rfcont}},
\code{\link{mice.impute.rfcat}},
\code{\link[mice]{mice.impute.rf}}
}


