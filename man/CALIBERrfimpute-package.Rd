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
statistical model to be specified. This package provides imputation functions for MICE
which use Random Forest instead of linear or logistic regression as the predictive model.

The \href{../doc/simstudy.pdf}{simulation study} shows the results obtained using different
Random Forest parameters.
}

\details{
\tabular{ll}{
Package: \tab CALIBERrfimpute\cr
Type: \tab Package\cr
Version: \tab 0.1-2\cr
Date: \tab 2013-10-30\cr
License: \tab GPL-3\cr
}


}
\author{
Anoop Shah

Maintainer: anoop@doctors.org.uk
}

\references{
[to be added]

}

\keyword{ package }

\seealso{
\code{\link{mice}}, \code{\link{randomForest}}
}


