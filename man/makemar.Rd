\name{makemar}
\alias{makemar}

\title{
Creates artificial missing at random missingness
}

\description{
Introduces missingness into x1 and x2 into a data.frame of the format produced by \code{\link{simdata}},
for use in the \href{../doc/simstudy.pdf}{simulation study}.
The probability of missingness depends on the logistic of the fully observed variables y and x3;
hence it is missing at random but not missing completely at random.
}

\usage{
makemar(simdata, prop = 0.2)
}

\arguments{
  \item{simdata}{
simulated dataset created by \code{\link{simdata}}.
}
  \item{prop}{
proportion of missing values to be introduced in x1 and x2.
}
}

\details{
This function is used for simulation and testing.
}

\value{
A data.frame with columns:

\item{y}{dependent variable, based on the model y = x1 + x2 + x3 + normal error}
\item{x1}{partially observed continuous variable}
\item{x2}{partially observed continuous or binary (factor) variable}
\item{x3}{fully observed continuous variable}
\item{x4}{variable not in the model to predict y, but associated with x1, x2 and x3;
used as an auxiliary variable in imputation}
}


\seealso{
\code{\link{simdata}}
}

\examples{
set.seed(1)
mydata <- simdata(n=100)
mymardata <- makemar(mydata, prop=0.1)
# Count the number of missing values
sapply(mymardata, function(x){sum(is.na(x))})
# y x1 x2 x3 x4 
# 0 11 10  0  0 
}
