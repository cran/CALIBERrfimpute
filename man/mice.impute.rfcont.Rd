\name{mice.impute.rfcont}
\alias{mice.impute.rfcont}

\title{
Impute continuous variables using Random Forest within MICE
}

\description{
This method can be used to impute continuous variables in MICE by specifying
method = 'rfcont'. It was developed independently from the
\code{\link[mice]{mice.impute.rf}} algorithm of Doove et al.,
and differs from it in drawing imputed values from a normal distribution.
}

\usage{
mice.impute.rfcont(y, ry, x, ntree_cont = NULL,
    nodesize_cont = NULL, maxnodes_cont = NULL, ntree = NULL, ...)
}

\arguments{
  \item{y}{
a vector of observed values and missing values of the variable to be imputed.
}
  \item{ry}{
a logical vector stating whether y is observed or not.
}
  \item{x}{
a matrix of predictors to impute y.
}
  \item{ntree_cont}{
number of trees, default = 10. A global option can be set by calling e.g. \code{setRFoptions(ntree_cont=10)}.
}
  \item{nodesize_cont}{
minimum size of nodes, default = 5.
A global option can be set by calling setRFoptions(nodesize_cont=XX).
Smaller values of nodesize create finer, more precise trees but makes the program slower.
}
  \item{maxnodes_cont}{
maximum number of nodes, default NULL. If NULL the number of nodes is determined by number of observations and nodesize_cont.
}
  \item{ntree}{
an alternative argument for specifying the number of trees, over-ridden by \code{ntree_cont}. This is for consistency with the \code{mice.impute.rf} function.
}
  \item{\dots}{
other arguments to pass to randomForest.
}
}

\details{
This Random Forest imputation algorithm has been developed as an alternative to normal-based
linear regression, and can accommodate non-linear relations and interactions among the
predictor variables without requiring them to be specified in the model. The algorithm takes
a bootstrap sample of the data to simulate sampling variability, fits a regression forest
trees and calculates the out-of-bag mean squared error. Each value is imputed as a random draw
from a normal distribution with mean defined by the Random Forest prediction and variance equal
to the out-of-bag mean squared error.

If only one tree is used (not recommended), a bootstrap sample is not taken in the first stage
because the Random Forest algorithm performs an internal bootstrap sample before fitting the tree.
}

\value{
A vector of imputed values of y.
}

\references{
Shah AD, Bartlett JW, Carpenter J, Nicholas O, Hemingway H. Comparison of Random Forest and parametric imputation models for imputing missing data using MICE: a CALIBER study. American Journal of Epidemiology 2014. doi: \href{http://dx.doi.org/10.1093/aje/kwt312}{10.1093/aje/kwt312}
}

\author{
Anoop Shah
}

\note{
This algorithm has been tested on simulated data with linear regression,
and in survival analysis of real data with artificially introduced missingness at random.
On the simulated data there was slight bias if the distribution of missing values was
very different from observed values, because imputed values were closer to the centre of
the data than the missing values. However in the survival analysis the hazard ratios
were unbiased and coverage of confidence intervals more conservative than normal-based MICE,
but the mean length of confidence intervals was shorter with mice.impute.rfcont.
}

\seealso{
\code{\link{setRFoptions}}, \code{\link{mice.impute.rfcat}},
\code{\link[mice]{mice}},
\code{\link[mice]{mice.impute.rf}},
\code{\link[mice]{mice.impute.cart}},
\code{\link[randomForest]{randomForest}}
}


\examples{
set.seed(1)
mydata <- simdata(100)

# Describe data
cat('\nSimulated multivariate normal data:\n')
print(data.frame(mean = colMeans(mydata), sd = sapply(mydata, sd)))

# Apply missingness pattern
mymardata <- makemar(mydata)
cat('\nNumber of missing values:\n')
print(sapply(mymardata, function(x){sum(is.na(x))}))

# Analyse data
cat('\nFull data analysis:\n')
print(summary(lm(y ~ x1 + x2 + x3, data=mydata)))

cat('\nMICE using normal-based linear regression:\n')
print(summary(pool(with(mice(mymardata,
    method = 'norm'), lm(y ~ x1 + x2 + x3)))))

# Set options for Random Forest
setRFoptions(ntree_cont = 10)

cat('\nMICE using Random Forest:\n')
print(summary(pool(with(mice(mymardata,
    method = 'rfcont'), lm(y ~ x1 + x2 + x3)))))
}

