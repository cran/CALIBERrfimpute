\name{mice.impute.rfcat}
\alias{mice.impute.rfcat}

\title{
Impute categorical variables using Random Forest within MICE
}

\description{
This method can be used to impute factor variables (binary or >2 levels)
in MICE by specifying
method = 'rfcat'. It was developed independently from the
\code{\link[mice]{mice.impute.rf}} algorithm of Doove et al.,
and differs from it in some respects.
}

\usage{
mice.impute.rfcat(y, ry, x, ntree_cat = NULL,
    nodesize_cat = NULL, maxnodes_cat = NULL, ntree = NULL, ...)
}

\arguments{
  \item{y}{
a factor vector of observed values and missing values of the variable to be imputed.
y must be a factor even if it has only 2 levels; it cannot be logical.
}
  \item{ry}{
a logical vector stating whether y is observed or not.
}
  \item{x}{
a matrix of predictors to impute y.
}
  \item{ntree_cat}{
number of trees, default = 10.

A global option can be set thus: \code{setRFoptions(ntree_cat=10)}.
}
  \item{nodesize_cat}{
minimum size of nodes, default = 1.

A global option can be set thus: \code{setRFoptions(nodesize_cat=1)}.
Smaller values of nodesize create finer, more precise trees but increase the computation time.
}
  \item{maxnodes_cat}{
maximum number of nodes, default NULL. If NULL the number of nodes is determined by number of observations and nodesize_cat.
}
  \item{ntree}{
an alternative argument for specifying the number of trees, over-ridden by \code{ntree_cat}. This is for consistency with the \code{mice.impute.rf} function.
}
  \item{\dots}{
other arguments to pass to randomForest.
}
}

\details{
This Random Forest imputation algorithm has been developed as an alternative to logistic or
polytomous regression, and can accommodate non-linear relations and interactions among the
predictor variables without requiring them to be specified in the model. The algorithm takes
a bootstrap sample of the data to simulate sampling variability, fits a set of classification
trees, and chooses each imputed value as the prediction of a randomly chosen tree.
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
This algorithm has been tested on simulated data and in survival analysis of real data
with artificially introduced missingness completely at random. There was slight bias in
hazard ratios compared to polytomous regression, but coverage of confidence intervals was
correct.
}

\seealso{
\code{\link{setRFoptions}},
\code{\link{mice.impute.rfcont}},
\code{\link[mice]{mice}},
\code{\link[mice]{mice.impute.rf}},
\code{\link[mice]{mice.impute.cart}},
\code{\link[randomForest]{randomForest}}
}

\examples{
set.seed(1)

# A small dataset with a single row to be imputed
mydata <- data.frame(x1 = as.factor(c('this', 'this', NA, 'that')),
	x2 = 1:4, x3 = as.factor(c('other', 'another', NA, 'another')))
mice(mydata, method = c('logreg', 'norm', 'logreg'), m = 2, maxit = 2)
mice(mydata[, 1:2], method = c('rfcat', 'rfcont'), m = 2, maxit = 2)
mice(mydata, method = c('rfcat', 'rfcont', 'rfcat'), m = 2, maxit = 2)

# A larger simulated dataset
mydata <- simdata(100, x2binary = TRUE)
mymardata <- makemar(mydata)

cat('\nNumber of missing values:\n')
print(sapply(mymardata, function(x){sum(is.na(x))}))

# Test imputation of a single column in a two-column dataset
cat('\nTest imputation of a simple dataset')
print(mice(mymardata[, c('y', 'x2')], method = 'rfcat'))

# Analyse data
cat('\nFull data analysis:\n')
print(summary(lm(y ~ x1 + x2 + x3, data = mydata)))

cat('\nMICE normal and logistic:\n')
print(summary(pool(with(mice(mymardata,
    method = c('', 'norm', 'logreg', '', '')), lm(y ~ x1 + x2 + x3)))))

# Set options for Random Forest
setRFoptions(ntree_cat = 10)

cat('\nMICE using Random Forest:\n')
print(summary(pool(with(mice(mymardata,
    method = c('', 'rfcont', 'rfcat', '', '')), lm(y ~ x1 + x2 + x3)))))

cat('\nDataset with unobserved levels of a factor\n')
data3 <- data.frame(x1 = 1:100, x2 = factor(c(rep('A', 25),
    rep('B', 25), rep('C', 25), rep('D', 25))))
data3$x2[data3$x2 == 'D'] <- NA
mice(data3, method = c('', 'rfcat'))
}
