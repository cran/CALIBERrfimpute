\name{setRFoptions}
\alias{setRFoptions}

\title{
Set Random Forest options for imputation using MICE
}

\description{
A convenience function to set global options for number of trees or number of nodes.
}

\usage{
setRFoptions(ntree_cat = NULL, ntree_cont = NULL,
    nodesize_cat = NULL, nodesize_cont = NULL,
    maxnodes_cat = NULL, maxnodes_cont = NULL)
}


\arguments{
  \item{ntree_cat}{
number of trees to be used for imputing categorical variables
(each imputed value is the prediction of a randomly chosen tree), default = 10.
}
  \item{ntree_cont}{
number of trees in the forest for imputing continuous variables, default = 10.
}
  \item{nodesize_cat}{
minimum node size for trees for imputing categorical variables, default = 1. 
A higher value can be used on larger datasets in order to save time.
}
  \item{nodesize_cont}{
minimum node size for trees for imputing continuous variables, default = 5.
A higher value can be used on larger datasets in order to save time.
}
  \item{maxnodes_cat}{
maximum number of nodes in trees for imputing categorical variables.
By default the size limit is set by the number of observations and nodesize_cat.
}
  \item{maxnodes_cont}{
maximum number of nodes in trees for imputing continuous variables.
By default the size limit is set by the number of observations and nodesize_cont.
}
}

\details{
This function sets the global options which have the prefix 'CALIBERrfimpute_'.
}

\value{
No return value. The function prints a message stating the new option setting.
}

\seealso{
\code{\link{mice.impute.rfcat}}, \code{\link{mice.impute.rfcont}}
}

\examples{
# Set option using setRFoptions
setRFoptions(ntree_cat=15)
options()$CALIBERrfimpute_ntree_cat

# Set option directly
options(CALIBERrfimpute_ntree_cat=20)
options()$CALIBERrfimpute_ntree_cat
}
