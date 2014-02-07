### R code from vignette source 'simstudy_survival.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: simstudy_survival.Rnw:13-36
###################################################
library(CALIBERrfimpute)
library(missForest)
library(survival)
library(xtable)
library(rpart)

kPmiss <- 0.2 # probability of missingness
kLogHR <- 0.5 # true log hazard ratio

# To analyse samples of more than 200 patients, (recommend about 2000,
# but this will slow down the program), set NPATS before running
# this vignette.
if (!exists('NPATS')){
	kSampleSize <- 200 # number of patients in simulated datasets
} else {
	kSampleSize <- NPATS
}
# e.g.
# NPATS <- 2000

# To analyse more than 3 samples, set N to a number greater than 3
# e.g.
# N <- 1000


###################################################
### code chunk number 2: simstudy_survival.Rnw:81-313
###################################################

#### DATA GENERATING FUNCTIONS ####

makeSurv <- function(n = 2000, loghr = kLogHR){
	# Creates a survival cohort of n patients. Assumes that censoring is
	# independent of all other variables
	
	# x1 and x2 are random normal variables
	data <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
	
	# Create the x3 variable
	data$x3 <- 0.5 * (data$x1 + data$x2 - data$x1 * data$x2) + rnorm(n)
	
	# Underlying log hazard ratio for all variables is the same
	data$y <- with(data, loghr * (x1 + x2 + x3))
	data$survtime <- rexp(n, exp(data$y))

	# Censoring - assume uniform distribution of observation times
	# up to a maximum
	obstime <- runif(nrow(data), min = 0,
		max = quantile(data$survtime, 0.5))
	data$event <- as.integer(data$survtime <= obstime)
	data$time <- pmin(data$survtime, obstime)
	data$cumhaz <- nelsonaalen(data, time, event)
	
	# True log hazard and survival time are not seen in the data
	# so remove them
	data$y <- NULL
	data$survtime <- NULL
	
	return(data)
}

makeMarSurv <- function(data, pmissing = kPmiss){
	# Introduces missing data dependent on event indicator
	# and cumulative hazard and x1 and x2
	
	logistic <- function(x){
		exp(x) / (1 + exp(x))
	}

	predictions <- function(lp, n){
		# uses the vector of linear predictions (lp) from a logistic model
		# and the expected number of positive responses (n) to generate
		# a set of predictions by modifying the baseline
	
		trialn <- function(lptrial){
			sum(logistic(lptrial))
		}
		stepsize <- 32
		lptrial <- lp
		while(abs(trialn(lptrial) - n) > 1){
			if (trialn(lptrial) > n){
				# trialn bigger than required 
				lptrial <- lptrial - stepsize
			} else {
				lptrial <- lptrial + stepsize
			}
			stepsize <- stepsize / 2
		}
		# Generate predictions from binomial distribution
		as.logical(rbinom(logical(length(lp)), 1, logistic(lptrial)))
	}

	data$x3[predictions(0.1 * data$x1 + 0.1 * data$x2 +
		0.1 * data$cumhaz + 0.1 * data$event, nrow(data) * pmissing)] <- NA
	return(data)	
}

#### IMPUTATION FUNCTIONS FROM DOOVE AND VAN BUUREN ####

mice.impute.rf <- function(y, ry, x, ntrees = 100,
	nodesize = 5, ...){
	# Use default mtry, i.e. one third the number of predictors for
	# categorical variables, square root of the number of predictors
	# for continuous dependent variables
	xobs <- x[ry,]
	xmis <- x[!ry,]
	yobs <- y[ry]
	# Function to create a single tree
	onetree <- function(xobs, xmis, yobs){
		fit <- randomForest(yobs ~ ., data = cbind(yobs, xobs),
			ntree = 1, replace = TRUE, type = regression,
			sampsize = length(yobs), nodesize = 5)
		leafnr <- predict(object = fit, newdata = xobs,
			nodes = T)
		nodes <- predict(object = fit, newdata = xmis,
			nodes = T)
		# Return a vector of observed y values that are
		# part of the same terminal leaf as the predicted
		# missing y value
		donor <- lapply(nodes, function(s){
			yobs[leafnr == s]
		})
		return(donor)
	}
	# Create a matrix of vectors of donors, from the number
	# of trees desired
	forest <- sapply(1:(ntrees = ntrees), FUN = function(s){
		onetree(xobs, xmis, yobs)
	})
	# For each missing value, randomly choose a donor value
	# from any of the possible donor values across all trees
	impute <- apply(forest, MARGIN = 1, FUN = function(s){
		sample(unlist(s), 1)
	})
	return(impute)
}

mice.impute.cart <- function(y, ry, x, minbucket = 5, cp = 1e-04,
	...){
	xobs <- x[ry,]
	xmis <- x[!ry,]
	yobs <- y[ry]
	if (is.factor(yobs)==F){
		fit <- rpart(yobs~., data = cbind(yobs,xobs), method = "anova",
			control = rpart.control(minbucket = minbucket, cp = cp), ...)
		leafnr  <- floor(as.numeric(row.names(fit$frame[fit$where,])))
		fit$frame$yval <- as.numeric(row.names(fit$frame))
		nodes <- predict(object = fit, newdata = xmis)
		donor <- lapply(nodes, function(s) yobs[leafnr == s])
		impute <- sapply(1:length(donor), function(s){
			sample(donor[[s]], 1)
		})
	} else {
		fit <- rpart(yobs~., data = cbind(yobs, xobs),
			method = "class", control = rpart.control(
			minbucket = minbucket, cp = cp), ...)
		nodes <- predict(object = fit, newdata = xmis)
		impute <- apply(nodes, MARGIN = 1, FUN = function(s){
			sample(colnames(nodes), size = 1, prob = s)
		})
	}
	return(impute)
}

#### OUR MICE RANDOM FOREST FUNCTIONS ####

mice.impute.rfcont5 <- function(y, ry, x, ...){
	mice.impute.rfcont(y = y, ry = ry, x = x, ntree_cont = 5)
}

mice.impute.rfcont10 <- function(y, ry, x, ...){
	mice.impute.rfcont(y = y, ry = ry, x = x, ntree_cont = 10)
}

mice.impute.rfcont20 <- function(y, ry, x, ...){
	mice.impute.rfcont(y = y, ry = ry, x = x, ntree_cont = 20)
}

mice.impute.rfcont50 <- function(y, ry, x, ...){
	mice.impute.rfcont(y = y, ry = ry, x = x, ntree_cont = 50)
}

mice.impute.rfcont100 <- function(y, ry, x, ...){
	mice.impute.rfcont(y = y, ry = ry, x = x, ntree_cont = 100)
}

#### FUNCTIONS TO DO THE ANALYSIS ####

coxfull <- function(data){
	# Full data analysis
	coefs <- summary(coxph(myformula, data = data))$coef
	# return a vector of coefficients (est), upper and lower 95% limits
	confint <- cbind(coefs[, 'coef'] - qnorm(0.975) * coefs[, 'se(coef)'],
		coefs[, 'coef'] + qnorm(0.975) * coefs[, 'se(coef)'])
	out <- cbind(coefs[, 'coef'], confint,
		kLogHR >= confint[,1] & kLogHR <= confint[,2])
	colnames(out) <- c('est', 'lo 95', 'hi 95', 'cover')
	out
}

coximpute <- function(imputed_datasets){
	# Analyses a list of imputed datasets
	docoxmodel <- function(data){
		coxph(myformula, data=data)
	}
	mirafits <- as.mira(lapply(imputed_datasets, docoxmodel))
	out <- summary(pool(mirafits))
	out <- cbind(out, kLogHR >= out[, 'lo 95'] & kLogHR <= out[, 'hi 95'])
	# Whether this confidence interval contains the true hazard ratio
	colnames(out)[length(colnames(out))] <- 'cover'
	out
}

domissf <- function(missdata, reps = NIMPS){
	# Imputation by missForest
	out <- list()
	for (i in 1:reps){
		invisible(capture.output(
			out[[i]] <- missForest(missdata)$ximp))
	}
	out
}

dorfimpute <- function(missdata, reps = NIMPS){
	# Cumulative hazard is the outcome, filling in missing data
	# using proximities, using a 300-tree Random Forest.
	out <- list()
	for (i in 1:reps){
		invisible(capture.output(
			out[[i]] <- rfImpute(cumhaz ~ ., missdata, iter = 10)))
	}
	out
}

domice <- function(missdata, functions, reps = NIMPS){
	mids <- mice(missdata, defaultMethod = functions,
		m = reps, visitSequence = 'monotone',
		printFlag = FALSE, maxit = 10)
	lapply(1:reps, function(x) complete(mids, x))
}

doanalysis <- function(x){
	# Creates a dataset, analyses it using different methods, and outputs
	# the result as a matrix of coefficients / SE and coverage 
	data <- makeSurv(kSampleSize)
	missdata <- makeMarSurv(data)
	out <- list()
	out$full <- coxfull(data)
	out$rfimpute <- coximpute(dorfimpute(missdata))
	out$missf <- coximpute(domissf(missdata))
	out$rf5 <- coximpute(domice(missdata, 'rfcont5'))
	out$rf10 <- coximpute(domice(missdata, 'rfcont10'))
	out$rf20 <- coximpute(domice(missdata, 'rfcont20'))
	out$rf50 <- coximpute(domice(missdata, 'rfcont50'))
	out$rf100 <- coximpute(domice(missdata, 'rfcont100'))
	out$rf <- coximpute(domice(missdata, 'rf'))
	out$cart <- coximpute(domice(missdata, 'cart'))
	out$mice <- coximpute(domice(missdata, 'norm'))
	out
}


###################################################
### code chunk number 3: simstudy_survival.Rnw:317-321
###################################################
mydata <- makeSurv(200)
plot(mydata[, c('x1', 'x2', 'x3')],
	main = "Associations between predictor variables in a sample dataset")
mydata <- makeSurv(20000)


###################################################
### code chunk number 4: simstudy_survival.Rnw:326-327
###################################################
summary(lm(x3 ~ x1*x2, data = mydata))


###################################################
### code chunk number 5: simstudy_survival.Rnw:330-341
###################################################
mydata <- makeSurv(2000)
mydata2 <- makeMarSurv(mydata)
# Plot non-missing data
plot(mydata$x1[!is.na(mydata2$x3)], mydata$x3[!is.na(mydata2$x3)],
	pch = 19, xlab = 'x1', ylab = 'x3')
# Plot missing data
points(mydata$x1[is.na(mydata2$x3)], mydata$x3[is.na(mydata2$x3)],
	col = 'red', pch = 19)
legend('bottomright', legend = c('x3 observed', 'x3 missing'),
	col = c('black', 'red'), pch = 19)
title('Association of predictor variables x1 and x3')


###################################################
### code chunk number 6: simstudy_survival.Rnw:346-370
###################################################
# Cox proportional hazards analysis
myformula <- as.formula(Surv(time, event) ~ x1 + x2 + x3)

# Analysis with 10,000 simulated patients (or more
# if the variable REFERENCE_SAMPLESIZE exists)
if (!exists('REFERENCE_SAMPLESIZE')){
	REFERENCE_SAMPLESIZE <- 10000
}

# Use parallel processing, if available, to create
# datasets more quickly.
if ('parallel' %in% loadedNamespaces() &
	!is.null(getOption('mc.cores')) &
	Sys.info()['sysname'] == 'Linux'){
	REFERENCE_SAMPLESIZE <- REFERENCE_SAMPLESIZE %/%
		getOption('mc.cores')
	simdata <- mclapply(1:getOption('mc.cores'),
		function(x) makeSurv(REFERENCE_SAMPLESIZE))
	simdata <- do.call('rbind', simdata)
} else {
	simdata <- makeSurv(REFERENCE_SAMPLESIZE)
}

summary(coxph(myformula, data = simdata))


###################################################
### code chunk number 7: simstudy_survival.Rnw:392-409
###################################################
# Setting analysis parameters: To analyse more than 3 samples,
# set N to the desired number before running this program
if (!exists('N')){
	N <- 3
}
# Number of imputations (set to at least 10 when
# running an actual simulation)
if (!exists('NIMPS')){
	NIMPS <- 3
}
# Use parallel processing if the 'parallel' package is loaded
if ('parallel' %in% loadedNamespaces() & Sys.info()['sysname'] == 'Linux'){
	cat('Using parallel processing\n')
	results <- mclapply(1:N, doanalysis)
} else {
	results <- lapply(1:N, doanalysis)
}


###################################################
### code chunk number 8: simstudy_survival.Rnw:439-475
###################################################
getParams <- function(coef, method){
	estimates <- sapply(results, function(x){
		x[[method]][coef, 'est']
	})
	bias <- mean(estimates) - kLogHR
	se_bias <- sd(estimates) / sqrt(length(estimates))
	z <- bias / se_bias
	ci_len <- mean(sapply(results, function(x){
		x[[method]][coef, 'hi 95'] - x[[method]][coef, 'lo 95']
	}))
	ci_cov <- mean(sapply(results, function(x){
		x[[method]][coef, 'cover']
	}))
	out <- c(bias, se_bias, z, sd(estimates), ci_len, ci_cov)
	names(out) <- c('bias', 'se_bias', 'z_bias', 'sd', 'ci_len', 'ci_cov')
	out
}

showTable <- function(coef){
	methods <- c('full', 'rfimpute', 'missf', 'cart', 'rf',
		'rf5', 'rf10', 'rf20', 'rf50', 'rf100', 'mice')
	methodnames <- c('Full data', 'rfImpute', 'missForest',
		'CART MICE', 'RF MICE (Doove)',
		paste('RFcont MICE, ', c(5, 10, 20, 50, 100), 'trees'),
		'Parametric MICE')
	out <- t(sapply(methods, function(x){
		getParams(coef, x)
	}))
	out <- formatC(out, digits = 3, format = 'fg')
	out <- rbind(c('', 'Standard', 'Z-score ', 'SD of', 'Mean 95%', '95% CI'),
		c('Bias', 'error of bias', 'for bias', 'estimate',
		'CI length', 'coverage'), out)
	out <- cbind(c('', '', methodnames), out)
	print(xtable(out), floating = FALSE, include.rownames = FALSE,
		include.colnames = FALSE, hline.after = c(0, 2, nrow(out)))
}


###################################################
### code chunk number 9: simstudy_survival.Rnw:490-491
###################################################
showTable('x1')


###################################################
### code chunk number 10: simstudy_survival.Rnw:500-501
###################################################
showTable('x2')


###################################################
### code chunk number 11: simstudy_survival.Rnw:511-512
###################################################
showTable('x3')


###################################################
### code chunk number 12: simstudy_survival.Rnw:522-542
###################################################
numtrees <- c(5, 10, 20, 50, 100)
bias <- sapply(numtrees, function(x){
	getParams('x3', paste('rf', x, sep=''))['bias']
})
se_bias <- sapply(numtrees, function(x){
	getParams('x3', paste('rf', x, sep=''))['se_bias']
})
lower_bias <- bias - 1.96*se_bias
upper_bias <- bias + 1.96*se_bias

# Blank plot
plot(-100, 0, type = 'p', pch = 15, cex = 1.3, ylab = 'Bias', 
	xlab = 'Number of trees', xlim = c(0,100), ylim = c(min(lower_bias), max(upper_bias)))
# Zero bias line
lines(c(0,100), c(0,0), lty = 2, color = 'gray')
# Confidence interval lines
for (i in 1:5){lines(rep(numtrees[i], 2), c(lower_bias[i], upper_bias[i]))}
# Points
points(numtrees, bias, pch = 15, cex = 1.3)
title('Bias in estimate of x3 coefficient after\nmultiple imputation using RFcont MICE')


###################################################
### code chunk number 13: simstudy_survival.Rnw:547-694
###################################################
# Comparing confidence interval coverage and bias between:
#    RF MICE 100 trees
#    RF MICE 10 trees
#    Parametric MICE

# Names of the variables in the comparison
variables <- c('x1', 'x2', 'x3')

pstar <- function(x){
	if (x < 0.001){
		'***'
	} else if (x < 0.01){
		'**'
	} else if (x < 0.05){
		'*'
	} else {
		''
	}
}

compareBias <- function(method1, method2){
	# Generates a table comparing bias
	# Comparison statistic is the difference in absolute bias
	# (negative means first method is better)
	
	compareBiasVar <- function(varname){
		# All coefficients should be kLogHR
		bias1 <- sapply(results, function(x){
			x[[method1]][varname, 'est']
		}) - kLogHR
		bias2 <- sapply(results, function(x){
			x[[method2]][varname, 'est']
		}) - kLogHR

		if (sign(mean(bias1)) == -1){
			bias1 <- -bias1
		}
		if (sign(mean(bias2)) == -1){
			bias2 <- -bias2
		}
		
		paste(formatC(mean(bias1) - mean(bias2), format = 'fg', digits = 3),
			pstar(t.test(bias1 - bias2)$p.value))
	}
	
	sapply(variables, compareBiasVar)
}

compareVariance <- function(method1, method2){
	# Generates a table comparing precision between two methods
	# Comparison statistic is ratio of variance
	# (smaller means first method is better)
	
	compareVarianceVar <- function(varname){
		e1 <- sapply(results, function(x){
			x[[method1]][varname, 'est']
		})
		e2 <- sapply(results, function(x){
			x[[method2]][varname, 'est']
		})
		paste(formatC(var(e1) / var(e2), format = 'fg', digits = 3),
			pstar(var.test(e1, e2)$p.value))
	}
	
	sapply(variables, compareVarianceVar)
}

compareCIlength <- function(method1, method2){
	# Generates a table comparing coverage percentage between two methods
	# Comparison statistic is the ratio of confidence interval lengths
	# (less than 1 = first better)
	
	compareCIlengthVar <- function(varname){
		# Paired t test for bias (difference in estimate)
		len1 <- sapply(results, function(x){
			x[[method1]][varname, 'hi 95'] -
				x[[method1]][varname, 'lo 95']
		})
		len2 <- sapply(results, function(x){
			x[[method2]][varname, 'hi 95'] -
				x[[method2]][varname, 'lo 95']
		})
		
		paste(formatC(mean(len1) / mean(len2),
			format = 'fg', digits = 4),
			pstar(t.test(len1 - len2)$p.value))
	}
	
	sapply(variables, compareCIlengthVar)	
}

compareCoverage <- function(method1, method2){
	# Generates a table comparing coverage percentage between two methods
	# Comparison statistic is the difference in coverage
	# (positive = first better)

	compareCoverageVar <- function(varname){
		# Paired t test for bias (difference in estimate)
		
		cov1 <- sapply(results, function(x){
			x[[method1]][varname, 'cover']
		})
		cov2 <- sapply(results, function(x){
			x[[method2]][varname, 'cover']
		})
				
		paste(formatC(100 * (mean(cov1) - mean(cov2)), format = 'f',
			digits = 1),
			pstar(binom.test(c(sum(cov1 == TRUE  & cov2 == FALSE),
			sum(cov1 == FALSE & cov2 == TRUE)))$p.value))
	}
	
	sapply(variables, compareCoverageVar)	
}

maketable <- function(comparison){
	# comparison is a function such as compareCoverage, compareBias
	compare <- cbind(comparison('rf10', 'mice'),
		comparison('rf100', 'mice'),
		comparison('rf10', 'rf100'))
	compare <- cbind(rownames(compare), compare)
	compare <- rbind(
		c('', 'RFcont MICE 10 vs', 'RFcont MICE 100 vs',
			'RFcont MICE 10 vs'),
		c('Coefficient', 'parametric MICE',
			'parametric MICE', 'RFcont MICE 100'),
		compare)
	print(xtable(compare), include.rownames = FALSE,
		include.colnames = FALSE, floating = FALSE,
		hline.after = c(0, 2, nrow(compare)))
	
	cat('\n\\vspace{1em}\n')
	
	compare <- cbind(comparison('rf10', 'rf'),
		comparison('rf10', 'cart'),
		comparison('rf', 'cart'))
	compare <- cbind(rownames(compare), compare)
	compare <- rbind(
		c('', 'RFcont MICE 10 vs', 'RFcont MICE 10 vs',
			'RF MICE (Doove) vs'),
		c('Coefficient', 'RF MICE (Doove)',
			'CART MICE', 'CART MICE'),
		compare)
	print(xtable(compare), include.rownames = FALSE,
		include.colnames = FALSE, floating = FALSE,
		hline.after = c(0, 2, nrow(compare)))
}


###################################################
### code chunk number 14: simstudy_survival.Rnw:703-704
###################################################
maketable(compareBias)


###################################################
### code chunk number 15: simstudy_survival.Rnw:713-714
###################################################
maketable(compareVariance)


###################################################
### code chunk number 16: simstudy_survival.Rnw:724-725
###################################################
maketable(compareCIlength)


###################################################
### code chunk number 17: simstudy_survival.Rnw:734-735
###################################################
maketable(compareCoverage)


###################################################
### code chunk number 18: simstudy_survival.Rnw:773-782
###################################################
showfunction <- function(functionname){
	cat(paste(functionname, '<-',
		paste(capture.output(print(get(functionname))),
		collapse = '\n')))
	cat('\n')
	invisible(NULL)
}
showfunction('makeSurv')
showfunction('makeMarSurv')


###################################################
### code chunk number 19: simstudy_survival.Rnw:787-800
###################################################
showfunction('coxfull')
showfunction('coximpute')
showfunction('domissf')
showfunction('dorfimpute')
showfunction('mice.impute.cart')
showfunction('mice.impute.rf')
showfunction('mice.impute.rfcont5')
showfunction('mice.impute.rfcont10')
showfunction('mice.impute.rfcont20')
showfunction('mice.impute.rfcont50')
showfunction('mice.impute.rfcont100')
showfunction('domice')
showfunction('doanalysis')


###################################################
### code chunk number 20: simstudy_survival.Rnw:805-810
###################################################
showfunction('pstar')
showfunction('compareBias')
showfunction('compareVariance')
showfunction('compareCIlength')
showfunction('compareCoverage')


###################################################
### code chunk number 21: simstudy_survival.Rnw:815-818
###################################################
showfunction('getParams')
showfunction('showTable')
showfunction('maketable')


