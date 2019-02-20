library(Rcpp)
sourceCpp("integrlogit.cpp")




IntegrLogit_fix_lambda2 = function(x, y, S, nlambda1, lambda_ratio, lambda2, xi){
	m <- ncol(S)
	p <- nrow(S)
	# standardize X and add 1_n for intercepts.
	X = vector("list", m)
	N = numeric(m)
	scales <- matrix(0, nrow = p, ncol = m)
	for (i in 1:m) {
		n <- nrow(x[[i]])
		N[i] <- n
		scales[, i] <- sqrt(colSums(x[[i]]^2)/n)
		X[[i]] = cbind(1, t(t(x[[i]])/scales[, i]))
	}


	fit <- integrlogit_fix_lambda2(X, y, S, nlambda1, lambda_ratio, lambda2, xi) 

	# scale back coefficients
	Beta <- array(0, dim = c(p+1, m, nlambda1))
	for (i in 1:m) {
		Beta[ 1, i, ] <- fit$beta[ 1 , i, ]
		Beta[-1, i, ] <- fit$beta[-1 , i, ]/scales[, i]
	}

	Beta[abs(Beta) < 10^(-4) ] = 0

	bic <- -2*c(fit$loglik) + apply(fit$beta[-1,,]!=0, 3, "sum")*log(sum(N))
	output = list()
	output$Beta    = Beta 
	output$bic     = bic
	output$lambda1 = c(fit$lambda1)
	output
}	


IntegrLogit = function(x, y, S, nlambda1, lambda_ratio, lambda2_seq, xi){
	M <- length(x)
	xtr <- xte <- vector("list", M)
	ytr <- yte <- vector("list", M)
	for (m in 1:M)	{
		n <- nrow(x[[m]])
		idx <- sample(1:n, round(n/2), replace = FALSE)
		xtr[[m]] <- x[[m]][ idx, ]
		xte[[m]] <- x[[m]][-idx, ]

		ytr[[m]] <- y[[m]][ idx]
		yte[[m]] <- y[[m]][-idx]
	}

	m <- ncol(S)
	p <- nrow(S)
	L <- length(lambda2_seq)
	BIC <- matrix(0, ncol = L, nrow = nlambda1)
	BETA <- array(0, dim = c(p+1, m, nlambda1, L))
	for (l in 1:L)	{
		lambda2 <- lambda2_seq[l]
		fit <- IntegrLogit_fix_lambda2(xtr, ytr, S, nlambda1, lambda_ratio, lambda2, xi)
		BETA[ , , , l] = fit$Beta
		BIC[, l] = fit$bic
	}

	opt_bic <- which(BIC == min(BIC), arr.ind = TRUE)

	CVL <- matrix(0, ncol = L, nrow = nlambda1)
	for (l in 1:L)
	{
		for (k in 1:nlambda1)	{
			for (m in 1:M)	{
				explp <- exp(c(BETA[1, m, k, l] + xte[[m]] %*% BETA[-1, m, k, l]))
				pr <- explp/(1 + explp)
				CVL[k, l] = CVL[k, l] +  sum(yte[[m]] * log(pr) + (1 - yte[[m]])*log(1 - pr))
			}
		}
	}

	opt_cv <- which(CVL == max(CVL), arr.ind = TRUE)

    output = list()
	output$BIC = BIC
	output$CVL = CVL
	output$BETA = BETA
	output$opt_bic= opt_bic
	output$opt_cv= opt_cv
	output$beta_bic = BETA[ , ,opt_bic[1], opt_bic[2]]
	output$beta_cv = BETA[ , ,opt_cv[1], opt_cv[2]]

	output
}	

