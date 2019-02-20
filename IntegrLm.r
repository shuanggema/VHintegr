library(Rcpp)
sourceCpp("integrlm.cpp")


standardize = function(x)
{
   p = ncol(x)
   n = nrow(x)
   center = colMeans(x)
   x.mean = x - matrix(rep(center, n), n, p, byrow=T)
   scale = sqrt(colSums(x.mean^2)/n)
   xx = t(t(x.mean)/scale)
   list(xx, center, scale)
}



IntegrLm_fix_lambda2 = function(x, y, S, nlambda1, lambda_ratio, lambda2, xi){
	m <- ncol(S)
	p <- nrow(S)
	# standardize X and add 1_n for intercepts.
	X = Y = vector("list", m)
	N = numeric(m)
	center <- matrix(0, nrow = p, ncol = m)
	scales <- matrix(0, nrow = p, ncol = m)
 
    ymeans <- numeric(m)
	for (i in 1:m) {
		n <- nrow(x[[i]])

		N[i] <- nrow(x[[i]])
		std <- standardize(x[[i]])
		X[[i]]      <- std[[1]]
		center[, i] <- std[[2]]
		scales[, i] <- std[[3]]

		ymeans[i] <- mean(y[[i]])
		Y[[i]] <- y[[i]] - ymeans[i]
	}


	fit <- integrlm_fix_lambda2(X, Y, S, nlambda1, lambda_ratio, lambda2, xi) 

	# scale back coefficients
	Beta <- array(0, dim = c(p+1, m, nlambda1))
	for (i in 1:m) {
		Beta[ 1, i, ] <- ymeans[i] - c( (center[, i]/scales[, i]) %*% fit$beta[, i, ] ) 
		Beta[-1, i, ] <- fit$beta[ , i, ]/scales[, i]
	}

	Beta[abs(Beta) < 10^-4 ] = 0

	bic <- -2*c(fit$loglik) + apply(fit$beta[-1,,]!=0, 3, "sum")*log(sum(N))
	output = list()
	output$Beta    = Beta 
	output$bic     = bic
	output$lambda1 = c(fit$lambda1)
	output
}	


IntegrLm_bic = function(x, y, S, nlambda1, lambda_ratio, lambda2_seq, xi){
	m <- ncol(S)
	p <- nrow(S)
	L <- length(lambda2_seq)
	BIC <- matrix(0, ncol = L, nrow = nlambda1)
	BETA <- array(0, dim = c(p+1, m, nlambda1, L))
	for (l in 1:L)	{
		lambda2 <- lambda2_seq[l]
		fit <- IntegrLm_fix_lambda2(x, y, S, nlambda1, lambda_ratio, lambda2, xi)
		BETA[ , , , l] = fit$Beta
		BIC[, l] = fit$bic
	}

	opt <- which(BIC == min(BIC), arr.ind = TRUE)

    output = list()
	output$BIC = BIC
	output$BETA = BETA
	output$opt = opt
	output$beta_hat = BETA[ , ,opt[1], opt[2]]
	output
}	



IntegrLm = function(x, y, S, nlambda1, lambda_ratio, lambda2_seq, xi){
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
		fit <- IntegrLm_fix_lambda2(xtr, ytr, S, nlambda1, lambda_ratio, lambda2, xi)
		BETA[ , , , l] = fit$Beta
		BIC[, l] = fit$bic
	}
	opt_bic <- which(BIC == min(BIC), arr.ind = TRUE)


	MSE <- matrix(0, ncol = L, nrow = nlambda1)
	for (l in 1:L)
	{
		for (k in 1:nlambda1)	{
			for (m in 1:M)	{
				MSE[k, l] = MSE[k, l] +  sum((yte[[m]] - 
											  BETA[1, m, k, l] - 
											   xte[[m]] %*% BETA[-1, m, k, l])^2)/length(yte[[m]])
			}
		}
	}


	opt_cv <- which(MSE == min(MSE), arr.ind = TRUE)
    output = list()
    output$MSE = MSE
	output$BIC = BIC
	output$BETA = BETA
	output$opt_cv = opt_cv
	output$opt_bic = opt_bic
	output$beta_cv = BETA[ , ,opt_cv[1], opt_cv[2]]
	output$beta_bic = BETA[ , ,opt_bic[1], opt_bic[2]]

	output
}	


