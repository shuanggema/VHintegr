# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;



// abc thresholding rule 
double threshold(double a, double b , double c)
{  
  double val;
  if (fabs(b) <= c) val = 0;
  else {
        if (b > c) val = (b - c)/a;
        else val = (b + c)/a;
  }
  return val;
}



// [[Rcpp::export]]
double loglik_m(List X, List Y, mat beta) 
{
    double loglik = 0;
    int m = beta.n_cols;
    for (int k = 0; k < m; k++) {
        mat x = as<mat>(X[k]);
        vec y = as<vec>(Y[k]);
        int n = x.n_rows;
        vec res = y - x * beta.col(k);

        loglik += -n/2.0 * (1 + log(sum(res % res)/n));
    }
    return loglik;
}


// [[Rcpp::export]]
void cdLM(List& X, List& Y, mat& S, mat& beta, double& lambda1, double& lambda2, double& loglik, double& xi)
{
    int it = 0, maxit = 200;
    int m = beta.n_cols;
    int p = beta.n_rows;
    double Lold = -INFINITY;
    double Lnew = 0;
    double eps = 10;
    double a = 0, b = 0, der1 = 0;

    while ((it < maxit) & (eps > 1e-4)) {

        for (int k = 0; k < m; k++) {
            mat x = as<mat>(X[k]);
            vec y = as<vec>(Y[k]);

           
            // slopes
            for (int j = 0; j < p; j++) {
                vec xb = x * beta.col(k);
                der1 = sum((y - xb) % x.col(j));

                if (xi == 0) {
                    a = x.n_rows*2.0 + lambda2;
                    double b3 = 0;
                    for (int l = 0; l < m; l++) {
                       if((l != k) & (S(j, k)*S(j, l) > 0)) b3 += beta(j, l); 
                    }
                    b = x.n_rows*2.0 * beta(j, k) + 2.0 * der1 + lambda2*b3;
                } else {    
                    double a3 = 0;
                    for (int l = 0; l < m; l++) {
                      if(l != k ) a3 += 1/pow(abs(beta(j, k)) + xi, 2.0); 
                    }
                    a = x.n_rows*2.0 + lambda2*a3;
                    double b3 = 0;
                    for (int l = 0; l < m; l++) {
                       if(l != k ) b3 += beta(j, l)/(abs(beta(j, l)) + xi); 
                    }
                    b = x.n_rows*2.0 * beta(j, k) + 2.0 * der1 + lambda2*b3/(abs(beta(j, k)) + xi);
                }    
                
                beta(j, k) = threshold(a, b, lambda1);
            }
        }
        Lnew = loglik_m(X, Y, beta);

        eps = fabs(Lnew - Lold);
        Lold = Lnew;  
        it++;     
    }
    loglik = Lnew;

    /*List output = List::create(Rcpp::Named("loglik") = Lnew,      
                               Rcpp::Named("beta") = beta,
                               Rcpp::Named("eps") = eps,
                               Rcpp::Named("iter") = it);
    return output;*/

}

// [[Rcpp::export]]
double calLambdaLM(List X, List Y, int m, int p){
    vec lambda_max(m);
    for (int k=0; k < m; k++) {
        mat x = as<mat>(X[k]);
        vec y = as<vec>(Y[k]);
        vec der(p, fill::zeros); 
        // omit the first covariates
        for  (int j=0; j < p; j++) {
            der(j) = sum(y % x.col(j));
        }
        lambda_max(k) = 2*max(abs(der));
    }
    return max(lambda_max);
}


// [[Rcpp::export]]
List integrlm_fix_lambda2(List X, List Y, mat S, int nlambda1, 
                            double lambda_min_ratio, double lambda2, double xi) 
//  lasso solver, for a sequence of lambda1, fixed lambda2.
//  each x need contain 1_n, to include the intercept.
//  x columns need to be standardized, except the intercept column.
//  output:
//         loglik - log likelihood of all data
//         beta   -  a matrix, each column is a coeffcient vector for each data set.
{   
    int m = S.n_cols;
    // intercept
    int p = S.n_rows;
    int i = 0;
    vec loglik(nlambda1, fill::zeros);

    // lambda1 sequence
    double lambdamax = calLambdaLM(X, Y, m, p);
    vec lambda1 = exp(linspace(log(lambdamax), log(lambdamax*lambda_min_ratio), nlambda1));
    cube beta(p, m, nlambda1);
    beta.fill(0.0);

    for (i = 0; i < nlambda1; i++)  {
         // solve lasso with fixed lambda1, lambda2, with warmstart.
        if (i != 0) beta.slice(i) = beta.slice(i-1);
        cdLM(X, Y, S, beta.slice(i), lambda1(i), lambda2, loglik(i), xi);
    }

    List output = List::create(Rcpp::Named("loglik") = loglik,
                               Rcpp::Named("lambda1") = lambda1,
                               Rcpp::Named("beta") = beta);
    return output;
}
