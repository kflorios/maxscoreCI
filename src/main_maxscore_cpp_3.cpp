#include <Rcpp.h>

using namespace Rcpp;

//[[Rcpp::export]]
double boot_score_cpp(NumericVector y, NumericMatrix X, 
                       NumericVector bootY, NumericMatrix bootX, 
                       NumericVector beta, NumericMatrix H, NumericVector beta_hat, double x1_coef = 1.0) {
    
    int n = y.size();
    int k = beta.size();
    
    // Compute bootstrap empirical objective: M_n*(beta)
    double boot_obj = 0.0;
    for (int i = 0; i < n; i++) {
        double idx = x1_coef * bootX(i, 0);
        for (int j = 0; j < k; j++) {
            idx += bootX(i, j + 1) * beta[j];
        }
        if (idx >= 0.0) boot_obj += (2.0 * bootY[i] - 1.0);
    }
    boot_obj /= n;
    
    // Compute sample empirical objective: M_n(beta)
    double sample_obj = 0.0;
    for (int i = 0; i < n; i++) {
        double idx = x1_coef * X(i, 0);
        for (int j = 0; j < k; j++) {
            idx += X(i, j + 1) * beta[j];
        }
        if (idx >= 0.0) sample_obj += (2.0 * y[i] - 1.0);
    }
    sample_obj /= n;
    
    // Quadratic centering term: 0.5 * (beta - beta_hat)' * H * (beta - beta_hat)
    NumericVector db(k);
    for (int j = 0; j < k; j++) db[j] = beta[j] - beta_hat[j];
    
    double quad = 0.0;
    for (int r = 0; r < k; r++) {
        for (int c = 0; c < k; c++) {
            quad += db[r] * H(r, c) * db[c];
        }
    }
    quad *= 0.5;
    
    return boot_obj - sample_obj + quad;
}