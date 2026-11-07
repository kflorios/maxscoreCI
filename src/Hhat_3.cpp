#include <Rcpp.h>

using namespace Rcpp;

double ker(double x){
    double PI = 3.1415926535897932;
    double anss = -x / sqrt(2.0 * PI) * exp(-pow(x, 2) / 2.0);
    return anss;
}

//[[Rcpp::export]]
double score_fnc_cpp(NumericVector y, NumericMatrix X, NumericVector beta, double x1_coef = 1.0) {
    int n = y.size();
    int k = beta.size(); 
    double obj = 0.0;
    
    for (int i = 0; i < n; i++) {
        // x1_coef * X(i, 0) + sum_{j=1}^{k} (X(i, j+1) * beta[j])
        double idx = x1_coef * X(i, 0);
        for (int j = 0; j < k; j++) {
            idx += X(i, j + 1) * beta[j];
        }
        
        if (idx >= 0.0) {
            obj += (2.0 * y[i] - 1.0);
        }
    }
    return obj / n;
}

//[[Rcpp::export]]
NumericMatrix hat_H_ms(NumericVector y, NumericMatrix X, double bw, NumericVector beta, double x1_coef = 1.0){
    int n = y.size();
    int k = beta.size(); 
    NumericMatrix H(k, k); 
    
    for (int i = 0; i < n; i++){
        double index = x1_coef * X(i, 0);
        for (int j = 0; j < k; j++) {
            index += X(i, j + 1) * beta[j];
        }
        index /= bw;
        
        double k_val = ker(index) / (bw * bw); 
        double weight = (2.0 * y[i] - 1.0) * k_val;
        
        for (int r = 0; r < k; r++) {
            for (int c = 0; c < k; c++) {
                H(r, c) += weight * X(i, r + 1) * X(i, c + 1);
            }
        }
    }
    
    for (int r = 0; r < k; r++) {
        for (int c = 0; c < k; c++) {
            H(r, c) /= n;
        }
    }
    return H;
}

//[[Rcpp::export]]
NumericMatrix hat_H_nd(NumericVector y, NumericMatrix X, NumericVector beta, double eps, double x1_coef = 1.0){
    int k = beta.size();
    NumericMatrix H(k, k);
    double f0 = score_fnc_cpp(y, X, beta, x1_coef);
    
    // Diagonal elements
    for (int i = 0; i < k; i++) {
        NumericVector b_plus = clone(beta);
        NumericVector b_minus = clone(beta);
        b_plus[i] += eps;
        b_minus[i] -= eps;
        
        double f_p = score_fnc_cpp(y, X, b_plus, x1_coef);
        double f_m = score_fnc_cpp(y, X, b_minus, x1_coef);
        
        H(i, i) = (f_p + f_m - 2.0 * f0) / (eps * eps);
    }
    
    // Off-diagonal elements (Cross-derivatives)
    for (int i = 0; i < k; i++) {
        for (int j = i + 1; j < k; j++) {
            NumericVector b_pp = clone(beta);
            NumericVector b_pm = clone(beta);
            NumericVector b_mp = clone(beta);
            NumericVector b_mm = clone(beta);
            
            b_pp[i] += eps; b_pp[j] += eps;
            b_pm[i] += eps; b_pm[j] -= eps;
            b_mp[i] -= eps; b_mp[j] += eps;
            b_mm[i] -= eps; b_mm[j] -= eps;
            
            double f_pp = score_fnc_cpp(y, X, b_pp, x1_coef);
            double f_pm = score_fnc_cpp(y, X, b_pm, x1_coef);
            double f_mp = score_fnc_cpp(y, X, b_mp, x1_coef);
            double f_mm = score_fnc_cpp(y, X, b_mm, x1_coef);
            
            double cross = (f_pp - f_pm - f_mp + f_mm) / (4.0 * eps * eps);
            H(i, j) = cross;
            H(j, i) = cross;
        }
    }
    
    return H;
}