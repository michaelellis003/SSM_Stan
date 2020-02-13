/*
Implement Kalman filter for multivariate Dynamic Linear Model with
known and constant covariance matices for the observation equation 
and the state transition equation. - See Chapter 2.7 of Dynamic Linear Models 
with R - Petris et. al.

Uses singular value decomposition to avoid numerical in instability when
updating filtering distribution covariance matrix - Appendix B of Dynamic Linear 
Models  with R - Petris et. al.
*/
functions {
#include SVD.stan
}
data {
    int<lower=1> N;         // Number of observations
    int<lower=1> M;         // Dimension of observation
    int<lower=1> P;         // Dimension of state vector
    
    // Observation Stuff
    vector[M] y[N];                             // Observations
    vector<lower=0, upper=1>[M] y_missing[N];   // Locations of missing observations. 1 if missing, 0 otherwise
    int<lower=0, upper=M> num_missing[N];       // Numer of missing observations in row
    matrix[M, P] FF;                            // Observation equation matrix
    cov_matrix[M] V;                            // Known constant covariance matix for observation equation
    
    // State Stuff
    matrix[P, P] GG;        // State/transition equation matrix
    cov_matrix[P] W;        // Known constant covariance matix for state transition equation
    
    vector[P] m0;           // Initial state vector prior mean
    matrix[P, P] C0;        // Initial state vector prior variance-covariance matix
}
parameters {
}
transformed parameters {
}
model {
}
generated quantities {
    // mean and variance-covariance matix for one-step-ahead (Gaussian) predictive
    // distribution of state t given all observations through time t-1
    vector[P] a[N];
    cov_matrix[P] R[N];
    matrix[P, P] U_R[N];
    matrix[P, P] D_R[N];
    
    // mean and variance-covariance matix for one-step-ahead (Gaussian) predictive
    // distribution of observation t given all observations through time t-1
    vector[M] f[N];
    cov_matrix[M] Q[N];
    
    // mean and variance-covariance matix for filtering distribution of
    // state t given all observations through time t
    vector[P] m[N+1];
    cov_matrix[P] C[N+1];
    matrix[P, P] U_C[N+1];
    matrix[P, P] D_C[N+1];
    
    // Store square root of W and V^-1
    matrix[P, P] N_W ;
    matrix[M, M] V_inv;
    matrix[M, M] N_V;
    
    // Get square root of W using SVD
    N_W = sqrt_svd(W);
    
    // Get square root of V^-1
    V_inv = inverse(V);
    N_V = sqrt_svd(V_inv);
    
    // Kalman filter
    // Intialize
    m[1] = m0;
    C[1] = C0;
    D_C[1] = D_svd(C0, 0);
    U_C[1] = V_svd(C0)';

    for(n in 1:N) {
        // local variables
        matrix[2*P, P] sqrt_D_R;
        matrix[2*P, P] M_R;      // M_R[n]' * M_R[n] = R[n]
        matrix[P+M, P] M_C;      // U_R[N] * M_C[n]' * M_C[n] * U_R[N]' = C[n]^-1
        matrix[P, P] V_M_C;
        matrix[P+M, P] Dinv_M_C;
        matrix[M, M] Q_inv;
        int num_miss = num_missing[n];
        
        a[n] = GG * m[n];
        M_R = make_MR(D_C[n], U_C[n], GG, N_W);
        R[n] = M_R' * M_R;
        U_R[n] = V_svd(M_R);
        sqrt_D_R = D_svd(M_R, 0);
        D_R[n] = sqrt_D_R' * sqrt_D_R;

        f[n] = FF * a[n];
        Q[n] = FF * R[n] * FF' + V;
        
        if(num_miss == M) {
            // All observations in a row are missing
            m[n+1] = a[n];
            U_C[n+1] = U_R[n];
            D_C[n+1] = D_R[n];
            C[n+1] = U_C[n+1] * D_C[n+1] * U_C[n+1]';
            
        } else if(num_miss > 0 && num_miss < M) {
            // Some but not all observations in a row are missing
            int M_tilde = M - num_missing[n];  // let M_tilde be number of non-missing observations
            matrix[M_tilde, M] Mt = rep_matrix(0, M_tilde, M);
            vector[M_tilde] f_tilde;
            matrix[M_tilde, M_tilde] Q_tilde;
            matrix[M_tilde, M_tilde] Qinv_tilde;
            matrix[P+M_tilde, P] M_C_tilde;
            matrix[P+M_tilde, P] Dinv_M_C_tilde;
            matrix[M_tilde, M_tilde] V_tilde;
            matrix[M_tilde, M_tilde] Vinv_tilde;
            matrix[M_tilde, M_tilde] N_V_tilde;
            matrix[M_tilde, P] FF_tilde;
            vector[M_tilde] y_tilde;
            int i_tilde = 1;
            
            for(i in 1:M){
                if(y_missing[n][i] == 0){
                   Mt[i_tilde, i] = 1;
                   y_tilde[i_tilde] = y[n][i];
                   i_tilde = i_tilde+1;
                }
            }
            
            FF_tilde = Mt * FF;
            V_tilde = Mt * V * Mt';
            
            f_tilde = FF_tilde * a[n];
            Q_tilde = FF_tilde * R[n] * FF_tilde' + V_tilde;
            Qinv_tilde = inverse(Q_tilde);
            
            m[n+1] = a[n] + R[n] * FF_tilde' * Qinv_tilde * (y_tilde - f_tilde);
            
            // get inverse of V_tilde
            Vinv_tilde = inverse(V_tilde);
            N_V_tilde = sqrt_svd(Vinv_tilde);
            
            M_C_tilde = make_MC(N_V_tilde, FF_tilde, U_R[n], D_R[n]);
            V_M_C = V_svd(M_C_tilde);
            Dinv_M_C_tilde = D_svd(M_C_tilde, 1);
        
            U_C[n+1] = U_R[n] * V_M_C;
            D_C[n+1] = Dinv_M_C_tilde' * Dinv_M_C_tilde;
            C[n+1] = U_C[n+1] * D_C[n+1] * U_C[n+1]';
            
        } else {
            // No observations in a row are missing
            Q_inv = inverse(Q[n]);
            m[n+1] = a[n] + R[n] * FF' * Q_inv * (y[n] - f[n]);
        
            M_C = make_MC(N_V, FF, U_R[n], D_R[n]);
            V_M_C = V_svd(M_C);
            Dinv_M_C = D_svd(M_C, 1);
        
            U_C[n+1] = U_R[n] * V_M_C;
            D_C[n+1] = Dinv_M_C' * Dinv_M_C;
            C[n+1] = U_C[n+1] * D_C[n+1] * U_C[n+1]';
        }
    }
}
