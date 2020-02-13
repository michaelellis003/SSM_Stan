/*
Estimate the unknown observation and state covariance matrices for a DLM.

For a DLM we can write the joint density of the observations as

        p(y_1:T | ...) = \prod_{t=1}^T p(y_t | y_{1:t-1}, ...)

The kalman filter gives us p(y_t | y_{1:t-1}, ...) at each time t. So we can use
the output of kalman filter to define the loglikelihood needed by Stan. That is,

        likelihood = \prod_{t=1}^T p(y_t | y_{1:t-1}, ...)
                   = \prod_{t=1}^T N(y_t | ft, Qt)
                
where, ft and Qt are given by the kalman filter.
See, Chapter 4 - Dynamic Linear Models with R (Petris et. al) and
Chapter 7 - Time Series Analysis by State Space Methods (Durbin and Koopman)
*/
functions {
#include SVD.stan
}
data {
    int<lower=1> N;         // Number of observations
    int<lower=1> M;         // Dimension of observation
    int<lower=1> P;         // Dimension of state vector
    
    vector[M] y[N];         // Observations
    vector<lower=0, upper=1>[M] y_missing[N];   // Locations of missing observations. 1 if missing, 0 otherwise
    int<lower=0, upper=M> num_missing[N]; 
    matrix[M, P] FF;        // Observation equation matrix
    matrix[P, P] GG;        // State/transition equation matrix
    
    vector[P] m0;           // Initial state vector prior mean
    matrix[P, P] C0;        // Initial state vector prior variance-covariance matix
}
parameters {
    cholesky_factor_corr[M] V_L_Omega;
    vector<lower=0,upper=pi()/2>[M] V_tau_unif;
    
    cholesky_factor_corr[P] W_L_Omega;
    vector<lower=0,upper=pi()/2>[P] W_tau_unif;
}
transformed parameters {
    vector<lower=0>[M] V_tau;
    vector<lower=0>[P] W_tau;
    cov_matrix[M] V;
    cov_matrix[P] W;
    
    for(i in 1:M){
        V_tau[i] = 2.5 * tan(V_tau_unif[i]);
    }
    for(p in 1:P){
        W_tau[p] = 2.5 * tan(W_tau_unif[p]);
    }
    
    V = diag_pre_multiply(V_tau, V_L_Omega) * diag_pre_multiply(V_tau, V_L_Omega)';
    W = diag_pre_multiply(W_tau, W_L_Omega) * diag_pre_multiply(W_tau, W_L_Omega)';
}
model {
    real log_lik_obs[N];
    real log_lik;
    V_L_Omega ~ lkj_corr_cholesky(1);
    W_L_Omega ~ lkj_corr_cholesky(1);
    
    {
        vector[P] a;
        matrix[P, P] R;
        vector[M] f;
        matrix[M, M] Q;
        matrix[M, M] Qinv;
        vector[P] m[N+1];
        matrix[P, P] C[N+1];
        matrix[P, P] U_C[N+1];
        matrix[P, P] D_C[N+1];
        vector[M] mean_zero = rep_vector(0, M);
        
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
        
        for(t in 1:N){
            matrix[2*P, P] M_R;                     // M_R[n]' * M_R[n] = R[n]
            matrix[P, P] U_R;
            matrix[P, P] D_R;
            matrix[2*P, P] sqrt_D_R;
            
            vector[M] err;
            matrix[P+M, P] M_C;                     // U_R[N] * M_C[n]' * M_C[n] * U_R[N]' = C[n]^-1
            matrix[P, P] V_M_C;
            matrix[P+M, P] Dinv_M_C;
            
            a = GG * m[t];
            M_R = make_MR(D_C[t], U_C[t], GG, N_W);
            R = M_R' * M_R;                         //R = GG * C[t] * GG' + W;
            U_R = V_svd(M_R);
            sqrt_D_R = D_svd(M_R, 0);
            D_R = sqrt_D_R' * sqrt_D_R;
            
            f = FF * a;
            Q = FF * R * FF' + V;
            log_lik_obs[t] = multi_normal_lpdf(y[t] | f, Q);
            
            err = y[t] - f;
            Qinv = inverse(Q);
            m[t+1] = a + R * FF' * Qinv * err;
            
            M_C = make_MC(N_V, FF, U_R, D_R);
            V_M_C = V_svd(M_C);
            Dinv_M_C = D_svd(M_C, 1);
            
            U_C[t+1] = U_R * V_M_C;
            D_C[t+1] = Dinv_M_C' * Dinv_M_C;
            C[t+1] = U_C[t+1] * D_C[t+1] * U_C[t+1]';
        }
    }
        
    log_lik = sum(log_lik_obs);
    target += log_lik;
    //y ~ gaussian_dlm_obs(FF, GG, V, W, m0, C0);
}
generated quantities {
}
