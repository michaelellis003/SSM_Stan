/*
Implement Kalman filter for multivariate Dynamic Linear Model with
known and constant covariance matices for the observation equation 
and the state transition equation. - See Chapter 2.7 of Dynamic Linear Models 
with R - Petris et. al.

Uses singular value decomposition to avoid numerical in instability when
updating filtering distribution covariance matrix - Appendix B of Dynamic Linear 
Models  with R - Petris et. al.

TO DO: Add handling of missing values 
*/
functions {
        
    matrix D_svd(matrix M, int nrow, int ncol) {
        // Singular Value Decomposition M = UDV'
        // D is a nrow × ncol diagonal matrix where the first r diagonal 
        // entries are the square roots of the eigenvalues of M'M in
        // decreasing order and all other entires are zero
    
        matrix[ncol, ncol] MT_M;
        vector[ncol] MT_M_eigenvalues;
        matrix[nrow, ncol] D = rep_matrix(0, nrow, ncol);
    
        // The eigenvalues_sym funciton returns a vector of eigenvalues 
        // of a symmetric matrix A in ascending order.
        MT_M = M' * M;
        MT_M_eigenvalues = eigenvalues_sym(MT_M);
    
        // Insert the eigenvalues in decreasing order.
        // largest eigenvalue in ncol
        for(i in 1:nrow) {
            for(j in 1:ncol) {
                if(i == j){
                    D[i, j] = sqrt(MT_M_eigenvalues[ncol-j+1]);
                }
            }
        }
    
        return D;
    }
    
    matrix U_svd(matrix M, int nrow, int ncol) {
        // Singular Value Decomposition M = UDV'
        // U is an nrow × nrow orthogonal matrix whose columns are the 
        // eigenvectors of MM'
        
        matrix[nrow, nrow] M_MT;
        matrix[nrow, nrow] backward_U;
        matrix[nrow, nrow] U;
        
        // The eigenvectors_sym function returns a matrix with the (column) 
        // eigenvectors of a symmetric matrix in the same order as the eigenvalues
        // returned by eigenvalues_sym
        M_MT = M * M';
        backward_U = eigenvectors_sym(M_MT);
        
        // Reverse the order of the  eigenvectors to match that the order of
        // the eigenvalues returned by eigenvalues_sym
        for(i in 1:nrow){
            U[, i] = backward_U[, nrow-i+1];
        }
        
        return U;
    }
    
    matrix V_svd(matrix M, int nrow, int ncol) {
        // Singular Value Decomposition M = UDV'
        // V is an ncol × ncol orthogonal matrix whose columns are the 
        // eigenvectors of M'M
        
        matrix[ncol, ncol] MT_M;
        matrix[ncol, ncol] backward_V;
        matrix[ncol, ncol] V;
        
        // The eigenvectors_sym function returns a matrix with the (column) 
        // eigenvectors of a symmetric matrix in the same order as the eigenvalues
        // returned by eigenvalues_sym
        MT_M = M' * M;
        backward_V = eigenvectors_sym(MT_M);
        
        // Reverse the order of the  eigenvectors to match that the order of
        // the eigenvalues returned by eigenvalues_sym
        for(i in 1:ncol){
            V[, i] = backward_V[, ncol-i+1];
        }
        
        return V;
    }
    
    matrix make_MR(matrix D_C, matrix U_C, matrix GG, matrix N_W, int P) {
        
        matrix[2*P, P] MR;
        matrix[P, P] S_C;
        matrix[P, P] MR_top;
        matrix[P, P] MR_bottom;
        
        S_C = sqrt(D_C);
        MR_top = S_C * U_C' * GG';
        MR_bottom = N_W;
        
        MR = append_row(MR_top, MR_bottom);
        
        return MR;
    }
    
    matrix make_MC(matrix N_V, matrix FF, matrix U_R, matrix D_R, int M, int P) {
        
        matrix[P+M, P] MC;
        matrix[P, P] S_R_inv;
        matrix[M, P] MC_top;
        matrix[P, P] MC_bottom;
        
        print(D_R);
        print(sqrt(D_R));
        S_R_inv = 1 ./ sqrt(D_R);
        print(S_R_inv);
        MC_top = N_V * FF * U_R;
        MC_bottom = S_R_inv;
        
        MC = append_row(MC_top, MC_bottom);
    
        return MC;
    }
    
}
data {
    int<lower=1> N;         // Number of observations
    int<lower=1> M;         // Dimension of observation
    int<lower=1> P;         // Dimension of state vector
    
    // Observation Stuff
    vector[M] y[N];         // Observations
    matrix[M, P] FF;        // Observation equation matrix
    cov_matrix[M] V;        // Known constant covariance matix for observation equation
    
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
    matrix[2*M, M] M_R[N];          // M_R[n]' * M_R[n] = R[n]
    
    // mean and variance-covariance matix for one-step-ahead (Gaussian) predictive
    // distribution of observation t given all observations through time t-1
    vector[M] f[N];
    cov_matrix[M] Q[N];
    matrix[M, M] Q_inv[N];
    
    // mean and variance-covariance matix for filtering distribution of
    // state t given all observations through time t
    vector[M] err[N];
    vector[P] m[N+1];
    cov_matrix[P] C[N+1];
    matrix[P, P] U_C[N+1];
    matrix[P, P] D_C[N+1];
    matrix[P+M, M] M_C[N];          // U_R[N] * M_C[n]' * M_C[n] * U_R[N]' = C[n]^-1
    matrix[P, P] U_M_C;
    matrix[M, M] V_M_C;
    matrix[P+M, M] D_M_C;
    
    // Stuff for square root of W
    matrix[M, M] N_W;
    matrix[M, M] U_W = U_svd(W, M, M);
    matrix[M, M] S_W = sqrt(D_svd(W, M, M));
    
    // Stuff for square root of V^-1
    matrix[P, P] N_V;
    matrix[P, P] U_V = U_svd(V, P, P);
    matrix[P, P] S_V = sqrt(D_svd(V, P, P));
    
    // Get square root of W and V^-1
    N_W = S_W * U_W';
    N_V = S_V * U_V';
    
    // Kalman filter
    // Intialize
    m[1] = m0;
    U_C[1] = U_svd(C0, P, P);
    D_C[1] = D_svd(C0, P, P);
    
    for(n in 1:N) {
        a[n] = GG * m[n];
        M_R[n] = make_MR(D_C[n], U_C[n], GG, N_W, P);
        R[n] = M_R[n]' * M_R[n];
        U_R[n] = V_svd(M_R[n], 2*P, P);
        D_R[n] = D_svd(M_R[n], 2*P, P);
        
        f[n] = FF * a[n];
        Q[n] = FF * R[n] * FF' + V;
        
        Q_inv[n] = inverse(Q[n]);
        err[n] = y[n] - f[n];
        m[n+1] = a[n] + R[n] * FF' * Q_inv[n] * err[n];
        
        M_C[n] = make_MC(N_V, FF, U_R[n], D_R[n], M, P);
        V_M_C = V_svd(M_C[n], P+M, P);
        D_M_C = D_svd(M_C[n], P+M, P);
        
        U_C[n+1] = U_R[n] * V_M_C;
        D_C[n+1] = 1 ./ (D_M_C .* D_M_C);
        
        C[n+1] = U_C[n+1] * D_C[n+1] * U_C[n+1]';
    }
}
