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
    
        int r = min(nrow, ncol);
        matrix[ncol, ncol] MT_M;
        vector[r] MT_M_eigenvalues;
        matrix[nrow, ncol] D = rep_matrix(0, nrow, ncol);
    
        // The eigenvalues_sym funciton returns a vector of eigenvalues 
        // of a symmetric matrix A in ascending order.
        MT_M = M' * M;
        MT_M_eigenvalues = eigenvalues_sym(MT_M);
    
        // Insert the eigenvalues in decreasing order.
        // largest eigenvalue in r
        for(i in 1:r) {
            for(j in 1:r) {
                if(i == j){
                    D[i, j] = sqrt(MT_M_eigenvalues[r-i+1]);
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
            U[, i] = backward_U[, ncol-i+1];
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
    cov_matrix[M] W;        // Known constant covariance matix for state transition equation
    
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
    cov_matrix[M] Q_inv[N];
    
    // mean and variance-covariance matix for filtering distribution of
    // state t given all observations through time t
    vector[M] err[N];
    vector[P] m[N+1];
    cov_matrix[P] C[N+1];
    matrix[P, P] U_C[N+1];
    matrix[P, P] D_C[N+1];
    
    // Kalman filter
    
    // Intialize
    m[1] = m0;
    
    
    C[1] = C0;
    
    
    for(n in 1:N) {
        a[n] = GG * m[n];
        R[n] = GG * C[n] * GG' + W;
        
        f[n] = FF * a[n];
        Q[n] = FF * R[n] * FF' + V;
        
        Q_inv[n] = inverse(Q[n]);
        err[n] = y[n] - f[n];
        m[n+1] = a[n] + R[n] * FF' * Q_inv[n] * err[n];
        C[n+1] = R[n] - R[n] * FF' * Q_inv[n] * FF * R[n];
    }
    }
