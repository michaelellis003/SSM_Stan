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
    
    matrix D_svd(matrix M) {
            // M is an p x q martix
            // Singular Value Decomposition M = UDV'
            // D is a p × q diagonal matrix where the first r diagonal 
            // entries are the square roots of the eigenvalues of M'M in
            // decreasing order and all other entires are zero
            
            int p = rows(M);
            int q = cols(M);
            int r = min(p, q);
            matrix[p, q] D = rep_matrix(0, p, q);
            vector[q] MT_M_eigenvalues;
            
            // The eigenvalues_sym funciton returns a vector of eigenvalues
            // of a symmetric matrix A in ascending order.
            MT_M_eigenvalues = eigenvalues_sym( M' * M);
                print(MT_M_eigenvalues);
                                                    
            // Insert the q eigenvalues in decreasing order.
            // largest eigenvalue in q
            for(i in 1:r) {
                real value = MT_M_eigenvalues[q-i+1];
                if(value >= 0) {
                    D[i, i] = sqrt(value);
                }
            }
                                                    
            return D;
    }
    
    matrix V_svd(matrix M) {
        // Singular Value Decomposition M = UDV'
        // M is an p x q martix
        // V is an q × q orthogonal matrix whose columns are the 
        // eigenvectors of M'M 
            
        int p = rows(M);
        int q = cols(M);
        matrix[q, q] MT_M;
        matrix[q, q] backward_V;
        matrix[q, q] V;
            
        // The eigenvectors_sym function returns a matrix with the (column)
        // eigenvectors of a symmetric matrix in the same order as the
        // eigenvalues returned by eigenvalues_sym
        MT_M = M' * M;
        backward_V = eigenvectors_sym(MT_M);
        
        
        // Reverse the order of the  eigenvectors to match that the order of
        // the eigenvalues returned by eigenvalues_sym
        for(i in 1:q){
            V[, i] = backward_V[, q-i+1];
        }
        
        return V;
    }
    
    matrix U_svd(matrix M, matrix V, matrix D) {
        // Singular Value Decomposition M = UDV'
        // M is an p x q martix
        // U is an p × p orthogonal matrix whose columns are the 
        // eigenvectors of M'M 
            
        int p = rows(M);
        int q = cols(M);
        int r = min(p, q);
        matrix[p, p] U = diag_matrix(rep_vector(1, p));
        real sing_value;
        
        for(i in 1:p) {
            if(i <= r) {
                sing_value = D[i,i];
                if(sing_value > 0){
                    U[, i] = (1/sing_value) * M * V[, i];
                }
            }
        }
        
        return U;
    }
    
    matrix sqrt_svd(matrix M) {
        // Function to calculate sqrt of covariance matrix using SVD
        // where the sqrt of a matrix is M^(1/2) such that
        // M = M^(1/2) * M^(1/2)
        
        int p = rows(M);
        int q = cols(M);
        matrix[q, q] V;
        matrix[p, q] D;
        matrix[p, q] S;
        matrix[p, p] U;
        matrix[p, p] N; // square root of M;
        
        if(p != q) {
            reject("Attempting to take square root of ", 
            "non-square matix using SVD");
        }
        
        V = V_svd(M);
        D = D_svd(M);
        U = U_svd(M, V, D);
        S = sqrt(D);        
        
        N = S * V';
        
        return(N);
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
        
        
        S_R_inv =  diag_matrix(1 ./ sqrt(diagonal(D_R)));
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
    matrix[P, P] C[N+1];
    matrix[P, P] U_C[N+1];
    matrix[P, P] D_C[N+1];
    matrix[P+M, M] M_C[N];          // U_R[N] * M_C[n]' * M_C[n] * U_R[N]' = C[n]^-1
    matrix[P, P] U_M_C;
    matrix[M, M] V_M_C;
    matrix[P+M, M] D_M_C;
    
    // Stuff for square root of W
    matrix[M, M] N_W;
    matrix[M, M] U_W = U_svd(W);
    matrix[M, M] S_W = sqrt(D_svd(W));
    
    // Stuff for square root of V^-1
    matrix[P, P] N_V;
    matrix[P, P] U_V = U_svd(V);
    matrix[P, P] S_V = sqrt(D_svd(V));
    
    // Get square root of W and V^-1
    N_W = S_W * U_W';
    N_V = S_V * U_V';
    
    // Kalman filter
    // Intialize
    m[1] = m0;
    U_C[1] = U_svd(C0);
    D_C[1] = D_svd(C0);
    
    for(n in 1:N) {
        a[n] = GG * m[n];
        M_R[n] = make_MR(D_C[n], U_C[n], GG, N_W, P);
        R[n] = M_R[n]' * M_R[n];
        U_R[n] = V_svd(M_R[n]);
        D_R[n] = D_svd(M_R[n]);
        
        f[n] = FF * a[n];
        Q[n] = FF * R[n] * FF' + V;
        
        Q_inv[n] = inverse(Q[n]);
        err[n] = y[n] - f[n];
        m[n+1] = a[n] + R[n] * FF' * Q_inv[n] * err[n];
        
        M_C[n] = make_MC(N_V, FF, U_R[n], D_R[n], M, P);
        V_M_C = V_svd(M_C[n]);
        D_M_C = D_svd(M_C[n]);
        
        U_C[n+1] = U_R[n] * V_M_C;
        D_C[n+1] = diag_matrix(1 ./ diagonal(D_M_C .* D_M_C));
        
        C[n+1] = U_C[n+1] * D_C[n+1] * U_C[n+1]';
        print(n+1);
        print(C[n+1]);
    }
}
