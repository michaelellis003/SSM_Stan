/*
Supporting functions for singular value decomposition based Kalman filter and
smooth

Resources
- Appendix B of Dynamic Linear Models with R  (Petris et. al.)
- Fixed-Interval Smoothing Algorithm Based on singular value decomposition (Zhang and Li)
*/
matrix D_svd(matrix M, int invert) {
    /* 
    INPUT:
    - M is an p x q martix. M = UDV'
    - invert = 0 or invert = 1
    OUTPUT:
    - Function returns "D" matrix of singular value decomposition (SVD)
    - If invert = 1 then return D_inv where D_inv[i,i] = 1/D[i,i]
    - D is a p × q diagonal matrix where the first r diagonal 
        entries are the square roots of the eigenvalues of M'M in
        decreasing order and all other entires are zero 
    */
    
    int p = rows(M);
    int q = cols(M);
    int r = min(p, q);
    matrix[p, q] D = rep_matrix(0, p, q);
    vector[q] MT_M_eigenvalues;
    
    if(invert != 0 && invert != 1) {
        reject("In D_svd function 'invert' variable can only equal 1 or 0");
    }
    
    // The eigenvalues_sym funciton returns a vector of eigenvalues
    // of a symmetric matrix A in ascending order.
    MT_M_eigenvalues = eigenvalues_sym( M' * M);
                                            
    // Insert the q eigenvalues in decreasing order.
    // largest eigenvalue in q
    for(i in 1:r) {
        real value = MT_M_eigenvalues[q-i+1];
        if(value >= 0) {
            if(invert == 0){
                D[i, i] = sqrt(value);
            }
            if(invert == 1) {
                D[i, i] = 1/sqrt(value);
            }
        }
    }
                                            
    return D;
}

matrix Dinv_svd(matrix D) {
    // M is an p x q martix
    // Singular Value Decomposition M = UDV'
    // return D_inv where D_inv[i,i] = 1/D[i,i]
    
    int p = rows(D);
    int q = cols(D);
    int r = min(p, q);
    matrix[p, q] D_inv = rep_matrix(0, p, q);
    
    for(i in 1:r) {
        real value = D[i, i];
        if(value >= 0) {
            D_inv[i, i] = 1/value;
        }
    }
                                            
    return D_inv;
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
    D = D_svd(M, 0);
    U = U_svd(M, V, D);
    S = sqrt(D);        
    
    N = S * V';
    
    return(N);
}

matrix make_MR(matrix D_C, matrix U_C, matrix GG, matrix N_W) {
    
    int p = rows(N_W);
    matrix[2*p, p] MR;
    matrix[p, p] S_C;
    matrix[p, p] MR_top;
    matrix[p, p] MR_bottom;
    
    S_C = sqrt(D_C);
    MR_top = S_C * U_C' * GG';
    MR_bottom = N_W;
    
    MR = append_row(MR_top, MR_bottom);

    return MR;
}

matrix make_MC(matrix N_V, matrix FF, matrix U_R, matrix D_R) {
    
    int m = rows(N_V);
    int p = rows(U_R);
    matrix[m+p, p] MC;
    matrix[p, p] S_R_inv = rep_matrix(0, p, p);
    matrix[m, p] MC_top;
    matrix[p, p] MC_bottom;
    
    S_R_inv = sqrt(Dinv_svd(D_R));
    MC_top = N_V * FF * U_R;
    MC_bottom = S_R_inv;
    
    MC = append_row(MC_top, MC_bottom);

    return MC;
}
