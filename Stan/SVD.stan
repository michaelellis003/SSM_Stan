/*
Functions to compute singular value decomposition

See Appendix B of Dynamic Linear Models with R - Petris et. al.

Using this script to mainly to test functions
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
        print(ncol);
        print(eigenvalues_sym(MT_M));
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
    
}
data {
    int<lower=1> p;
    int<lower=1> q;
    matrix[p, q] M;       // Matrix to decompose
}
parameters {
}
transformed parameters {
}
model {
}
generated quantities {
    matrix[p, p] U = U_svd(M, p, q);
    matrix[p, q] D = D_svd(M, p, q);
    matrix[q, q] V = V_svd(M, p, q);
}
