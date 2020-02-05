/*
    Functions to compute singular value decomposition

See Appendix B of Dynamic Linear Models with R - Petris et. al.

Using this script to mainly to test functions
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
    matrix[p, q] D = D_svd(M);
    matrix[q, q] V = V_svd(M);
    matrix[p, p] U = U_svd(M, V, D);
    matrix[p, q] M_svd = U * D * V';
}
