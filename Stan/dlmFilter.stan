/*
Implement Kalman filter for multivariate Dynamic Linear Model with
known and constant covariance matices for the observation equation 
and the state transition equation.

See Chapter 2.7 of Dynamic Linear Models with R - Petris et. al.

TO DO: Add handling of missing values 
*/
data {
    int<lower=1> N;         // Number of observations
    int<lower=1> M;         // Dimension of observation
    int<lower=1> P;         // Dimension of state vector
    
     // Observation Stuff
    vector[M] y[N];         // Observations
    matrix[M, P] FF;         // Observation equation matrix
    cov_matrix[M] V;        // Known constant covariance matix for observation equation
    
     // State Stuff
    matrix[P, P] GG;         // State/transition equation matrix
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
    
    // Intialize
    m[1] = m0;
    C[1] = C0;
    
    // Kalman filter
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
