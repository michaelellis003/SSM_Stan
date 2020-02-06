library(data.table)
library(dlm)
library(rstan)

dat <- data.table( 
    y = as.vector(Nile),
    Year = 1871:1970
)

ggplot(dat, aes(x = Year, y = y)) + 
    geom_point() + 
    geom_line() + 
    theme_minimal() + 
    ylab("Nile")

mod1 <- dlmModPoly( 
    order = 1,
    dV = 15100,
    dW = 755
)

# build and fit model using stan -----------------------------------------------
stan_args <- list( 
    N = nrow(dat),          # Number of observations
    M = 1,            # Dimension of observation
    P = 1,            # Dimension of state vector
    
    # Observation Stuff
    y = matrix(dat$y, nrow = nrow(dat), ncol = 1),         # Observations
    FF = mod1$FF,      # Observation equation matrix
    V = mod1$V,        # Known constant covariance matix for observation equation
    
    # State Stuff
    GG = mod1$GG,       # State/transition equation matrix
    W = mod1$W,         # Known constant covariance matix for state transition equation
    
    m0 = as.array(mod1$m0),       # Initial state vector prior mean
    C0 = mod1$C0        # Initial state vector prior variance-covariance matix
)

mod1_Stan <- stan( 
    file = "Stan/dlmFilterSVD.stan",
    data = stan_args,
    chains = 1,
    iter = 1,
    warmup = 0,
    algorithm = "Fixed_param",
    cores = 1
)

dat[, mod1_stan_filtered := as.vector(unlist(extract(mod1_Stan, pars = "m")))[-1]]

# build and model using dlm ----------------------------------------------------

NileFilt <- dlmFilter( 
    Nile,
    mod1
)

dat[, mod1_filtered := NileFilt$m[-1]]

# plot and compare -------------------------------------------------------------

ggplot(dat) + 
    geom_point(aes(x = Year, y = y)) + 
    geom_line(aes(x = Year, y = mod1_filtered, color = "Model 1 - DLM")) + 
    geom_line(aes(x = Year, y = mod1_stan_filtered, color = "Model 1 - Stan"), linetype="dashed") + 
    theme_minimal() + 
    ylab("Nile")
