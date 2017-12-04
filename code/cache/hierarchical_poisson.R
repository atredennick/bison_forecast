
#### HIERARCICAL POISSON
x_model <- "  
model{

#### Variance Priors
#tau_proc ~ dgamma(0.0001, 0.0001)
#sigma_proc <- 1/sqrt(tau_proc)
shape_p ~ dunif(0, 100)

#### Fixed Effects Priors
b0 ~ dnorm(0, 0.001)
b1 ~ dnorm(0, 0.001)
b2 ~ dnorm(0, 0.001)

#### Initial Conditions
N0      ~ dunif(0,10)
Nmed[1] <- b0 + b1*N0 + b2*x[1]
N[1]    ~ dgamma(shape_p, shape_p / exp(Nmed[1])) 

#### Process Model
for(t in 2:n){
Nmed[t] <- b0 + b1*N[t-1] + b2*x[t]
N[t]    ~ dgamma(shape_p, shape_p / exp(Nmed[t])) 
}

#### Data Model
for(t in 1:n){
var_obs[t] <- sd_obs[t]*sd_obs[t]
shape[t]   <- N[t]*N[t]/var_obs[t]
rate[t]    <- N[t]/var_obs[t]
lambda[t]  ~ dgamma(shape[t], rate[t])
Nobs[t]    ~ dpois(lambda[t])
}

}"
