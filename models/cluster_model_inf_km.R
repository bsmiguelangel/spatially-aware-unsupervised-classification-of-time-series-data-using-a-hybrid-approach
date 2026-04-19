cluster_model <- nimbleCode({
  
  for(Station in 1:NStations) {
    # Likelihood
    y[Station, 1:NBasis] ~ dmnorm(mu_c[1:NBasis, z[Station]], 
                                  w_tau[1:NBasis, 1:NBasis, Station])
    # Precision matrix w_tau
    w_tau[1:NBasis, 1:NBasis, Station] <- w[Station] * tau[1:NBasis, 1:NBasis]
    # Scaling parameter w
    w[Station] ~ dgamma(2, 2)
    # Allocation variable z
    z[Station] ~ dcat(omega[Station, 1:NClusters])
    
    for (Cluster in 1:NClusters) {
      omega[Station, Cluster] <- phi[Station, Cluster]/sum(phi[Station, 1:NClusters])
      # Linear predictor for the allocation probabilities of each station
      log(phi[Station, Cluster]) <- psi * (beta_0[Cluster] + 
        sd.alpha[Cluster] * alpha[Station, Cluster] + 
        sd.theta[Cluster] * theta[Station, Cluster]) + 
        (1 - psi) * (gamma[Cluster] * z_0[Station, Cluster])
      
      # Prior distribution: Non-structured spatial random effect alpha
      alpha[Station, Cluster] ~ dnorm(0, sd = 1)
    }
  }
  
  # Prior distribution: Weighting parameter psi
  # # Non-informative prior
  # psi ~ dbeta(0.5, 0.5) # Jeffreys prior
  
  # # Informative priors for k-means
  psi ~ dbeta(1, 14) # Moderately informative
  
  beta_0[1] <- 0
  for (Cluster in 2:NClusters) {
    # Prior distribution: Intercept
    beta_0[Cluster] ~ dnorm(0, sd = 1)
  }
  
  for (Cluster in 1:NClusters) {
    # Prior distribution: Cluster means mu_c
    mu_c[1:NBasis, Cluster] ~ dmnorm(mu_0[1:NBasis], tau_0[1:NBasis, 1:NBasis])
    
    # Structured spatial random effect: ICAR prior distribution theta
    theta[1:NStations, Cluster] ~ dcar_normal(adj = adj[1:Nadj], weights = weights[1:Nadj], 
                                              num = num[1:NStations], tau = 1, zero_mean = 1)
    
    # Prior distribution: 1st method (kmeans) gamma
    gamma[Cluster] ~ T(dnorm(0, sd = 1), 0, )
    
    # Prior distribution: Hyperparameters sd.alpha and sd.theta
    sd.alpha[Cluster] ~ T(dnorm(0, sd = 1), 0, )
    sd.theta[Cluster] ~ T(dnorm(0, sd = 1), 0, )
    
    # Zero-mean constraint for alpha
    zero.alpha[Cluster] ~ dnorm(mean.alphas[Cluster], 10000)
    mean.alphas[Cluster] <- mean(alpha[1:NStations, Cluster])
  }
  
  # Prior distribution: Precision matrix tau
  tau[1:NBasis, 1:NBasis] ~ dwish(wish_V[1:NBasis, 1:NBasis], NBasis)
})
