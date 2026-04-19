### Required packages ###

pacman::p_load(sf, spdep, ggplot2, RColorBrewer, nimble, MCMCvis, ggpubr, 
               install = FALSE)

### Data loading ###

load(file.path("outputs", "data.rda"))
load(file.path("outputs", "df_coefs.rda"))
estaciones_coords <- read.csv(file = file.path("data", "estaciones_coords.csv"))

unique(df_coefs$station[df_coefs$pollutant == "O3"])
estaciones_data <- match(unique(df_coefs$code[df_coefs$pollutant == "O3"]), unique(estaciones_coords$codigo))
estaciones_coords$nombre[estaciones_data]
# Number of basis functions
num_basis <- 50
y <- matrix(df_coefs$spline_coefficients[df_coefs$pollutant == "O3"], 
            ncol = num_basis, byrow = TRUE)

NStations <- nrow(y)
NBasis <- ncol(y)

### Models ###

NClusters <- 4

load(file = file.path("outputs", "dist_km.rda"))
unique(df_coefs$station[df_coefs$pollutant == "O3"])
estaciones_coords$nombre[estaciones_data]

# Number of neighbors of each station
num <- rep(NStations - 1, NStations)
# Neighbors of each station
adj <- c()
for (Station in 1:NStations) {
  adj <- c(adj, (1:NStations)[-Station])
}
# Matrix of weights
W <- 1/dist_km
# Vector of weights
weights <- c()
for (Station in 1:NStations) {
  weights <- c(weights, (W[, Station])[-Station])
}
weights <- as.numeric(weights)
# Sum of all the neighbor numbers of all stations
Nadj <- length(adj) # sum(num)

# k-means assignment
kmeans_assigned <- readRDS(file = file.path("outputs", paste0("kmeans_", NClusters, ".rds")))
kmeans_assigned <- as.numeric(kmeans_assigned)
table(kmeans_assigned)
kmeans_assigned

# Re-assignment of clusters
freqs <- sort(table(kmeans_assigned), decreasing = TRUE)
freqs

order_clusters <- names(freqs)
order_clusters

relabel_clusters <- setNames(seq_along(order_clusters), order_clusters)
relabel_clusters

kmeans_assigned <- as.numeric(relabel_clusters[as.character(kmeans_assigned)])
table(kmeans_assigned)
kmeans_assigned

z_0 <- matrix(nrow = NStations, ncol = NClusters)
for (Station in 1:NStations) {
  for (Cluster in 1:NClusters) {
    z_0[Station, Cluster] <- ifelse(kmeans_assigned[Station] == Cluster, 1, 0)
  }
}

### Model code ###

source(file = file.path("models", "cluster_model_non_inf.R"))
# source(file = file.path("models", "cluster_model_inf_km.R"))

modelCode <- cluster_model

### Data to be loaded ###

modelData <- list(y = y, zero.alpha = rep(0, NClusters))

modelConstants <- list(NBasis = NBasis, NStations = NStations, z_0 = z_0, 
                       NClusters = NClusters, mu_0 = rep(0, NBasis),
                       wish_V = diag(0.01, NBasis), tau_0 = diag(0.01, NBasis),
                       adj = adj, weights = weights, num = num, Nadj = Nadj)

### Parameters to be saved ###

modelParameters <-  c("mu_c", "tau", "z", "w", "omega", "beta_0",
                      "alpha", "theta", "sd.theta", "sd.alpha",
                      "gamma", "psi")

### Nimble call ###

# Let’s create the NIMBLE model, creates the nodes (inits should be passed now)
# calculate = FALSE helps speed
model <- nimbleModel(code = modelCode, 
                     constants = modelConstants, 
                     data = modelData, 
                     inits = list(tau = diag(1, NBasis),
                                  z = sample(1:NClusters, NStations, replace = TRUE),
                                  mu_c = matrix(rep(0, NBasis * NClusters), ncol = NClusters),
                                  w = rep(1, NStations),
                                  alpha = matrix(rnorm(NStations * NClusters, sd = 0.1),
                                                 nrow = NStations, ncol = NClusters),
                                  beta_0 = c(NA, rnorm(NClusters - 1)),
                                  gamma = runif(NClusters),
                                  sd.alpha = runif(NClusters),
                                  sd.theta = runif(NClusters),
                                  psi = rbeta(1, 1/2, 1/2),
                                  theta = matrix(rnorm(NStations * NClusters, sd = 0.1),
                                                 nrow = NStations, ncol = NClusters))
                     , calculate = FALSE)

# Compile the model, which means generating C++ code, compiling that code, and loading it back into R
Cmodel <- compileNimble(model)
# Let's change the default configuration
# The conjugacy checking is a slow part of MCMC configuration, skipping it 
# (useConjugacy = FALSE) helps speed
modelMCMCconfiguration <- configureMCMC(model, useConjugacy = TRUE,
                                        enableWAIC = TRUE)

# Add new monitors
modelMCMCconfiguration$monitors <- c()
modelMCMCconfiguration$addMonitors(modelParameters)
# Build MCMC object
modelMCMC <- buildMCMC(modelMCMCconfiguration)
# Need to reset the nimbleFunctions in order to add the new MCMC
CmodelMCMC <- compileNimble(modelMCMC, project = model,
                            resetFunctions = TRUE)
# # Results
# results <- runMCMC(CmodelMCMC, niter = 150000, nburnin = 50000, thin = 100,
#                    nchains = 1, summary = TRUE, WAIC = TRUE)

# saveRDS(results, file = file.path("results", paste0("cluster_model_non_inf_", NClusters, ".rds")))
# saveRDS(results, file = file.path("results", paste0("cluster_model_inf_km_", NClusters, ".rds")))

### Figure 7 ###

results_non_inf <- readRDS(file = file.path("results", paste0("cluster_model_non_inf_", NClusters, ".rds")))
results_inf_km <- readRDS(file = file.path("results", paste0("cluster_model_inf_km_", NClusters, ".rds")))

# Spatial model with non-informative k-means

omega <- matrix(nrow = NStations, ncol = NClusters)

index_omegas <- grep("omega", rownames(results_non_inf$summary))
aux_omegas <- results_non_inf$summary[index_omegas, 1]

omega <- matrix(aux_omegas, byrow = FALSE, nrow = NStations)

cluster_assigned_non_inf <- apply(omega, 1, which.max)
table(cluster_assigned_non_inf)

# Spatial model with informative k-means

omega <- matrix(nrow = NStations, ncol = NClusters)

index_omegas <- grep("omega", rownames(results_inf_km$summary))
aux_omegas <- results_inf_km$summary[index_omegas, 1]

omega <- matrix(aux_omegas, byrow = FALSE, nrow = NStations)

cluster_assigned_inf_km <- apply(omega, 1, which.max)
table(cluster_assigned_inf_km)

waic_spat_non_inf <- results_non_inf$WAIC$WAIC
waic_spat_inf_km <- results_inf_km$WAIC$WAIC

p_kmeans <- ggplot() +
  geom_sf(data = carto_muni_sf, fill = "#FFF9F0", color = "grey50") +
  geom_sf(data = cbind(locations, kmeans_assigned = factor(kmeans_assigned)), 
          aes(color = factor(kmeans_assigned)), size = 2) +
  scale_color_manual(values = c("1" = "#E41A1C", "2" = "#377EB8",
                                "3" = "#FF7F00", "4" = "#4DAF4A", 
                                "5" = "#984EA3")) +
  theme_minimal() + labs(color = "Cluster") + ggtitle("Stage 1: k-means")

p_non_inf <- ggplot() +
  geom_sf(data = carto_muni_sf, fill = "#FFF9F0", color = "grey50") +
  geom_sf(data = cbind(locations, cluster_assigned_non_inf = factor(cluster_assigned_non_inf)), 
          aes(color = factor(cluster_assigned_non_inf)), size = 2) +
  scale_color_manual(values = c("1" = "#E41A1C", "2" = "#377EB8",
                                "3" = "#FF7F00", "4" = "#4DAF4A", 
                                "5" = "#984EA3")) +
  theme_minimal() + labs(color = "Cluster", x = NULL, y = NULL) +
  ggtitle(expression(paste("Stage 2: Spatial model with ", psi, "  ~ Beta(1/2,1/2)"))) + 
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
    label = paste0("WAIC = ", round(waic_spat_non_inf, 2)), size = 3, 
    color = "black")

p_inf_km <- ggplot() +
  geom_sf(data = carto_muni_sf, fill = "#FFF9F0", color = "grey50") +
  geom_sf(data = cbind(locations, cluster_assigned_inf_km = factor(cluster_assigned_inf_km)), 
          aes(color = factor(cluster_assigned_inf_km)), size = 2) +
  scale_color_manual(values = c("1" = "#E41A1C", "2" = "#377EB8",
                                "3" = "#FF7F00", "4" = "#4DAF4A", 
                                "5" = "#984EA3")) +
  theme_minimal() + labs(color = "Cluster", x = NULL, y = NULL) +
  ggtitle(expression(paste("Stage 2: Spatial model with ", psi, "  ~ Beta(1,14)"))) + 
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
    label = paste0("WAIC = ", round(waic_spat_inf_km, 2)), size = 3,
    color = "black")

ggarrange(p_kmeans, p_non_inf, p_inf_km, ncol = 3,
          common.legend = TRUE, legend = "bottom")

# ggsave(paste0(file.path("Figures", "Figure7.jpg")), 
#        device = "jpg", width = 15, height = 7, dpi = 1200)

### Convergence assessment ###

NClusters <- 4
results <- readRDS(file = file.path("results", paste0("cluster_model_non_inf_", NClusters, ".rds")))
# results <- readRDS(file = file.path("results", paste0("cluster_model_inf_km_", NClusters, ".rds")))

MCMCsummary(object = results$samples, params = "psi",
            # exact = TRUE,
            # ISB = FALSE,
            round = 4)

MCMCtrace(object = results$samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          Rhat = TRUE,
          n.eff = TRUE,
          params = "alpha")
