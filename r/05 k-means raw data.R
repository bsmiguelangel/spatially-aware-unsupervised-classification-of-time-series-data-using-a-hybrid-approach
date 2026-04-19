### Required packages ###

# install.packages("pacman")
pacman::p_load(fda, zoo, ggplot2, install = FALSE)

### Figure 4 ###

Sys.setlocale("LC_ALL", "English")

# Beta priors
beta_params <- matrix(c(14, 1, 0.5, 0.5, 1, 14), ncol = 2, byrow = TRUE)
data_plot <- c()
for (i in 1:nrow(beta_params)) {
  data_plot <- rbind(data_plot, 
                     data.frame(x = seq(0, 1, 0.001),
                                y = dbeta(seq(0, 1, 0.001),
                                          beta_params[i, 1], beta_params[i, 2]),
                                Prior = paste0("Beta(", beta_params[i, 1], ",", beta_params[i, 2],")")))
}

ggplot(aes(x = x, y = y, col = Prior), data = data_plot) +
  geom_line(size = 1) + theme_bw() + 
  xlab(expression(psi)) + ylab(expression(italic(p)(psi))) +
  scale_color_manual(values = brewer.pal(3, "Set2")) +
  theme(text = element_text(size = 18),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# ggsave(paste0(file.path("figures", "Figure4.jpg")),
#        width = 7, height = 6, dpi = 1200)

### Data loading ###

load(file = file.path("outputs", "data.rda"))
load(file = file.path("outputs", "df_missing.rda"))

### Figure 5 ###

df_raw <- c()
pollutants <- c("NO2", "SO2", "O3")
for (pollutant in pollutants) {
  stations <- unique(df_missing$station[df_missing$pollutant == pollutant & df_missing$Prop_20 == 1])
  for (station in stations) {
    print(paste0(pollutant, " ", station))
    t_points <- seq(0.5, 366, by = 1)
    time_series <- data[data$NOM_ESTACION == station, pollutant]
    time_series <- na.approx(time_series, na.rm = FALSE, rule = 2)
    df_raw <- rbind(df_raw, data.frame(pollutant,
                                       code = unique(data$COD_ESTACION[data$NOM_ESTACION == station]),
                                       station, 
                                       time = 1:366, 
                                       value = as.numeric(time_series)))
  }
}

# K-means on the raw time series
df_raw_filt <- df_raw[df_raw$pollutant == "O3", ]
values <- matrix(df_raw_filt$value, nrow = length(unique(df_raw_filt$station)))
nstart <- 50
wss <- function(k) {
  kmeans(values, centers = k, nstart = nstart)$tot.withinss
}

k.values <- 1:40
wss_values <- sapply(k.values, wss)

df <- data.frame(k = k.values, wss = wss_values)

ggplot(df, aes(x = k, y = wss)) + geom_line() + geom_point(size = 2) +
  labs(x = "Number of clusters", y = "Within-cluster Sum of Squares (WSS)") +
  ylim(0, 3000000) + theme_minimal()

# ggsave(paste0(file.path("figures", "Figure5.jpg")), 
#        device = "jpg", width = 6.5, height = 4.5, dpi = 1200)

centers <- 4
aux <- kmeans(values, centers = centers, nstart = nstart)
table(aux$cluster)
cluster_assigned <- aux$cluster

names(cluster_assigned) <- unique(df_raw_filt$station)
df_raw_filt$cluster_assigned <- cluster_assigned[df_raw_filt$station]

kmeans_assigned <- as.numeric(cluster_assigned)
table(kmeans_assigned)
kmeans_assigned

# saveRDS(cluster_assigned, file = file.path("outputs", paste0("kmeans_", centers, ".rds")))

### Figure 6 ###

NClusters <- 4

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

names(kmeans_assigned) <- unique(df_raw_filt$station)
df_raw_filt$kmeans_assigned <- kmeans_assigned[df_raw_filt$station]
label_map <- c("1" = "Cluster 1", "2" = "Cluster 2",
               "3" = "Cluster 3", "4" = "Cluster 4")
ggplot(aes(x = time, y = value, group = station, col = as.factor(kmeans_assigned)), 
       data = df_raw_filt) + geom_line(linewidth = 0.1) +
  scale_colour_manual(name = "Cluster", 
                      values = c("1" = "#E41A1C", "2" = "#377EB8",
                                 "3" = "#FF7F00", "4" = "#4DAF4A")) +
  labs(y = "Ozone level (µg/m³)", x = NULL) + theme_bw() +
  facet_wrap(~kmeans_assigned, labeller = as_labeller(label_map)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
  theme(legend.position = "none", text = element_text(size = 18))

# ggsave(paste0(file.path("figures", "Figure6.jpg")), 
#        device = "jpg", width = 9.5, height = 6.5, dpi = 1200)
