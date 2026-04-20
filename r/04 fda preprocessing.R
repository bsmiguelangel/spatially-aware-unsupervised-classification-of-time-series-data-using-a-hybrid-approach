### Required packages ###

# install.packages("pacman")
pacman::p_load(fda, zoo, ggplot2, gridExtra, install = FALSE)

### Data loading ###

load(file = file.path("outputs", "data.rda"))
load(file = file.path("outputs", "df_missing.rda"))

pollutants <- c("NO2", "SO2", "O3")
df_coefs <- c()
df_smoothed <- c()
# Number of basis functions
num_basis <- 50 
for (pollutant in pollutants) {
  stations <- unique(df_missing$station[df_missing$pollutant == pollutant & df_missing$Prop_20 == 1])
  for (station in stations) {
    print(paste0(pollutant, " ", station))
    t_points <- seq(0.5, 366, by = 1)
    time_series <- data[data$NOM_ESTACION == station, pollutant]
    time_series <- na.approx(time_series, na.rm = FALSE, rule = 2)
    time_range <- c(0, 366)
    # Order 4 corresponds to cubic B-splines
    spline_order <- 4 
    # Create the B-spline basis
    basis_obj <- create.bspline.basis(rangeval = time_range, nbasis = num_basis, norder = spline_order)
    # This gives us the coefficients of the basis functions, which represent our curve in a low-dimensional space
    smoothed_fd <- smooth.basis(argvals = t_points, y = time_series, fdParobj = basis_obj)$fd
    # Extract the coefficients
    spline_coefficients <- smoothed_fd$coefs
    # Visualization of the results
    fine_t_points <- seq(0, 366, length.out = 500)
    smoothed_curve <- eval.fd(fine_t_points, smoothed_fd)
    df_coefs <- rbind(df_coefs, data.frame(pollutant, 
                                           code = unique(data$COD_ESTACION[data$NOM_ESTACION == station]),
                                           station, 
                                           spline_coefficients = as.numeric(spline_coefficients)))
    df_smoothed <- rbind(df_smoothed, data.frame(pollutant, 
                                                 code = unique(data$COD_ESTACION[data$NOM_ESTACION == station]),
                                                 station, fine_t_points, 
                                                 smoothed_curve = as.numeric(smoothed_curve)))
  }
}
# save(df_coefs, file = file.path("Outputs", "df_coefs.rda"))
# save(df_smoothed, file = file.path("Outputs", "df_smoothed.rda"))

load(file = file.path("outputs", "df_coefs.rda"))
load(file = file.path("outputs", "df_smoothed.rda"))

df_smoothed_filt <- df_smoothed[df_smoothed$pollutant == "O3", ]
Sys.setlocale("LC_ALL", "English")
stations <- unique(df_smoothed_filt$station)
aux <- 1:length(stations)
names(aux) <- stations
df_smoothed_filt$num_station <- aux[df_smoothed_filt$station]
p1 <- ggplot(aes(x = fine_t_points, y = smoothed_curve, group = station, col = num_station), 
             data = df_smoothed_filt) + geom_line() + theme_bw() +
        xlab(NULL) + ylab(expression("Ozone level" ~ (mu*g/m^3)))+
        scale_y_continuous(limits=c(0, 140)) +
        theme(text = element_text(size = 18), 
              plot.title = element_text(face = "bold"),
              legend.position = "none")

# Raw vs. smoothed Station 1
load(file.path("outputs", "data_plot.rda"))
df_smoothed_filtb <- df_smoothed_filt[df_smoothed_filt$num_station==1, ]
df_smoothed_filtb$Type <- "Smoothed"
df_smoothed_filtb <- rbind(df_smoothed_filtb,data.frame(pollutant = "O3",
                                                        code = "3009006",
                                                        station = "ALCOI - VERGE DELS LLIRIS",
                                                        fine_t_points = seq(0.5, 365.5, 1),
                                                        smoothed_curve = data_plot$Value[data_plot$Station == "ALCOI - VERGE DELS LLIRIS"],
                                                        num_station = 1, Type = "Raw"))

p2 <- ggplot(aes(x = fine_t_points, y = smoothed_curve, lty = Type, col = Type), 
             data = df_smoothed_filtb) + geom_line() + theme_bw() +
        xlab(NULL) + ylab(expression("Ozone level" ~ (mu*g/m^3))) +
        scale_y_continuous(limits = c(0, 140)) +
        scale_color_manual(values = c("gray60", "black")) +
        theme(text = element_text(size = 18),
              plot.title = element_text(face = "bold"),
              legend.position = "none")

p <- grid.arrange(p2, p1, nrow = 1)
# ggsave(plot = p, filename = paste0(file.path("figures", "Figure2.jpg")),
#        width = 14, height = 6, dpi = 1200)

# Spline coefficients
num_basis <- 50
df_coefs_filt <- df_coefs[df_coefs$pollutant == "O3", ]
df_coefs_filt$NumFun <- rep(1:num_basis, length(unique(df_coefs_filt$station)))
stations <- unique(df_coefs_filt$station)
aux <- 1:length(stations)
names(aux) <- stations
df_coefs_filt$num_station <- aux[df_coefs_filt$station]
ggplot(aes(x = num_station, y = NumFun, fill = spline_coefficients), data = df_coefs_filt) +
  geom_tile() + coord_equal() + theme_bw() +
  xlab(expression("Station ("*italic("i")*")")) +
  ylab(expression("Spline basis function ("*italic("k")*")")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_distiller(palette = "Spectral", name = expression(italic(b)[italic(i)*italic(k)])) +
  theme(text=element_text(size = 18), plot.title = element_text(face = "bold"),
        legend.position = "right")

# ggsave(paste0(file.path("figures", "Figure3.jpg")), 
#        width = 7, height = 6, dpi = 1200)
