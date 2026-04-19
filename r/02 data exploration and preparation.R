### Required packages ###

# install.packages("pacman")
pacman::p_load(ggplot2, install = FALSE)

### Data loading ###

data <- read.csv(file = file.path("data", "contaminacion-atmosferica-y-ozono-promedios-diarios_202401.csv"), sep = ";")
head(data)
length(unique(data$COD_ESTACION))
length(unique(data$NOM_ESTACION))
# plot.ts(data$NO2[data$COD_ESTACION == "46250047"])

files <- dir("data")[2:13]
data <- c()
for (file in files) {
  print(file)
  aux <- read.csv(file.path("data", file), sep = ";")
  data <- rbind(data, aux)
}

### Transform to numeric ###

for (Var in 4:ncol(data)) {
  data[, Var] <- as.numeric(gsub(",", "\\.", data[, Var]))
}
# save(data, file = file.path("outputs", "data.rda"))

### Missing data proportion for several pollutants ###

# pollutants <- c("NO2", "SO2", "O3")
pollutants <- colnames(data)[4:ncol(data)]
df_missing <- c()
for (pollutant in pollutants) {
  for (station in unique(data$NOM_ESTACION)) {
    df_missing <- rbind(df_missing, data.frame(pollutant,
                                               code = unique(data$COD_ESTACION[data$NOM_ESTACION == station]),
                                               station,
                                               Prop = sum(is.na(data[which(data$NOM_ESTACION == station), pollutant]))/366))
  }
}
# plot.ts(data$NO2[data$NOM_ESTACION == "BENICASSIM"])
df_missing$Prop_10 <- as.numeric(df_missing$Prop <= 0.10)
df_missing$Prop_20 <- as.numeric(df_missing$Prop <= 0.20)
table(df_missing$pollutant, df_missing$Prop_10)
table(df_missing$pollutant, df_missing$Prop_20)
# save(df_missing,  file = file.path("outputs", "df_missing.rda"))

### Figure 1 ###

data_plot <- c()
sum_station <- 0
for (station in unique(data$NOM_ESTACION)) {
  print(station)
  if (df_missing$Prop[df_missing$station == station & df_missing$pollutant == "O3"] < 0.12) {
    sum_station <- sum_station + 1
    for (date in as.character(seq.Date(as.Date("2024-01-01"), as.Date("2024-12-31"), 1))) {
      data_plot <- rbind(data_plot, data.frame(Num_station = sum_station,
                                               Station = station,
                                               Day = date,
                                               Value = data$O3[data$NOM_ESTACION == station & data$FECHA == date]))
    }
  }
}
# save(data_plot, file = file.path("outputs", "data_plot.rda"))

Sys.setlocale("LC_ALL", "English")
ggplot(aes(x = as.Date(Day), y = Value, group = Station, col = Num_station),
       data = data_plot) + geom_line() + theme_bw() + xlab(NULL) +
  ylab(expression("Ozone level" ~ (mu*g/m^3))) +
  scale_x_continuous(breaks = c(as.Date("2024-01-01"), as.Date("2024-03-01"),
                                as.Date("2024-05-01"), as.Date("2024-07-01"),
                                as.Date("2024-09-01"), as.Date("2024-11-01"))) +
  scale_y_continuous(limits = c(0, 140)) +
  theme(text = element_text(size = 18),
        plot.title = element_text(face = "bold"),
        legend.position = "none")

# ggsave(paste0(file.path("figures", "Figure1.jpg")), 
#        width = 10, height = 6, dpi = 1200)
