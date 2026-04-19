### Required packages ###

# install.packages("pacman")
pacman::p_load(rvest, dplyr, stringr, install = FALSE)

### Extract coordinates ###

file_path <- file.path("data", "estacionesweb.txt")
html_content <- read_html(file_path)

stations <- html_nodes(html_content, xpath = "//*[@data-lat and @data-lng]")
coords_df <- data.frame(lat = as.numeric(html_attr(stations, "data-lat")),
                        lng = as.numeric(html_attr(stations, "data-lng")),
                        stringsAsFactors = FALSE)

### Extract names with code ###

titles <- html_nodes(html_content, "h2.location-title[property='name']") %>% 
  html_text(trim = TRUE)
titles_df <- data.frame(codigo = str_extract(titles, "^[0-9]+"),
                        nombre = str_trim(str_replace(titles, "^[0-9]+ - ", "")),
                        stringsAsFactors = FALSE)

### Combine ###

estaciones_df <- bind_cols(titles_df, coords_df)

### Remove duplicates ###

estaciones_df <- estaciones_df %>% distinct(codigo, nombre, .keep_all = TRUE)

### Save as CSV ###

# write.csv(estaciones_df, file.path("data", "estaciones_coords.csv"), 
#           row.names = FALSE)
