# Cluster Analysis - simplying using broom and map

packages <- c('knitr','tidyverse','readxl','stringr','factoextra','doParallel','dbscan','mclust','extrafont', 'broom')

# Install packages if they don't exist
packages <- lapply(packages, FUN = function(x) {
  if(!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)}})

set_urban_defaults(style = "print")


# Run DBSCAN ----
# 1. Set up grid and function that reads cluster results into the grid
# 2. Add distance measure (silhouette score) into grid
# 3. Plot silhouette score

# read in weighted data
wgtdata <- read_csv('./Data/5. DatasetForcluster.csv') %>% 
  mutate(fips = as.character(fips)) %>% 
  mutate(fips = str_pad(fips, 5, pad = "0"))

wgtdata_clean <- wgtdata %>% 
  select_if(is.numeric) %>% 
  filter(complete.cases)

# Kmeans
# See here for helpful kmeans code: https://uc-r.github.io/kmeans_clustering#kmeans
fit_kmeans <- function(data, centers, nstart, ...) {
  # drop non-numeric vectors
  input_data <- select_if(.tbl = data, is.numeric) %>%  
  # only keeps complete observations - this shouldn't drop anything since we imputed
  filter(complete.cases(data))
  # fits kmeans model with specficied data, centers, and nstart
  kmeans(x = input_data, centers = centers, nstart = nstart)
}


# test the function with a single specification
fit_kmeans(data = wgtdata, centers = 3, nstart = 1) %>% 
  glance()

# create a grid of parameters for the clustering model
# notice the nested data frame
tuning_grid <- expand.grid(data = list(wgtdata), nstart = 1, centers = 2:20) %>%
  mutate(model_number = row_number())

model_output <- tuning_grid %>%
  # iterate fit_kmeans over every row in the tuning grid
  mutate(fits = pmap(list(data, centers, nstart), fit_kmeans)) %>%
  # extract tot.withinss
  mutate(tot.withinss = map_dbl(fits, ~glance(.) %>% pull(tot.withinss)))

#View(model_output)

# Plot tot.withinss over number of clusters
ggplot(data = model_output) +
  geom_path(mapping = aes(x = centers, y = tot.withinss))



# Optimal number of clusters ----
# Total within sum of squares
fviz_nbclust(wgtdata_clean, kmeans, method = "wss")
# Silhouette score
fviz_nbclust(wgtdata_clean, kmeans, method = "silhouette")








