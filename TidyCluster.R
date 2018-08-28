# Tidy Cluster Analysis - simplifying using broom and map

packages <- c('knitr','tidyverse','readxl','stringr','factoextra','doParallel','dbscan','mclust','extrafont', 'broom', 'fpc', 'NbClust')

# Install packages if they don't exist
packages <- lapply(packages, FUN = function(x) {
  if(!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)}})

devtools::install_github("hadley/multidplyr")

# set_urban_defaults(style = "print")


# read in weighted data
wgtdata <- read_csv('./Data/5. DatasetForcluster.csv') %>% 
  mutate(fips = as.character(fips)) %>% 
  mutate(fips = str_pad(fips, 5, pad = "0"))

wgtdata_clean <- select_if(.tbl = wgtdata, is.numeric) %>% 
  filter(complete.cases(wgtdata))

# read in unweighted data
unwgtdata <- read_csv('./Data/2. Unweighted data.csv') %>% 
  mutate(fips = as.character(fips)) %>% 
  mutate(fips = str_pad(fips, 5, pad = "0"))

unwgtdata_clean <- select_if(.tbl = unwgtdata, is.numeric) %>% 
  filter(complete.cases(wgtdata))

# K-means ----

# Kmeans
# See here for helpful kmeans code: https://uc-r.github.io/kmeans_clustering#kmeans

set.seed(123)
# Set up function that runs the kmeans model
fit_kmeans <- function(data, centers, nstart, ...) {
  # drop non-numeric vectors
  input_data <- select_if(.tbl = data, is.numeric) %>%  
  # only keeps complete observations - this shouldn't drop anything since we imputed
  filter(complete.cases(data))
  # fits kmeans model with specficied data, centers, and nstart
  kmeans(x = input_data, centers = centers, nstart = nstart)
}

# test the function with a single specification
fit_kmeans(data = wgtdata, centers = 3, nstart = 5) %>% 
  glance()

# create a grid of parameters for the clustering model
# notice the nested data frame
tuning_grid <- expand.grid(data = list(wgtdata), nstart = 25, centers = 2:20) %>%
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


# Optimal number of clusters
# Total within sum of squares
fviz_nbclust(wgtdata_clean, kmeans, method = "wss", k.max = 20)
# Silhouette score
fviz_nbclust(wgtdata_clean, kmeans, method = "silhouette", k.max = 20)



# Viewing clusters using the fviz_cluster package

# 3 clusters on Kmeans
km.res <- kmeans(wgtdata_clean, 3, nstart = 25)
fviz_cluster(km.res, wgtdata_clean, frame = FALSE, geom = "point")

# 6 clusters on Kmeans
km.res <- kmeans(wgtdata_clean, 6, nstart = 25)
fviz_cluster(km.res, wgtdata_clean, frame = FALSE, geom = "point")

# 12 clusters on kmeans
km.res <- kmeans(wgtdata_clean, 10, nstart = 25)
fviz_cluster(km.res, wgtdata_clean, frame = FALSE, geom = "point")

#Hierarchical
# 3 clusters on Kmeans
hc.cut <- hcut(wgtdata_clean, 3, hc_method = "complete")
fviz_cluster(hc.cut, wgtdata_clean, frame = FALSE, geom = "point")

# 6 clusters on Kmeans
hc.cut <- hcut(wgtdata_clean, 6, hc_method = "complete")
fviz_cluster(hc.cut, wgtdata_clean, frame = FALSE, geom = "point")

# 12 clusters on kmeans
hc.cut <- hcut(wgtdata_clean, 12, hc_method = "complete")
fviz_cluster(hc.cut, wgtdata_clean, frame = FALSE, geom = "point")



# DBSCAN ----
# See here for helpful DBSCAN documentation: http://www.sthda.com/english/wiki/wiki.php?id_contents=7940

# 1. Set up grid and function that reads cluster results into the grid
# 2. Add distance measure (silhouette score) into grid
# 3. Plot silhouette score

# From Graham and Sybil ----
# Goal: get this to output tidy dfs
ncores <- detectCores() - 2
cl <- makeCluster(ncores)
registerDoParallel(cl)

results <- foreach(pt = 2:10, .combine = 'c') %:% 
  foreach(ep = 1:75) %dopar% {
    library(dbscan)
    ep1 <- ep + 75
    dbscan(wgtdata_clean, eps = ep1, minPts = pt)
  }

# Filter out ones with less than 2 or more than 20, choose only the ones with the least noise
num_results <- 400
results_dbscan <- vector("list", 19)
noise_count <- numeric(19)
for (i in 1:19){ noise_count[[i]] <- 201 }
for (i in 1:num_results){
  num_clusters <- max(results[[i]]$cluster) + 1
  if (num_clusters < 21 && num_clusters > 1){
    noise <- sum(results[[i]]$cluster == 0)
    if (noise < noise_count[[num_clusters - 1]]){
      noise_count[[num_clusters - 1]] <- noise
      results_dbscan[[num_clusters - 1]] <- results[[i]]$cluster
    }
  }
}

library(extrafont)
library(cluster)
dist_m <- dist(wgtdata_clean, method = "euclidean")
silhouette_score <- foreach (i = 1:19) %dopar% {
  library(cluster)
  #print(results_dbscan[[i]])
  s <- silhouette(results_dbscan[[i]],dist_m)
  suppressWarnings(if (is.na(s)){ score <- 0 } else { score <- mean(s[,3]) })
  score
}
n_groups <- 2:20
silhouette_score <- unlist(silhouette_score)
sil_score <- tibble(n_groups,silhouette_score)
result <- ggplot(sil_score, mapping = aes(x = n_groups, y = silhouette_score)) + 
  geom_line() + 
  expand_limits(y=0) + 
  labs(title = "Goodness of fit",
       subtitle = "Average Silhouette Score",
       caption = "Urban Institute",
       x = "Number of Groups",
       y = "Average Score")

# Tidy version ----

set.seed(123)
# Set up function that runs the kmeans model
fit_db <- function(data, eps, MinPts, ...) {
  # drop non-numeric vectors
  input_data <- select_if(.tbl = data, is.numeric) %>%  
    # only keeps complete observations - this shouldn't drop anything since we imputed
    filter(complete.cases(data))
  # fits kmeans model with specficied data, centers, and nstart
  dbscan(input_data, eps, MinPts)
}

# test the function with a single specification
fit_db(data = wgtdata, eps = 100, MinPts = 5)

tuning_grid <- expand.grid(data = list(wgtdata), eps = 95:105, MinPts = 2:5) %>%
  mutate(model_number = row_number())

model_output <- tuning_grid %>%
  # iterate fit_db over every row in the tuning grid
  mutate(fits = pmap(list(data, eps, MinPts), fit_db)) %>%
  # extract tot.withinss
  mutate(tot.withinss = map_dbl(fits, ~glance(.) %>% pull(tot.withinss)))



# Determining the optimal number of clusters
# With 6 clusters
dbscan::kNNdistplot(wgtdata_clean, k =  6)
abline(h = 0.15, lty = 2)

# With 12 clusters
dbscan::kNNdistplot(wgtdata_clean, k =  12)
abline(h = 0.15, lty = 2)


db100 <- dbscan(wgtdata_clean, eps=100, MinPts=3)
hullplot(wgtdata_clean, db100)
db100
str(db100)



# Looks like optimal epsilon is around 100

# Run DBSCAN cluster
set.seed(123)
# fpc package
res.fpc <- fpc::dbscan(wgtdata_clean, eps = 99, MinPts = 5)
# Visualize results of dbscan
fviz_cluster(res.fpc, wgtdata_clean, geom = "point")

# dbscan package
res.db <- dbscan::dbscan(wgtdata_clean, 100, 4)

# Make sure FPC and DBSCAN packages produce the same results
all(res.fpc$cluster == res.db) # ERROR






NbClust(wgtdata_clean, distance = "euclidean", 
        min.nc = 2, max.nc = 20, method = "complete", index = "all")

# Hierarchical clustering
# Compute pairewise distance matrices
dist.res <- dist(wgtdata_clean, method = "euclidean")
# Hierarchical clustering results
hc <- hclust(dist.res, method = "complete")
# Visualization of hclust
plot(hc, labels = FALSE, hang = -1)
# Add rectangle around 3 groups
rect.hclust(hc, k = 3, border = 2:4)
rect.hclust(hc, k = 6, border = 2:4)
rect.hclust(hc, k = 12, border = 2:4)

# Cut into 3 groups
hc.cut <- cutree(hc, k = 3)
head(hc.cut, 20)




# Misc. ----
# Use the scale function to scale our data for us?
scaleddata <- scale(wgtdata_clean, center = TRUE, scale = TRUE)
glimpse(scaleddata)


hc.cut <- hcut(wgtdata_clean, k = 3, hc_method = "complete")
# Visualize dendrogram
fviz_dend(hc.cut, show_labels = FALSE, rect = TRUE)
# Visualize silhouhette information
fviz_silhouette(hc.cut)
