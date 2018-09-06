# Tidy Cluster Analysis - simplifying using broom and map

packages <- c('knitr','tidyverse','readxl','stringr','factoextra','doParallel','dbscan','mclust','extrafont', 'broom', 'fpc', 'NbClust')

# Install packages if they don't exist
packages <- lapply(packages, FUN = function(x) {
  if(!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)}})

# devtools::install_github("hadley/multidplyr")

library(urbnthemes)
set_urban_defaults(style = "print")

# Steps
# 1. Weight the input variables
# 2. Run cluster analysis
# 3. Combine analysis with ensemble methodology - We aren't sure if we'll do this step, but included anyway
# 4. Select the optimal number of groups
# 5. Validate cluster stability
# 6. Assign a "central" score

# 1. Weight the input variables ----
# Analyze Variable Correlations
site17 <- suppressMessages(suppressWarnings(read_csv('./Data/2. Unweighted Data.csv')))

site_num0 <- site17 %>%
  select_if( is.numeric)

site_num <- site_num0 %>%
  select( -fips)

## Use this to dump the cor.prob output to a 4 column matrix
## with row/column indices, correlation, and p-value.
## See StackOverflow question: http://goo.gl/fCUcQ
#https://www.r-bloggers.com/more-on-exploring-correlations-in-r/
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(rown = rownames(m)[row(m)[ut]],
             coln = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])
}


# Obtain correlations and "flatten" the table

cor_33 <- cor(as.matrix(site_num))
cor_333 <- flattenSquareMatrix(cor_33)
cor_4 <- cor_333 %>%
  mutate(abs_cor =abs(cor))

#writing it out as a CSV
# Write CSV in R
write.csv(cor_4, file = './Data/3. Correlations.csv', row.names=FALSE)

# Eliminate highly correlated vars
cor_4 <- suppressMessages(suppressWarnings(read_csv('./Data/3. Correlations.csv')))

#high correlation pairs
high.cor <- cor_4 %>%
  select(rown, coln, cor, abs_cor) %>%
  filter(abs_cor  > 0.9)

#non-high correlation pairs
non.high.cor <- cor_4 %>%
  select(rown, coln, cor, abs_cor) %>%
  filter(abs_cor< 0.9)
#View(non.high.cor)

#creating a list of highly correlated columns
#The rown values are to be eliminated
high.cor.var.r <- high.cor %>%
  mutate(rown = as.character(rown)) %>%
  pull(rown)

high.cor.to.drop <- high.cor.var.r

#loop through the non high correlation data frame and delete any pair that contains a value in high.cor.to.drop
pairs.with.high.var.rown <- non.high.cor[non.high.cor$rown %in% high.cor.to.drop,]
pairs.with.high.var.coln <- non.high.cor[non.high.cor$coln %in% high.cor.to.drop,]

pairs.with.high.var <- rbind(pairs.with.high.var.rown, pairs.with.high.var.coln)

pairs.under.consideration <-  anti_join(non.high.cor, pairs.with.high.var, by = c("rown","coln","cor","abs_cor"))

# Duplicate Variables and Down-Weight
# Vars 0.75 and above
cor.75 <- pairs.under.consideration %>%
  mutate(rown = as.character(rown), coln = as.character(coln) ) %>%
  filter(abs_cor  >=  0.75)

#View(cor.75)
#we subtract 2 from column rown and 1 from column coln
#rown processing
#t1 is a data frame that contains value -2 for each row. 
t1 <- data.frame("reps" = rep.int(-2, nrow(cor.75)))
analy.var <- data.frame( "analysis.var" = cor.75$rown) 
d1 <-  bind_cols(analy.var, t1)

#coln
#t2 is a data frame that contains value -1 for each row. 
t2 <- data.frame("reps" = rep.int(-1, nrow(cor.75)))
analy.var.1 <- data.frame( "analysis.var" = cor.75$coln) 
d2 <-  bind_cols(analy.var.1, t2)
#View(d2)

#appending the row and col
wt.1 <- suppressWarnings(suppressMessages(bind_rows(d2, d1)))

# 0.5 <= x < 0.75 - eliminate 1 copy from each variable in the pair
cor.50 <- pairs.under.consideration %>%
  filter(abs_cor>= 0.5 & abs_cor< 0.75)

nrep.2 <- nrow(cor.50)
t3 <- data.frame("reps" = rep.int(-1, nrep.2*2))

analy.var.2 <- data.frame( "analysis.var" = cor.50$coln)
analy.var.3 <- data.frame( "analysis.var" = cor.50$rown)
analy.var.23 <- suppressWarnings(suppressMessages(bind_rows(analy.var.2, analy.var.3)))
wt.2 <- bind_cols(analy.var.23, t3)

# 0.25 <= x < 0.5 - eilminate 1 copy from one variable in the pair
cor.25 <- pairs.under.consideration %>%
  filter(abs_cor< 0.50 & abs_cor >= 0.25)

#eliminate 1 from the coln
t4 <- data.frame("reps" = rep.int(-1, nrow(cor.25)))
analy.var.4 <- data.frame( "analysis.var" = cor.25$coln)
analy.44 <- cbind(analy.var.4, t4)

#eliminate 0 from rown
t5 <- data.frame("reps" = rep.int(0, nrow(cor.25)))
analy.var.5 <- data.frame( "analysis.var" = cor.25$rown)
analy.55 <- cbind(analy.var.5, t5)

wt.3 <- suppressWarnings(suppressMessages(bind_rows(analy.44, analy.55)))
#View(wt.3)

#putting the 3 weight data frames together
#aggregating the number of reps we need to subtract out

wt.all.0 <- suppressWarnings(suppressMessages(bind_rows(wt.1, wt.2, wt.3)))
#View(wt.all.0)
#View(wt.1)

wt.all.agg <- aggregate(reps~analysis.var, wt.all.0 ,sum)
#View(wt.all.agg)

# Up-weight key variables
# We "up-weight" certain variables considered high-priority for their research agenda by creating additional copies of those variables in the dataset.
# We create n x (the number of existing copies) in the dataset, where n is determined according to the list below. 
# For example, if after the previous step a variable has three copies, and n = 2 (our pick), it will have 6 copies after this step.
# HH: We pick n. Number of rows with suggested weight in comment following var
imp <- data.frame("analysis.var" =c("fdinsec",
                                    "fdinsec_cd",
                                    "diabetic",
                                    "r_lbw", 
                                    "fpl200",
                                    "snap", 
                                    "shcb",
                                    "ruralpop_pct" 
),
"mult.fact"=c(6,3,2,3,3,3,4,4))

#View(imp)
wt.all.1 <- merge(x = imp, y = wt.all.agg, by = "analysis.var", all = TRUE)
#View(wt.all.1)

#Putting together all the weight data frames to get the final weights
site17 <- suppressMessages(suppressWarnings(read_csv('./Data/2. Unweighted Data.csv')))
# site17 <- site17 %>% select (-MemberID, -MemberName, -fips)
site17 <- site17 %>% select ( -fips)

#adding column orig.wt that has value 4. All columns have an initial weight of 4
t6 <- data.frame("orig.wt" = rep.int(4, nrow(wt.all.1)))
wt.all.1.2 <- cbind(wt.all.1, t6)
#View(wt.all.1)

#determining the weight for each col
#replacing NA's with 0. 
wt.all.1.2[c("reps")][is.na(wt.all.1.2[c("reps")])] <- 0
wt.all.2 <- wt.all.1.2 %>%
  mutate(Diff.origwt.reps = reps+orig.wt)

#for all the weights that are negative we change final weight to 1
wt.all.3 <- wt.all.2 %>%
  mutate(final.wt0 = ifelse(Diff.origwt.reps < 1 , 1 ,Diff.origwt.reps)) 

#multiplying the mult.fact
#replacing the NA's in final.wt and mult.fact by 1
wt.all.3[c("mult.fact")][is.na(wt.all.3[c("mult.fact")])] <- 1

wt.all.4 <- wt.all.3 %>%
  mutate(final.wt=final.wt0 * mult.fact)

#adding in the variables that are not in the wt.all.4
col.pre.corr <- colnames(site17)
wt.all.col.names <- unique(wt.all.4$analysis.var)

vals.to.add <-data.frame("analysis.var"= col.pre.corr[!(col.pre.corr %in% wt.all.col.names)])
t7 <- data.frame("final.wt0" = rep.int(NA, nrow(vals.to.add)),
                 "mult.fact" =NA,
                 "reps" =NA,
                 "orig.wt"=NA,
                 "final.wt" =4,
                 "Diff.origwt.reps" =NA
)

vals.to.add2 <-  bind_cols(vals.to.add, t7)

#Member name gets a weight of 1
vals.to.add3 <- vals.to.add2 %>%
  mutate(final.wt =ifelse(analysis.var =="county", 1, final.wt))
#View(vals.to.add3)

#for the variables with high correlation that were eliminated at the start we assign a weight of 1
vals.to.add4 <- vals.to.add3 %>%
  mutate(final.wt =ifelse(analysis.var %in% high.cor.to.drop, 1, final.wt))
#View(vals.to.add4)

#---
wt.all.5 <- rbind(vals.to.add4, wt.all.4)
#View(wt.all.5)

#write out the weights file
write.csv(wt.all.5 %>% arrange(desc(final.wt)), file = './Data/4. Weights.csv', row.names = FALSE)

# Replicating columns based on the weight file
col.replicates <- data.frame(row.names = 1:nrow(site17))
for(i in 1:nrow(wt.all.5)) {
  rowlop <- wt.all.5 %>%
    slice(i)
  
  var.yui <- rowlop %>%
    mutate(analysis.var = as.character(analysis.var)) %>%
    pull(analysis.var)
  
  n <- rowlop %>%
    pull(final.wt)
  
  col.replicates <- cbind(col.replicates, as.data.frame(replicate(n, site17[var.yui])))
}
#View(col.replicates)
#colnames(col.replicates)

# writing it out as a CSV
site17.1 <- suppressMessages(suppressWarnings(read_csv('./Data/2. Unweighted Data.csv')))
col.replicates <- col.replicates %>% mutate(fips = pull(site17.1 %>% select(fips)), county = pull(site17.1 %>% select(county)))

# HH line below:
# col.replicates <- col.replicates %>% mutate(fips = pull(site17.1%>% select(fips)))
write.csv(col.replicates, file = "./Data/5. DatasetForcluster.csv", row.names=FALSE)
glimpse(col.replicates)
#mutate(group = groups[[num_groups - 1]], fips = pull(tc %>% select(fips)), County = pull(tc %>% select(County)), MemberID = pull(tc %>% select(MemberID))) %>% arrange(group)






# read in weighted data
wgtdata <- read_csv('./Data/5. DatasetForcluster.csv') %>% 
  mutate(fips = as.character(fips)) %>% 
  mutate(fips = str_pad(fips, 5, pad = "0"))

wgtdata_clean <- select_if(.tbl = wgtdata, is.numeric) %>% 
  filter(complete.cases(wgtdata))

# read in unweighted data
unwgtdata <- read_csv('./Data/2. Unweighted data.csv') %>% 
  mutate(fips = as.character(fips)) %>% 
  mutate(fips = str_pad(fips, 5, pad = "0")) %>% 
  subset(select = -county)

unwgtdata_clean <- select_if(.tbl = unwgtdata, is.numeric) %>% 
  filter(complete.cases(wgtdata))


# DBSCAN ----
# See here for helpful DBSCAN documentation: http://www.sthda.com/english/wiki/wiki.php?id_contents=7940

# 1. Set up grid and function that reads cluster results into the grid
# 2. Add distance measure (silhouette score) into grid
# 3. Plot silhouette score

# Set up DBSCAN to look at range of epsilons and minimum points ----

set.seed(665544)
# Set up function that runs the dbscan model
fit_db <- function(data, eps, MinPts, ...) {
  # drop non-numeric vectors
  input_data <- select_if(.tbl = data, is.numeric) %>%  
    # only keeps complete observations - this shouldn't drop anything since we imputed
    filter(complete.cases(data))
  # fits dbscan model with specficied data, centers, and nstart
  res <- dbscan(input_data, eps, MinPts)
  return(unlist(max(res$cluster)))
  #fit <- list(eps = res$eps, clusters = max(res$cluster))
  #res$clusters <- max(res$cluster)
  #return(fit$clusters)
}

# test the function with a single specification
test <- fit_db(data = wgtdata, eps = 100, MinPts = 5)


tuning_grid <- expand.grid(data = list(wgtdata), eps = 1:60, MinPts = 4:5) %>%
  mutate(model_number = row_number())

model_output <- tuning_grid %>%
  # iterate fit_db over every row in the tuning grid
  mutate(clusters = pmap(list(data, eps, MinPts), fit_db)) %>% 
  mutate(clusters = unlist(clusters))

ggplot(model_output, mapping = aes(eps, clusters)) + 
  geom_line() + 
  facet_wrap(~ MinPts) +
  expand_limits(y=0) + 
  theme_light() +
  labs(title = "Epsilon vs. Number of clusters",
       subtitle = "Faceted by minimum number of points",
       caption = "Urban Institute",
       x = "Epsilon",
       y = "Number of groups")
#ggsave("./Plots/eps_facet.png", width = 10, height = 5)


# Run DBSCAN and get into proper format for silhouette and stability analysis ----
ncores <- detectCores() - 2
cl <- makeCluster(ncores)
registerDoParallel(cl)

num_results <- 80
results <- foreach(pt = 4:5, .combine = 'c') %:% 
  foreach(ep = 20:60) %dopar% {
    library(dbscan)
    ep1 <- ep
    dbscan(wgtdata_clean, eps = ep1, minPts = pt)
  }



# Filter out ones with less than 2 or more than 20, choose only the ones with the least noise
results_dbscan <- vector("list", 19)
noise_count <- numeric(19)
for (i in 1:19){ noise_count[[i]] <- 3500 }
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
silhouette_score <- foreach(i = 1:19) %dopar% {
  library(cluster)
  print(i)
  if (is.null(results_dbscan[[i]])) { score <- NA } else {
    s <- silhouette(results_dbscan[[i]],dist_m)
    score <- mean(s[,3])
    }
  score
}

n_groups <- 2:20
silhouette_score <- unlist(silhouette_score)
sil_score <- tibble(n_groups,silhouette_score)
#sil_score <- tibble(n_groups,score)

ggplot(sil_score, mapping = aes(x = n_groups, y = silhouette_score)) + 
  geom_line() + 
  expand_limits(y=0) + 
  labs(title = "Goodness of fit",
       subtitle = "Silhouette Score",
       caption = "Urban Institute",
       x = "Number of Groups",
       y = "Score")


num_groups = 12
#sub1 <- sub %>% mutate(group = groups[[num_groups - 1]], fips = pull(tc %>% select(fips)), County = pull(tc %>% select(County)), MemberID = pull(tc %>% select(MemberID)), MemberName = pull(tc %>% select(MemberName))) %>% arrange(group)
clusterdata <- wgtdata_clean %>%
  mutate(group = results_dbscan[[num_groups-1]], fips = pull(wgtdata %>% select(fips))) %>%
  arrange(group)

table(clusterdata$group)



# Map of cluster ----
library(urbnthemes)
library(urbnmapr)
library(stringr)
# set_urban_defaults(style = "map")
setnames(counties, 'county_fips', 'fips')

peergroups <- sub1 %>% 
  select(group, fips) %>% 
  mutate(fips = as.character(fips)) %>% 
  mutate(fips = str_pad(fips, 5, pad = "0")) %>% 
  rename(clustergroup = group) %>% 
  mutate(clustergroup = as.numeric(clustergroup))

peermap <- left_join(peergroups, counties, by = "fips")
glimpse(peermap)

map <- function(fill_var, label_nm) {
  peergroups %>% 
    fill_cut <- factor(fill_var[1:10]) %>% 
      ggplot(aes_string("long", "lat", group = "group", fill = fill_cut)) +
      geom_polygon(color = "#ffffff", size = 0.05) +
      coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
      scale_fill_discrete() +
      theme(legend.position = "right",
            legend.direction = "vertical",
            legend.title = element_text(face = "bold", size = 11),
            legend.key.height = unit(.2, "in")) +
      labs(fill = label_nm)
    ggsave(paste0("./Maps/", fill_var, ".png"), width = 12, height = 9, units = "in")
}
