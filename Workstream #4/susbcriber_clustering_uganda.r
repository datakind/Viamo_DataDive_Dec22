#clustering Uganda data 
#by Jahnavi G

#packages----
library(arrow) #for reading parquet
library(ggplot2)
library(PerformanceAnalytics) #visualize corr matrix
library(KernSmooth) #for pairwise density estimation
library(ResourceSelection) #for kdepairs
library(purrr) #for map function to repeat calculation for k times - useful for elbow/scree plot, silhouette model
library(cluster) #pam model for silhouette analysis
library(dplyr)
library(psych)
library(factoextra) #create beautiful graphs of clusters

#data----
raw_data <- read_parquet("Viamo_user_data_Uganda.parquet")
data <- subset(raw_data, select = c(n_calls, 
                                    n_themes,
                                    n_topics,
                                    median_call_duration,
                                    gender_numeric,
                                    age_numeric,
                                    fav_theme_numeric, 
                                    second_fav_theme_numeric,
                                    fav_topic_numeric,
                                    second_fav_topic_numeric,
                                    subscriber_id))

data.scaled <- scale(data[1:10]) #scaling all data save subscriber id

#Distribution checks----
summary(data)
ggplot(data, aes(n_calls)) + 
  geom_histogram(binwidth = 3) + scale_x_continuous(limits = c(0,100))
ggplot(data, aes(n_themes)) + geom_histogram() 

#Correlation checks----
#taking a sample of n subscriber data
data.sample <- data[data$subscriber_id %in% 
                      sample(x = data$subscriber_id, size = 100000, replace = FALSE), ]
data.sample.60k <- data[data$subscriber_id %in% 
                      sample(x = data$subscriber_id, size = 60000, replace = FALSE), ]
data.sample.30k <- data[data$subscriber_id %in% 
                          sample(x = data$subscriber_id, size = 30000, replace = FALSE), ]
data.sample.10k <- data[data$subscriber_id %in% 
                          sample(x = data$subscriber_id, size = 10000, replace = FALSE), ]

#visualizing corr matrix
chart.Correlation(scale(data.sample[1:10]), 
                  method = "pearson", 
                  histogram = TRUE, 
                  pch = 20) #solid circle

#Density checks----
pairwise_df <- data[,c(1,2)]
df_density <- bkde2D(pairwise_df, sapply(pairwise_df, dpik))
plot(pairwise_df, pch = 19)
contour(x = df_density$x1, y = df_density$x2, z = df_density$fhat, add = TRUE)

data.for.density <- data.sample[1:10]
data.for.density[data.for.density < 0] <- 0 #assigning 0 to negatives for density pairs

data.for.denisty.log <- log10(data.for.density)
data.for.denisty.log[is.na(data.for.denisty.log)] <- 0
data.for.denisty.log[is.infinite(data.for.denisty.log)] <- 0

kdepairs(x = data.for.denisty.log, 
         density=TRUE, 
         contour=TRUE)

#selecting relevant fields
#after all the distribution, corr, density checks
data.scaled.log <- log10(data.scaled)
# data.scaled.pruned <- subset(data.scaled, drop = x)

#applying K-means clustering-----------------------------------------------------------------
#converting into matrix for processing through k-means clustering functions.
matrix_scaled <- as.matrix(data.scaled)

##Elbow scree plot----
#calculate total-within cluster sum of squares
tot_withinss <- map_dbl(1:10, function(k){
  model <- kmeans(matrix_scaled, centers = k)
  model$tot.withinss
})
#append tot_withinss vector to correspond to its k value
elbow_df <- data.frame(k = 1:10, tot_withinss = tot_withinss)
#elbow plot
elbow_plot <- ggplot(elbow_df, aes(k, tot_withinss)) + 
  geom_line() +
  scale_x_continuous(breaks = 1:10)
#RESULT: elbow plot shows an optimal k of 2, 3 or 6 <- let's try for 6 clusters

##Silhouette analysis----
#checking validity of results from elbow plot against those from silhouette analysis
#calculating avg silhouette width across multiple values of k
data.sample.60k.scaled <- scale(data.sample.60k[,1:10])
matrix.sample.60k.scaled <- as.matrix(data.sample.60k.scaled)

data.sample.30k.scaled <- scale(data.sample.30k[,1:10])
matrix.sample.30k.scaled <- as.matrix(data.sample.30k.scaled)

data.sample.10k.scaled <- scale(data.sample.10k[,1:10])
matrix.sample.10k.scaled <- as.matrix(data.sample.10k.scaled)

start_time <- Sys.time()
sil_width <- map_dbl(2:10, function(k){ ##!!CAUTION: 2:10 took over 1.5 hrs to run!!
  pam_model <- pam(matrix.sample.30k.scaled, k = k)
  pam_model$silinfo$avg.width
})
end_time <- Sys.time()
end_time - start_time

#append sil_width vector to correspond to its k value
sil_df <- data.frame(k = 2:10, sil_width = sil_width)
#visualize diff btw k & avg silhouette width
sil_plot <- ggplot(sil_df, aes(k, sil_width)) +
  geom_line() +
  scale_x_continuous(breaks = 2:10)

#using function from factoextra package to calculate & visualize sil_width
start_time <- Sys.time()
fviz_nbclust(matrix.sample.10k.scaled, #reducing the sample from 30 to 10k massively reduced the time.
             FUNcluster = pam, 
             method = "silhouette",
             k.max = 10)
end_time <- Sys.time()
end_time - start_time
#RESULT: sil_plot shows an optimal k of 2 clusters - interesting. 

##k-means clustering----
### k=2----
model <- kmeans(matrix.sample.10k.scaled, centers = 6)
#appending assigned cluster back to original df
df_clustered_k2 <- mutate(data.sample.10k, cluster = model$cluster)
#calculating count of users assigned per cluster
count(df_clustered_k2, cluster) 

#plotting clusters - bivariate cluster plot
par(mfrow = c(1,2))
clusplot(df_clustered_k2[1:10], df_clustered_k2$cluster,
         color = TRUE, labels = 1, shade = TRUE, lines = 0, pch = 20)
clusplot(df_clustered_k6[1:10], df_clustered_k6$cluster,
         color = TRUE, labels = 1, shade = TRUE, lines = 0, pch = 20)

#explaining the clusters
numeric_fields <- df_clustered_k2[, c(1:10, 12)]
summary_cluster <- describeBy(numeric_fields, group = numeric_fields$cluster)
summary_cluster <- data.frame(t(sapply(summary_cluster,c)))

summary(df_clustered_k2[df_clustered_k2$cluster == 2, ])
df_clustered_k2[, c(2:10, 12)] %>% 
  group_by(cluster) %>% 
  summarise(median_calc = quantile(n_calls, c(0.25, 0.75)), prob = c(0.25, 0.75))

#comparing the distribution btw clusters
ggplot(df_clustered_k2, aes(n_calls)) +
  geom_histogram() +
  facet_wrap(~cluster)

#field -- cluster 1 --    cluster 2
#n------- 3k callers --   7k callers
#n_calls--median=15 calls--median = 2 calls
#n_themse--median=7 --    median=2 
#n_topics--median=3 --    median = 1
#median_call_duration-med=12--med=9 
#gender_numeric--med=males-med = NAs (others)
#age_numeric- med=18-24     med=NAs
#fav_theme_numeric-med=NA   med=NA
#second_fav_theme_numeric-  med=2(health) med=1(weather)
#fav_topic_numeric-med=NA   med=NA
#second_fav_topic_numeric-med=2-srh med=NA

#using vioplot to share differences. 
library(vioplot)
with(df_clustered_k2, 
     vioplot(n_calls[cluster==1], n_calls[cluster==2], 
             names = c("1", "2"),
             col = NA,
             main = "Number of Calls",
             pchMed = 16, colMed = "black", rectCol = "white",
             areaEqual = TRUE, 
             cex = 1, cex.main = 1.5,
             plotCentre = "points"))

median(subset(df_clustered_k2$n_calls, df_clustered_k2$cluster == 1))

### k=6----
df_clustered_k6 <- mutate(data.sample.10k, cluster = model$cluster)
df_clustered_k6$cluster <- as.numeric(df_clustered_k6$cluster)
#calculating count of users assigned per cluster
count(df_clustered_k6, cluster) 

clusplot(df_clustered_k6[1:10], df_clustered_k6$cluster,
         color = TRUE, labels = 1, shade = TRUE, lines = 0, pch = 20)

#explaining the clusters
ggplot(df_clustered_k6, aes(n_calls)) +
  geom_histogram() +
  facet_wrap(~cluster)
names(df_clustered_k6)

summary(df_clustered_k6[df_clustered_k6$cluster == 6, ])

