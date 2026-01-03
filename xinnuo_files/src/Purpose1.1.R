## ------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(glue)
library(tidyverse)
library(scatterplot3d)
library(synthpop)
library(factoextra)
library(parallel)
library(mvtnorm)
library(NbClust)
library(clue)
library(DescTools)
library(clusterGeneration) # Added from later chunk
library(MixSim)            # Added from later chunk
library(tidyverse)         # Added from later chunk (this provides dplyr and %>%)
library(patchwork)         # Added from later chunk

synt.data.generation <- function(real.data, m, methods){
    syn.obj <- syn(real.data,
                       method = rep(methods, ncol(real.data)),
                       m = m,
                       print.flag = FALSE, 
                       proper = TRUE)
    synt.data <- syn.obj$syn
    specks.val <- numeric()
    if(m==1){
      names(synt.data) <- names(real.data)
      data <- rbind(
            cbind(real.data, label = 0),
            cbind(synt.data, label = 1)
          )
      specks.val <- c(specks.val, specks.correct(data = data)$Specks)
      specks.value <- data.frame(SPECKS = specks.val)
    } else {
      for(k in 1:m){
          names(synt.data[[k]]) <- names(real.data)
          
          data <- rbind(
            cbind(real.data, label = 0),
            cbind(synt.data[[k]], label = 1)
          )
          
          specks.val <- c(specks.val, specks.correct(data = data)$Specks)
          specks.value <- data.frame(SPECKS = specks.val)
      }
    }
    results <- list("Specks" = specks.value, "Syn data" = synt.data)
    return(results)
}


label.switching <- function(km.real, km.synt){ # Label switching function
  if(!is.null(km.synt)){
      centroids <- rbind(km.real$centers, km.synt$centers)
      dist.global <- as.matrix(dist(centroids, method = "euclidean"))
      n.real <- nrow(km.real$centers) #numero de cluster en real
      n.synt <- nrow(km.synt$centers) #numero de cluster en synt
      dist <- dist.global[1:n.real, (n.real + 1) : (n.real + n.synt)]
      label.ass <- solve_LSAP(dist, maximum = FALSE) 
      
      # Synthetic data reetiquetation
      km.synt$recluster <- km.real$cluster # nova variable
      
      for (i in 1:length(label.ass)) {
        km.synt$recluster[km.synt$cluster == label.ass[i]] <- i
      }
  }
  return(km.synt)
}


## ------------------------------------------------------------------------------------------------------------------------------------
#### Manual generation ####
library(mvtnorm)

## For p = 2
mu1 <- c(160, 50)
mu2 <- c(170, 65)

sd <- abs(mu1 - mu2)/2    # 2 sd -> separation between clusters
Sigma <- diag(sd^2)

n  <- 200
clust1 <- rmvnorm(n, mean = mu1, sigma = Sigma)
clust2 <- rmvnorm(n, mean = mu2, sigma = Sigma)

data2 <- rbind(
  data.frame(var1 = clust1[,1], var2 = clust1[,2], group = 1),
  data.frame(var1 = clust2[,1], var2 = clust2[,2], group = 2)
)

ggplot(data2, aes(var1, var2, colour=factor(group))) +
  geom_point(alpha = 1, size = 2) +
  scale_colour_brewer(palette = "Set1", name = "Group") +
  theme_minimal()


## For p = 3
set.seed(123)

# Parameters
n <- 100
p <- 3
sigma <- 1              # standard deviation
dist.sep <- 2           # separation between centroids

# Spherical covariance matrix: σ² * I
sigma_mat <- diag(sigma^2, p)

# Define the centroids: separation = 2
# mu1 at (0,0,0), mu2 along a diagonal distance of 2
# Choose a uniform direction for simplicity
direction <- c(1, 1, 1)
direction <- direction / sqrt(sum(direction^2))  # normalize unit vector

mu1 <- c(0, 0, 0)
mu2 <- direction * dist.sep                      # 2SD distance in 3D

dist.centroids <- sqrt(sum((mu2 - mu1)^2))       # Check if distance=2
print(dist.centroids)

cluster1 <- rmvnorm(n, mean = mu1, sigma = sigma_mat)  # data generation
cluster2 <- rmvnorm(n, mean = mu2, sigma = sigma_mat)

data <- rbind(cluster1, cluster2)
labels <- factor(rep(1:2, each = n))

# 3D visualization
scatterplot3d(data, color = as.numeric(labels), pch = 19,
              main = "Two clusters separated by 2σ in 3D")


## ------------------------------------------------------------------------------------------------------------------------------------
# function 
centroid.generation <- function(k = k, p = p, separation = separation) {
  centroids <- matrix(NA, nrow = k, ncol = p)
  centroids[1, ] <- rep(0, p)                 # first centroid at the origin 
  for (i in 2:k) {
    tries <- 0
    repeat {
      vector <- rnorm(p)                      # a random point
      vector <- vector / sqrt(sum(vector^2))  # normalize the point vector
      vector <- vector * separation           # distance of number of separations from the origin
      
      dists <- sqrt(rowSums((centroids[1:(i - 1), , drop = FALSE] -  # distances from the new point to each of the fixed points
                             matrix(vector, nrow = i - 1, ncol = p, byrow = TRUE))^2))
      if (all(dists >= separation)) {         # accept the new point if this complete the separation
        centroids[i, ] <- vector
        break
      }
      tries <- tries + 1                      # avoid infinite loop
      if (tries > 1000) stop("It has not been possible to generate sufficiently separate centroids.")
    }
  }
  return(centroids)
}

# test
centroids <- centroid.generation(k = 4, p = 2, separation = 2)
dist(centroids)


## ------------------------------------------------------------------------------------------------------------------------------------
#### Real data generation function ####
data.generation <- function(N, p, rho, separation, k){
  # centroid position
  centroids <- centroid.generation(k = k, p = p, separation = separation)
  
  # correlation matrix
  R <- matrix(rho, nrow = p, ncol = p)
  diag(R) <- 1
  
  # covariance matrix
  sd <- rep(1, p)
  D <- diag(sd)
  Sigma <- D %*% R %*% D  #
  
  # data generation for each cluster
  n <- round(N/k, digits = 0)
  clusters <- lapply(1:k, function(i) {
    rmvnorm(n,
            mean  = centroids[i, ],
            sigma = Sigma)
  })
  data <- do.call(rbind, clusters)
  data <- data.frame(data, group = factor(rep(1:k, each = n)))
  
  # graphic representation
  Graphic <- ggplot(data, aes(x = X1, y = X2, color = group)) +
                geom_point(size = 2) +
                scale_color_brewer(palette = "Set1") +
                labs(
                  title = "First dimention representation",
                  color = "Cluster"
                ) +
                theme_minimal(base_size = 14)
  result <- list("Data" = data[,-ncol(data)], "Graphic"=Graphic, "Clusters" = clusters)
}


synt.data.generation <- function(real.data, m, methods){
    syn.obj <- syn(real.data,
                       method = rep(methods, ncol(real.data)),
                       m = m,
                       print.flag = FALSE, 
                       proper = TRUE)
    synt.data <- syn.obj$syn
    specks.val <- numeric()
    if(m==1){
      names(synt.data) <- names(real.data)
      data <- rbind(
            cbind(real.data, label = 0),
            cbind(synt.data, label = 1)
          )
      specks.val <- c(specks.val, specks.correct(data = data)$Specks)
      specks.value <- data.frame(SPECKS = specks.val)
    } else {
      for(k in 1:m){
          names(synt.data[[k]]) <- names(real.data)
          
          data <- rbind(
            cbind(real.data, label = 0),
            cbind(synt.data[[k]], label = 1)
          )
          
          specks.val <- c(specks.val, specks.correct(data = data)$Specks)
          specks.value <- data.frame(SPECKS = specks.val)
      }
    }
    results <- list("Specks" = specks.value, "Syn data" = synt.data)
    return(results)
}


# Test
data <- data.generation(N = 1000, p = 2, rho = 0, separation = 6, k = 2)
data$Data
data$Graphic
data$Clusters


## ------------------------------------------------------------------------------------------------------------------------------------
# Generation using R function
library(clusterGeneration)

## Cluster separation depends on the sepVal parameter, where:
  # -0.99 -> null separation
  # -0.32 -> 2 standard deviations
  #  0.01 -> 4 standard deviations
  #  0.21 -> 6 standard deviations
set.seed(12345)
sim <- genRandomClust(numClust      = 3,             # number of clusters
                      sepVal        = -0.32,         # separation between clusters
                      numNonNoisy   = 2,             # number of informative variables
                      numNoisy      = 0,             # number of noisy variables
                      clustszind    = 1,             # equal cluster sizes
                      clustSizeEq   = 25,            # size of each cluster
                      numReplicate  = 1)             # number of replicates

data3  <- as.data.frame(sim$datList[[1]])            # generated data matrix
data3$group  <- sim$memList[[1]]                     # cluster assignment vector for each observation

ggplot(data3, aes(x1, x2, colour=factor(group))) +   # graphical representation of the data
  geom_point(alpha = 1, size = 2) +
  scale_colour_brewer(palette = "Set1", name = "Group") +
  theme_minimal()

## correlation between variables:
# Depends on covMethod: "eigen" or "unifcorrmax"

# "eigen"
R <- matrix(0.7, nrow=2, ncol=2)
diag(R) <- rep(1, 2)

eigen.values <- eigen(R)$values

sim4 <- genRandomClust(numClust      = 2,            # number of clusters
                      sepVal        = -0.32,         # separation between clusters
                      numNonNoisy   = 2,             # number of informative variables
                      numNoisy      = 0,             # number of noisy variables
                      clustszind    = 1,             # equal cluster sizes
                      clustSizeEq   = 25,            # size of each cluster
                      numReplicate  = 1,             # number of replicates
                      rotateind = FALSE)
                      
data4  <- as.data.frame(sim4$datList[[1]])            # generated data matrix
data4$group  <- sim4$memList[[1]]                     # cluster assignment vector for each observation

ggplot(data4, aes(x1, x2, colour=factor(group))) +   # graphical representation of the data
  geom_point(alpha = 1, size = 2) +
  scale_colour_brewer(palette = "Set1", name = "Group") +
  theme_minimal()

cor(data4$x1, data4$x2)

# "unifcorrmax" – setting alphad value
rcorrmatrix(3, alphad = 0.35)

f <- function(a) gamma(a)^2 / gamma(2*a) - 0.7

# find the root: choose an interval where the signs are opposite
a <- uniroot(f, lower = 0.5, upper = 2)$root

(gamma(a)*gamma(a))/gamma(2*a)
d <- 2
alphad <- a - ((d-2)/2)
rcorrmatrix(d, alphad = alphad)


## ------------------------------------------------------------------------------------------------------------------------------------
library(MixSim)

## Cluster separation
  ## MaxOmega = 0.01 → at most 1% of the generated points can be misclassified. Well-separated clusters.
  ## ecc = 0 → no correlation between variables
set.seed(1234)
sd <- 6
omega <- 2 * pnorm(-(1/2) * sd)
Q <- MixSim(BarOmega = omega, K = 3, p = 3, hom = TRUE, sph = TRUE, ecc = 0.1)

# Build dataset using simdataset
sim <- simdataset(n = 50, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
data5 <- as.data.frame(sim$X)
data5$group <- sim$id

ggplot(data5, aes(V1, V2, colour = factor(group))) +                      # graphical representation of the data
  geom_point(alpha = 1, size = 2) +
  scale_colour_brewer(palette = "Set1", name = "Group") +
  theme_minimal()

cov(data5$V1, data5$V2)

data5 %>%
  group_by(group) %>%
  summarise(corr = cor(V1, V2))

# Data generation function
data.generation.test <- function(n, p, separation, k, rho = FALSE){
  set.seed(1234)
  omega <- 2 * pnorm(-(1/2) * separation)
  if (rho == FALSE){
    Q <- MixSim(MaxOmega = omega, K = k, p = p, hom = TRUE, sph = TRUE)   # spherical → no correlation
  } else {
    Q <- MixSim(MaxOmega = omega, K = k, p = p, hom = TRUE, sph = FALSE)  # non-spherical → correlation allowed
  }
  
  sim <- simdataset(n = n, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
  data <- as.data.frame(sim$X)
  data$group <- sim$id
  return(data)
}

# Test
set.seed(1234)
data <- data.generation.test(50, 2, rho = TRUE, separation = 6, k = 2)
ggplot(data, aes(V1, V2, colour = factor(group))) +
  geom_point(alpha = 1, size = 2) +
  scale_colour_brewer(palette = "Set1", name = "Group") +
  theme_minimal()

cor(data$V1, data$V2)

data %>%
  group_by(group) %>%
  summarise(corr = cor(V1, V2))


## ------------------------------------------------------------------------------------------------------------------------------------
# SPECKS function overfitting bias correction
specks.correct <- function(data, sample = FALSE){   # With p-value
  if(sample == TRUE){
    data[,3] <- sample(data[,3])
  }
  
  N <- nrow(data)
  idx <- sample(1:N, round(N * 0.3, digits = 0))    # 30% test split
  
  data.train <- data[-idx,]   # Training set
  data.test <- data[idx,]     # Test set
  
  logit_model <- glm(label ~ ., data = data.train, family = "binomial")      # Logistic regression
  
  pred.prob <- predict(logit_model, newdata = data.test, type = "response")  # Predict probabilities
  
  label.vector <- data.test$label                    # Extract true labels
  
  prop.score.real  <- pred.prob[label.vector == 0]  # Scores for real data
  prop.score.synt  <- pred.prob[label.vector == 1]  # Scores for synthetic data
  
  ecdf.real  <- ecdf(prop.score.real)
  ecdf.synth <- ecdf(prop.score.synt)
  
  all.scores <- sort(unique(c(prop.score.real, prop.score.synt)))
  
  specks.value <- max(abs(ecdf.real(all.scores) - ecdf.synth(all.scores)))        # Kolmogorov-Smirnov distance
  p.value <- suppressWarnings(ks.test(prop.score.real, prop.score.synt)$p.value)  # KS p-value
  
  return(data.frame(
    "Specks" = specks.value, 
    "P-value" = p.value))
}

# Test
my.seed <- 12345

real.data.1 <- data.generation(N = 500, p = 8, rho = 0.5, separation = 2, k = 2)$Data

sds.default.cart1 <- syn(real.data.1, method = "cart", seed = my.seed, m = 1, print.flag = FALSE)
synt.data.cart1 <- sds.default.cart1$syn

combined.data <- rbind(
  cbind(real.data.1, label = 0),    # Original data (label = 0)
  cbind(synt.data.cart1, label = 1) # Synthetic data (label = 1)
)

# Corrected SPECKS
specks.correct(combined.data)


## ------------------------------------------------------------------------------------------------------------------------------------
# Old function to determine the optimal number of clusters
k.decision <- function(data, k.real = 0, synthetic = TRUE){
  if(synthetic == FALSE){
    kopt.obj <- NbClust(data = data, diss = NULL, distance = "euclidean", 
                        method = "kmeans", index = "all", alphaBeale = 0.1)
    kopt <- as.numeric(names(table(as.vector(kopt.obj$Best.nc[1,])))[
              which.max(table(as.vector(kopt.obj$Best.nc[1,])))])

    sdata <- scale(data)
    clust <- kmeans(sdata, kopt, nstart = 25)
    
    clustering.result <- list("K optim" = kopt, "Clustering" = clust)
    return(clustering.result)
    
  } else {
    kopt.obj <- NbClust(data = data, diss = NULL, distance = "euclidean", 
                        method = "kmeans", index = "all", alphaBeale = 0.1)
    kopt <- as.numeric(names(table(as.vector(kopt.obj$Best.nc[1,])))[
              which.max(table(as.vector(kopt.obj$Best.nc[1,])))])

    if(kopt == k.real){
      sdata <- scale(data)
      clust <- kmeans(sdata, kopt, nstart = 25)
      
      clustering.result <- list("K optim" = kopt, "Clustering" = clust)
      return(clustering.result)
    }
  }
}

# Get optimal k using silhouette method
sil <- fviz_nbclust(real.data.1, kmeans, method = "silhouette", k.max = 10, nstart = 25)$data
k <- as.numeric(sil$clusters[which.max(sil$y)])

# Improved version using silhouette method
k.decision <- function(data, k.real = 0, synthetic = TRUE){ 
  if(synthetic == FALSE){
    sil <- fviz_nbclust(data, kmeans, method = "silhouette")$data
    kopt <- as.numeric(sil$clusters[which.max(sil$y)])
    
    sdata <- scale(data)
    clust <- kmeans(sdata, kopt, nstart = 25)
    
    clustering.result <- list("K optim" = kopt, "Clustering" = clust)
    return(clustering.result)
    
  } else {
    sil <- fviz_nbclust(data, kmeans, method = "silhouette")$data
    kopt <- as.numeric(sil$clusters[which.max(sil$y)])
    
    if(kopt == k.real){
      sdata <- scale(data)
      clust <- kmeans(sdata, kopt, nstart = 25)  # clustering result
      
      clustering.result <- list("K optim" = kopt, "Clustering" = clust)
      return(clustering.result)
    }
  }
}

# Test
k.decision(real.data.1, synthetic = FALSE)$`K optim`


## ------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

simulate_one <- function(factor.combination) # Computation function for each combination
{
  N <- factor.combination$N
  p <- factor.combination$p
  k <- factor.combination$k
  rho <- factor.combination$rho
  method <- factor.combination$method
  separation <- factor.combination$separation
  total.m <- total.m
  
  real.data <- data.generation(N, p, rho, separation, k)$Data
  
  kd        <- k.decision(real.data, synthetic = FALSE)
  k.real    <- kd$`K optim`
  km.real   <- kd$Clustering
  vars <- names(real.data)
  
  gini.real <- Gini(km.real$size)
  
  lreal.data <- real.data
  lreal.data$label <- as.factor(km.real$cluster)
  
  out <- data.frame()
  
  # Synthetic data generation
  count <- 0
  m <- 0
  synt.data <- list()
  km.synt <- list()
  k.synt <- numeric()
  speck <- numeric()
  while(m <= total.m){
    set.seed(m + count)
    count <- count + 1
    syn   <- synt.data.generation(real.data, m = 1, methods = method)
    clust <- k.decision(syn$`Syn data`, k.real = k.real)
    if(!is.null(clust)){
      m <- m + 1
      synt.data[[m]] <- syn$`Syn data`
      km.synt[[m]] <- clust$Clustering
      k.synt[m] <- clust$`K optim`
      speck[m] <- syn$Specks$SPECKS
    }
  }
  
  ls.km.synt <- lapply(km.synt, label.switching, km.real) # apply label switching to synthetic clusters
  
  # Metrics
  gini.syn <- sapply(km.synt, function(km) Gini(km$size))
        
  mean.distance <- numeric(m)
  for (i in seq_along(km.synt)){
    lsynt.data <- synt.data[[i]]
    lsynt.data$label <- factor(ls.km.synt[[i]]$recluster)
    real.centroid <- aggregate(. ~ label,
                                data = lreal.data[, c(vars, "label")],
                                FUN = base::mean)
    synt.centroid <- aggregate(. ~ label,
                                data = lsynt.data[, c(vars, "label")],
                                FUN = base::mean)
    real.centroid <- real.centroid[order(real.centroid$label), ]
    synt.centroid <- synt.centroid[order(synt.centroid$label), ]
        
    mat.real <- as.matrix(real.centroid[, vars])
    mat.synt <- as.matrix(synt.centroid[, vars])
    distances <- sqrt(rowSums((mat.real - mat.synt)^2))
    mean.distance[i] <- mean(distances)
  }      
  out <- rbind(out,
    data.frame(N = rep(N, m),
               p = rep(p, m),
               k = rep(k, m),
               rho = rep(rho, m),
               method = rep(method, m),
               separation = rep(separation, m),
               k.real = rep(k.real, m),
               Speck = speck,
               count = rep(count, m),
               Gini.real = rep(gini.real, m),
               Gini.synt = gini.syn,
               Mean.distance = mean.distance))
  return(out)
}

n.cores <- detectCores() - 1
cl <- makeCluster(n.cores)
total.m <- 50

clusterEvalQ(cl, { # load packages on each core
  library(mvtnorm)
  library(NbClust)
  library(synthpop)
  library(clue)
  library(DescTools)
  library(factoextra)
})

clusterExport(cl, # export functions and variables to cluster
  varlist = c("simulate_one",
              "data.generation",
              "centroid.generation",
              "synt.data.generation",
              "k.decision",
              "label.switching",
              "specks.correct",
              "total.m"),
  envir = .GlobalEnv)

clusterSetRNGStream(cl, iseed = 123)

# Experimental setup
method <- c("cart")
N <- c(250)
p <- c(2, 5)
rho <- c(0.4)
k <- c(2, 3)
separation <- c(2, 10)

param.grid <- expand.grid(N = N,            # all possible combinations
                          p = p,
                          k = k,
                          rho = rho,
                          method = method,
                          separation = separation,
                          KEEP.OUT.ATTRS = FALSE)

# Run simulations in parallel
result.list <- parLapply(
  cl,
  X   = split(param.grid, seq_len(nrow(param.grid))),
  fun = simulate_one
)

result_cart_250 <- do.call(rbind, result.list)

# Compute Gini difference
result <- result_cart_250 %>%
  mutate(Diff.gini = Gini.real - Gini.synt)

# Gini plot setup
x_min <- min(result$Speck, na.rm = TRUE)
x_max <- max(result$Speck, na.rm = TRUE)
y_min <- min(result$Diff.gini, na.rm = TRUE)
y_max <- max(result$Diff.gini, na.rm = TRUE)

# Plot for Gini difference vs SPECKS
result %>% 
  mutate(index = 50/count) %>%
  group_by(N, p, k, method, separation, rho) %>% 
  group_walk(~ {
    idx <- round(unique(.x$index)*100, 2)       ## % of matching clusters
    
    model <- lm(Diff.gini ~ Speck, data = .x)   ## slope
    slope <- round(coef(model)[2], 2)          
    
    ci <- confint(model, "Speck", level = 0.95) ## confidence interval for slope
    ci_low  <- round(ci[1], 2)
    ci_high <- round(ci[2], 2)
    
    p <- ggplot(.x, aes(Speck, Diff.gini)) +
      geom_point(color = "#1f77b4", alpha = 0.7, size = 2.5) +
      geom_smooth(method = "lm", se = FALSE, color = "#d62728") +
      
      annotate("text", x = x_max, y = y_max, hjust = 1, vjust = 1,
               label = glue("% Equal num of cluster = {idx}"), size = 3, fontface = "plain") +
      
      annotate("text", x = x_max, y = y_max - 0.05*(y_max - y_min), hjust = 1, vjust = 1,
               label = glue("Slope = {slope} [{ci_low}, {ci_high}]"), size = 3, fontface = "plain") +
      
      coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
      labs(
        title = glue( "N={.y$N}  p={.y$p} k={.y$k} rho={.y$rho} ",
                      "Method={.y$method}  Separation={.y$separation}"),
        x = "Speck", y = "Diff Gini", size = "Index"
      ) +
      theme_light() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none",
        panel.grid.minor = element_blank()
      )
    file <- glue("corplot_{.y$N}_{.y$p}_{.y$k}_{.y$rho}",
                 "{.y$method}_{.y$separation}.png")
    ggsave(file, p, width = 6, height = 4, dpi = 300)
  })

# Plot for Mean Distance vs SPECKS
x_min <- min(result$Speck, na.rm = TRUE)
x_max <- max(result$Speck, na.rm = TRUE)

result %>% 
  mutate(index = 50 / count) %>%
  group_by(N, p, k, method, separation, rho) %>% 
  group_walk(~ {
    idx <- round(unique(.x$index) * 100, 2)       ## % of matching clusters
    
    model <- lm(Mean.distance ~ Speck, data = .x) ## slope
    slope <- round(coef(model)[2], 2)
    
    ci <- confint(model, "Speck", level = 0.95)   ## confidence interval
    ci_low  <- round(ci[1], 2)
    ci_high <- round(ci[2], 2)
    
    y_max <- max(.x$Mean.distance, na.rm = TRUE)
    y_min <- min(.x$Mean.distance, na.rm = TRUE)

    p <- ggplot(.x, aes(Speck, Mean.distance)) +
      geom_point(color = "#1f77b4", alpha = 0.7, size = 2.5) +
      geom_smooth(method = "lm", se = FALSE, color = "#d62728") +
      
      annotate("text", x = x_max, y = y_max, hjust = 1, vjust = 1,
               label = glue("% Equal num of cluster = {idx}"), 
               size = 3, fontface = "plain") +
      
      annotate("text", x = x_max, 
               y = y_max - 0.05 * (y_max - y_min), hjust = 1, vjust = 1,
               label = glue("Slope = {slope} [{ci_low}, {ci_high}]"),
               size = 3, fontface = "plain") +
      
      labs(title = glue("N={.y$N} p={.y$p} k={.y$k} rho={.y$rho} ",
                        "Method={.y$method} Separation={.y$separation}σ"),
           x = "Speck", y = "Mean distance") +
      theme_light() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none",
        panel.grid.minor = element_blank()
      ) +
      coord_cartesian(xlim = c(x_min, x_max))  # restrict X-axis only

    file <- glue("corplot_{.y$N}_{.y$p}_{.y$k}_{.y$rho}",
                 "{.y$method}_{.y$separation}.png")
    ggsave(file, p, width = 6, height = 4, dpi = 300)
  })

stopCluster(cl)  # stop the parallel cluster


## ------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

distance <- result %>%
  mutate(separation = factor(separation,
                             levels = sort(unique(separation)),
                             labels = paste0(sort(unique(separation)), "σ")))

# Diff.gini
distance.gini <- ggplot(distance, aes(x = separation, y = Diff.gini)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Separació",
    y     = "Diff.gini"
  ) +
  theme_minimal(base_size = 14)

# Mean.distance
distance.mean <- ggplot(distance, aes(x = separation, y = Mean.distance)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Separació",
    y     = "Mean.distance"
  ) +
  theme_minimal(base_size = 14)

distance.speck <- ggplot(distance, aes(x = separation, y = Speck)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Separació",
    y     = "Speck"
  ) +
  theme_minimal(base_size = 14)

print(distance.gini)
print(distance.mean)
print(distance.speck)


## ------------------------------------------------------------------------------------------------------------------------------------
clusternum <- result %>%
  mutate(separation = factor(k,
                             levels = sort(unique(k)),
                             labels = sort(unique(k))))

# Diff.gini
clusternum.gini <- ggplot(clusternum, aes(x = factor(k), y = Diff.gini)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Clúster"
  ) +
  theme_minimal(base_size = 14)

# Mean.distance
clusternum.mean <- ggplot(clusternum, aes(x = factor(k), y = Mean.distance)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Clúster"
  ) +
  theme_minimal(base_size = 14)

clusternum.speck <- ggplot(clusternum, aes(x = factor(k), y = Speck)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Clúster",
    y     = "Speck"
  ) +
  theme_minimal(base_size = 14)


print(clusternum.gini)
print(clusternum.mean)
print(clusternum.speck)


## ------------------------------------------------------------------------------------------------------------------------------------
variable <- result %>%
  mutate(separation = factor(p,
                             levels = sort(unique(p)),
                             labels = sort(unique(p))))

# Diff.gini
variable.gini <- ggplot(variable, aes(x = factor(p), y = Diff.gini)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Variables"
  ) +
  theme_minimal(base_size = 14)

# Mean.distance
variable.mean <- ggplot(variable, aes(x = factor(p), y = Mean.distance)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Variables"
  ) +
  theme_minimal(base_size = 14)

variable.speck <- ggplot(variable, aes(x = factor(p), y = Speck)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Clúster",
    y     = "Speck"
  ) +
  theme_minimal(base_size = 14)


print(variable.gini)
print(variable.mean)
print(variable.speck)


## ------------------------------------------------------------------------------------------------------------------------------------
correlation <- result %>%
  mutate(separation = factor(rho,
                             levels = sort(unique(rho)),
                             labels = sort(unique(rho))))

# Diff.gini
correlation.gini <- ggplot(clusternum, aes(x = factor(rho), y = Diff.gini)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Correlation"
  ) +
  theme_minimal(base_size = 14)

# Mean.distance
correlation.mean <- ggplot(clusternum, aes(x = factor(rho), y = Mean.distance)) +
  geom_boxplot(fill = "#6E9FC6", outlier.size = 1) +
  labs(
    x     = "Correlation"
  ) +
  theme_minimal(base_size = 14)

print(correlation.gini)
print(correlation.mean)


## ------------------------------------------------------------------------------------------------------------------------------------
library(patchwork)
          
layout <- ((distance.gini  | clusternum.gini  | variable.gini  | correlation.gini) /
          (distance.mean  | clusternum.mean  | variable.mean  | correlation.mean)) + 
  plot_annotation(title = "N = 250")

print(layout)

ggsave(
  filename = "clustering_summary_n250.png",
  plot     = layout,
  width    = 12,
  height   = 6, 
  dpi      = 300
)


## ------------------------------------------------------------------------------------------------------------------------------------
# 1. Load necessary libraries and functions first (if not already loaded)
library(tidyverse)
library(synthpop)
library(mvtnorm)

# 2. Generate ONE specific dataset (e.g., N=250, p=2, k=2, sep=10)
set.seed(123) # Set seed so it's reproducible
rdata <- data.generation(N = 250, p = 2, rho = 0, separation = 10, k = 2)
real.data <- rdata$Data

# 3. Generate the Synthetic version (using CART)
sdata <- synt.data.generation(real.data, m = 1, methods = "cart")
synt.data <- sdata$`Syn data` # Extract the dataframe

# 4. Export to CSV
# Make sure the 'data' folder exists, or just save to root: "original_data.csv"
write.csv(real.data, file = "../data/original_data.csv", row.names = FALSE)
write.csv(synt.data, file = "../data/synthetic_data.csv", row.names = FALSE)

message("Files exported successfully to the data/ folder.")

