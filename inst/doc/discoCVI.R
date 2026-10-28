## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 6, fig.height = 4.2,
                      eval = requireNamespace("dbscan", quietly = TRUE))

## ----data---------------------------------------------------------------------
library(discoCVI)
library(dbscan)

set.seed(42)
make_moons <- function(n, sd = 0.08) {
  n1 <- n %/% 2
  t1 <- runif(n1, 0, pi)
  t2 <- runif(n - n1, 0, pi)
  moons <- rbind(cbind(cos(t1), sin(t1)),
                 cbind(1 - cos(t2), 0.5 - sin(t2)))
  moons + matrix(rnorm(2 * n, sd = sd), ncol = 2)
}
X <- rbind(make_moons(600),
           cbind(runif(40, -1.5, 2.5), runif(40, -1, 1.5)))
dim(X)
plot(X, pch = 19, cex = 0.5, asp = 1, xlab = "x1", ylab = "x2")

## ----hdbscan------------------------------------------------------------------
hdb <- hdbscan(X, minPts = 10)
labels_hdb <- hdb$cluster
labels_hdb[labels_hdb == 0] <- -1L
table(labels_hdb)

disco_score(X, labels_hdb)

## ----kmeans-------------------------------------------------------------------
labels_km <- kmeans(X, centers = 2, nstart = 10)$cluster
disco_score(X, labels_km)

## ----tuning-------------------------------------------------------------------
score_minpts <- function(m) {
  lab <- hdbscan(X, minPts = m)$cluster
  lab[lab == 0] <- -1L
  disco_score(X, lab)
}
grid <- c(3, 5, 10, 15, 20, 30, 50)
scores <- sapply(grid, score_minpts)
data.frame(minPts = grid, DISCO = round(scores, 3))
plot(grid, scores, type = "b", pch = 19, xlab = "minPts",
     ylab = "DISCO score")

## ----pointwise----------------------------------------------------------------
s <- disco_samples(X, labels_hdb)
summary(s)

pal <- hcl.colors(100, "RdYlBu")
col <- pal[cut(s, breaks = seq(-1, 1, length.out = 101), include.lowest = TRUE)]
plot(X, col = col, pch = ifelse(labels_hdb == -1, 4, 19), cex = 0.7,
     asp = 1, xlab = "x1", ylab = "x2",
     main = "Pointwise DISCO (blue = good, red = poor; x = noise)")

## ----noise--------------------------------------------------------------------
nz <- p_noise(X, labels_hdb)
noise_scores <- data.frame(point    = which(labels_hdb == -1),
                           p_sparse = nz$p_sparse,
                           p_far    = nz$p_far)
head(noise_scores[order(pmin(noise_scores$p_sparse, noise_scores$p_far)), ])

## ----reuse--------------------------------------------------------------------
D <- compute_dc_distances(X, min_points = 5)
dim(D)
nn <- labels_hdb != -1
summary(p_cluster(D[nn, nn], labels_hdb[nn], precomputed_dc_dists = TRUE))

