# -------------------------------------------------------------------------------------------------------------
#' k clusters tuning function
#' @export
#' @param dat A proximity or similarity matrix that use to cluster
#' @param return_cluster A logical parameter that determine whether the clustering results to be returned
#' @param plot_k A logical parameter that determine whether the plot of tuning should be showed
#' @param input Input for clustering in PAM. The default is using the dissimilarity matrix 1-dat. The alternative option is
#' using the first k + 1 eigenvectors from laplacian matrix (embed).
#'
#' @return optimal number of k
#' @inheritParams cluster::pam

tune_k_clusters <- function(x,...){
  UseMethod("tune_k_clusters")
}

#' @rdname tune_k_clusters
#' @export
tune_k_clusters.default <- function(x, return_cluster = F, plot_k = F,
                                    method = "PAM", tune_method = "silhouette", gap_w = "uniform", prox = F,...){

  k_tune <- seq(2,12,by = 1)

  if(method == "Spectral"){
    cl <- spectral_cl(x, k_tune = k_tune, gap_w = gap_w,  ...)
  }
  
  if(method == "PAM"){
    if(prox){
      x <- 1-x
      diss <- T
    } else {diss <- F}
    cl <- pam_cl(x, k_tune = k_tune, diss = diss, tune_method = tune_method, ...)
  }

  if(plot_k){
    if(method == "PAM"){
      if(tune_method == "silhouette"){
        m <- 2:9
        ylab <- "Silhouette"
      } 
      if(tune_method == "ratio"){
        ylab <- "Ratio"
        m <- 2:10
        par(mfrow = c(1,2))
        col <- rep(1,8)
        col[cl$best_k - 1] <- 2
        plot(2:9, cl$diff, pch = 19, ylab = "Ratio Diff", xlab = "k", col = col)
      } 
      col <- rep(1,length(m))
      col[cl$best_k - 1] <- 2
      plot(m, cl$sil, pch = 19, ylab = ylab, xlab = "k", col = col)
    } else {
      par(mfrow = c(1,2))
      col <- rep(1,8)
      col[cl$best_k - 1] <- 2
      plot(2:9, cl$diff_e, pch = 19, ylab = "Diff eigenvalue", xlab = "k", col = col)
      col <- rep(1,9)
      col[cl$best_k - 1] <- 2
      plot(2:10, cl$ev, pch = 19, ylab = "Eigenvalues", xlab = "k", col = col)
    } 
      
  
  }


  if(return_cluster){
    # cl <- cluster::pam(p, k = k, diss = T, cluster.only = T, ...)
    k <- cl
  } else {k <- cl$best_k}

  return(k)
}

#' @rdname tune_k_clusters
#' @export
tune_k_clusters.mrf3 <- function(x, return_cluster = F, plot_k = F, method = "PAM", tune_method = "silhouette", gap_w = "uniform", prox = F, ...){

  if(class(x)[2] %in% "cl"){
    x <- x$dat
  } else {
    stop("Must be a mrf cl object")
  }

  k_tune <- seq(2,9,by = 1)
  
  if(method == "Spectral"){
    cl <- spectral_cl(x, k_tune = k_tune, gap_w = gap_w, ...)
  }
  
  if(method == "PAM"){
    if(prox){
      x <- 1-x
      diss <- T
    } else {diss <- F}
    cl <- pam_cl(x, k_tune = k_tune, diss = diss, tune_method = tune_method, ...)
  }
  
  if(plot_k){
    if(method == "PAM"){
      if(tune_method == "silhouette"){
        m <- 2:9
        ylab <- "Silhouette"
      } 
      if(tune_method == "ratio"){
        ylab <- "Ratio"
        m <- 2:10
        par(mfrow = c(1,2))
        col <- rep(1,8)
        col[cl$best_k - 1] <- 2
        plot(2:9, cl$diff, pch = 19, ylab = "Ratio Diff", xlab = "k", col = col)
      } 
      col <- rep(1,length(m))
      col[cl$best_k - 1] <- 2
      plot(m, cl$sil, pch = 19, ylab = ylab, xlab = "k", col = col)
    } else {
      par(mfrow = c(1,2))
      col <- rep(1,8)
      col[cl$best_k - 1] <- 2
      plot(2:9, cl$diff_e, pch = 19, ylab = "Diff eigenvalue", xlab = "k", col = col)
      col <- rep(1,9)
      col[cl$best_k - 1] <- 2
      plot(2:10, cl$ev, pch = 19, ylab = "Eigenvalues", xlab = "k", col = col)
    } 
    
    
  }
  
  
  if(return_cluster){
    # cl <- cluster::pam(p, k = k, diss = T, cluster.only = T, ...)
    k <- cl
  } else {k <- cl$best_k}
  return(k)
}

#' @rdname tune_k_clusters
#' @export
spectral_cl <- function(x, k_tune = seq(2,12,by = 1), gap_w = "uniform", d = NULL, ...){
  
  if(is.null(d)){
    d <- rowSums(x)
  }
  l <- diag(1, nrow(x)) - diag(d^(-1/2)) %*% x %*% diag(d^(-1/2))
  
  e <- eigen(l, symmetric = T)
  eigenvectors <- e$vectors
        
  if(length(k_tune) > 1){
    
    eigenvalues <- rev(e$values)[2:(length(k_tune) + 2)]
    
    if(gap_w == "log"){
      w <- log(k_tune)
    }
    if(gap_w == "uniform"){
      w <- rep(1, length(k_tune))
    }
    
    diff_e <- diff(eigenvalues) * w
    k <- which.max(diff_e) + 1

  } else  {
    k <- k_tune

    eigenvalues <- rev(e$values)[2:(k + 1)]
    diff_e <- diff(eigenvalues)[k-1]
  }

  n <- nrow(x)
  
  mat <- as.matrix(eigenvectors[,(n-k+1):(n)])
  mat_norm <- mat/sqrt(rowSums(mat^2))
 
  cl <- pam(mat, k = k, ...)

  embed <- as.matrix(eigenvectors[,(n-k):(n-1)])
  ei <- rev(e$values)[2:(k+1)]
  embed <- embed[,2:k]
  
  return(list(
    best_k = k, cl = cl$cluster, diff_e = diff_e, ev = eigenvalues, embed = embed, obj = cl$objective[2], sil = cl$silinfo$avg.width
  ))
}

#' @rdname tune_k_clusters
#' @export
pam_cl <- function(x, k_tune = seq(2,9,by = 1), diss = T, tune_method = "silhouette", ...){
  
  sil <- c()
  if(length(k_tune) > 1){
    if(!diss){
      x <- as.matrix(dist(x))
    } 
    
    if(tune_method == "ratio"){
      k_tune <- c(k_tune, max(k_tune) + 1)
      m <- 1-x
      dist_all <- mean(m)
    } 
    for (k in k_tune){
      cl <- pam(x, k = k, diss = T, ...)
      cl_cluster <- cl$cluster
      if(tune_method == "silhouette"){
        sil <- c(sil, cl$objective[1] - cl$objective[2])
      }
      if(tune_method == "ratio"){
        cl_unique <- unique(cl$cluster)
        
        class_dist <- sapply(cl_unique, function(i){
          mean(m[cl$cluster == i,cl$cluster != i])
        })
        
        sil <- c(sil,(mean(class_dist))/dist_all)
      }
     
    }
    if(tune_method == "silhouette"){
      k <- which.max(sil) + 1
      diff_S <- NULL
    }
    if(tune_method == "ratio"){

      diff_S <- abs(diff(sil))*log(k_tune[1:(max(k_tune) - 2)])

      k <- which.max(diff_S) + 1
    }
  } else {
    k <- k_tune
    diff_S <- NULL
  }
  
  cl <- pam(x, k = k, diss = diss, ...)
  
  return(
    list(best_k = k, cl = cl$cluster, sil = sil, diff = diff_S, obj = cl$objective[2])
  )
  
}
