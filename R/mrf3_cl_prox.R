#' MRF unsupervised clustering -- Proximity method
#' @param rfit A model list of random forest models
#' @param k Pre-defined number of clusters. The default is selecting the optimal k by tuning method
#' @param enhanced A logical parameter that determines whether enhanced proximity is calculated when selecting Proximity as clustering method.
#' The default is True.
#'
#' @return mrf3 clustering object
#' @export
#' @inheritParams randomForestSRC::rfsrc
#' @import randomForestSRC
#'
# -------------------------------------------------------------------------------------------------------------
# Proximity clustering method
# -------------------------------------------------------------------------------------------------------------

mrf3_cl_prox <- function(rfit, k = NULL,
                         enhanced = T, thres = NULL,
                         size_min = 5, use = "X", symm = T, calc_imp = F,
                         parallel = T,
                         sparse = F,
                         method_cl = "PAM",
                         cores = max(detectCores() - 2,20),
                         ...){

  if(length(rfit) == 1){
    if(is.null(rfit[[1]]$yvar)) symm <- F
  }
  if(enhanced){
    cl_mod <- plyr::llply(
      rfit,
      .fun = function(r){
        cl_forest(r,
                  size_min = size_min,
                  parallel = parallel,
                  use = use,
                  symm = symm,
                  calc_imp = calc_imp,
                  cores = cores,
                  ...)
      }
    )

    prox <- purrr::map(cl_mod, "prox")
    prox <- Reduce("+", prox)
  } else {

    cl_mod <- NULL
    prox <- purrr::map(rfit, "proximity")

    prox <- Reduce("+", prox)

  }

  num_dim <- diag(prox)[1]
  prox <- prox/num_dim
  
  # Make sparse proximity
  if(enhanced & sparse){
    diag(prox) <- 0
    prox[prox < estimate_mode(prox)] <- 0
    diag(prox) <- 1
  }
  
  
  rownames(prox) <- colnames(prox) <- rownames(rfit[[1]]$xvar)

  if(method_cl == "PAM") {
    p <- 1 - prox 
  } 
  if(method_cl == "Spectral"){
    p <- prox
    diag(p) <- 0
  }
  if(is.null(k)){
    message("Start tuning k step..")
    
    k <- tune_k_clusters(p, return_cluster = T, method = method_cl)
    cl <- k$cl

  } else {
    if(method_cl == "PAM"){
      cl <- pam_cl(p, k_tune = k, diss = F)
    } 
    if(method_cl == "Spectral"){
      cl <- spectral_cl(p, k_tune = k)
    }
    cl <- cl$cl
  }

  message("Done!")
  out <- list(dat = prox,
              cl = cl,
              cl_mod = cl_mod,
              enhanced = enhanced,
              method = "Proximity")
  class(out) <- "prox"

  return(out)
}

get_prox <-  function(class_mem){

  if(length(unique(class_mem)) == 1){

    class_new <- matrix(0, nrow = length(class_mem), ncol = length(class_mem))

  } else {

    class_new <- data.frame(class = as.factor(class_mem))
    one_hot <-  model.matrix(~class + 0, class_new)
    class_new <- one_hot %*% t(one_hot)

  }

  return(class_new)
}

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

