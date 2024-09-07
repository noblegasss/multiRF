# -------------------------------------------------------------------------------------------------------------
# Reconstruction clustering method
# -------------------------------------------------------------------------------------------------------------

#' MRF unsupervised clustering -- Reconstruction method
#' @param rfit A model list of random forest models
#' @param k Pre-defined number of clusters. The default is selecting the optimal k by tuning method
#' @param dat_use_to_cluster If clustering method selected as Reconstr, the user can choose the data to conduct clustering.
#' Default is the data that has most connections
#'
#' @return mrf3 clustering object
#' @export mrf3_cl_reconstr
#' @inheritParams randomForestSRC::rfsrc
#' @import randomForestSRC

mrf3_cl_reconstr <- function(rfit, 
                             k = NULL,
                             weights_cluster = T,
                             dat_use_to_cluster = "ALL",
                             parallel = T,
                             nearest_coef = 10,
                             ...){

  cl_mod <- get_reconstr_matrix(rfit, nearest_coef = nearest_coef)
  embed <- NULL
  if(weights_cluster) {
    dat_ls <- cl_mod$W
    if(toupper(dat_use_to_cluster) == "ALL"){
      dat <- dat_ls$W_all
    } else {
      dat <- (dat_ls$W_s)[[dat_use_to_cluster]]
    }
    dat <- dat %*% t(dat)
    diag(dat) <- 0
  } else {
    dat_ls <- cl_mod$fused_mat
    if(toupper(dat_use_to_cluster) == "ALL"){
      dat <- Reduce(cbind, dat_ls)
    } else {
      dat <- dat_ls[[dat_use_to_cluster]]
    }
  }
  rownames(dat) <- rownames(rfit[[1]]$xvar)
  if(weights_cluster){
    method <- "Spectral"
  } else {
    method = "PAM"
  }
  if(is.null(k)){
    message("Start tuning k step..")
    clm <- tune_k_clusters(dat, return_cluster = T, method = method)
    cl <- clm$cl
    if(method == "Spectral") embed <- clm$embed

  } else {
    if(method == "PAM"){
      clm <- pam_cl(dat, k_tune = k, diss = F)
    }
    if(method == "Spectral"){
      clm <- spectral_cl(dat, k_tune = k)
      embed <- clm$embed
    }
    cl <- clm$cl
  }

  message("Done!")
  out <- list(dat = dat,
              cl = cl,
              cl_mod = clm,
              embed = embed,
              fused_mat = cl_mod$fused_mat,
              W_mat = cl_mod$W,
              dat_used = dat_use_to_cluster,
              method = "Reconstr")

  class(out) <- "reconstr"
  return(out)

}


#' @export
#' @param rfit A model list of random forest models
#' @rdname mrf3_cl_reconstr

get_reconstr_matrix <- function(rfit, nearest_coef = 10){

  fused_mat <- plyr::llply(
    names(rfit),
    .fun = function(m){

      mod_names <- unlist(purrr::map(m, ~stringr::str_split(.,"_")))
      mod <- rfit[[m]]

      forest.wt <- mod$forest.wt
      W <- forest.wt/(1 - diag(forest.wt))
      
      diag(W) <- 0
      if(!is.null(nearest_coef)){
        W <- apply(W, 1, function(i) {

          o <- sort(i, decreasing =T)[1:nearest_coef] 
          i[i < min(o)] <- 0
          i})
      }
      W <- t(W)
      xvar <- mod$xvar
      yvar <- mod$yvar
      #pred_X <- forest.wt %*% as.matrix(xvar)

      if(!is.null(yvar)){
        out <- list(W %*% as.matrix(xvar),
                    W %*% as.matrix(yvar))
        #pred_Y <- forest.wt %*% as.matrix(yvar)

        #oob <- mean(c(colMeans((pred_Y - W %*% as.matrix(yvar))^2) , colMeans((pred_X - W %*% as.matrix(xvar))^2)))
          
      } else {
        out <- list(W %*% as.matrix(xvar))
        #oob <- mean(colMeans((pred_X - W %*% as.matrix(xvar))^2))
      }

      names(out) <- rev(mod_names)
      out <- list(mat = out, W = W)#, oob = oob)
      return(out)
    }
  )

  mat <- purrr::map(fused_mat, "mat")
  mat <- Reduce(c,mat)
  W <- purrr::map(fused_mat, "W")
  names(W) <- names(rfit)
  #oob <- purrr::map(fused_mat, "oob")
  #oob <- Reduce("+",oob)
  
  # All omics-weights
  W_all <- Reduce("+", W)/length(W)

  full_fused_ls <- plyr::llply(
    unique(names(mat)),
    .fun = function(i){

      f_l <- mat[i]
      df <- Reduce("+", f_l)/length(f_l)
      
      # Omics-specific weights
      W_s <- NULL
      k <- 0
      for(m in names(rfit)){
        mod_names <- stringr::str_split(m, "_")[[1]]
        W_s <- W_s + W[[m]]
        k <- k + 1
      }
      W_s <- W_s/k

      list(m = df,
           W_s = W_s)
    }
  )
  full_fused_mat <- purrr::map(full_fused_ls, "m")
  W_s <- purrr::map(full_fused_ls, "W_s")
  names(full_fused_mat) <- names(W_s) <- unique(names(mat))
  
  list(fused_mat = full_fused_mat,
       W = list(W_all = W_all,
                W_s = W_s))
       #oob = oob)

}
