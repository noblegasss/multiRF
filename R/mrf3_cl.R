#' MRF unsupervised clustering
#'
#' @param dat.list A list that contains multi-omics datasets with samples in columns and features in rows. Samples should be matched in each dataset.
#' @param mod A pre-defined model list from random forest. If the list is provided, the clustering will based on the models in the list
#' @param k Pre-defined number of clusters. The default is selecting the optimal k by tuning method
#' @param method Clustering methods. You can choose from "Proximity" and "Reconstr". The default is "Proximity".
#' @param tune_k_method Method for tuning number of clusters. The default is "Ratio", which is calculated by separation/scaled diameter.
#' The other option is using "Silhouette" information.
#' @param dat_use_to_cluster If clustering method selected as Reconstr, the user can choose the data to conduct clustering.
#' Default is the data that has most connections
#' @param enhanced A logical parameter that determines whether enhanced proximity is calculated when selecting Proximity as clustering method. The default is True.
#' @param feature_selection A logical parameter that determines whether feature selection is conducted before clustering. The default is True.
#' @param connect_list A pre-defined connection list between datasets. If provided, variable selection will be conducted based on this connection list.
#' If not provided, the algorithm will find the optimal connection between datasets.
#' @param dim Embedding dimensions of datasets when finding optimal connections. The default is 20.
#' @param direct A logical parameter that determines whether to keep both directions in the connection list for optimal connections.
#' @param thres Penalize threshold matrix. If thres is NULL, the function will derive the optimal thresholds from a pre-selected threshold list.
#' @param ntree1 Number of trees for variable selection. Default is 100
#' @param ntree2 Number of trees for refitting the MRF model. The default is 300.
#' @param screening A logical parameter that determines whether to conduct a screening of noise variables before model fitting. The default is False.
#' @param parallel A logical parameter that determines whether to use parallel computation in the calculation of weights.
#'
#' @inheritParams randomForestSRC::rfsrc
#'
#' @return mrf3 clustering object
#' @export
#' @import randomForestSRC
#' @examples
#' library(dplyr)
#' dat.list <- InterSIM::InterSIM(n.sample = 150, cluster.sample.prop=c(0.20,0.30,0.40,0.10),
#' delta.methyl = 1, delta.expr = 1, delta.protein = 1 )
#' dat.list$dat.protein <- dat.list$dat.protein %>% janitor::clean_names()
#'
#' # use proximity method
#' mod <- mrf3_cl(dat.list[1:3], method = "Proximity")
#' mclust::adjustedRandIndex(mod$cl, dat.list$clustering.assignment$cluster.id)
#'
#' # use reconstr method
#' mod <- mrf3_cl(dat.list[1:3], method = "Reconstr")
#' mclust::adjustedRandIndex(mod$cl, dat.list$clustering.assignment$cluster.id)

# -------------------------------------------------------------------------------------------------------------
# Wrapper function for MRF unsupervised clustering
# -------------------------------------------------------------------------------------------------------------
mrf3_cl <- function(mod, dat.list, k = NULL,
                    method = "Reconstr",
                    dat_use_to_cluster = NULL,
                    weights_cluster = T,
                    enhanced = F,
                    parallel = T,
                    method_cl = "Spectral", 
                    t_max = 10,
                    cores = detectCores() - 2, ...){

  if(length(dat.list) == 1) type = "unsupervised"
  if(length(dat.list) > 1) type = "regression"
  
  rfit <- mod$mod
  connect_list <- mod$connection

  message("Start clustering step..")
  if(method == "Proximity"){
    out0 <- mrf3_cl_prox(
      rfit = rfit,
      type = type,
      k = k,
      enhanced = enhanced,
      method_cl = method_cl,
      cores = cores,
      ...
    )
  }

  if(method == "Reconstr"){
    
    if(weights_cluster){
      if(is.null(dat_use_to_cluster)){
        dat_use_to_cluster <- "ALL"
      }
    } else {
      if(is.null(dat_use_to_cluster)){
        
        connect_freq <- table(unlist(connect_list))
        dat_use_to_cluster <- names(connect_freq)[which.max(connect_freq)]
        
      }
    }
    
    out0 <- mrf3_cl_reconstr(
      rfit = rfit,
      weights_cluster = weights_cluster,
      dat_use_to_cluster = dat_use_to_cluster,
      type = type,
      k = k,
      t_max = t_max,
      ...
    )
  }
  out <- c(out0, list(mod = mod$mod))

  class(out) <- c("mrf3", "cl", class(out0))

  return(out)

}

