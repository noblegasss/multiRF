
#' Find pairwise IMD between variables across datasets
#' @param x A mrf3 object
#' @param all_var A logical parameter that determines whether to compute the connections of all variables or selected variables.
#' The default is FALSE for memory saving.
#' @return A pairwise adjacency matrix between variables and a pairwise connection matrix between datasets
#' @export pairwise_imd
pairwise_imd <- function(x, ...) {
  UseMethod("pairwise_imd")
}

#' @rdname pairwise_imd
#' @method pairwise_imd mrf3
#' @export

pairwise_imd.mrf3 <- function(x, all_var = F, normalized = F){
  
  if(!is.na(class(x)[2])) {
    if(class(x)[2] == "cl"){
      mod <- x$mod
    } else {
      mod <- x
    }
  } else {
    if(is.na(class(x)[2])) {
      if(is.na(class(x)[1])) stop("Must be a mrf3 object") else {
        mod <- x
      }
    }
  }
  
 
  net <- mod$net
  connection <- mod$connection
  dat_names <- names(mod$weights)
  
  nt <- mod$ntree
  
  var_names <- purrr::map(mod$weights, ~names(.))
  
  num_dat <- length(mod$weights)
  am <- matrix(0, nrow = num_dat, ncol = num_dat,
               dimnames = list(dat_names, dat_names))
  
  d <- do.call(rbind, connection) 
  length_count <- table(d)[dat_names]
  
  for (i in 1:nrow(d)) {
    am[d[i, 1], d[i, 2]] <- 1
    am[d[i, 2], d[i, 1]] <- 1
  }
  diag(am) <- length_count
  
  large_mat <- matrix(0, nrow = length(unlist(var_names)), ncol = length(unlist(var_names)),
                      dimnames = list(unlist(var_names), unlist(var_names)))
  
  l_ply(
    1:length(net),
    .fun = function(i){
      ne <- net[[i]]
      vn <- var_names[connection[[i]]]
      names(vn) <- c("yvar", "xvar")

      mat <- pairwise_imd_mod(ne, vn)
      mat <- mat[rownames(mat) %in% rownames(large_mat), colnames(mat) %in% colnames(large_mat)]
      large_mat[rownames(mat),colnames(mat)] <<- large_mat[rownames(mat),colnames(mat)] + mat
    }
  )
  
  for (i in 1:ncol(am)) {
    for(j in (i:ncol(am))){
      m1 <- colnames(am)[i]
      m2 <- colnames(am)[j]
      if(am[i,j] != 0){
        large_mat[var_names[[m1]],var_names[[m2]]] <- large_mat[var_names[[m1]],var_names[[m2]]]/am[i,j]
        if(i != j){
          large_mat[var_names[[m2]],var_names[[m1]]] <- large_mat[var_names[[m2]],var_names[[m1]]]/am[i,j]
        }
      }
    }
  }
  
  if(!all_var) {
    large_mat <- large_mat[names(Reduce(c,mod$weights))[Reduce(c,mod$weights) != 0],
                           names(Reduce(c,mod$weights))[Reduce(c,mod$weights) != 0]]
    var_names <- purrr::map(mod$weights, ~names(.)[. > 0])
    
  }

  large_mat <- large_mat/nt
  if(normalized) {
    d <- rowSums(large_mat)
    l <- diag(d^(-1/2)) %*% large_mat %*% diag(d^(-1/2))
    dimnames(l) <- dimnames(large_mat)
    large_mat <- l
  }
  
  out <- list(
    adj_var_mat = large_mat,
    adj_dat_mat = am,
    var_use = var_names
  )
  
  return(out)
  
}


pairwise_imd_mod <- function(net, var_name){
  
  var_mat <- matrix(0, nrow = length(var_name$xvar) + length(var_name$yvar),
                    ncol = length(var_name$xvar) + length(var_name$yvar))
  
  dimnames(var_mat) <- list(c(var_name$xvar, var_name$yvar),
                            c(var_name$xvar, var_name$yvar))
  
  plyr::l_ply(
    net,
    .fun = function(i){
      netdf <- data.frame(y = i$Y_id, x = gsub("^(.*)_.*", "\\1",i$from), imd = i$inv_d)
      netdf <- na.omit(netdf)
      netdf <- netdf %>% filter(
        y %in% var_name$yvar & x %in% var_name$xvar
      )

      plyr::l_ply(1:nrow(netdf),
                  .fun = function(j) {

                    var_mat[netdf[j,"x"],netdf[j,"y"]] <<- var_mat[netdf[j,"x"],netdf[j,"y"]] + netdf[j,"imd"]
                    var_mat[netdf[j,"y"],netdf[j,"x"]] <<- var_mat[netdf[j,"y"],netdf[j,"x"]] + netdf[j,"imd"]
                  })
    }
  )
  
  return(var_mat)
  
}
