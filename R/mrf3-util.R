## Utility functions for mrf3.R

# Get tree net from random forest
get_tree_net <- function(mod, tree.id){

  xvar.names <- mod$xvar.names
  xvar.factor <- mod$xvar.factor
  native.array <- mod$forest$nativeArray
  native.f.array <- mod$forest$nativeFactorArray

  node.stat <- mod$node.stats
  native.array <- cbind(native.array, node.stat) %>% data.frame


  tree.df <- native.array %>% dplyr::filter(treeID == tree.id)

  #converted.tree <- display.tree
  vars.id <- data.frame(var = c("<leaf>", xvar.names), parmID = 0:length(xvar.names), stringsAsFactors = FALSE)
  tree.df$var <- vars.id$var[match(tree.df$parmID, vars.id$parmID)]

  num.node <- tree.df$nodeID %>% unique %>% length

  var.count <- integer(nrow(tree.df))
  for (v in unique(tree.df$var)) {
    pt <- tree.df$var == v
    var.count[which(pt)] <- seq_len(sum(pt))
  }

  tree.df$var_count <- var.count
  tree.df$var_conc <- paste0(tree.df$var, "_", tree.df$var_count)

  var.id <- c((num.node+1): nrow(tree.df))
  tree.df$var.tip.id <- 1:nrow(tree.df)
  tree.df$var.tip.id[tree.df$var != "<leaf>"] <- var.id
  tree.df$var.tip.id[tree.df$var == "<leaf>"] <- tree.df$var_count[tree.df$var == "<leaf>"]

  from_node <- ""
  edge_n <- max(0L, nrow(tree.df) - 1L)
  network <- data.frame(
    from = character(edge_n),
    to = character(edge_n),
    from_id = numeric(edge_n),
    to_id = numeric(edge_n),
    inv_d = numeric(edge_n),
    edge = numeric(edge_n),
    nodesize = numeric(edge_n),
    stringsAsFactors = FALSE
  )
  edge_idx <- 0L
  num.children <- data.frame(tree.df, children = 0)
  num.children <- num.children[num.children$var != "<leaf>",, drop = FALSE]
  num.children <- num.children[!duplicated(num.children$var_conc),, drop = FALSE]
  num_children <- as.list(rep(0, nrow(num.children)))
  names(num_children) <- num.children$var_conc

  for (i in seq_len(nrow(tree.df))) {
    rowi <- tree.df[i, ]
    if (i == 1) {
      from_node <- rowi$var_conc
      from_id <- rowi$var.tip.id
      ns_all <- rowi$nodeSZ
      next
    }

    to_node <- rowi$var_conc
    to_id <- rowi$var.tip.id
    ns_part <- rowi$nodeSZ
    dpthST <- 1/(rowi$dpthST)
    edge_idx <- edge_idx + 1L
    network$from[edge_idx] <- from_node
    network$to[edge_idx] <- to_node
    network$from_id[edge_idx] <- from_id
    network$to_id[edge_idx] <- to_id
    network$inv_d[edge_idx] <- dpthST
    network$edge[edge_idx] <- ns_all/ns_part
    network$nodesize[edge_idx] <- ns_part
    num_children[[from_node]] <- num_children[[from_node]] + 1

    if (rowi$var != "<leaf>") {
      from_node <- to_node
      from_id <- to_id
      ns_all <- ns_part
    } else if (i != nrow(tree.df)) {
      while (num_children[[from_node]] == 2) {
        from_node <- network$from[network$to == from_node]
        from_id <- network$from_id[network$to_id == from_id]
        ns_all <- ns_part
      }
    }
  }


  return(network)

}

# Get leaf node information
get_leaf_ds <- function(mod, tree.membership, net){
  
  # Find upper level node
  prev_leaf <- dplyr::filter(.data = net, is_leaf == 1)
  
  # Update mem_id and membership
  net$mem_id_old <- net$mem_id
  mem_update <- slice_max(.data = prev_leaf, order_by = mem_id, by = from, with_ties = F)
  mem_old <- slice_min(.data = prev_leaf, order_by = mem_id, by = from, with_ties = F)
  mem_new <- setNames(mem_update$mem_id, mem_update$from)
  mem_org_id <- mem_old$mem_id
  names(mem_org_id) <- mem_new
  sapply(1:nrow(mem_old), function(id) net$mem_id[net$mem_id == mem_org_id[id]] <<- as.numeric(names(mem_org_id)[id]))
  
  find_node <- unique(prev_leaf$from)
  
  leaf_select <- filter(net, to %in% find_node)
  map_id <- unique(leaf_select$from)
  
  leaf_select <- filter(net, from %in% map_id)
  leaf_select$mem_id[match(names(mem_new), leaf_select$to)] <- mem_new
  drop_leaf <- leaf_select$from[leaf_select$mem_id == 0]
  leaf_select <- leaf_select[!leaf_select$from %in% drop_leaf,]
  
  # Update leaf information
  
  net$is_leaf[net$to %in% prev_leaf$to] <- 0
  net$is_leaf[net$from %in% leaf_select$from] <- 1
  
  net$mem_id[match(names(mem_new), net$to)] <- mem_new
  
  tree_mem <- tree.membership
  
  sapply(1:length(mem_org_id),
         function(id){
           tree_mem[tree_mem == mem_org_id[id]] <<-
             as.numeric(names(mem_org_id)[id])
         })
  
  
  if(any(leaf_select$mem_id %in% 0)) {
    from_id <- leaf_select[leaf_select$mem_id %in% 0,"from"]
    leaf_select[leaf_select$from %in% from_id,"is_leaf"] <- 0
    leaf_drop <- leaf_select[leaf_select$is_leaf == 0,]
    net[net$to %in% leaf_drop$to,"is_leaf"] <- 0
    net[net$to %in% leaf_drop$to,"leaf_stack"] <- 1
  }
  
  return(
    list(net = net,
         tree.mem = tree_mem)
  )
}

# Get response outcome splitting scores for splits in a tree
get_Y_imp <- function(net, tree.membership, dat, 
                      robust = F, w = NULL, yprob = 1, seed = -5,
                      permute = T,
                      nperm     = 10,   # number of permutations per node
                      alpha     = 0.05, 
                      top = 5){ # p‐value cutoff)
                      
  
  node_id <- unique(net$from)
  var_imp_ls <- llply(
    node_id,
    .fun = function(id){
      children <- filter(.data = net, from %in% id)
      mem_id <- children$mem_id
      idx_node<- which(tree.membership %in% mem_id)
      mem_selected <- tree.membership[idx_node]
      Yt <- scale(dat[idx_node,])
      
      idx_L <- which(mem_selected %in% mem_id[1])
      idx_R <- which(mem_selected %in% mem_id[2])
      nL <- length(idx_L); nR <- length(idx_R)
      q <- ncol(Yt)
      
      split_stat <- sapply(seq_len(q), function(j){
        sum(Yt[idx_L, j])^2/nL +
          sum(Yt[idx_R, j])^2/nR
      })
      split_stat_raw <- split_stat

      #split_stat_sign <- sign(split_stat[[1]] - split_stat[[2]])
      #split_stat <- Reduce("+", split_stat)
      if (robust) {
        if (permute) {
          # null distribution via parallel plyr on permutation index
          G_null <- matrix(0, nrow=q, ncol=nperm)
          for(b in seq_len(nperm)){
            set.seed(b)
            permuted <- sample(nrow(Yt))
            Yp       <- Yt[permuted, , drop = FALSE]
            G_null[,b] <- sapply(seq_len(q), function(j){
              sum(Yp[idx_L, j])^2/nL   +
                sum(Yp[idx_R, j])^2/nR 
            })
          }
          
          pvals <- (rowSums(G_null >= split_stat) + 1)/(nperm + 1)
          split_stat[pvals > alpha] <- 0
        } else {
          split_stat[order(split_stat, decreasing = T)[-c(1:top)]] <- 0
        }
        # Fallback: if robust filtering removes all signal, use unfiltered statistic.
        if (!any(is.finite(split_stat) & split_stat != 0)) {
          split_stat <- split_stat_raw
        }
      } else {

        samp <- rep(1, length(split_stat))
        if(!is.null(w)) {
          ns <- min(ceiling(length(split_stat)*yprob), length(w[w != 0]))
          w <- w/sum(w)
          set.seed(seed)
          samp0 <- sample.int(length(split_stat), ns, prob = w)
          samp0 <- (1:length(split_stat))[-samp0]
          samp[samp0] <- 0
        } else {
          set.seed(seed)
          samp0 <- sample.int(length(split_stat), ceiling(length(split_stat)/3))
          samp[samp0] <- 0
        }
        
      }
      
      
      idx <- which.max(split_stat)
      varY <- mean(split_stat)
      varY_all <- split_stat#/(nL + nR)
      names(varY_all) <- colnames(dat)
      names(varY) <- colnames(dat)[idx]
      
      list(varY = varY, varY_all = varY_all)
    }, .parallel = F)

  var_imp <- unlist(var_imp_ls %>% purrr::map("varY"))
  var_imp_all <- Reduce(rbind, var_imp_ls %>% purrr::map("varY_all"))
  
  list(var_imp = var_imp, var_imp_all = var_imp_all)
  
}

# Get importance for leaves in a tree
get_tree_imp <- function(mod, dat = NULL, robust = F, tree.membership, net, calc = "Both", M = NULL, w = NULL, yprob = 1, weighted = F, seed = -5,
                         permute = T,
                         nperm     = 10,   # number of permutations per node
                         alpha     = 0.05 ){

  if(is.null(dat)){
    dat <- mod$xvar
  }

  if(calc %in% c("Y", "Both")){
    datY <- mod$yvar[rownames(dat),]
  }

  updated_net <- net

  top_node_info <- filter(.data = updated_net, is_leaf == 1)

  node_ds <- setNames(top_node_info$inv_d, top_node_info$from)

  old_net <- updated_net[updated_net$from %in% top_node_info$from,]
  scores_imp <- NULL

  old_net_corr <- group_by(.data = old_net, from)
  old_net_corr <- dplyr::summarise_at(old_net_corr, .vars = "inv_d", .funs = mean)

  node_ds <- node_ds[old_net_corr$from]

  use_sample <- old_net_corr$from

  match_old_net <- old_net[old_net$from %in% use_sample,]

  # Get important variable
  # Get Y

  if(calc %in% c("Both", "Y")){

    impY_ls <- get_Y_imp(net = match_old_net, tree.membership = tree.membership, robust = robust, dat = datY, yprob = yprob, seed = seed, 
                         permute = permute,
                         nperm     = nperm,   # number of permutations per node
                         alpha     = alpha )
    impY <- impY_ls$var_imp
    updated_net$Y_id[unique(match(match_old_net$from, updated_net$from))] <- names(impY)
    scores_impY <- updated_net$inv_d[match(unique(match_old_net$from), updated_net$from)] 
    if(weighted) {
      scores_impY <- scores_impY * updated_net$edge[match(unique(match_old_net$from), updated_net$from)]
    }
    names(scores_impY) <- names(impY)

    impY_mat <- impY_ls$var_imp_all * impY
    # scores_impY <- scores_impY
  }

  # Get X
  if(calc %in% c("Both", "X")){

    imp_var <- top_node_info[match(use_sample, top_node_info$from),"from"]
    if(robust) {
      impX_ls <- get_Y_imp(net = match_old_net, tree.membership = tree.membership, robust = robust, dat = dat, yprob = yprob, seed = seed, 
                           permute = permute,
                           nperm     = nperm,   # number of permutations per node
                           alpha     = alpha )
      impX <- impX_ls$var_imp
    }
   
    sub_net <- old_net_corr[match(unique(match_old_net$from), old_net_corr$from),]
    scores_imp <- (sub_net$inv_d) 
    if(weighted) {
      scores_imp <- scores_imp * old_net_corr$edge
    }

    # scores_imp <- drop_case[match(imp_var,drop_case$to),] %>% pull(corr)
    names(scores_imp) <- gsub("^(.*)_.*", "\\1",imp_var[match(sub_net$from, imp_var)])
    if (robust) impX_mat <- impX_ls$var_imp_all * (sub_net$inv_d) * impX
  }
  if (robust) {
    if(!is.null(impY) | !is.null(impX)) {
      
      if(is.null(dim(impY_mat))) {
        M <- M + sub_net$inv_d * tcrossprod(impX_mat, impY_mat)
      } else {
        M_s <- lapply(1:nrow(impY_mat), function(k) {
          x <- impX_mat[k,]
          y <- impY_mat[k,]
          sub_net$inv_d[k] * tcrossprod(x,y)
        })
        M <- M + Reduce("+", M_s)
      }
    }
    #M[names(scores_imp),] <- M[names(scores_imp),] + impY_mat
    #M[,names(scores_impY)] <- M[,names(scores_impY)] + t(impX_mat)
  }  else { M <- NULL }
  
  if(calc == "Both"){
    scores_imp <- list(X = scores_imp, Y = scores_impY)
  } else if(calc == "X"){
    scores_imp <- list(X = scores_imp)
  } else {
    scores_imp <- list(Y = scores_impY)
  }


  return(
    list(
      net = updated_net,
      dat = dat,
      imp_var = scores_imp,
      M = M
    )
  )
}

# Update tree importance from bottom to top
update_iter_imp <- function(mod, tree.id, calc = "Both", robust = F, w = NULL, yprob = 1, weighted = F, seed = -5,
                            permute = T,
                            nperm     = 10,   # number of permutations per node
                            alpha     = 0.05 ) {

  # Get the tree structure for the specified tree.id
  net <- get_tree_net(mod = mod, tree.id = tree.id)

  # Get the membership information for the specified tree.id
  mem <- mod$membership[, tree.id]

  # Initialize an empty data frame 'dat'
  dat <- NULL

  # Initialize an index 'i' and k
  i <- 1
  k <- 1

  # Initialize an empty list 'imp'
  imp <- list()

  net$Y_id <- NA
  if(robust) {
    mat <- matrix(0, ncol = ncol(mod$yvar), nrow = ncol(mod$xvar))
    
    colnames(mat) <- colnames(mod$yvar)
    rownames(mat) <- colnames(mod$xvar)
  }
  else mat <- NULL 


  # Iterate until all nodes in the tree are leaves or there are only 2 or fewer nodes left
  while (any(is.null(net$is_leaf), sum(net$used_leaf == 0) > 0, is.null(net$used_leaf))) {

    if(k > 1){
      if(is.null(net$used_leaf)) {
        net$used_leaf <- ifelse(net$is_leaf == 1, 1, 0)
      } else{
        net$used_leaf[net$is_leaf == 1] <- 1
      }
      # Get initial leaf information
      update_net <- get_leaf_ds(mod = mod,
                        tree.membership = mem,
                        net = net)
      mem <- update_net$tree.mem
      net <- update_net$net
    } else {
      # Create new columns
      if(is.null(net$is_leaf)) {
        net$is_leaf <- ifelse(grepl("leaf" ,net$to), 1, 0)
      }
      if(is.null(net$mem_id)) {
        net$mem_id <- ifelse(net$is_leaf == 1, net$to_id, 0)
      }
      find_node <- net[net$is_leaf == 1,]
      find_node <- dplyr::summarise(.data = find_node, n = n(), .by = from)
      find_node <- filter(.data = find_node, n == 1)
      net[net$is_leaf == 1,"is_leaf"][net[net$is_leaf == 1,]$from %in% find_node$from] <- 0
    }

    # Get the updated tree importance statistics
    update_ls <- get_tree_imp(
      mod = mod,
      dat = dat,
      robust = robust,
      tree.membership = mem,
      net = net,
      calc = calc,
      M = mat,
      w = w,
      yprob = yprob,
      weighted = weighted,
      seed = seed,
      permute = permute,
      nperm     = nperm,   # number of permutations per node
      alpha     = alpha 
    )

    # Update 'net', 'mem', and 'dat' with the results from 'update_ls'
    net <- update_ls$net
    dat <- update_ls$dat
    
    mat <- update_ls$M

    # If importance for a variable is available, add it to the 'imp' list
    if (!is.null(update_ls$imp_var)) {
      imp[[i]] <- update_ls$imp_var
      # if(i == 1) {
      #   varls <- update_ls$imp_var
      #   varls <- purrr::map(varls, ~.-.)
      #   imp[[i]] <- varls
      # }
      i <- i + 1
    }

    k <- k + 1

  }


  # Filter out the root node(s)
  root_net <- filter(net, is_leaf == 1)

  # Initialize variables for variable importance for X and Y
  # if (calc %in% c("X", "Both")) {
  #
  #   imp_root <- 1/2
  #   names(imp_root) <- gsub("^(.*)_.*", "\\1", unique(root_net$from))
    impX <- unlist(purrr::map(imp, "X"))
  #   impX <- c(impX, imp_root)
  # }
  #
  if (calc %in% c("Y", "Both")) {
    impY <- unlist( purrr::map(imp, "Y") )
  }

  # Create lists of variable importance and variable names
  if (calc == "Both") {
    imp_ls <- list(impX, impY)
    var_name <- list(colnames(mod$xvar), colnames(mod$yvar))
  } else if (calc == "X") {
    imp_ls <- list(impX)
    var_name <- list(colnames(mod$xvar))
  } else {
    imp_ls <- list(impY)
    var_name <- list(colnames(mod$yvar))
  }

  # Calculate variable importance scores for each variable
  imp_col <- plyr::llply(
    1:length(imp_ls),
    .fun = function(l) {
      get_iv(var_name[[l]], imp_ls[[l]])
    }
  )

  # Set appropriate names for the variable importance scores
  if (calc == "Both") {
    names(imp_col) <- c("X", "Y")
  } else if (calc == "X") {
    names(imp_col) <- "X"
  } else {
    names(imp_col) <- "Y"
  }

  # add penalty
  # x_freq <- mod$var.used[tree.id,]
  # x_freq <- x_freq[x_freq != 0]
  # imp_col <- add_lambda(imp_col, net, x_freq, lambda)

  out <- list(imp_ls = imp_col, net = net, M = mat)
  return(out)
}

# Add penalty to weights
add_lambda <- function(imp_ls, net, x_freq, lambda){

  if(all(is.na(net$Y_id))){
    freq <- list(X = x_freq)
  } else {
    freq <- list(X = x_freq, Y = table(net$Y_id))
  }

  new_ww <- plyr::llply(
    names(imp_ls),
    .fun = function(i){

      var <- freq[[i]]
      var_imp_all <- imp_ls[[i]]
      var_imp <- var_imp_all[var_imp_all != 0]
      var <- var[names(var_imp)]

      if(any(var == 1)){
        var_imp[var == 1] <- var_imp[var == 1] * lambda
        var_imp_all[names(var_imp)] <- var_imp
      }

      var_imp_all
    }
  )

  names(new_ww) <- names(imp_ls)

  new_ww

}

#' Get forest importance
#' @param mod A fitted randomForestSRC model.
#' @param parallel Logical; whether to parallelize across trees.
#' @param robust Logical; whether to use robust matrix-based aggregation.
#' @param calc Which importance side to compute: `"X"`, `"Y"`, or `"Both"`.
#' @param weighted Logical; whether to use weighted importance updates.
#' @param use_depth Logical; whether to aggregate non-zero depths instead of simple mean.
#' @param normalized Logical; whether to l2-normalize returned importance.
#' @param w Optional case weights.
#' @param yprob Response sampling proportion used in node-level updates.
#' @param cores Number of CPU cores used when `parallel = TRUE`.
#' @param seed Random seed passed to stochastic components.
#' @param permute Logical; whether to use permutation testing in robust updates.
#' @param nperm Number of permutations per node when `permute = TRUE`.
#' @param alpha Significance level for permutation filtering.
#' @rdname get_imp_forest
get_imp_forest <- function(mod, parallel = T, robust = F, calc = "Both", weighted = F, use_depth = F, normalized = F, 
                           w = NULL, yprob = 1, cores = max(1, detectCores() - 2), seed = -5,
                           permute = T,
                           nperm     = 10,   # number of permutations per node
                           alpha     = 0.05 ){

  nt <- mod$ntree

  if(parallel){
    cores <- max(1L, min(as.integer(cores), as.integer(parallel::detectCores())))
    if(Sys.info()["sysname"] == "Windows"){
      cluster <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cluster)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
    } else {
      doParallel::registerDoParallel(cores)
    }
  }

    if(calc == "Both") {
      cc <- "Both"
      idx <- c("X", "Y")
    }
    if(calc == "X") {cc <- "X";idx <- c("X")}
    if(calc == "Y") {cc <- "Y";idx <- c("Y")}

    results <- plyr::llply(1:nt,
                       .fun = function(t){
                         update_iter_imp(mod,
                                         tree.id = t,
                                         calc = cc,
                                         robust = robust,
                                         yprob = yprob,
                                         w = w,
                                         weighted = weighted,
                                         seed = seed,
                                         permute = permute,
                                         nperm     = nperm,   # number of permutations per node
                                         alpha     = alpha )
                       }, .parallel = parallel)

    imp <- purrr::map(results, "imp_ls")
    net <- purrr::map(results, "net")

    if (robust) {

      M <- purrr::map(results, "M")
      imp <- Reduce("+", M)/nt
      impY <- apply(imp, 2, function(x) sqrt(sum(x^2))/nrow(imp)); impX <- apply(imp, 1, function(x) sqrt(sum(x^2))/ncol(imp))
      # M1 <- M1/max(M1); M2 <- M2/max(M2)
      # impX <- diag(M1); impY <- diag(M2)
      # impX <- rowMeans(M); impY <- colMeans(M)
      if(normalized){
        imp_ls <- list(X = impX/sqrt(sum(impX^2)), Y = impY/sum(sqrt(impY^2)))
      } else {
        imp_ls <- list(X = impX, Y = impY)
      }
      
    } else {
      imp_ls <- plyr::llply(idx,
                            .fun = function(im){
                              im_df <- Reduce(cbind, purrr::map(imp, im))
                              if(!is.null(ncol(im_df))){
                                if(use_depth) {
                                  rs <- rowSums(im_df != 0)
                                  rs[rs == 0] <- -999
                                  imp <- rowSums(im_df)/rs
                                } else {
                                  imp <- rowMeans(im_df)
                                }
                              } else imp <- im_df
                              if(normalized){
                                imp <- imp/sqrt(sum(imp^2))
                              }
                              return(imp)
                            })
    }
    
    
  
    names(imp_ls) <- idx


  out <- list(
    imp_ls = imp_ls,
    imp_ls_init = imp,
    net = net
  )
  return(out)
}

# Match importance name to column name of the data frame
get_iv <- function(var_name, imp){

  var_v <- rep(0, length(var_name))
  names(var_v) <- var_name

  if (length(imp) == 0L || is.null(names(imp))) {
    return(var_v)
  }

  nm <- names(imp)
  ok <- !is.na(nm) & nzchar(nm) & is.finite(imp)
  if (!any(ok)) {
    return(var_v)
  }

  imp <- as.numeric(imp[ok])
  nm <- nm[ok]
  keep <- nm %in% var_name
  if (!any(keep)) {
    return(var_v)
  }

  imp_max <- tapply(imp[keep], nm[keep], max)
  var_v[names(imp_max)] <- as.numeric(imp_max)

  return(var_v)

}

#' Get multi-omics weights
#' @param mod_list A named list of fitted random forest models.
#' @param dat.list A named list of omics datasets used to align weights.
#' @param y Optional response vector (reserved for extensions).
#' @param weighted Logical; whether to use weighted importance updates.
#' @param use_depth Logical; whether to average depth-aware importances.
#' @param robust Logical; whether to use robust matrix-based aggregation.
#' @param parallel Logical; whether to parallelize across models.
#' @param normalized Logical; whether to normalize the merged weights.
#' @param calc Which importance side to compute: `"X"`, `"Y"`, or `"Both"`.
#' @param yprob Response sampling proportion used in node-level updates.
#' @param w Optional case weights.
#' @param cores Number of CPU cores used when `parallel = TRUE`.
#' @param seed Random seed passed to stochastic components.
#' @param permute Logical; whether to use permutation testing in robust updates.
#' @param nperm Number of permutations per node when `permute = TRUE`.
#' @param alpha Significance level for permutation filtering.
#' @param ... Additional arguments for downstream helper functions.
#' @rdname get_multi_weights
get_multi_weights <- function(mod_list, dat.list, y = NULL, weighted = F,  use_depth = F, robust = F,
                              parallel = T, normalized = T, calc = "Both", yprob = 1,
                              w = NULL, cores = max(1, detectCores() - 2), seed = -5, 
                              permute = T,
                              nperm     = 10,   # number of permutations per node
                              alpha     = 0.05, ...){

  mod_names <- names(mod_list)
  # if(length(lambda) == 1) {
  #   lambda <- rep(lambda, length(mod_list))
  #   names(lambda) <- mod_names
  # }

  results <- get_results(mod_list = mod_list, parallel = parallel, robust = robust, weighted = weighted, normalized = F, use_depth = use_depth,
                         calc = calc, w = w, cores = cores, yprob = yprob, seed = seed, permute = permute,
                         nperm     = nperm,   # number of permutations per node
                         alpha     = alpha)

  net <- purrr::map(results, "net")
  weight_l <- purrr::map(results, "wl")
  weight_l_init <- purrr::map(results, "wl_init")
  # M_list <- purrr::map(results, "M")
  
  if(length(weight_l) > 1){
    weight_list <- plyr::llply(
      names(dat.list),
      .fun = function(i){
        w <- purrr::map(weight_l, i)
        w <- purrr::compact(w)
        w <- (Reduce("+", w))/length(w)
        if(normalized) {
          denom <- sqrt(sum(w^2))
          if (is.finite(denom) && denom > 0) {
            w <- w/denom
          }
        }
        w
      }
    )
    names(weight_list) <- names(dat.list)
  } else {
    weight_list <- weight_l[[1]]
    wl <- plyr::llply(
      names(weight_list),
      .fun = function(i){
        w <- weight_list[[i]]
        if(normalized) {
          denom <- sqrt(sum(w^2))
          if (is.finite(denom) && denom > 0) {
            w <- w/denom
          }
        }
        w
      }
    )
    names(wl) <- names(weight_list)
    weight_list <- wl
    weight_list <- weight_list[names(dat.list)]
  }

  out <- list(weight_list = weight_list,
              weight_list_init = weight_l_init,
              net = net)

  return(out)
}


cal_freq <- function(mod, net){

  x_freq <- colMeans(mod$var.used != 0)
  freq <- list(X = x_freq)

  if(is.null(mod$yvar) | class(mod)[3] == "class+"){
    return(freq)
  } else {

    f <- rep(0, ncol(mod$yvar))
    names(f) <- colnames(mod$yvar)
    y_freq <- purrr::map(
      net,
      ~{
        id <- .[["Y_id"]]
        tb <- table(id)
        names(tb)
      }
    )
    y_freq <- table(unlist(y_freq))/mod$ntree
    f[names(y_freq)] <- y_freq

    freq <- c(freq,list(Y = f))

    return(freq)
  }


}

# step_two_weight <- function(two_step, rm_noise, normalized, weight_l, dat.list, freq_ls, s){
# 
#   if(two_step){
#     if(length(freq_ls) > 1){
#       freq <- plyr::llply(
#         names(dat.list),
#         .fun = function(i){
# 
#           fr <- purrr::map(freq_ls, i)
#           fr <- purrr::compact(fr)
#           Reduce("+", fr)/length(fr)
#         }
#       )
#       names(freq) <- names(dat.list)
#     } else {
#       freq <- freq_ls[[1]]
#     }
#   }
# 
#   if(length(weight_l) > 1){
#     weight_list <- plyr::llply(
#       names(dat.list),
#       .fun = function(i){
#         w <- purrr::map(weight_l, i)
#         w <- purrr::compact(w)
#         w <- (Reduce("+", w))/length(w)
#         if(two_step){
#           x <- freq[[i]]
#           # o <- sort(x, decreasing = T)[1:ceiling(a[i] * length(x))]
#           # x[x < min(o)] <- 0
#           w <- w * x
#         }
#         if(rm_noise){
#           x <- s * sd(w)
#           w[w < x] <- 0
#           w[w > x] <- w[w > x] - x
#           
#         } 
#         if(normalized) {
#           w <- w/sqrt(sum(w^2))
#         }
#         
#         w
#       }
#     )
#     names(weight_list) <- names(dat.list)
#   } else {
#     weight_list <- weight_l[[1]]
#     wl <- plyr::llply(
#       names(weight_list),
#       .fun = function(i){
#         w <- weight_list[[i]]
#         if(two_step){
# 
#           x <- freq[[i]]
#           # o <- sort(x, decreasing = T)[1:ceiling(a[i] * length(x))]
#           # x[x < min(o)] <- 0
#           w <- w * x
#         }
#         if(rm_noise){
#           x <- s * sd(w)
#           w[w < x] <- 0
#           w[w > x] <- w[w > x] - x
# 
#         } 
#         if(normalized) {
#           w <- w/sqrt(sum(w^2))
#         }
#         
#         w
#       }
#     )
#     names(wl) <- names(weight_list)
#     weight_list <- wl
#     weight_list <- weight_list[names(dat.list)]
#   }
# 
#   weight_list
# }

get_results <- function(mod_list, parallel, 
                        normalized = F, weighted = F, robust = F,
                        use_depth = F, calc, w = NULL, yprob = 1, cores = max(1, detectCores() - 2), seed = -5, 
                        permute = T,
                        nperm     = 10,   # number of permutations per node
                        alpha     = 0.05){

  mod_names <- names(mod_list)
  plyr::llply(
    mod_names,
    .fun = function(m_name){

      mod <- mod_list[[m_name]]
      # l <- lambda[m_name]
      if(!is.null(w)) {
        w0 <- w[[gsub("_.*", "", m_name)]]
      } else {w0 <- NULL}
      results <- get_imp_forest(mod, parallel = parallel, robust = robust, normalized = normalized, 
                                weighted = weighted, calc = calc,  w = w0, 
                                yprob = yprob, cores = cores, use_depth = use_depth, seed = seed, permute = permute,
                                nperm = nperm, alpha = alpha)

      wl <- results$imp_ls
      wl_init <- results$imp_ls_init
      
      
      net <- results$net

      m_name_sep <- unlist(stringr::str_split(m_name, "_"))

      names(wl) <- rev(m_name_sep)
      # freq <- cal_freq(mod, net)
      # names(freq) <- rev(m_name_sep)
      return(list(
        wl = wl,
        wl_init = wl_init,
        net = net
        # freq = freq
      ))
    },.parallel = F
  )


}
