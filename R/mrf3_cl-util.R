## Utility functions for mrf3_cl.R

# Get Spearman correlation of two memberships
get_corr <- function(tree.membership, weights, net, map_id, symm){

  node_id <- unique(tree.membership)
  mem_select <- filter(net, from %in% map_id)
  mem_select <- mem_select$mem_id

  mem_idx <- which(node_id %in% mem_select)

  if(symm){
    cor_w <- plyr::llply(
      weights,
      .fun = function(w){
        cor(w[,mem_idx[1]], w[,mem_idx[2]], method = "spearman")
      }
    )
    cor_w <- unlist(cor_w)
    cor_w <- mean(cor_w)
  } else {
    cor_w <- cor(weights[,mem_idx[1]], weights[,mem_idx[2]])
  }


  return(cor_w)

}

# Get weights of variable based on correlation
get_tree_weight <- function(dat, tree.membership){

  dat.mean <- colMeans(dat, na.rm = T)
  node.index <- unique(tree.membership)

  class.info <- plyr::llply(
    node.index,
    .fun = function(ni){

      dat.ni <- as.matrix(dat[tree.membership == ni,])
      class.mean <- mean(dat.ni)
      dat.class.mean <- colMeans(dat.ni, na.rm = T)

      d0 <- dat.class.mean - dat.mean

      d <- d0/sqrt(sum(d0^2))

      return(list(
        class.mean = class.mean,
        d = d
      ))
    }
  )

  class_mean <- purrr::map(class.info, "class.mean")
  class_mean <- unlist(class_mean)

  weights <- purrr::map(class.info, "d")
  weights <- purrr::reduce(weights, cbind)

  return(list(w = weights, class_mean = class_mean))


}

# Get correlation in terminal nodes
get_leaf_corr <- function(mod, tree.membership, net, symm = T, use = "X"){

  # Create new columns
  if(is.null(net$corr)) {
    net$corr <- 0
  }
  if(is.null(net$is_leaf)) {
    net$is_leaf <- ifelse(grepl("leaf" ,net$to), 1, 0)
  }
  if(is.null(net$leaf_used)) net$leaf_used <- 0
  if(is.null(net$mem_id)) {
    net$mem_id <- ifelse(net$is_leaf == 1, net$to_id, 0)
  }


  # Get weights
  if(symm){
    dat <- list(X = mod$xvar, Y = mod$yvar)
    w <- llply(
      dat,
      .fun = function(d){

        get_tree_weight(d, tree.membership)
      }
    )
    weights <- purrr::map(w, "w")
  } else {
    if(use == "X"){
      dat <- mod$xvar
    }

    if(use == "Y"){
      dat <- mod$yvar
    }

    w <- get_tree_weight(dat, tree.membership)
    weights <- w$w
  }

  # Find upper level node
  find_node <- dplyr::filter(.data = net, is_leaf == 1)
  find_node <- dplyr::summarise(.data = find_node, n = n(), .by = from)
  find_node <- filter(.data = find_node, n == 2)
  find_node <- find_node$from

  leaf_select <- filter(net, from %in% find_node)

  map_id <- unique(leaf_select$from)

  # Calculate correlation
  corr_leaf <- map_id %>%
    purrr::map(~get_corr(map_id = .,
                         tree.membership = tree.membership,
                         weights = weights,
                         net = leaf_select,
                         symm = symm)) %>%
    unlist()


  # Map correlation to node
  net$leaf_used[match(leaf_select$to, net$to)] <- 1
  net$corr[match(map_id, net$to)] <- corr_leaf

  # Update leaf information

  net$is_leaf[match(map_id, net$to)] <- 1
  net$is_leaf[match(leaf_select$to, net$to)] <- 0

  # Update mem_id and membership
  net$mem_id_old <- net$mem_id

  mem_update <- slice_max(.data = leaf_select, order_by = mem_id, by = from)
  mem_new <- setNames(mem_update$mem_id, mem_update$from)
  mem_new <- mem_new[leaf_select$from]

  mem_org <- filter(leaf_select, !mem_id %in% mem_new)
  mem_org_id <- setNames(mem_org$mem_id, mem_org$from)

  mem_org_id <- mem_org_id[unique(names(mem_new))]
  names(mem_org_id) <- unique(mem_new)

  net$mem_id[match(mem_update$from, net$to)] <- mem_update$mem_id
  net$mem_id[which(net$from %in% names(mem_new))] <- mem_new

  sapply(1:length(mem_org_id), function(id) net$mem_id[net$mem_id == mem_org_id[id]] <<- as.numeric(names(mem_org_id)[id]))
  net$mem_id <- as.numeric(net$mem_id)

  sapply(1:length(mem_org_id),
         function(id){
           tree.membership[tree.membership == mem_org_id[id]] <<-
             as.numeric(names(mem_org_id)[id])
         })


  return(
    list(net = net,
         tree.mem = tree.membership)
  )
}

# Update tree leaf weights
update_tree_leaf <- function(mod, dat = NULL, tree.membership, net, size_min = 10, symm = T,
                             use = "X", calc_imp = T){

  if(is.null(dat)){
    if(use == "X") dat <- mod$xvar
    if(use == "Y") dat <- mod$yvar
  }

  datY <- mod$yvar
  if(!is.null(datY)){
    datY <- datY[rownames(dat),]
  }


  # Get initial correlation
  update <- get_leaf_corr(mod = mod,
                          tree.membership = tree.membership,
                          net = net,
                          symm = symm,
                          use = use)

  updated_net <- update$net

  top_node_info <- filter(.data = updated_net, is_leaf == 1 & !grepl("<leaf>", to))

  node_corr <- setNames(top_node_info$corr, top_node_info$to)

  old_net <- updated_net[updated_net$from %in% top_node_info$to,]
  classified <- NULL
  scores_imp <- NULL

  if(!all(old_net$corr == 0)){

    old_net_corr <- filter(old_net, corr != 0)
    old_net_corr <- group_by(old_net_corr, from)
    old_net_corr <- summarise_at(old_net_corr, c("corr", "nodesize"), list(max = max, min = min))

    node_corr <- node_corr[old_net_corr$from]

    old_net_corr <- mutate(old_net_corr,
                           new_corr = node_corr,
                           mem_drop = ifelse((new_corr < corr_max|corr_max < 0) & nodesize_min > size_min, 1, 0))

    if(any(old_net_corr$mem_drop == 1)){

      drop_df <- old_net_corr[old_net_corr$mem_drop == 1,]
      drop_sample <- drop_df$from

      match_old_net <- old_net[old_net$from %in% drop_sample,]
      drop_case <- slice_max(match_old_net, order_by = corr, by = from)

      if(calc_imp){
        # Get important variable
        # Get X

        imp_var <- top_node_info[match(drop_sample, top_node_info$to),"from"]
        scores_imp <- abs(drop_df$new_corr - drop_df$corr_max)

        # scores_imp <- drop_case[match(imp_var,drop_case$to),] %>% pull(corr)
        names(scores_imp) <- gsub("_.*", "",imp_var)

        # Get Y
        if(!is.null(datY)){
          scores_impY <- get_Y_imp(net = match_old_net, tree.membership = tree.membership, dat = datY)
          scores_imp <- list(X = scores_imp, Y = scores_impY)
        } else {scores_imp <- list(X = scores_imp)}

      }

      new_case <- updated_net[updated_net$to %in%  drop_case$from,]
      new_case <- new_case[match(drop_case$from, new_case$to),]

      updated_net2 <- filter(updated_net, !mem_id_old %in% drop_case$mem_id_old)

      updated_net2[match(drop_case$from, updated_net2$to), "is_leaf"] <- 0
      updated_net2[match(new_case$to, updated_net2$from), c("from", "is_leaf")] <- cbind(new_case$from, rep(1,nrow(drop_case)))
      updated_net2 <-  filter(updated_net2, !to %in% drop_case$from)

      dat_update <- dat[!tree.membership %in% drop_case$mem_id_old,]

      classified <- list(
        dat = dat[tree.membership %in% drop_case$mem_id_old,],
        mem = tree.membership[tree.membership %in% drop_case$mem_id_old]
      )
      m_idx <- which(!tree.membership %in% drop_case$mem_id_old)
      names(tree.membership) <- update$tree.mem
      m_case <- update$tree.mem
      sapply(1:length(unique(drop_case$mem_id_old)), function(id){
        m_case[names(tree.membership) %in% unique(drop_case$mem_id_old)[id]] <<- tree.membership[names(tree.membership) %in% unique(drop_case$mem_id_old)[id]]
      })
      mem_new <- m_case[m_idx]

      updated_net <- updated_net2
      dat <- dat_update
      updated_net$is_leaf <- as.numeric(updated_net$is_leaf)


      updated_net$mem_id[updated_net$mem_id %in% drop_case$mem_id_old] <- updated_net$mem_id_old[updated_net$mem_id %in% drop_case$mem_id_old]

    } else {
      mem_new <- update$tree.mem

    }

  } else {

    mem_new <- update$tree.mem
  }



  return(
    list(
      net = updated_net,
      dat = dat,
      classified = classified,
      mem = mem_new,
      imp_var = scores_imp
    )
  )
}

# Update terminal nodes from bottom to top
update_iter_cl <- function(mod, tree.id, size_min = 10, symm = T, use = "X", calc_imp = T){

  net <- get_tree_net(mod = mod, tree.id = tree.id)
  mem <- mod$membership[,tree.id]
  dat <- NULL
  i <- 1
  imp <- list()
  class_ls <- list()

  while(any(is.null(net$is_leaf), sum(net$is_leaf) > 2)){

    update_ls <- update_tree_leaf(
      mod = mod,
      dat = dat,
      tree.membership = mem,
      net = net,
      size_min = size_min,
      symm = symm,
      use = use,
      calc_imp = calc_imp
    )

    net <- update_ls$net
    mem <- update_ls$mem
    dat <- update_ls$dat

    if(!is.null(update_ls$classified)){

      class_ls[[i]] <- update_ls$classified
      i <- i + 1
      if(calc_imp){
        imp[[i]] <- update_ls$imp_var
      }


    }
  }

  if(calc_imp){
    # Get root variable imp
    root_net <- net %>% filter(is_leaf == 1)
    imp_root <- abs(-1 - mean(root_net$corr))
    names(imp_root) <- gsub("^(.*)_.*", "\\1", unique(root_net$from))
    impX <- c(unlist(purrr::map(imp, "X") ), imp_root)

    if(!is.null(mod$yvar)){
      # Get root variable imp for Y
      root_net$mem_id_old <- root_net$mem_id
      datY <- mod$yvar[rownames(dat),]
      imp_rootY <- get_Y_imp(root_net, mem, datY)
      impY <- c(unlist(purrr::map(imp, "Y")), imp_rootY)
      imp_ls <- list(X = impX, Y = impY)
      var_name <- list(X = colnames(mod$xvar), Y = colnames(mod$yvar))
    } else {
      imp_ls <- list(X = impX)
      var_name <- list(X = colnames(mod$xvar))
    }

    if(!is.null(mod$yvar)){idx <- c("X","Y")} else {idx <- c("X")}
    # create var imp vector
    imp_col <- plyr::llply(
      idx,
      .fun = function(l){
        get_iv(var_name[[l]],imp_ls[[l]])
      }
    )
    names(imp_col) <- idx
  } else imp_col <- NULL

  if(any(table(mem)) < size_min){

    unique_mem <- unique(mem)
    mem[mem == unique_mem[1]] <- unique_mem[2]

  }

  dat_new <- purrr::map(class_ls, "dat")
  dat_new <- Reduce(rbind, dat_new)
  dat_new <-  rbind(dat_new, dat)

  class_d <- purrr::map(class_ls, "mem")
  class_d <- unlist(class_d)
  class_d <- c(class_d, mem)

  class_d <- class_d[match(rownames(mod$xvar), rownames(dat_new))]
  prox <- get_prox(class_d)

  return(
    list(class_mem = class_d,
         prox = prox,
         imp = imp_col)
  )
}

# Create new membership and new proximity matrix of forest
cl_forest <- function(mod, size_min = 5, use = "X", symm = F, calc_imp = T, parallel = T, cores = max(detectCores() - 2,20), ...){

  nt <- mod$ntree

  if(parallel){
    if(Sys.info()["sysname"] == "Windows"){
      cluster <- parallel::makeCluster(cores)
    } else {cluster <- cores}
    doParallel::registerDoParallel(cluster)
  }
  `%myinfix%` <- ifelse(parallel, `%dopar%`, `%do%`)

  forest_ls <- foreach(t = 1:nt, .errorhandling = 'remove') %myinfix% {
    update_iter_cl(mod,
                   tree.id = t,
                   size_min = size_min,
                   use = use,
                   symm = symm,
                   calc_imp = calc_imp)
  }
  nt <- length(forest_ls)

  if(is.null(mod$yvar)) {idx <- "X"} else {idx <- c("X", "Y")}
  imp_var <- plyr::llply(
    idx,
    .fun = function(im){
      df <- purrr::map(forest_ls, "imp")
      df <- purrr::map(df, im)
      Reduce("+", df)/nt
    }
  )

  prox <- purrr::map(forest_ls, "prox")
  prox <- Reduce("+", prox)/nt

  new_mem <- purrr::map(forest_ls, "class_mem")
  new_mem <- Reduce(cbind, new_mem)


  return(list(
    prox = prox,
    imp = imp_var,
    membership.new = new_mem
  ))

}

# Transfer proximity to cosine distance
trans_cos_mat <- function(mat){
  p <- mat/sqrt(rowSums(mat^2))
  prox <- p %*% t(p)
}

# Screening function
mrf_XY.screening <- function(mod,...){

  fw <- mod$forest.wt

  ex <- colMeans((fw %*% as.matrix(mod$xvar) - as.matrix(mod$xvar))^2)
  if(is.null(mod$yvar)){
    out <- list(X = colnames(mod$xvar[,ex < 1]))
  } else {
    if(class(mod)[3] == "class+"){
      out <- list(X = colnames(mod$xvar[,ex < 1]))
    } else {
      ey <- colMeans((fw %*% as.matrix(mod$yvar) - as.matrix(mod$yvar))^2)
      out <- list(
        X = colnames(mod$xvar[,ex < 1]),
        Y = colnames(mod$yvar[,ey < 1])
      )
    }
  }
  

  out

}
