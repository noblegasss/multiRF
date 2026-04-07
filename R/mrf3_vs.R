#' MRF variable selection
#'
#' @param mod A fitted `mrf3_fit` object, or an `mrf3`-like object that already
#' contains IMD weights in `mod$imd`.
#' @param dat.list A list of omics matrices used for feature selection/refit.
#' @param method Feature-selection rule: `"filter"`, `"test"`, `"mixture"`, or `"thres"`.
#' @param signal Which signal component to select on: `"shared"` (cross-modal),
#'   `"specific"` (per-block residual), or `"all"` (both, returned as a list).
#' @param se Multiplier on standard deviation for thresholding in `"thres"` mode.
#' @param c1 Distribution family for component 1 in mixture mode.
#' @param c2 Distribution family for component 2 in mixture mode.
#' @param level Significance level for test/mixture selection. Can be:
#'   `"auto"` — adaptive Bayesian FDR control (recommended for mixture);
#'   a single numeric (applied to all blocks); or a named vector/list for
#'   per-block thresholds, e.g. `c(mri = 0.2, dnam = 0.01)`.
#' @param tscore Logical; whether to use test scores in mixture mode.
#' @param use_distribution Logical; whether to use distributional approximation in thresholding.
#' @param re_weights Logical; whether to pass selected weights into refitting.
#' @param re_fit Logical; whether to refit forests after feature selection.
#' @param ntree Number of trees used in optional refit.
#' @param scale Logical; whether to standardize data before refit.
#' @param k Number of folds/repeats used by filtering heuristics.
#' @param tol Improvement tolerance in adaptive threshold search.
#' @param iter Maximum EM iterations for mixture mode.
#' @param eps EM convergence tolerance for mixture mode.
#' @param normalized Logical; whether to renormalize selected weights.
#' @param select Data block(s) to apply selection to; default `"ALL"`.
#' @param ... Additional arguments passed to underlying fitting functions.
#' @export
mrf3_vs <- function(mod,
                    dat.list = NULL,
                    method = "filter",
                    signal = c("shared", "specific", "all"),
                    se = 1,
                    c1 = "normal",
                    c2 = "normal",
                    level = 0.05,
                    tscore = FALSE,
                    use_distribution = TRUE,
                    re_weights = FALSE,
                    re_fit = TRUE,
                    ntree = 300,
                    scale = FALSE,
                    k = 3,
                    tol = 0.01,
                    iter = 1000,
                    eps = 1e-05,
                    normalized = TRUE,
                    select = "ALL",
                    ...){

  signal <- match.arg(signal)

  ## ---- signal = "all": run both shared and specific, return combined ----
  if (signal == "all") {
    message("Running shared + specific variable selection...")
    cl <- match.call()
    cl$signal <- "shared"
    shared_res <- eval(cl, parent.frame())
    cl$signal <- "specific"
    specific_res <- eval(cl, parent.frame())
    return(list(shared = shared_res, specific = specific_res))
  }

  if (inherits(mod, "mrf3_fit")) {
    wf <- mod
    if (is.null(dat.list)) {
      dat.list <- wf$data
    }
    if (is.null(dat.list)) {
      stop(
        "`mrf3_vs()` with `mrf3_fit` requires data in `wf$data`. ",
        "Run `mrf3_fit(..., return_data = TRUE)` or provide `dat.list` explicitly."
      )
    }

    if (signal == "specific") {
      ## ---- Specific signal: use residual-RF IMD ----
      spec <- wf$specific
      if (is.null(spec) || is.null(spec$imd)) {
        stop(
          "`mrf3_vs(signal = 'specific')` requires specific IMD weights ",
          "in `wf$specific$imd`. Ensure shared-specific decomposition was run."
        )
      }
      wf_weights <- spec$imd
      ## Build imd_ls from specific imd_per_tree (p x ntree per block)
      ## Format needed by test_fn: list of 1 pseudo-connection, each = list
      ## of ntree elements, each = named list(block = vec)
      wf_weights_ls <- NULL
      spec_pt <- spec$imd_per_tree
      if (!is.null(spec_pt) && length(spec_pt) > 0L) {
        block_names <- names(spec_pt)
        ## Each block's unsupervised RF is independent, so wrap each as a
        ## separate pseudo-connection so test_fn processes them per-block
        wf_weights_ls <- lapply(block_names, function(bn) {
          pt_mat <- spec_pt[[bn]]  # p x ntree
          if (is.null(pt_mat)) return(NULL)
          nt <- ncol(pt_mat)
          lapply(seq_len(nt), function(t) {
            out <- list(pt_mat[, t])
            names(out) <- bn
            out
          })
        })
        wf_weights_ls <- Filter(Negate(is.null), wf_weights_ls)
      }
    } else {
      ## ---- Shared signal (original path) ----
      wf_weights <- wf$imd
      wf_weights_ls <- wf$imd_init
    }

    if (is.null(wf_weights)) {
      stop(
        "`mrf3_vs()` with `mrf3_fit` requires IMD weights. ",
        "Run `mrf3_fit(..., run_imd = TRUE)` first."
      )
    }

    ## Build connection list for specific: one pseudo-connection per block
    if (signal == "specific") {
      connect_list_vs <- lapply(names(wf_weights), function(bn) c(bn, bn))
    } else {
      connect_list_vs <- wf$connection
    }

    mod <- list(
      imd = wf_weights,
      imd_ls = wf_weights_ls,
      connection = connect_list_vs,
      ytry = wf$config$ytry,
      ntree = wf$config$ntree,
      type = wf$type,
      oob_err = wf$oob_err,
      mod = wf$models
    )
    class(mod) <- "mrf3"
  }

  if (is.null(dat.list)) {
    stop("`dat.list` must be provided.")
  }
  
  weights <- mod$imd
  if (is.null(weights)) {
    stop(
      "`mrf3_vs()` requires IMD weights in `mod$imd`. ",
      "`mrf3_init()` no longer computes IMD weights; run `mrf3_fit(..., run_imd = TRUE)` ",
      "or construct `mod` with weights from `get_multi_weights()`."
    )
  }
  new_dat <- dat.list
  dat_names <- names(new_dat)
  connect_list <- mod$connection

  is_all_selected <- identical(select, "ALL")
  normalize_weight_vec <- function(w) {
    if (!isTRUE(normalized)) {
      return(w)
    }
    norm_w <- sqrt(sum(w^2))
    if (!is.finite(norm_w) || norm_w <= 0) {
      return(w)
    }
    w / norm_w
  }
  keep_feature_names <- function(dat_block, weight_vec, block_name, allow_empty = FALSE) {
    cols <- colnames(dat_block)
    if (is.null(cols)) {
      keep_idx <- which(weight_vec > 0)
      if (!length(keep_idx) && !allow_empty && length(weight_vec)) {
        keep_idx <- which.max(weight_vec)
      }
      keep_idx <- keep_idx[keep_idx <= ncol(dat_block)]
      return(cols[keep_idx])
    }

    if (is.null(names(weight_vec))) {
      keep_idx <- which(weight_vec > 0)
      if (!length(keep_idx) && !allow_empty && length(weight_vec)) {
        keep_idx <- which.max(weight_vec)
      }
      keep_idx <- keep_idx[keep_idx <= length(cols)]
      return(cols[keep_idx])
    }

    keep_names <- intersect(names(weight_vec)[weight_vec > 0], cols)
    if (!length(keep_names) && !allow_empty) {
      best_name <- names(weight_vec)[which.max(weight_vec)]
      if (!is.na(best_name) && nzchar(best_name) && best_name %in% cols) {
        keep_names <- best_name
      } else if (length(cols)) {
        keep_names <- cols[[1]]
      }
    }
    keep_names
  }
  subset_block <- function(dat_block, weight_vec, block_name, allow_empty = FALSE) {
    if (!is_all_selected && !block_name %in% select) {
      return(dat_block)
    }
    keep_names <- keep_feature_names(dat_block, weight_vec, block_name, allow_empty = allow_empty)
    if (!length(keep_names)) {
      return(dat_block[, 0, drop = FALSE])
    }
    dat_block[, keep_names, drop = FALSE]
  }

  message("Variable selection..")

  if(method == "thres") {

    thres <- chooss_thres3(
      mod$imd,
      se = se
    )
    
  }
  
  
  if(method == 'test'|tscore) {

    if (signal == "specific") {
      ## Direct per-block z-test for specific signal (no cross-modal connections)
      imd_ls <- mod$imd_ls
      if (is.null(imd_ls)) {
        warning(
          "`method = 'test'` for specific signal requires per-tree IMD. ",
          "Falling back to `method = 'mixture'`.",
          call. = FALSE
        )
        method <- "mixture"
      } else {
        block_names <- names(mod$imd)
        thres_list <- lapply(seq_along(block_names), function(i) {
          bn <- block_names[i]
          w <- imd_ls[[i]]
          raw <- purrr::map(w, bn)
          raw <- raw[vapply(raw, length, integer(1)) > 0L]
          if (length(raw) == 0L) return(list(keep_idx = NULL, pval = NULL, ts = NULL))
          mat <- do.call(cbind, raw)  # p x ntree
          mu <- mean(mat)
          se <- apply(mat, 1, function(k) {
            if (sum(k != 0) > 1) sd(k) / sqrt(length(k) - 1) else -1
          })
          z <- ifelse(se != -1, (rowMeans(mat) - mu) / se, 0)
          p <- 1 - pnorm(z)
          sig <- if (is.numeric(level)) level[1] else 0.05
          keep <- ifelse(p < sig, 1, 0)
          list(keep_idx = keep, pval = p, ts = z)
        })
        names(thres_list) <- block_names
        thres <- list(
          keep_idx = lapply(thres_list, `[[`, "keep_idx"),
          pval = lapply(thres_list, `[[`, "pval"),
          ts = lapply(thres_list, `[[`, "ts")
        )
      }
    } else if (is.null(mod$imd_ls)) {
      warning(
        "`method = 'test'` requires per-tree weight distributions (unavailable ",
        "with native engine pre-computed IMD). Falling back to `method = 'mixture'`.",
        call. = FALSE
      )
      method <- "mixture"
    } else {
      thres <- test_fn(
        wl = mod$imd_ls,
        connection = connect_list,
        dat_names = dat_names,
        sig.thres = level)
    }
  }
  
  if(method == "filter" && signal == "specific") {
    warning(
      "`method = 'filter'` is not supported for specific signal. ",
      "Falling back to `method = 'mixture'`.",
      call. = FALSE
    )
    method <- "mixture"
  }

  if(method == "mixture") {
    if(tscore) ts <- thres$ts

    thres <- plyr::llply(
      names(weights),
      .fun = function(t) {
        if(tscore) {
          ww0 <- ts[[t]]
        } else {
          ww0 <- weights[[t]]
        }
        
        param <- get_param(ww0)
        e <- em(ww0, p = param$p, mu = param$mu, sigma = param$sigma, iter = iter, eps = eps, c1 = c1, c2 = c2)
        cbind(e$par$post[,1]/rowSums(e$par$post),
              e$par$post[,2]/rowSums(e$par$post),
              e$par$post[,3]/rowSums(e$par$post))
        
      }
    )
    
    names(thres) <- names(weights)

  }
  
  if(method == "filter") {
    weights_norm <- purrr::map(weights, ~./sqrt(sum(.^2)))
    thres <- choose_thres2(
      weights,
      connection = connect_list,
      new_dat = new_dat,
      ytry = mod$ytry,
      ntree = mod$ntree,
      type = mod$type,
      oob_init = mod$oob_err,
      k = k,
      select = select,
      tol = tol,
      ...)

  }

  
  weights_new <- plyr::llply(
    dat_names,
    .fun = function(i) {
      w <- weights[[i]]
      if(!method %in% c("mixture", "test") ) {
        if(select != 'ALL') {
          
          if(i %in% select) {
            w[w < thres[i]] <- 0
            w[w > thres[i]] <- w[w > thres[i]] - thres[i]
          } 
        } else {
          w[w < thres[i]] <- 0
          w[w > thres[i]] <- w[w > thres[i]] - thres[i]
        }
      } else {
        if(method == "mixture") {
          lfdr <- thres[[i]][,1]  # local FDR = P(noise | data)
          if (identical(level, "auto")) {
            ## Bayesian FDR control: sort by local FDR, find largest
            ## set where cumulative mean lfdr <= target_fdr
            target_fdr <- 0.05
            nonzero <- which(lfdr < 1 - 1e-8)  # skip exact-zero features
            if (length(nonzero) > 0) {
              ord <- order(lfdr[nonzero])
              cum_fdr <- cumsum(lfdr[nonzero][ord]) / seq_along(ord)
              n_sel <- max(0L, sum(cum_fdr <= target_fdr))
              lvl_i <- if (n_sel > 0) lfdr[nonzero][ord][n_sel] + 1e-12 else 0
            } else {
              n_sel <- 0L
              lvl_i <- 0
            }
            if (n_sel == 0L) {
              message(sprintf("  [mixture auto] %s: no features at FDR<=%.2f; ",
                              i, target_fdr),
                      "selecting top feature by lowest local FDR")
              best <- which.min(lfdr)
              t <- rep(0L, length(lfdr))
              t[best] <- 1L
            } else {
              message(sprintf("  [mixture auto] %s: %d features at FDR<=%.2f (lfdr cutoff=%.4f)",
                              i, n_sel, target_fdr, lvl_i))
              t <- ifelse(lfdr < lvl_i, 1, 0)
            }
          } else {
            lvl_i <- if (is.null(names(level))) level[1] else {
              if (i %in% names(level)) level[[i]] else level[1]
            }
            t <- ifelse(lfdr < lvl_i, 1, 0)
          }
        } else {
          t <- thres$keep_idx[[i]]
          
        }

        if(select != 'ALL') {
          
          if(i %in% select) {
            
            w <- w * t
          } 
        } else {

          w <- w * t
        }
        
      }
      
      
      normalize_weight_vec(w)
    }
  )

  names(weights_new) <- dat_names
  mod$imd <- weights_new
  if (isTRUE(re_fit)) {
    message("Refit model..")
  }
    
  new_dat <- plyr::llply(
    dat_names,
    .fun = function(i) subset_block(new_dat[[i]], weights_new[[i]], i)
  )
  names(new_dat) <- dat_names
  m <- purrr::map(weights_new, ~.[.>0])
  
  if(scale) {
    new_dat2 <- purrr::map(new_dat, ~scale(.))
  } else {
    new_dat2 <- new_dat
  }
    
  if(!re_weights) {
    m <- NULL
  } 
  if(re_fit && signal == "specific") {
    message("Skipping refit for specific signal (self-connections not meaningful for supervised RF).")
    re_fit <- FALSE
  }
  if(re_fit) {
    refit <- fit_multi_forest(new_dat2, connect_list = connect_list, ntree = ntree,
                             type = mod$type, var.wt = m, ytry = mod$ytry, ...)

    oob_err <- purrr::map(refit, ~get_r_sq(.))
    oob_err <- Reduce("+", oob_err)

    mod$oob_err <- oob_err
    mod$mod <- refit
  }
  
  mod$dat.list <- new_dat2
  mod$thres <- thres
  mod$ntree <- ntree
  
  class(mod) <- c(class(mod), "vs")
  
  mod
  
}

choose_thres2 <- function(weights, connection, new_dat, ytry, ntree, type, oob_init, k = 3, tol = 0.01, select = "ALL", ...) {
  
  num <- seq(.8,3.1,by = 0.1)
  oob <- c()
  i <- 1
  is_all_selected <- identical(select, "ALL")
  apply_threshold <- function(block_name, w, threshold) {
    if (!is_all_selected && !block_name %in% select) {
      return(w)
    }
    w[w < threshold] <- 0
    w[w > threshold] <- w[w > threshold] - threshold
    w
  }
  subset_block <- function(dat_block, weight_vec, block_name) {
    if (!is_all_selected && !block_name %in% select) {
      return(dat_block)
    }
    cols <- colnames(dat_block)
    if (is.null(names(weight_vec))) {
      keep_idx <- which(weight_vec > 0)
      if (!length(keep_idx) && length(weight_vec)) {
        keep_idx <- which.max(weight_vec)
      }
      keep_idx <- keep_idx[keep_idx <= ncol(dat_block)]
      return(dat_block[, keep_idx, drop = FALSE])
    }
    keep_names <- intersect(names(weight_vec)[weight_vec > 0], cols)
    if (!length(keep_names) && length(cols)) {
      best_name <- names(weight_vec)[which.max(weight_vec)]
      keep_names <- if (!is.na(best_name) && best_name %in% cols) best_name else cols[[1]]
    }
    dat_block[, keep_names, drop = FALSE]
  }

  while(i <= length(num)) {
    t <- num[i]

    weights_new <- purrr::imap(weights, ~apply_threshold(.y, .x, t * stats::sd(.x)))
    
    new_dat2 <- plyr::llply(
      names(new_dat),
      .fun = function(i) subset_block(new_dat[[i]], weights_new[[i]], i)
    )
    names(new_dat2) <- names(new_dat)

    oob_err0 <- plyr::laply(1:k,
                            .fun = function(k) {
                              refit <- fit_multi_forest(new_dat2, connect_list = connection,
                                                       ntree = ntree, type = type, ytry = ytry,
                                                       var.wt = NULL,
                                                       forest.wt = "oob")

                              oob_err <- purrr::map(refit, ~get_r_sq(.))
                              Reduce("+", oob_err)/length(new_dat2)
                            })
    
    oob_err <- mean(unlist(oob_err0))
    oob <- c(oob, oob_err)
    oob_init <- oob_err
    i <- i + 1
  }

  diff <- abs(diff(oob[-1]))
  t <- sort(num[-1][which(diff < tol)])[1]
  
  message("Choose ", t, " times sd")
  thres <- plyr::llply(
    weights,
    .fun = function(w) {
      
      t * sd(w)
    }
  )
  
  unlist(thres)
  
}

chooss_thres3 <- function(weights, se) {
  t <- plyr::laply(
    weights,
    .fun = function(w) {
      mean(w) + se * sd(w)
    }
  )
  names(t) <- names(weights)
  t
}

test_fn <- function(wl, connection, dat_names, sig.thres = 0.05) {

  res <- plyr::llply(1:length(wl),
              .fun = function(i) {
                w <- wl[[i]]
                conn <- rev(connection[[i]])
                pls <- plyr::llply(
                  names(w[[1]]),
                  .fun = function(j) {
                    raw_j <- purrr::map(w, j)
                    raw_j <- raw_j[vapply(raw_j, length, integer(1)) > 0L]
                    if (length(raw_j) == 0L) return(list(keep_idx = NULL, pval = NULL, ts = NULL))
                    mat <- Reduce(cbind, raw_j)

                    mu <- mean(mat)

                    se <- apply(mat, 1, function(k){
                      if(length(k[k!=0]) > 1) {
                        sd(k)/sqrt(length(k) - 1)
                      } else {
                        -1
                      }

                    } )

                    z <- ifelse(se != -1, (rowMeans(mat) - mu)/se, 0)

                    # p <- pt(z,df = nrow(mat) - 1 ,lower.tail = F)
                    p <- 1-pnorm(z)

                    keep <- ifelse(p < sig.thres, 1, 0)

                    list(keep_idx = keep,
                         pval = p,
                         ts = z)
                  }
                )

                p <- lapply(pls, `[[`, "pval")
                keep_idx <- lapply(pls, `[[`, "keep_idx")
                ts <- lapply(pls, `[[`, "ts")
                names(p) <- names(keep_idx) <- names(ts) <- conn


                list(keep_idx = keep_idx,
                     pval = p, 
                     ts = ts)
              })

  keep_idx <- purrr::map(res, "keep_idx")
  p <- purrr::map(res, "pval")
  ts <- purrr::map(res, "ts")

  keep_res <- plyr::llply(dat_names,
              .fun = function(d) {
                keep_raw <- purrr::map(keep_idx, d)
                keep_raw <- keep_raw[vapply(keep_raw, function(x) length(x) > 0L, logical(1))]
                if (length(keep_raw) == 0L) {
                  # d is a block name; get feature count from first non-empty keep_idx
                  n_feat <- length(purrr::compact(purrr::map(keep_idx, d))[[1]])
                  if (n_feat == 0L) n_feat <- 1L
                  return(list(keep_idx = rep(0L, n_feat), pval = NULL, ts = NULL))
                }
                keep <- Reduce(cbind, keep_raw)
                ts_raw <- purrr::map(ts, d)
                ts_raw <- ts_raw[vapply(ts_raw, function(x) length(x) > 0L, logical(1))]
                ts <- if (length(ts_raw) > 0L) Reduce(cbind, ts_raw) else NULL
                if(is.null(ncol(ts))) ts <- ts else ts <- colMeans(ts)


                if(!is.null(ncol(keep))) {
                  keep <- rowMeans(keep) >= 0.5
                } else {
                  keep <- keep == 1
                }
                keep_idx <- ifelse(keep, 1, 0)
                p_raw <- purrr::map(p, d)
                p_raw <- p_raw[vapply(p_raw, function(x) length(x) > 0L, logical(1))]
                p <- if (length(p_raw) > 0L) Reduce(cbind, p_raw) else NULL
                list(keep_idx = keep_idx,
                     pval = p,
                     ts = ts)

              })

  names(keep_res) <- dat_names
  keep_idx <- purrr::map(keep_res, "keep_idx")
  p <- purrr::map(keep_res, "pval")
  ts <- purrr::map(keep_res, "ts")

  list(keep_idx = keep_idx,
       pval = p,
       ts = ts)

}
