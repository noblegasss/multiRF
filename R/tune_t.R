#' Tune reconstruction parameter t
#' @rdname tune_nearest_coef
#' @export

tune_nearest_coef <- function(dat.list, mod, tmin = 10, by = 5, k = NULL, object = "diss"){

  mod_init <- mrf3_cl(
    mod = mod,
    dat.list,
    nearest_coef = nrow(dat.list[[1]]),
    k = k
  )
  obj_init <- mod_init$cl_mod$obj
  sil_init <- mod_init$cl_mod$sil
  diff_init <- max(mod_init$cl_mod$diff_e)
  
  
  t_max <- seq(tmin, nrow(dat.list[[1]]), by = by)

  k0 <- c()

  obj <- sil <- diff <- c()
  i <- 1

  while(i <= length(t_max)){
    
    suppressMessages(
      m <- mrf3_cl(
        dat.list,
        mod = mod,
        nearest_coef = t_max[i],
        k = k
      )
    )
    
    obj <- c(obj, m$cl_mod$obj)
    sil <- c(sil, m$cl_mod$sil)
    k0 <- c(k0, length(unique(m$cl)))
    diff <- c(diff, max(m$cl_mod$diff_e))
   
    i <- i + 1
  } 
  
  diff <- c(diff, diff_init)
  k_all <- c(k0, length(unique( mod_init$cl)))
  df <- data.frame(obj = c(obj,obj_init), 
                   sil = c(sil, sil_init),
                   diffe = diff,
                   nearest_coef = c(t_max, nrow(dat.list[[1]])),
                   k = k_all)
  
  if(object == "diss") {
    dfsumm <- df %>% 
      group_by(k) %>%
      slice_min(obj, with_ties = F)
  } 
  if(object == "silhouette"){
    dfsumm <- df %>% 
      group_by(k) %>%
      slice_max(sil, with_ties = F)
  }
  if(object == "eigen") {
    dfsumm <- df %>% 
      group_by(k) %>%
      slice_max(diffe, with_ties = F)
  }
  
  
  list(tmax_tb = dfsumm,
       object = df)
}
