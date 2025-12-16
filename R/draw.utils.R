#' Extract Posterior Samples
#'
#' Extracts posterior samples from a fitted model object (Stan or similar).
#'
#' @param fit_obj A fitted model object. Can be of class `stanreg`, `CmdStanFit`, or similar.
#' @param chain.subset An integer indicating the chain to extract. Use 0 for all chains.
#'
#' @return A list of posterior samples, one element per variable, with class `"mcmc.draws"`.
#' @export
extract_samples <- function(fit_obj, chain.subset = 0) {
  es <- NULL
  r <- NULL

  if("stanreg" %in% class(fit_obj)){
    return( extract_samples_from_rstan( fit_obj, chain.subset = chain.subset ) )
  }

  if(typeof(fit_obj) == "list"){
    vars <- fit_obj$metadata$stan_variables
    draws <- posterior::as_draws_rvars(fit_obj$post_warmup_draws)
    es <- lapply(vars, \(var_name){ posterior::draws_of(draws[[var_name]], with_chains = FALSE) }) |> setNames(vars)
  }else{
    vars <- fit_obj$metadata()$stan_variables
    draws <- posterior::as_draws_rvars(fit_obj$draws())
    es <- lapply(vars, \(var_name){ posterior::draws_of(draws[[var_name]], with_chains = FALSE) }) |> setNames(vars)
  }

  if(chain.subset == 0) r <- es
  else{
    nc <- fit_obj$num_chains()
    stopifnot(chain.subset <= nc)

    ess <- list()
    k <- 1
    whi <- chain.idx(fit_obj) == chain.subset
    for(s in es){
      ess[[k]] <- s[whi, ]
      k <- k+1
    }
    names(ess) <- names( es )
    r <- ess
  }

  class(r) <- "mcmc.draws"
  return( r )
}

#' Extract Samples from rstan Object
#'
#' Extracts posterior samples from a `stanreg` object from the rstan package.
#'
#' @param fit_obj An object of class `stanreg`.
#' @param chain.subset An integer indicating the chain to extract. Use 0 for all chains.
#'
#' @return A list of samples with class `"mcmc.draws"`.
#' @export
extract_samples_from_rstan <- function(fit_obj, chain.subset = 0) {
  es <- NULL
  r <- NULL

  m <- as.matrix( fit_obj )
  es <- list()
  for(j in 1:ncol(m)) es[[j]] <- m[,j]
  names( es ) <- colnames( m )

  if(chain.subset == 0) r <- es
  else{
    nc <- fit_obj$num_chains()
    stopifnot(chain.subset <= nc)

    ess <- list()
    k <- 1
    whi <- chain.idx(fit_obj) == chain.subset
    for(s in es){
      ess[[k]] <- s[whi, ]
      k <- k+1
    }
    names(ess) <- names( es )
    r <- ess
  }

  class(r) <- "mcmc.draws"
  return( r )
}

#' Create Density Object with Optional Normalization
#'
#' Computes the kernel density estimate for a numeric vector with optional normalization.
#'
#' @param x A numeric vector.
#' @param normalize Logical or numeric. If `TRUE`, normalize to peak of 1. If numeric, compute and store but do not apply.
#' @param ... Additional arguments passed to [stats::density()].
#'
#' @return A list of class `"dens"` representing the density object.
#' @export
dens <- function(x,
                 normalize = NULL,
                 ...){



  d <- density(x, ...)
  d$yorig <- d$y
  if(is.null(normalize)) d$y <- d$yorig
  else if(is.logical(normalize) && normalize == T) d$y <- d$yorig / max(d$yorig)
  else if(is.logical(normalize) && normalize == F) d$y <- d$yorig
  else if(is.numeric(normalize)) d$yn <- d$yorig / max(d$yorig)

  class(d) <- "dens"
  d
}

#' Evaluate Density at Specific Locations
#'
#' Returns the approximate density at specified x values.
#'
#' @param d A density object created by [dens()].
#' @param x.at Numeric vector of x-locations at which to evaluate the density.
#'
#' @return A numeric vector of density values.
#' @export
dens.at <- function(d, x.at){
  if(missing(d)) stop("No density object specified.")
  if(missing(x.at)) stop("No location to approximate density specified.")

  dens_appr <- approxfun(d$x, d$y)
  dens.at <- dens_appr(x.at)

  return( dens.at )
}

#' Log-Sum-Exp Trick
#'
#' Computes the log of the sum of exponentials in a numerically stable way.
#'
#' @param x A numeric vector.
#'
#' @return A single numeric value: \eqn{\log(\sum \exp(x_i))}.
#' @export
log.sum.exp <- function(x){
  is.ign <- which( is.nan(x) | is.infinite(x)  )
  if( length(is.ign) > 0){
    warning("NaNs or Inf ignored.")
    x <- x[ -is.ign ]
  }

  xmax <- max(x)
  xsum <- sum(exp(x - xmax))
  xmax + log(xsum)
}

#' Normalize Log Posterior Values
#'
#' Normalize a log posterior vector using the log-sum-exp trick.
#'
#' @param x A numeric vector of log posterior values.
#' @param scale.to.prob Logical; if `TRUE`, exponentiate to return normalized probabilities.
#'
#' @return A normalized log posterior vector or probability vector.
#' @export
normalize.log.post <- function(x, scale.to.prob = FALSE){
  xs <- x - log.sum.exp(x)
  if(scale.to.prob) exp( xs )
  else return(xs)
}

#' Extract Variable Names from Draws
#'
#' Retrieves a flat list of variable names including array indices.
#'
#' @param draws A list of posterior samples.
#' @param ... Additional draw lists to include.
#' @param unique.names Logical; if `TRUE`, returns unique variable names.
#'
#' @return A character vector of variable names.
#' @export
draws.varnames <- function(draws, ..., unique.names = TRUE){

  n <- character(0)
  for(v in names(draws)){
    dv <- dim(draws[[v]])
    if(is.null(dv)){
      n <- append(n, v)
    }
    if(length(dv) == 2){
      if(dv[2] == 1) n <- append(n, v)
      if(dv[2] > 1){
        for(i in 1:dv[2]) n <- append(n, paste0(v, "[", i ,"]"))
      }
    }

    if(length(dv) == 3){
      for(i in 1:dv[2]){
        for(j in 1:dv[3]){
          n <- append(n, paste0(v, "[", i, ",", j ,"]"))
        }
      }
    }

  }

  dadds <- list(...)
  if(length( dadds ) > 0){
    for(d in dadds) n <- c(n, draws.varnames( d ))
  }

  n <- if(unique.names){
    unique(n)
  }else{
    n
  }

  return(n)
}

numdraws <- function(fit_obj){
  # if x is a fitobj returned by CmdStanR
  if("CmdStanFit" %in% class(fit_obj)){
    nd <- posterior::ndraws(fit_obj$draws())
  }else if("stanreg" %in% class(fit_obj)){
    nd <-  nrow( as.data.frame(fit_obj))
  }
  # if of S3 class mcmc.draws
  else if("mcmc.draws" %in% class(fit_obj)){
    nd <- dim(draws[[1]])[1]
  }
  nd
}

numchains <- function(fit_obj){
  # if x is a fitobj returned by CmdStanR
  if("CmdStanFit" %in% class(fit_obj)){
    nchains <- posterior::nchains(fit_obj$draws())
  }else if("stanreg" %in% class(fit_obj)){
    dd <-posterior::as_draws_df( fit_obj$stanfit )
    nchains <-  length( unique( dd$.chain ))
  }
  nchains
}


#' Generate Chain Indices
#'
#' Returns the chain index for each draw in a posterior sample.
#'
#' @param fit_obj A CmdStanFit object or similar.
#'
#' @return An integer vector indicating the chain each draw belongs to.
#' @export
chain.idx <- function(fit_obj){

  nd <- numdraws( fit_obj )
  nc <- numchains( fit_obj)

  dpc <- nd / nc
  rep(1:nc, each = dpc)
}

#' Internal: Gather MCMC Output; not meant to be called directly.
#'
#' Gathers and summarizes posterior samples from a fitted object or draws list.
#'
#' @param fit_obj A fitted model object.
#' @param draws Optional. A list of draws (if not using fit_obj).
#' @param vars Character vector of variables to include. If `NULL`, includes all.
#' @param gather Character vector; one or more of `"matrix"`, `"summary"`, `"densities"`.
#' @param ignore.lp Logical; if `TRUE`, drops `lp__` from draws.
#' @param PI.lvls Numeric vector of probability intervals (e.g. `c(0.5, 0.9)`).
#'
#' @return A list containing matrices and/or summaries depending on `gather`.
#' @keywords internal
..mcmc.gather <- function(fit_obj, draws = NULL, vars = NULL, gather = c("matrix", "summary", "densities"), ignore.lp = FALSE, PI.lvls = c(.5, .9)){
  if(missing(fit_obj) && is.null(draws)) stop("No CmStanFit or list of samples provided.")
  if(is.null(draws))  draws <- extract_samples( fit_obj )


  if(ignore.lp) draws$lp__ <- NULL

  varsp <- NULL
  if(is.null(vars)){
    varsp <- draws.varnames(draws)
  }else if(!is.null(vars) && is.character(vars)){
    varsp <- vars
  }

  nvars <- length( varsp )
  nsams <- dim(draws[[1]])[1]

  gather_mat <- "matrix" %in% gather
  gather_dens <- "densities" %in% gather
  gather_sum <- "summary" %in% gather

  PI.lvls.or <- PI.lvls

  if(gather_mat){
    X <- matrix(nr = nsams, nc = nvars)
    colnames(X) <- varsp
  }
  if(gather_sum){
    dfs <- data.frame(
      mean   = rep(NA, nvars ),
      sd     = rep(NA, nvars ),
      median = rep(NA,nvars)
    )


    if(is.null(PI.lvls)) PI.lvls <- 0.5

    PInames <- sort( apply(expand.grid(100*PI.lvls, c(".lwr", ".upr")),1, \(e) paste0("PI", e[1], e[2])))
    PImat <- matrix(nc = 2*length(PI.lvls), nr = nvars)
    colnames(PImat) <- PInames
    dfs <- cbind(dfs, PImat)

    rownames(dfs) <- varsp
  }
  #TODO gather densities or consider ditching this idea


  for(v in names(draws)){

    dv <- dim(draws[[v]])

    if(is.null(dv)){
      if(v %in% varsp){
        s <- draws[[v]]
        k <- which(v == varsp)
        if(gather_mat){
          X[,k] <- s
        }

        if(gather_sum){
          dfs[k,]$mean <- mean( s )
          dfs[k,]$sd   <- sdN( s )
          dfs[k,]$median   <- median( s )

          for(p in PI.lvls){
            np <- c( paste0("PI", 100 * p, ".lwr"),
                     paste0("PI", 100 * p, ".upr")
            )

            pii <- PI(s, p)
            dfs[k, np[1] ] <- pii[1]
            dfs[k, np[2] ] <- pii[2]
          }
        }
      }
    }

    if(length(dv) == 2){
      if(dv[2] == 1){
        if(v %in% varsp){

          s <- draws[[v]]
          k <- which(v == varsp)
          if(gather_mat){
            X[,k] <- s
          }

          if(gather_sum){

            dfs[k,]$mean <- mean( s )
            dfs[k,]$sd   <- sdN( s )
            dfs[k,]$median   <- median( s )

            for(p in PI.lvls){
              np <- c( paste0("PI", 100 * p, ".lwr"),
                       paste0("PI", 100 * p, ".upr")
              )

              pii <- PI(s, p)
              dfs[k, np[1] ] <- pii[1]
              dfs[k, np[2] ] <- pii[2]
            }
          }

        }

      }
      if(dv[2] > 1){
        for(i in 1:dv[2]){
          n <-  paste0(v, "[", i ,"]")
          if(n %in% varsp){
            s <- draws[[v]][,i]
            k <- which(n == varsp)
            if(gather_mat){
              X[,k] <- s
            }

            if(gather_sum){
              dfs[k,]$mean <- mean( s )
              dfs[k,]$sd   <- sdN( s )
              dfs[k,]$median   <- median( s )

              for(p in PI.lvls){
                np <- c( paste0("PI", 100 * p, ".lwr"),
                         paste0("PI", 100 * p, ".upr")
                )


                pii <- PI(s, p)
                dfs[k, np[1] ] <- pii[1]
                dfs[k, np[2] ] <- pii[2]
              }

            }

          }
        }
      }
    }

    if(length(dv) == 3){
      for(i in 1:dv[2]){
        for(j in 1:dv[3]){
          n <- paste0(v, "[", i, ",", j ,"]")
          if(n %in% varsp){
            s <- draws[[v]][,i,j]
            k <- which(n == varsp)

            if(gather_mat){
              X[,k] <- s
            }

            if(gather_sum){
              dfs[k,]$mean <- mean( s )
              dfs[k,]$sd   <- sdN( s )
              dfs[k,]$median   <- median( s )

              for(p in PI.lvls){
                np <- c( paste0("PI", 100 * p, ".lwr"),
                         paste0("PI", 100 * p, ".upr")
                )

                pii <- PI(s, p)
                dfs[k, np[1] ] <- pii[1]
                dfs[k, np[2] ] <- pii[2]
              }
            }

          }
        }
      }
    }
  }

  gather <- list()
  if(gather_mat) gather[["draws_mat"]] <- X
  if(gather_sum){
    if(is.null(PI.lvls.or)){ dfs$PI50.lwr <- NULL;dfs$PI50.upr <- NULL}
    gather[["summary"]] <- dfs
  }

  return( gather )
}

#' Summarize Posterior Draws
#'
#' Provides summary statistics (mean, median, SD, intervals) of posterior samples.
#'
#' @param fit_obj A fitted model object.
#' @param draws Optional. A list of posterior draws.
#' @param vars Optional character vector of variables to include.
#' @param ignore.lp Logical; if `TRUE`, removes `lp__` from summary.
#' @param digits Integer; number of decimal places to round.
#' @param PI.lvls Numeric vector of probability intervals.
#'
#' @return A data frame of summary statistics.
#' @export
mcmc.summary <- function(fit_obj, draws = NULL, vars = NULL, ignore.lp = FALSE, digits = 3, PI.lvls = c(.5, .9)){
  round( ..mcmc.gather(fit_obj = fit_obj,
                     draws = draws,
              vars = vars,
              ignore.lp = ignore.lp,
              PI.lvls = PI.lvls,
              gather = "summary")$summary, digits)
}

#' Extract Draws as Matrix
#'
#' Returns posterior samples as a matrix for selected variables.
#'
#' @param draws A list of posterior samples.
#' @param vars Optional character vector of variables to extract.
#' @param ignore.lp Logical; if `TRUE`, removes `lp__`.
#'
#' @return A numeric matrix of posterior draws.
#' @export
mcmc.drawsmat <- function(draws, vars = NULL, ignore.lp = FALSE){
  ..mcmc.gather(fit_obj = draws,
              vars = vars,
              ignore.lp = ignore.lp,
              gather = "matrix")$draws_mat
}
