#' Plotting of empirical quantiles against quantiles of replicated data
#'
#' Given some data vector x and a function to simulate fake data, qqsim_plot
#' produces a set of fake data replicates and plots empirical data quantiles
#' against the mean fake data quantiles. A 90 \% percentile interval is added to the plot to
#' account for sampling variation. Alternatively, a matrix of pre-generated fake
#' data can be provided.
#'
#' @param x A vector of empirical data
#' @param sim_fun The function to simulate a vector of fake data. Defaults to sampling normal fake data with data mean and data std as sampling parameters.
#' @param Nsims The number of data simulation replicates.
#' @param Xhat A matrix of simulated data (Nsims times length(x)). When provided, the sim_fun argument is ignored.
#' @param produce_plot Should a plot actually be generated?
#' @param return_data If TRUE, underlying plotting data is returned.
#'
#' @returns NULL or a list of plotting data.
#' @export
#'
#' @examples
qqsim_plot <- function(x,
                       sim_fun = function(x) rnorm( length(x), mean(x), sd(x)),
                       Nsims = 100,
                       Xhat = NULL, produce_plot = TRUE, return_data = FALSE){

  if ( is.null( Xhat ) ) {
    Xhat <- t( replicate( Nsims, sim_fun(x) ))
  }

  sx <- sort( x )
  nhat <- ncol( Xhat )

  Qhat <- t( apply(Xhat, 1, function(xhat){
    qxhat <- sapply(sx, function( sxi) sum(xhat <= sxi) /nhat)
    return( qxhat )
  }))

  qx <- 1:length( sx ) / length( sx )
  mean_qhat <- apply(Qhat, 2, mean )
  pi_qhat <- apply(Qhat, 2, PI)

  if ( produce_plot ) {
    plot(1, 1, type = "n", xlim = c(0,1), ylim = c(0, 1), xlab = "emp. quantile", ylab = "sim. data quantile")
    polygon(x = c(qx, rev(qx)), y = c(pi_qhat[1,], rev( pi_qhat[2, ])),
            col = "skyblue",
            border = F)
    lines( qx, mean_qhat, col = "darkblue", lwd = 2)
    segments(x0 = 0, x1 = 1, y0 = 0, y1 = 1, col = "red")
  }

  if ( return_data ) {
    dr <- list(
      qx = qx,
      Xhat = Xhat,
      Qhat = Qhat,
      qhat_mean = mean_qhat,
      qhat_pi = pi_qhat,
      sim_fun = sim_fun
    )
    return( dr )
  }
}
