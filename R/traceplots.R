#' Function for producing a simple traceplot
#'
#' @param x A vector of draws
#' @param lab Change y-axis label
#' @param col Color of trace line
#'
#' @returns A ggplot2 line plot object
#' @export
#'
#' @examples
simple_traceplot <- function( x, lab = "x", col = bu.color( 1 )) {

  p <- ggplot2::ggplot( data.frame(iteration = 1:length(x),
                     x = x),
  ggplot2::aes(x = iteration, y = x )) +
  ggplot2::geom_line( col = col ) +
  ggplot2::theme_minimal(  ) +
  ggplot2::ylab( lab )

  return( p )
}


if( F ){
  d <- data.frame(x = iris$Sepal.Length, y = iris$Sepal.Width)
  f <- rstanarm::stan_glm(y ~ x, data = d, chains = 1)

  dr <- extract_samples_from_rstan( f )
  chain.idx( f )
  simple_traceplot( dr$x )
}
