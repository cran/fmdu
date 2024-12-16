#' Ultra Fast Multidimensional Unfolding Function
#'
#' \code{ultrafastmds} performs simple (weighted) metric multidimensional unfolding.
#' The function follows algorithms given de Leeuw (1977), Agrafiotis (2003), Rajawat and Kumar (2017), and Busing (submitted).
#' The memory footprint is delta, x, and y, and w, fixed x and fixed y if present, all provided as input parameter.
#'
#' @param data an n by m dissimilarity matrix
#' @param x an n by p (p < m) initial row coordinates matrix (required).
#' @param y an m by p (p < m) initial column coordinates matrix (required).
#' @param w (optional) an n by m weights matrix
#' @param fx (optional) an n by p (p < m) fixed row coordinates indicator matrix (0=free;1=fixed) (optional).
#' @param fy (optional) an m by p (p < m) fixed column coordinates indicator matrix (0=free;1=fixed) (optional).
#' @param NSTEPS (optional) minimum number of learning rate steps (default = 4096).
#' @param RCRIT (optional) relative convergence criterion, i.e., lowest learning rate (default = 0.00000001).
#' @param seed (optional) seed passed to the C functions.
#'
#' @return x final n by p matrix with row coordinates.
#' @return y final m by p matrix with column coordinates.
#'
#' @references de Leeuw (1977).
#'             Application of convex analysis to multidimensional scaling.
#'             Recent developments in statistics, pp. 133-145, North-Holland.
#'
#'             Agrafiotis, D.K. (2003). Stochastic Proximity Embedding.
#'             Journal of computational chemistry, volume 24, number 10, pages 1215-1221. Wiley Online Library.
#'
#'             Rajawat, K. and Kumar, S. (2017). Stochastic Multidimensional Scaling.
#'             IEEE Transactions on Signal and Information Processing over Networks, Vol. 3, No. 2, June 2017.
#'
#'             Busing (submitted). Node Localization by Multidimensional Scaling with Iterative Majorization:
#'             An Overview of a Comprehensive Algorithm.
#'
#' @examples
#' \dontrun{
#' library(fmdu)
#'
#' n <- 1000
#' m <- 10
#' data <- matrix( runif( n * m ), n, m )
#' p <- 2
#' x <- matrix( rnorm( n * p ), n, p )
#' y <- matrix( rnorm( m * p ), m, p )
#' r <- ultrafastmdu( data, x, y )
#' }
#'
#' @importFrom stats runif
#' @export
#' @useDynLib fmdu, .registration=TRUE

ultrafastmdu <- function( data, x, y, w = NULL, fx = NULL, fy = NULL, NSTEPS = 4096, RCRIT = 0.00000001, seed = runif( 1, 1, as.integer( .Machine$integer.max ) ) )
{
  # parameter handling
  data <- as.matrix( data )
  n <- nrow( data )
  m <- ncol( data )
  x <- as.matrix( x )
  y <- as.matrix( y )
  p <- ncol( x )
  if ( NSTEPS <= 0 ) NSTEPS <- 1024
  if ( RCRIT <= 0.0 ) RCRIT <- 0.00000001

  # .C execution
  if ( is.null( w ) ) {
    if ( is.null( fx ) & is.null( fy ) ) result <- ( .C( "CRultrafastmdu", n=as.integer(n), m=as.integer(m), data=as.double(t(data)), p=as.integer(p), x=as.double(t(x)), y=as.double(t(y)), NSTEPS=as.integer(NSTEPS), RCRIT=as.double(RCRIT), seed=as.integer( seed ), PACKAGE= "fmdu" ) )
    else {
      fx <- ifelse( is.null( fx ), matrix( 0, n, p ), as.matrix( fx ) )
      fy <- ifelse( is.null( fy ), matrix( 0, m, p ), as.matrix( fy ) )
      result <- ( .C( "CRultrafastmdufxd", n=as.integer(n), m=as.integer(m), data=as.double(t(data)), p=as.integer(p), x=as.double(t(x)), fx=as.integer(t(fx)), y=as.double(t(y)), fy=as.integer(t(fy)), NSTEPS=as.integer(NSTEPS), RCRIT=as.double(RCRIT), seed=as.integer( seed ), PACKAGE= "fmdu" ) )
    }
  }
  else {
    w <- as.matrix( w )
    if ( is.integer( w ) ) {
      if ( is.null( fx ) & is.null( fy ) ) result <- ( .C( "CRultrafastwgtmdu", n=as.integer(n), m=as.integer(m), data=as.double(t(data)), w=as.integer(t(w)), p=as.integer(p), x=as.double(t(x)), y=as.double(t(y)), NSTEPS=as.integer(NSTEPS), RCRIT=as.double(RCRIT), seed=as.integer( seed ), PACKAGE= "fmdu" ) )
      else {
        fx <- ifelse( is.null( fx ), matrix( 0, n, p ), as.matrix( fx ) )
        fy <- ifelse( is.null( fy ), matrix( 0, m, p ), as.matrix( fy ) )
        result <- ( .C( "CRultrafastwgtmdufxd", n=as.integer(n), m=as.integer(m), data=as.double(t(data)), w=as.integer(t(w)), p=as.integer(p), x=as.double(t(x)), fx=as.integer(t(fx)), y=as.double(t(y)), fy=as.integer(t(fy)), NSTEPS=as.integer(NSTEPS), RCRIT=as.double(RCRIT), seed=as.integer( seed ), PACKAGE= "fmdu" ) )
      }
    }
    else {
      if ( is.null( fx ) & is.null( fy ) ) result <- ( .C( "CRultrafastmdu2", n=as.integer(n), m=as.integer(m), data=as.double(t(data)), p=as.integer(p), x=as.double(t(x)), y=as.double(t(y)), NSTEPS=as.integer(NSTEPS), RCRIT=as.double(RCRIT), seed=as.integer( seed ), PACKAGE= "fmdu" ) )
      else {
        fx <- ifelse( is.null( fx ), matrix( 0, n, p ), as.matrix( fx ) )
        fy <- ifelse( is.null( fy ), matrix( 0, m, p ), as.matrix( fy ) )
        result <- ( .C( "CRultrafastwgtmdufxd", n=as.integer(n), m=as.integer(m), data=as.double(t(data)), w=as.double(t(w)), p=as.integer(p), x=as.double(t(x)), fx=as.integer(t(fx)), y=as.double(t(y)), fy=as.integer(t(fy)), NSTEPS=as.integer(NSTEPS), RCRIT=as.double(RCRIT), seed=as.integer( seed ), PACKAGE= "fmdu" ) )
      }
    }
  }

  # finalization
  x <- matrix( result$x, n, p, byrow = TRUE )
  y <- matrix( result$y, m, p, byrow = TRUE )

  r <- list( x = x, y = y )
  r

} # ultrafastmdu
