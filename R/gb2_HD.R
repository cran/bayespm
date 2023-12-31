
gb2_HD <- function( cover = NULL, a = NULL, c = NULL, d = NULL, plot= FALSE, xlab = "x",
                    ylab = "Density" ){

  # 'cover' (i) missing (ii) non-numeric (iii) out of the range (0,1)
  if ( is.null(cover) ) {
    stop("'cover' has not been defined")
  } else {
    if ( length(unlist(cover))>1 ) { message("More than one value for 'cover', the first one will only be used")
      if ( !is.numeric(cover[1]) | cover<=0 | cover>=1 ) { stop("Invalid 'cover' value") } else { cover <- cover[1] }
    } else { if ( !is.numeric(cover) | cover<=0 | cover>=1 ) { stop("Invalid 'cover' value") } }
  }

  # Likelihood shape1 input (i) more than one value for parameters (ii) non-numeric input (iii) non-positive
  if ( is.null(a) ) {
    stop("'a' has not been defined")
  } else {
    if ( length(unlist(a))>1 ) { message("More than one value for 'k', the first one will only be used")
      if ( !is.numeric(a) | a<=0 ) { stop("Invalid 'a' value") } else { a <- a[1] }
    } else { if ( !is.numeric(a) | a<=0 ) { stop("Invalid 'a' value") } }
  }

  # Likelihood shape2 input (i) more than one value for parameters (ii) non-numeric input (iii) non-positive
  if ( is.null(c) ) {
    stop("'c' has not been defined")
  } else {
    if ( length(unlist(c))>1 ) { message("More than one value for 'c', the first one will only be used")
      if ( !is.numeric(c) | c<=0 ) { stop("Invalid 'c' value") } else { c <- c[1] }
    } else { if ( !is.numeric(c) | c<=0 ) { stop("Invalid 'c' value") } }
  }

  # Likelihood scale input (i) more than one value for parameters (ii) non-numeric input (iii) non-positive
  if ( is.null(d) ) {
    stop("'d' has not been defined")
  } else {
    if ( length(unlist(d))>1 ) { message("More than one value for 'd', the first one will only be used")
      if ( !is.numeric(d) | d<=0 ) { stop("Invalid 'd' value") } else { d <- d[1] }
    } else { if ( !is.numeric(d) | d<=0 ) { stop("Invalid 'd' value") } }
  }


  # Calculation of the HD region
  far <- 1-cover
  f <- function(x){ 1/qbetapr( far-x, shape1 = a, shape2 = c, scale = d ) -
                    1/qbetapr( 1-x, shape1 = a, shape2 = c, scale = d ) }
  out <- optimize( f, c(0, far), tol = .Machine$double.eps )
  ed <- c( 1/qbetapr( 1-out$minimum, shape1 = a, shape2 = c, scale = d ),
           1/qbetapr( far-out$minimum, shape1 = a, shape2 = c, scale = d ) )

  if ( plot==T ) {

    # Graphical parameters for the range and the plotted region
    range <- ed[2] - ed[1]
    xi <- seq( max(0.001, ed[1] - 0.15*range), ed[2] + 0.15*range, length.out = 10^4 )
    yi <- dbetapr( 1/xi, shape1 = a, shape2 = c, scale = d )/(xi^2)

    # Graphical parameters for the main of the plot
    percov <- 100*cover
    ed1 <- round( ed[1], 2 )
    ed2 <- round( ed[2], 2 )

    plot( xi, yi, xlim = c( max(0, 0.978*min(xi)), 1.025*max(xi) ), ylim = c( 0, 1.05*max(yi) ), type = "l",
          xlab = xlab, ylab = ylab, main = bquote("Generalized Beta Prime: "~.(percov)*"% HD = ["*.(ed1)*", "~ .(ed2)*"]" ), axes = F,  lwd = 2 )

    # Green vertical segments on the bounds of the HD region
    segments( ed[1], 0, ed[1], dbetapr( 1/ed[1], shape1 = a, shape2 = c, scale = d )/(ed[1]^2), lwd = 3, col = "green" )
    segments( ed[2], 0, ed[2], dbetapr( 1/ed[2], shape1 = a, shape2 = c, scale = d )/(ed[2]^2), lwd = 3, col = "green" )

    axis(1) ; axis(2)

    # Adding the light green area in the graph
    xi2 <- seq( ed[1], ed[2], length.out = 10^4 )
    polygon( c(xi2, rev(xi2)), c( dbetapr( 1/xi2, shape1 = a, shape2 = c, scale = d )/(xi2^2), rep(0, 10^4) ), col = rgb(0, 1, 0, 0.3), border = NA )

  }

  # The data frame of the output
  RES <- data.frame( lower.bound = ed[1], upper.bound = ed[2], coverage = cover )

  return(RES)

}



