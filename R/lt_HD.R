
lt_HD <- function( cover = NULL, df = NULL, mulog = NULL, sdlog = NULL, plot = FALSE, xlab = "x",
                   ylab = "Density" ){

  # 'cover' (i) missing (ii) non-numeric (iii) out of the range (0,1)
  if ( is.null(cover) ) {
    stop("'cover' has not been defined")
  } else {
    if ( length(unlist(cover))>1 ) { message("More than one value for 'cover', the first one will only be used")
      if ( !is.numeric(cover[1]) | cover<=0 | cover>=1 ) { stop("Invalid 'cover' value") } else { cover <- cover[1] }
    } else { if ( !is.numeric(cover) | cover<=0 | cover>=1 ) { stop("Invalid 'cover' value") } }
  }

  # Likelihood degrees of freedom input (i) more than one value for parameters (ii) non-numeric input (iii) non-positive
  if ( is.null(df) | is.na(df) ) {
    message("'df' has not been defined")
    return( data.frame( lower.bound = NA, upper.bound = NA, coverage = NA ) )
  } else {
    if ( length(unlist(df))>1 ) { message("More than one value for 'df', the first one will only be used")
      if ( !is.numeric(df) | df<=0 ) { stop("Invalid 'df' value") } else { df <- df[1] }
    } else { if ( !is.numeric(df) | df<=0 ) { stop("Invalid 'df' value") } }
  }

  # Likelihood scale input (i) more than one value for parameters (ii) non-numeric input (iii) non-positive
  if ( is.null(mulog) ) {
    stop("'mulog' has not been defined")
  } else {
    if ( length(unlist(mulog))>1 ) { message("More than one value for 'mulog', the first one will only be used")
      if ( !is.numeric(mulog) ) { stop("Invalid 'mulog' value") } else { mulog <- mulog[1] }
    } else { if ( !is.numeric(mulog)  ) { stop("Invalid 'mulog' value") } }
  }

  # Likelihood shape input (i) more than one value for parameters (ii) non-numeric input (iii) non-positive
  if ( is.null(sdlog) ) {
    stop("'sdlog' has not been defined")
  } else {
    if ( length(unlist(sdlog))>1 ) { message("More than one value for 'sdlog', the first one will only be used")
      if ( !is.numeric(sdlog) | sdlog<=0 ) { stop("Invalid 'sdlog' value") } else { sdlog <- sdlog[1] }
    } else { if ( !is.numeric(sdlog) | sdlog<=0 ) { stop("Invalid 'sdlog' value") } }
  }

  far <- 1-cover
  f <- function(x){ exp(mulog + qt( 1-far+x, df = df )*sdlog) -
      exp(mulog + qt( x, df = df )*sdlog ) }
  out <- optimize(f, c(0, far), tol = .Machine$double.eps)
  ed <- c( exp(mulog + qt( out$minimum, df = df )*sdlog ),
           exp(mulog + qt( 1-far+out$minimum, df = df )*sdlog ) )


  if ( plot==T ) {

    # Graphical parameters for the range and the plotted region
    range <- ed[2] - ed[1]
    xi <- seq( max(0.001, ed[1] - 0.15*range), ed[2] + 0.15*range, length.out = 10^4 )
    yi <- dt( (log(xi)-mulog)/sdlog, df = df )/xi

    # Graphical parameters for the main of the plot
    percov <- 100*cover
    ed1 <- round( ed[1], 2 )
    ed2 <- round( ed[2], 2 )

    plot( xi, yi, xlim = c( max(0, 0.978*min(xi)), 1.025*max(xi) ), ylim = c( 0, 1.05*max(yi) ), type = "l",
          xlab = xlab, ylab = ylab, main = bquote("Logt: "~.(percov)*"% HD = ["*.(ed1)*", "~ .(ed2)*"]" ), axes = F, lwd = 2 )

    # Green vertical segments on the bounds of the HD region
    segments( ed[1], 0, ed[1], dt( (log(ed[1])-mulog)/sdlog, df = df )/ed[1], lwd = 3, col = "green" )
    segments( ed[2], 0, ed[2], dt( (log(ed[2])-mulog)/sdlog, df = df )/ed[2], lwd = 3, col = "green" )

    axis(1) ; axis(2)

    # Adding the light green area in the graph
    xi2 <- seq( ed[1], ed[2], length.out = 10^4 )
    polygon( c( xi2, rev(xi2) ), c( dt( (log(xi2)-mulog)/sdlog, df = df )/xi2, rep(0, 10^4) ), col = rgb(0, 1, 0, 0.3), border = NA )

  }

  # The data frame of the output
  RES <- data.frame( lower.bound = ed[1], upper.bound = ed[2], coverage = cover )

  return(RES)

}




