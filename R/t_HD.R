
t_HD <- function( cover = NULL, df = NULL, mu = NULL, sdv = NULL, plot = FALSE, xlab = "x",
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

  # mu input (i) more than one value for parameters (ii) non-numeric input
  if ( is.null(mu) ) {
    stop("'mu' has not been defined")
  } else {
    if ( length(unlist(mu))>1 ) { message("More than one value for 'mu', the first one will only be used")
      if ( !is.numeric(mu) ) { stop("Invalid 'mu' value") } else { mu <- mu[1] }
    } else { if ( !is.numeric(mu)  ) { stop("Invalid 'mu' value") } }
  }

  # Likelihood sd input (i) more than one value for parameters (ii) non-numeric input (iii) non-positive
  if ( is.null(sdv) ) {
    stop("'sdv' has not been defined")
  } else {
    if ( length(unlist(sdv))>1 ) { message("More than one value for 'sdv', the first one will only be used")
      if ( !is.numeric(sdv) | sdv<=0 ) { stop("Invalid 'sdv' value") } else { sdv <- sdv[1] }
    } else { if ( !is.numeric(sdv) | sdv<=0 ) { stop("Invalid 'sdv' value") } }
  }


  # Calculation of the HD region
  far <- 1-cover

  if ( is.na(df) ) {
    ed <- rep( NA, 2 )
  } else {
    ed <- c( mu + qt( far/2, df = df ) * sdv, mu + qt( 1-far/2, df = df ) * sdv )
  }

  if ( plot==T ) {

    # Graphical parameters for the range and the plotted region
    range <- ed[2] - ed[1]
    xi <- seq( ed[1]-0.15*range, ed[2]+0.15*range, length.out = 10^4 )
    yi= dt( (xi-mu)/sdv, df = df ) / sdv

    # Graphical parameters for the main of the plot
    percov <- 100*cover
    ed1 <- round(ed[1], 2)
    ed2 <- round(ed[2], 2)

    plot( xi, yi, xlim = c( 0.978*min(xi), 1.025*max(xi) ), ylim = c (0, 1.05*max(yi) ), type = "l",
          xlab = xlab, ylab = ylab, main = bquote("t-Student: "~.(percov)*"% HD = ["*.(ed1)*", "~ .(ed2)*"]" ), axes = F,  lwd = 2)

    # Green vertical segments on the bounds of the HD region
    segments( ed[1], 0, ed[1], dt( (ed[1]-mu)/sdv, df = df ) / sdv, lwd = 3, col = "green")
    segments( ed[2], 0, ed[2], dt( (ed[2]-mu)/sdv, df = df ) / sdv, lwd = 3, col = "green")

    axis(1) ; axis(2)

    # Adding the light green area in the graph
    xi2 <- seq( ed[1], ed[2], length.out = 10^4 )
    polygon( c(xi2, rev(xi2)), c( dt( (xi2-mu)/sdv, df = df ) / sdv, rep(0, 10^4) ), col = rgb(0, 1, 0, 0.3), border = NA )

  }

  # The data frame of the output
  RES <- data.frame( lower.bound = ed[1], upper.bound = ed[2], coverage = cover )

  return(RES)



}

