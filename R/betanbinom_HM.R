

betanbinom_HM <- function( cover = NULL, r = NULL, a = NULL, b = NULL, plot = FALSE, xlab = "x",
                           ylab = "Probability" ){

  # 'cover' (i) missing (ii) non-numeric (iii) out of the range (0,1)
  if ( is.null(cover) ) {
    stop("'cover' has not been defined")
  } else {
    if ( length(unlist(cover))>1 ) { message("More than one value for 'cover', the first one will only be used")
      if ( !is.numeric(cover[1]) | cover<=0 | cover>=1 ) { stop("Invalid 'cover' value") } else { cover <- cover[1] }
    } else { if ( !is.numeric(cover) | cover<=0 | cover>=1 ) { stop("Invalid 'cover' value") } }
  }

  # Likelihood r input (i) more than one value for parameters (ii) non-numeric input (iii) non-positive
  if ( is.null(r) ) {
    stop("'r' has not been defined")
  } else {
    if ( length(unlist(r))>1 ) { message("More than one value for 'r', the first one will only be used")
      if ( !is.numeric(r) | r<=0 ) { stop("Invalid 'r' value") } else { n <- n[1] }
    } else { if ( !is.numeric(r) | r<=0 ) { stop("Invalid 'r' value") } }
  }

  # Likelihood a input (i) more than one value for parameters (ii) non-numeric input (iii) non-positive
  if ( is.null(a) ) {
    stop("'a' has not been defined")
  } else {
    if ( length(unlist(a))>1 ) { message("More than one value for 'a', the first one will only be used")
      if ( !is.numeric(a) | a<=0 ) { stop("Invalid 'a' value") } else { a <- a[1] }
    } else { if ( !is.numeric(a) | a<=0 ) { stop("Invalid 'a' value") } }
  }

  # Likelihood b input (i) more than one value for parameters (ii) non-numeric input (iii) non-positive
  if ( is.null(b) ) {
    stop("'b' has not been defined")
  } else {
    if ( length(unlist(b))>1 ) { message("More than one value for 'b', the first one will only be used")
      if ( !is.numeric(b) | b<=0 ) { stop("Invalid 'b' value") } else { b <- b[1] }
    } else { if ( !is.numeric(b) | b<=0 ) { stop("Invalid 'b' value") } }
  }


  far <- 1-cover

  ##################################################################
  # Algorithm for the calculation of the HM region, using ordered probabilities

  mean_pr <- ifelse( a > 1, r*b/(a-1), b*10^2 )
  var_pr <- ifelse( a > 2, r*(a + r - 1)*b*(a + b - 1) / ((a - 2)*(a - 1)^2), (b^2)*10^4 )

  lb <- max( floor(mean_pr - sqrt((1/far)*var_pr)), 0 )
  ub <- ceiling( mean_pr + sqrt( (1/far)*var_pr ) )

  # Locations of the ordered probabilities
  Pi <- order( dbnbinom( lb:ub, size = r, alpha = a, beta = b ), decreasing = T ) + lb - 1
  nnn <- 1
  sumprob <- 0
  diff <- 1
  E <- c()
  stopp <- 0

  # Loop which ends when the absolute difference with the desired coverage is minimized
  while ( stopp==0 ) {
    sumprob <- sumprob + dbnbinom( Pi[nnn], size = r, alpha = a, beta = b )
    if ( abs(sumprob - (1-far)) < diff ) {
      E <- c( E, Pi[nnn] )
      diff <- abs( sumprob - (1-far) )
      nnn <- nnn + 1
    } else { stopp <- 1 }
  }

  # HM region
  ed <- c( min(E), max(E) )

  if ( plot == T) {

    # Graphical parameters for the range and the plotted region
    range <- ed[2] - ed[1]
    xi <- max( 0, floor( ed[1] - 0.15*range ) ):ceiling( ed[2] + 0.15*range )
    pi <- dbnbinom( xi, size = r, alpha = a, beta = b )

    # Graphical parameters for the colors
    ind0 <- which( xi<min(ed) | xi > max(ed) )
    cols <- rep( rgb(0, 1, 0, 0.3), times = length(xi) )
    cols[ind0] <- "white"

    # Graphical parameter for the main of the plot
    percov <- round( 100*sum( dbnbinom( ed[1]:ed[2], size = r, alpha = a, beta = b ) ), 2 )

    barplot( pi, xlab = xlab, ylab = ylab, main = bquote("Beta-Negative Binomial: "~.(percov)*"% HM = ["*.(ed[1])*", "~ .(ed[2])*"]" ), axes = F, names.arg = xi, col = cols )

    axis(2)

  }

  # The data frame of the output
  RES <- data.frame( lower.bound = ed[1], upper.bound = ed[2],
                     coverage = sum( dbnbinom( ed[1]:ed[2], size = r, alpha = a, beta = b ) ) )

  return(RES)

}




