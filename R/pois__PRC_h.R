# Function for obtaining h - Poisson with unknown rate

pois_PRC_h <- function( ARL_0 = 370.4, FAP = NULL, N = NULL, s = NULL, historical_data = NULL, historical_s = NULL,
                        c0 = NULL, d0 = NULL, alpha_0 = NULL, k = 2, it = 1e4, ARL0tol = 10/it )

{
  ### Initial checks before proceeding to the main body of function
  ### Mainly this chunk of code will correspond to invalid general input before running stuff

  # 'ARL_0' (i) non-numeric (ii) negative
  if( !is.null(ARL_0) ) {
    if ( length(unlist(ARL_0))>1 ) { message("More than one value for 'ARL_0', the first one will only be used")
      if ( !is.numeric(ARL_0[1]) | ARL_0<=0 ) { stop("Invalid 'ARL_0' value") } else { ARL_0 <- ARL_0[1] }
    } else { if ( !is.numeric(ARL_0) | ARL_0<=0 ) { stop("Invalid 'ARL_0' value") } }
  }

  # 'FAP' (i) non-numeric (ii) negative
  if (!missing(FAP)){
    if ( length(unlist(FAP))>1 ) { message("More than one value for 'FAP', the first one will only be used")
      if ( !is.numeric(FAP[1]) | FAP<=0 | FAP>=1 ) { stop("Invalid 'FAP' value") } else { FAP <- FAP[1] }
    } else { if ( !is.numeric(FAP) | FAP<=0 | FAP>=1 ) { stop("Invalid 'FAP' value") } }
  }

  # 'N' (i) non-numeric (ii) negative (iii) non integer
  if( !missing(N) ) {
    if ( length(unlist(N))>1 ) { message("More than one value for 'N', the first one will only be used")
      if ( !is.numeric(N[1]) | N<=0 | N%%1!=0) { stop("Invalid 'N' value") } else { N <- N[1] }
    } else { if ( !is.numeric(N) | N<=0 | N%%1!=0) { stop("Invalid 'N' value") } }
  }

  # 'ARL0tol' (i) non-numeric (ii) negative
  if( !missing(ARL0tol) ) {
    if ( length(unlist(ARL0tol))>1 ) { message("More than one value for 'ARL0tol', the first one will only be used")
      if ( !is.numeric(ARL0tol[1]) | ARL0tol<=0 ) { stop("Invalid 'ARL0tol' value") } else { ARL0tol <- ARL0tol[1] }
    } else { if ( !is.numeric(ARL0tol) | ARL0tol<=0 ) { stop("Invalid 'ARL0tol' value") } }
  }

  # 'it' (i) non-numeric (ii) negative (iii) non integer
  if( !missing(it) ) {
    if ( length(unlist(it))>1 ) { message("More than one value for 'it', the first one will only be used")
      if ( !is.numeric(it[1]) | it<=0 | it%%1!=0) { stop("Invalid 'it' value") } else { it <- it[1] }
    } else { if ( !is.numeric(it) | it<=0 | it%%1!=0) { stop("Invalid 'it' value") } }
  }

  # 'historical_data' (i) not in vector (ii) contain non-numeric value
  if ( !is.null(historical_data) ) {
    if ( any(!is.numeric((unlist(historical_data)))) ) stop("Invalid 'historical_data' input")
    if ( !is.vector(historical_data) ) stop("'historical data' must be in vector form")
  }

  # 'k' (i) non-numeric (ii) negative
  if( !missing(k) ) {
    if ( length(unlist(k))>1 ) { message("More than one value for 'k', the first one will only be used")
      if ( !is.numeric(k[1]) | k<=0 ) { stop("Invalid 'k' value") } else { k <- k[1] }
    } else { if ( !is.numeric(k) | k<=0 ) { stop("Invalid 'k' value") } }
  }


  ###############################################################
  ###############################################################
  ## START (1) Only this bit changes from function to function ##
  ## In some cases 'data' and 'historical_data' restrictions   ##
  ## change as well at the beginning of the function           ##
  ###############################################################
  ###############################################################

  if ( is.null(s) & !is.null(N)) { s <- rep(1, times=N)
  } else {
    if ( is.null(s) & is.null(N)) { s <- 1}
    if ( !is.vector(s) ) { stop("'rates - s' must be in vector form")
    } else {
        if ( any(!is.numeric((unlist(s)))) ) { stop("Invalid 'rate - s' input")
        } else { if( any((unlist(s)<=0)) ) stop("Invalid 'rate - s' input, s must be positive") }
    }
  }

  if ( is.null(s) & is.null(N)) { s <- 1}


  if( !missing(c0) ) {
    if ( length(unlist(c0))>1 ) { message("More than one value for 'c0', the first one will only be used")
      if ( !is.numeric(c0) | c0<=0 ) { stop("Invalid 'c0' value") } else { c0 <- c0[1] }
    } else { if ( !is.numeric(c0) | c0<=0 ) { stop("Invalid 'c0' value") } }
  }

  if( !missing(d0) ) {
    if ( length(unlist(d0))>1 ) { message("More than one value for 'd0', the first one will only be used")
      if ( !is.numeric(d0) | d0<=0 ) { stop("Invalid 'd0' value") } else { d0 <- d0[1] }
    } else { if ( !is.numeric(d0) | d0<=0 ) { stop("Invalid 'd0' value") } }
  }

  ### Main body of function - PCC illustration - USING FAR (or FAP equivelantly)
  ## Histotic data and processing
  if ( !is.null(historical_data) ){

    if ( is.null(historical_s) ) { historical_s <- rep(1, times=length(historical_data))
    } else {
      if ( !is.vector(historical_s) ) { stop("'historical rates - historical_s' must be in vector form")
      } else {
        if ( length(historical_s)!=length(historical_data) ) { stop("Vector of 'historical_s rates' must have the same length as 'historical_data'")
        } else {
          if ( any(!is.numeric((unlist(historical_s)))) ) { stop("Invalid 'historical rates - historical_s' input")
          } else { if( any((unlist(historical_s)<=0)) ) stop("Invalid 'historical rates - historical_s' input, historical_s must be positive") }
        }
      }
    }

    N_historicaldata <- length(historical_data)    # Check about alpha_0
    # If no chosen value for alpha_0 use default setting
    if (is.null(alpha_0)) { alpha_0 <-1/N_historicaldata
    } else {
      if ( length(unlist(alpha_0))>1 ) {
        message("More than one value for 'alpha_0', the first one will only be used")
        if ( !is.numeric(alpha_0) | alpha_0<0 | alpha_0>1) { stop("Invalid 'alpha_0' value")
        } else { if ( !is.numeric(alpha_0) | alpha_0<0 | alpha_0>1 ) { stop("Invalid 'alpha_0' value") } }
      }
    }
    # Process historical data
    # Power Prior parameters
    c_PowerP  <- c0 + alpha_0*sum(historical_data)
    d_PowerP  <- d0 + alpha_0*sum(historical_s)
    # Keep similar notation as input
    c0 <- c_PowerP ; d0 <- d_PowerP
  }

  s <- round( mean(s) )
  rr <- 1 - ( 1/(d0+1) )

  if ( rr < 0.9 ) {
    message("Warning message: The expected ratio of the likelihood over the marginal is quite small (<0.9), \nso the resulting desicion limit may be too conservative")
  }

  if ( !is.null(ARL_0) & !is.null(FAP) ) {
    message("Both ARL_0 and FAP are defined as input, so ARL_0 is used by default. \nIn order to use FAP instead, set ARL_0 = NULL")
    FAP <- NULL
    # If only FAP is chosen
  } else if ( is.null(ARL_0) & !is.null(FAP) ) {

    set.seed(47)
    ic <- matrix( rnbinom( it*N, size = c0,prob = d0/(s+d0) ), nrow = it, ncol = N )

    step <- 1
    vector_h <- c()

    while ( step<=length(ic[,1]) ) {

      PRCSTATS <- c()
      n <- 1
      xsum <- 0
      ssum <- 0
      xsum <- xsum + ic[step,n]
      ssum <- ssum + s

      A <- c0 + xsum
      B <- d0 + ssum
      Lplus <- 0
      PRCSTATS[1] <- 0
      n <- 2

      while ( n<=length(ic[1,]) ){

        Lu <- ( A + ic[step,n] ) * log( (( s + B)/(s + B/k) ) )  + A* log ( (1/k) )

        Lplus <- max( 0, Lplus + Lu )
        PRCSTATS[n] <- Lplus
        xsum <- sum + ic[step, n]
        ssum <- ssum + s

        A <- c0 + xsum
        B <- d0 + ssum

        n <- n + 1
      }
      vector_h <- c( vector_h, max(PRCSTATS) )
      step <- step + 1
    }

    h <- unname( quantile(vector_h, 1-FAP) )
    return(h)

      # If only ARL0 is chosen
  } else if ( !is.null(ARL_0) & is.null(FAP) ) {


    arlfun <- function( IC, S, c, d, hs ) {

      STOPS <- c()
      step <- 1


      for ( step in 1:it ) {

        ic <- IC[step, ]
        s <- S

        n <- 1
        xsum <- 0
        ssum <- 0
        xsum <- xsum + ic[n]
        ssum <- ssum + s

        A <- c + xsum
        B <- d + ssum

        Lplus <- 0
        n <- 2

        end <- 0

        while ( end==0 ) {

          Lu <- ( A + ic[n] ) * log( (( s + B)/(s + B/k) ) ) + A* log ( (1/k) )
          Lplus <- max( 0, Lplus + Lu )

          if ( Lplus>hs ) {

            STOPS <- c(STOPS, n) ; end <- 1

          } else {

            xsum <- xsum + ic[n]
            ssum <- ssum + s

            A <- c + xsum
            B <- d + ssum
            if ( n==length(ic) ) {

              set.seed( 3*step + 2 + n )
              icnew <- rnbinom( round(4*ARL_0), size = c, prob = d/(s+d) )
              ic <- c( ic, icnew )

            }
            n <- n + 1
          }
        }
      }
      return( mean(STOPS) )
    }


    set.seed(5)

    IC <- matrix( rnbinom( it*round(4*ARL_0), size = c0, prob = d0/(s+d0) ), nrow = it, ncol = round(4*ARL_0), byrow = T )


    ind <- 1
    h1 <- max( (log(ARL_0*k^2/2+1) - k^(1.1)-0.6)/rr, 1.5 )
    ARL1 <- arlfun( IC = IC, S = s, c = c0, d = d0, hs = h1 )
    print( paste0("It ", ind, ": Suggested h=", round(h1, 3), ", achieved ARL_0=", round(ARL1, 3)) )

    if ( abs(ARL1-ARL_0)/ARL_0 < ARL0tol ) { return(h1) }

    ind <- 2
    if ( ARL1>ARL_0 ) {
      h2 <- max( ( log(ARL_0*k^2/2+1) - k^(1.1)-0.6 )/rr - 0.4, 0.5 )
    } else {
      h2 <- max( ( log(ARL_0*k^2/2+1) - k^(1.1)+0.2 )/rr, 3 )
    }

    ARL2 <- arlfun( IC = IC, S = s, c = c0, d = d0, hs = h2 )
    print( paste0("It ", ind, ": Suggested h=", round(h2, 3), ", achieved ARL_0=", round(ARL2, 3)) )


    while( ARL2==ARL1 ) {

      h2 <- max(h2+rnorm(1, 0, 0.1), 0.001)
      ind <- ind + 1
      ARL2 <- arlfun( IC = IC, S = s, c = c0, d = d0, hs = h2 )
      print( paste0("It ", ind, ": Suggested h=", round(h2, 3), ", achieved ARL_0=", round(ARL2, 3)) )

    }

    if ( abs(ARL2-ARL_0)/ARL_0 < ARL0tol ) { return(h2) }

    hnew <- max( h2 + (ARL_0-ARL2)*(h2-h1)/(ARL2-ARL1), 0.001 )
    ind <- ind + 1
    ARLnew <- arlfun( IC = IC, S = s, c = c0, d = d0, hs = hnew )
    print( paste0("It ", ind, ": Suggested h=", round(hnew, 3), ", achieved ARL_0=", round(ARLnew, 3)) )

      vh <- c( round(h1, 3), round(h2, 3), round(hnew, 3) )

      if (length(unique(vh)) == 1) {

        message("Warning message: the desired accuracy cannot be achieved.")
        return(hnew)

      }

    while( abs(ARLnew-ARL_0)/ARL_0 > ARL0tol ) {
      h1 <- h2
      ARL1 <- ARL2
      h2 <- hnew
      ARL2 <- ARLnew

      while( ARL2==ARL1 ) {
        h2 <- max( h2 + rnorm(1, 0, 0.1), 0.001 )
        ind <- ind + 1
        ARL2 <- arlfun( IC = IC, S = s, c = c0, d = d0, hs = h2 )
        print( paste0("It ", ind, ": Suggested h=", round(hnew, 3), ", achieved ARL_0=", round(ARLnew, 3)) )

      }

      hnew <- max( h2 + (ARL_0-ARL2)*(h2-h1)/(ARL2-ARL1), 0.001 )
      ind <- ind + 1
      ARLnew <- arlfun( IC = IC, S = s, c = c0, d = d0, hs = hnew )
      print( paste0("It ", ind, ": Suggested h=", round(hnew, 3), ", achieved ARL_0=", round(ARLnew, 3)) )

      vh <- c( round(h1, 3), round(h2, 3), round(hnew, 3) )
      if (length(unique(vh)) == 1) {
        message( "Warning message: the desired accuracy cannot be achieved." )
        return(hnew)
      }

    }

    return(hnew)

  }

}








