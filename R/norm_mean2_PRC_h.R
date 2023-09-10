# Function for obtaining h - Normal with both parameters unknown - mean model

norm_mean2_PRC_h <- function( ARL_0 = 370.4, FAP = NULL, N=NULL, historical_data = NULL,
                              l0 = 0, a0 = -1/2, alpha_0 = NULL, k = 1, it = 1e4, ARL0tol = 10/it )

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

  # 'k' (i) non-numeric
  if( !missing(k) ) {
    if ( length(unlist(k))>1 ) { message("More than one value for 'k', the first one will only be used")
      if ( !is.numeric(k[1]) ) { stop("Invalid 'k' value") } else { k <- k[1] }
    } else { if ( !is.numeric(k)) { stop("Invalid 'k' value") } }
  }


  ###############################################################
  ###############################################################
  ## START (1) Only this bit changes from function to function ##
  ## In some cases 'data' and 'historical_data' restrictions   ##
  ## change as well at the beginning of the function           ##
  ###############################################################
  ###############################################################


  # Prior parameter input (i) more than one value for parameters (ii) non-numeric input

  if( !missing(l0) ) {
    if ( length(unlist(l0))>1 ) { message("More than one value for 'l0', the first one will only be used")
      if ( !is.numeric(l0) | l0<0 ) { stop("Invalid 'l0' value") } else { l0 <- l0[1] }
    } else { if ( !is.numeric(l0) | l0<0 ) { stop("Invalid 'l0' value") } }
  }

    if( !missing(a0) ) {
    if ( length(unlist(a0))>1 ) { message("More than one value for 'a0', the first one will only be used")
      if ( !is.numeric(a0) | a0<=0 ) { stop("Invalid 'a0' value") } else { a0 <- a0[1] }
    } else { if ( !is.numeric(a0) | a0<=0 ) { stop("Invalid 'a0' value") } }
  }

  ### Main body of function - PCC illustration - USING FAR (or FAP equivelantly)
  ## Historic data and processing
  if ( !is.null(historical_data) ){
    N_historicaldata <- length(historical_data)
    # Check about alpha_0
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
    l0_PowerP  <- l0 + alpha_0*N_historicaldata
    a0_PowerP  <- a0 + alpha_0*N_historicaldata/2
    # Keep similar notation as input
    l0 <- l0_PowerP ; a0 <- a0_PowerP
  }


  nn <- ifelse( a0==-1/2, 3, 2 )


  if ( !is.null(ARL_0) & !is.null(FAP) ) {
    message("Both ARL_0 and FAP are defined as input, so ARL_0 is used by default. \nIn order to use FAP instead, set ARL_0 = NULL")
    FAP <- NULL
    # If only FAP is chosen
  } else if ( is.null(ARL_0) & !is.null(FAP) ) {

    ic <- matrix( 0, it, N )

    iy <- rep( nn:N, each = it )
    set.seed(47)
    ic[, nn:N] <- rt( it*(N-nn+1), df = 2*(a0+iy/2) )


    step <- 1
    vector_h <- c()

    while ( step <= length(ic[, 1]) ){

      PRCSTATS <- c()
      n <- nn - 1

      an <- a0 + n/2
      ln <- l0 + n

      Lplus <- 0
      PRCSTATS[1:2] <- Lplus

      n <- nn

      while ( n <= length(ic[1, ]) ){

        Lu <- ( an + 1/2 )*log( (2*an + ic[step, n]^2) / (2*an + (ic[step, n] - k*ln/(ln+1)  )^2) )
        Lplus <- max( 0, Lplus + Lu )
        PRCSTATS[n] <- Lplus

        an <- a0 + n/2
        ln <- l0 + n

        n <- n + 1
      }

      vector_h <- c( vector_h, max(PRCSTATS) )
      step <- step + 1
    }

    h <- unname( quantile(vector_h, 1-FAP) )
    return(h)


      # If only ARL0 is chosen
  } else if ( !is.null(ARL_0) & is.null(FAP) ){

  arlfun <- function( IC, l, a, hl ) {

    STOPS <- c()
    step <- 1


    for ( step in 1:it ) {

      ic <- IC[step, ]

      n <- nn-1

      an <- a + n/2
      ln <- l + n

      Lplus <- 0
      n <- nn

      end <- 0

      while ( end==0 ) {


        Lu <- ( an + 1/2 )*log( (2*an + ic[n]^2) / (2*an + (ic[n] - k*ln/(ln+1)  )^2) )
        Lplus <- max( 0, Lplus + Lu )

        if ( Lplus > hl ) {

          STOPS <- c(STOPS,n) ; end <- 1

        } else {

            an <- a + n/2
            ln <- l + n

            if ( n==length(ic) ) {

              iii <- 1:(round(4*ARL_0))

              set.seed(3*step+2+n)
              icnew <- rt( (round(4*ARL_0)), df = 2*(an+iii/2) )
              ic <- c( ic, icnew )

            }

            n <- n + 1
        }
      }
    }

    return(mean(STOPS))

  }


  IC <- matrix( 0, it, round(4*ARL_0) )

  ii <- rep( nn:(round(4*ARL_0)), each = it )
  set.seed(5)
  IC[, nn:(round(4*ARL_0))] <- rt(it*(round(4*ARL_0)-nn+1), df = 2*(a0+ii/2))

  ind <- 1
  h1 <- max( log( ARL_0*k^2/2 + 1 ) - k^(1.1) - 0.2, 0.6 )
  ARL1 <- arlfun( IC = IC, l = l0, a = a0, hl = h1 )
  print( paste0("It ", ind, ": Suggested h=", round(h1, 3), ", achieved ARL_0=", round(ARL1, 3)) )
  if ( abs(ARL1-ARL_0)/ARL_0 < ARL0tol ) { return(h1) }

  ind <- 2
  if ( ARL1 > ARL_0 ) {

    h2 <- max( log( ARL_0*k^2/2 + 1 ) - k^(1.1) - 0.6, 0.4 )

  } else {

    h2 <- max(log(ARL_0*k^2/2+1)-k^(1.1)+0.2, 0.9) }
    ARL2 <- arlfun( IC=IC, l=l0, a=a0, hl=h2 )
    print( paste0("It ", ind, ": Suggested h=", round(h2, 3), ", achieved ARL_0=", round(ARL2, 3)) )

  while( ARL2 == ARL1 ) {

    h2 <- max( h2 + rnorm(1, 0, 0.1), 0.001 )
    ind <- ind + 1
    ARL2 <- arlfun( IC = IC, l = l0, a = a0, hl = h2 )
    print( paste0("It ", ind, ": Suggested h=", round(h2, 3), ", achieved ARL_0=", round(ARL2, 3)) )

  }

  if ( abs(ARL2-ARL_0)/ARL_0 < ARL0tol ) { return(h2) }

  hnew <- max( h2 + (ARL_0-ARL2)*(h2-h1)/(ARL2-ARL1), 0.001 )
  ind <- ind + 1
  ARLnew <- arlfun( IC = IC, l = l0, a = a0, hl = hnew )
  print( paste0("It ", ind,": Suggested h=", round(hnew, 3), ", achieved ARL_0=", round(ARLnew, 3)) )

      vh <- c( round(h1, 3), round(h2, 3), round(hnew, 3) )
      if ( length(unique(vh)) == 1 ) {

        message( "Warning message: the desired accuracy cannot be achieved." )
        return(hnew)

      }

  while( abs(ARLnew-ARL_0)/ARL_0 > ARL0tol ) {

    h1 <- h2
    ARL1 <- ARL2
    h2 <- hnew
    ARL2 <- ARLnew

    while( ARL2 == ARL1 ) {

      h2 <- max( h2 + rnorm(1, 0, 0.1), 0.001 )
      ind <- ind + 1
      ARL2 <- arlfun( IC = IC, l = l0, a = a0, hl = h2 )
      print( paste0("It ", ind, ": Suggested h=", round(h2, 3), ", achieved ARL_0=", round(ARL2, 3)) )

    }

    hnew <- max (h2 + (ARL_0-ARL2)*(h2-h1)/(ARL2-ARL1), 0.001 )
    ind <- ind + 1
    ARLnew <- arlfun( IC = IC, l = l0, a = a0, hl = hnew )
    print( paste0("It ", ind, ": Suggested h=", round(hnew, 3), ", achieved ARL_0=", round(ARLnew, 3)) )

      vh <- c( round(h1, 3), round(h2, 3), round(hnew, 3) )
      if ( length(unique(vh)) == 1 ) {

        message( "Warning message: the desired accuracy cannot be achieved." )
        return(hnew)

      }

  }

  return(hnew)

  }

}



