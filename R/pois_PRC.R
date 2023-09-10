# Function for running PRC Poisson with rate unknown

pois_PRC <- function( data = NULL, s = NULL, historical_data = NULL, historical_s = NULL,
                      c0 = 1/2, d0 = 0, alpha_0 = NULL, k = 2, two.sided=FALSE,
                      h = log(100), FIR = FALSE, fFIR = 1/2, dFIR = 3/4,
                      summary_list = TRUE, PRC_PLOT = TRUE, pdf_report = FALSE, path_pdf_report = tempdir(),
                      xlab = "Observation Order", ylab = "PRC cumulative statistics",
                      main = "PRC Poisson with unknown rate" )
{
  ### Initial checks before procceeding to the main body of function
  ### Mainly this chunk of code will correspond to invalid general input before running stuff
  # 'data' (i) not defined (ii) not in vector (iii) contain non-numeric value
  if ( is.null(data) ) {
    stop("'data' have not been defined")
  } else { if ( any(!is.numeric((unlist(data)))) ) stop("Invalid 'data' input")
    if ( !is.vector(data) ) stop("'data' must be in vector form")
  }
  # 'historical_data' (i) not in vector (ii) contain non-numeric value
  if ( !is.null(historical_data) ) {
    if ( any(!is.numeric((unlist(historical_data)))) ) stop("Invalid 'historical_data' input")
    if ( !is.vector(data) ) stop("'historical data' must be in vector form")
  }

  # 'k' (i) non-numeric (ii) negative
  if( !missing(k) ) {
    if ( length(unlist(k))>1 ) { message("More than one value for 'k', the first one will only be used")
      if ( !is.numeric(k[1]) | k<=0 ) { stop("Invalid 'k' value") } else { k <- k[1] }
    } else { if ( !is.numeric(k) | k<=0 ) { stop("Invalid 'k' value") } }
  }


  # 'h' (i) non-numeric (ii) negative
  if( !missing(h) ) {
    if ( length(unlist(h))>1 ) { message("More than one value for 'h', the first one will only be used")
      if ( !is.numeric(h[1]) | h<=0 ) { stop("Invalid 'h' value") } else { h <- h[1] }
    } else { if ( !is.numeric(h) | h<=0 ) { stop("Invalid 'h' value") } }
  }


  # 'FIR' (i) logical (ii) fFIR - dFIR conditions
  # fFIR - dFIR conditions if  FIR
  if ( FIR ) {
    if ( !missing(dFIR) ) {
      if ( length(unlist(dFIR))>1 ) {
        message("More than one value for 'dFIR', the first one will only be used")
        if ( !is.numeric(dFIR[1]) | dFIR[1]<=0 | dFIR[1]>=1 ) {
          stop("Invalid 'dFIR' value")
        } else { dFIR <- dFIR[1] }
      } else {
        if ( !is.numeric(dFIR) | dFIR<=0 | dFIR>=1 ) {
          stop("Invalid 'dFIR' value")
        }
      }
    }

    if ( !missing(fFIR) ) {
      if ( length(unlist(fFIR))>1 ) {
        message("More than one value for 'fFIR', the first one will only be used")
        if ( !is.numeric(fFIR[1]) | fFIR[1]<=0 ) {
          stop("Invalid 'fFIR' value")
        } else { fFIR <- fFIR[1] }
      } else {
        if ( !is.numeric(fFIR) | fFIR<=0 ) {
          stop("Invalid 'fFIR' value")
        }
      }
    }
  }


  # data length
  N <- length(data)

  # If FIR PRC is chosen
  if ( FIR ) {
    tf <- 1:N
    fir_index <- c((  1 + fFIR * dFIR^(tf-1) ) )
  }

  ###############################################################
  ###############################################################
  ## START (1) Only this bit changes from function to function ##
  ## In some cases 'data' and 'historical_data' restrictions   ##
  ## change as well at the beginning of the function           ##
  ###############################################################
  ###############################################################

  # Prior parameter input (i) more than one value for parameters (ii) non-numeric input
  if ( is.null(s) ) { s <- rep(1, times=length(data))
  } else {
    if ( !is.vector(s) ) { stop("'rates - s' must be in vector form")
    } else {
      if ( length(s)!=length(data) ) { stop("Vector of 'rates - s' must have the same length as 'data'")
      } else {
        if ( any(!is.numeric((unlist(s)))) ) { stop("Invalid 'rate - s' input")
        } else { if( any((unlist(s)<=0)) ) stop("Invalid 'rate - s' input, s must be positive") }
      }
    }
  }

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

  ### PRC implementation
  # Sum of observations
  dataSum <- cumsum(data)[seq(1, N-1)]
  # Sum of rates
  dataRates <- cumsum(s)[seq(1, N-1)]

  # Posterior parameters
  c0_Post <- c0 + dataSum
  d0_Post <- d0 + dataRates

  # PRC Statistics
  Lu <- ( c0_Post[1:(N-1)]+data[2:N] ) * log(   (( s[2:N] + d0_Post[1:(N-1)])/(s[2:N] + d0_Post[1:(N-1)]/k) )  ) -
        ( c0_Post[1:(N-1)] ) * log (  k  )


  # FIR option
  if (FIR) { Lu <- Lu*fir_index[1:(N-1)] }

  Splus <- c() ; Splus[1] <- 0
  for (i in 2:N) { Splus[i] <- max( 0, Splus[i-1] + Lu[i-1] ) }

  # Two-sided PRC

  if (two.sided) {
    Ld <- ( c0_Post[1:(N-1)]+data[2:N] ) * log( (( s[2:N] + d0_Post[1:(N-1)])/(s[2:N] + k*d0_Post[1:(N-1)]) ) ) +
            c0_Post[1:(N-1)]* log (k)

    # FIR option
    if (FIR) { Ld <- Ld*fir_index[1:(N-1)] }

    Sminus <- c() ; Sminus[1] <- 0
    for (i in 2:N) { Sminus[i] <- min( 0, Sminus[i-1] - Ld[i-1] ) }
  }


  ####################################################################
  ####################################################################
  ## END (1) Only the above bit changes from function to function   ##
  ####################################################################
  ####################################################################

  ## Output
  if (!two.sided)  { # Construction of 'In' and 'Out' of control column for return results
    States <- rep( "", times = N )
    States[ifelse(Splus > h, TRUE, FALSE)] <- "Alarm"
    # Return results
    PRC_summary <- data.frame( data = data, Sn = Splus, Alarms = States )   } else {

      States <- rep("", times=N)
      States[ifelse(Splus > h, TRUE, FALSE)] <- "Alarm (U)" ; U_alarms <- ifelse(Splus > h, TRUE, FALSE)
      States[ifelse(Sminus < -h, TRUE, FALSE)] <- "Alarm (D)" ; D_alarms <- ifelse(Sminus < -h, TRUE, FALSE)
      States[ifelse(Splus > h & Sminus < -h, TRUE, FALSE)] <- "Alarm (Both)"
      # Return results
      PRC_summary <- data.frame( data = data, Snplus = Splus, Snminus = Sminus, Alarms = States )

    }


  ## Dynamic recalculation of PRC plot's y axis


  ### Output of function
  ## PRC plot
  if ( PRC_PLOT | pdf_report ) {
    # Creation of PRC plot
    PRC_PlotSummary <- cbind( Indices = 1:N, PRC_summary )

    if (!two.sided){

      PRC <-ggplot( PRC_PlotSummary, aes(PRC_PlotSummary[, "Indices"], PRC_PlotSummary[, "Sn"]) ) +
        geom_line( aes(x = PRC_PlotSummary[, "Indices"], y = PRC_PlotSummary[, "Sn"]), na.rm = TRUE ) +
        geom_line( aes(x = PRC_PlotSummary[, "Indices"], y = h), color = "red", linetype = "solid", size = 1 ) +
        geom_ribbon( aes(x = PRC_PlotSummary[, "Indices"], ymin = 0, ymax = h, fill = TRUE), alpha = 0.25, show.legend = FALSE ) +
        scale_fill_manual( values = c("TRUE"="green")) +
        geom_point( aes(group = PRC_PlotSummary[, "Indices"], color = as.factor(PRC_PlotSummary[, "Alarms"]), stroke = 1.5), show.legend = FALSE, na.rm = TRUE ) +
        scale_color_manual( values = c("black", "red", "red"), na.value = "black" ) +
        labs( title = main, x = xlab, y = ylab ) +
        theme( legend.position = "top",
               legend.title = element_blank(),
               axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
               panel.background = element_blank(),
               plot.title = element_text(size = 18, hjust = 0.5),
               text = element_text(size = 15),
               axis.text.x = element_text(colour="black", size = 12),
               axis.text.y = element_text(colour="black", size = 12) )
      } else{

      PRC <- ggplot( PRC_PlotSummary, aes(PRC_PlotSummary[, "Indices"], PRC_PlotSummary[, "Snplus"], PRC_PlotSummary[, "Snminus"]) ) +
        geom_line( aes(x = PRC_PlotSummary[, "Indices"], y = PRC_PlotSummary[, "Snplus"]), na.rm = TRUE ) +
        geom_line( aes(x = PRC_PlotSummary[, "Indices"], y = PRC_PlotSummary[, "Snminus"]), na.rm = TRUE ) +
        geom_line( aes(x = PRC_PlotSummary[, "Indices"], y = h), color = "red", linetype="solid", size = 1 ) +
        geom_line( aes(x = PRC_PlotSummary[, "Indices"], y = -h), color = "red", linetype="solid", size = 1 ) +
        geom_line( aes(x = PRC_PlotSummary[, "Indices"], y = 0), color = "gray75", linetype = "dashed", size = 0.5 ) +
        geom_ribbon( aes(x = PRC_PlotSummary[, "Indices"], ymin = -h, ymax = h, fill = TRUE), alpha = 0.25, show.legend = FALSE ) +
        scale_fill_manual(values = c("TRUE"="green")) +
        geom_point( aes(x = PRC_PlotSummary[, "Indices"], y = PRC_PlotSummary[, "Snplus"], color = as.factor(U_alarms), stroke = 1.5), show.legend = FALSE, na.rm = TRUE)  +
        geom_point( aes(x = PRC_PlotSummary[, "Indices"], y = PRC_PlotSummary[, "Snminus"], color = as.factor(D_alarms), stroke = 1.5), show.legend = FALSE, na.rm = TRUE ) +
        scale_color_manual( values = c("black", "red", "black", "red"), na.value = "black") +
        labs( title = main, x = xlab, y = ylab ) +
        theme( legend.position = "top",
               legend.title = element_blank(),
               axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
               panel.background = element_blank(),
               plot.title = element_text(size = 18, hjust = 0.5),
               text = element_text(size = 15),
               axis.text.x = element_text(colour="black", size = 12),
               axis.text.y = element_text(colour="black", size = 12) )
      }

    if ( PRC_PLOT) { print(PRC) }
  }
  # List of results
  if ( summary_list ) { print(PRC_summary) }

  # List of results return in pdf
  if ( pdf_report ) {

    # save pdf
    pdf(
      paste0( path_pdf_report, "\\", "PRC_results_", paste0( unlist(strsplit(date(), " "))[c(1,2,3,5)], collapse = "_" ), "_",
              paste0( unlist(strsplit( unlist(strsplit(date(), " "))[4], ":" )), collapse = "." ),
              ".pdf" ),
      height = 8.264, width = 11.694)

    # PRC plot on pdf
    print(PRC)


    # Results matrix on pdf
    # Chunk of code to split results matrix to different pages - Set a default number based on pdf height/width
    NRowsPerPage <- 25
    if(NRowsPerPage > nrow(PRC_summary)){ FloatingRow <- nrow(PRC_summary) } else { FloatingRow <- NRowsPerPage }
    sapply(1:ceiling(nrow(PRC_summary)/NRowsPerPage), function(index) {
      if (index==1) { StartingRow <- 1 }
      grid.newpage()
      grid.table(PRC_summary[StartingRow:FloatingRow, ])
      StartingRow <<- FloatingRow + 1
      if( sum(NRowsPerPage, FloatingRow) < nrow(PRC_summary)){ FloatingRow <<-  NRowsPerPage + FloatingRow } else { FloatingRow <<- nrow(PRC_summary) }
    })

    dev.off()

  }


}





