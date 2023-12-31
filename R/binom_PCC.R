# Function for running PCC Binomial with probability unknown

binom_PCC <- function( data = NULL, n = NULL, historical_data = NULL,
                       historical_n = NULL, a0 = 1/2, b0 = 1/2, alpha_0 = NULL,
                       ARL_0 = 370.4, FAP = NULL, FIR = FALSE, fFIR = .99,
                       aFIR = 1/8, summary_list = TRUE, PCC_PLOT = TRUE, PriorPosterior_PLOT = FALSE,
                       historical_data_PLOT = FALSE, pdf_report = FALSE, path_pdf_report = tempdir(),
                       xlab = "Observations", ylab = "Quality characteristic Values",
                       main = "PCC Binomial with unknown probability" )
{
  ### Initial checks before procceeding to the main body of function
  ### Mainly this chunk of code will correspond to invalid general input before running stuff
  # 'data' (i) not defined (ii) not in vector (iii) contain non-numeric value
  if ( is.null(data) ) {
    stop("'data' have not been defined")
  } else { if ( any(!is.numeric((unlist(data)))) | any(unlist(data)<0) | any(unlist(data)%%1!=0) ) stop("Invalid 'data' input")
    if ( !is.vector(data) ) stop("'data' must be in vector form")
  }
  # 'historical_data' (i) not in vector (ii) contain non-numeric value
  if ( !is.null(historical_data) ) {
    if ( any(!is.numeric((unlist(historical_data)))) | any(unlist(historical_data)<0) | any(unlist(historical_data)%%1!=0) ) stop("Invalid 'historical_data' input")
    if ( !is.vector(historical_data) ) stop("'historical data' must be in vector form")
  }


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

  # 'FIR' (i) logical (ii) fFIR - aFIR conditions
  if ( length(unlist(FIR))>1 ) {
    message("More than one value for 'FIR', the first one will only be used")
    if ( !is.logical(FIR[1]) ) { stop("Invalid 'FIR' value ; 'FIR' must be logical") } else { FIR <- FIR[1] }
  } else {
    if ( !is.logical(FIR) ) { stop("Invalid 'FIR' value ; 'FIR' must be logical") }
  }

  # fFIR - aFIR conditions if  FIR
  if ( FIR ) {
    if ( !missing(fFIR) ) {
      if ( length(unlist(fFIR))>1 ) {
        message("More than one value for 'fFIR', the first one will only be used")
        if ( !is.numeric(fFIR[1]) | fFIR[1]<=0 | fFIR[1]>=1 ) {
          stop("Invalid 'fFIR' value")
        } else { fFIR <- fFIR[1] }
      } else {
        if ( !is.numeric(fFIR) | fFIR<=0 | fFIR>=1 ) {
          stop("Invalid 'fFIR' value")
        }
      }
    }

    if ( !missing(aFIR) ) {
      if ( length(unlist(aFIR))>1 ) {
        message("More than one value for 'aFIR', the first one will only be used")
        if ( !is.numeric(aFIR[1]) | aFIR[1]<=0 ) {
          stop("Invalid 'aFIR' value")
        } else { aFIR <- aFIR[1] }
      } else {
        if ( !is.numeric(aFIR) | aFIR<=0 ) {
          stop("Invalid 'aFIR' value")
        }
      }
    }
  }

  ### Setting the False Alarm Probability & False Alarm Rate based on the Sidak correction
  # data length
  N <- length(data)
  # If both ARL_0 and FAP chosen
  if ( !is.null(ARL_0) & !is.null(FAP) ) {
    message("Both ARL_0 and FAP are defined as input, so ARL_0 is used by default. \nIn order to use FAP instead, set ARL_0 = NULL")
    FAR <- 1/ARL_0
    # If only FAP is chosen
  } else if ( is.null(ARL_0) & !is.null(FAP) ) {
    FAR <- 1-(1-FAP)^(1/(N-1))
    # If only ARL0 is chosen
  } else if ( !is.null(ARL_0) & is.null(FAP) ){
    FAR <- 1/ARL_0
  }
  # If FIR PCC is chosen - default value for f=0.99
  if ( FIR ) {
    tf <- 1:N
    Afir <- c((  1- (1-fFIR)^(1+aFIR*(tf-1)) ) )
    FAR <- 1-(1-FAR)*Afir
  }


  ###############################################################
  ###############################################################
  ## START (1) Only this bit changes from function to function ##
  ## In some cases 'data' and 'historical_data' restrictions     ##
  ## change as well at the beginning of the function           ##
  ###############################################################
  ###############################################################

  # Prior parameter input (i) more than one value for parameters (ii) non-numeric input
  if ( is.null(n) ) { stop("'number of trials - n' have not been defined")
  } else {
    if ( !is.vector(n) ) { stop("'number of trials - n' must be in vector form")
    } else {
      if ( length(n)!=length(data) ) { stop("Vector of 'number of trials - n' must have the same length as 'data'")
      } else {
        if ( any(!is.numeric((unlist(n)))) ) { stop("Invalid 'number of trials - n' input")
        } else { if( any((unlist(n)<=0)) ) stop("Invalid 'number of trials - n' input, n must be positive")
                 if( any((unlist(data)>unlist(n))) ) stop("Invalid 'number of trials - n' input, some 'data' larger than 'n'") }
      }
    }
  }

  if( !missing(a0) ) {
    if ( length(unlist(a0))>1 ) { message("More than one value for 'a0', the first one will only be used")
      if ( !is.numeric(a0) | a0<=0 ) { stop("Invalid 'a0' value") } else { a0 <- a0[1] }
    } else { if ( !is.numeric(a0) | a0<=0 ) { stop("Invalid 'a0' value") } }
  }

  if( !missing(b0) ) {
    if ( length(unlist(b0))>1 ) { message("More than one value for 'b0', the first one will only be used")
      if ( !is.numeric(b0) | b0<=0 ) { stop("Invalid 'b0' value") } else { b0 <- b0[1] }
    } else { if ( !is.numeric(b0) | b0<=0 ) { stop("Invalid 'b0' value") } }
  }

  ### Main body of function - PCC illustration - USING FAR (or FAP equivelantly)
  ## Histotic data and processing
  if ( !is.null(historical_data) ){

    if ( is.null(historical_n) ) { historical_n <- rep(1, times=length(historical_data))
    } else {
      if ( !is.vector(historical_n) ) { stop("'historical_n - number of trials' must be in vector form")
      } else {
        if ( length(historical_n)!=length(historical_data) ) { stop("Vector of 'historical_n - number of trials' must have the same length as 'historical_data'")
        } else {
          if ( any(!is.numeric((unlist(historical_n)))) ) { stop("Invalid 'historical_n - number of trials' input")
          } else { if( any((unlist(historical_n)<=0)) ) stop("Invalid 'historical_n - number of trials' input, historical_n must be positive")
                   if( any((unlist(historical_data)>unlist(historical_n))) ) stop("Invalid 'historical_n - number of trials' input, some 'historical_data' larger than 'historical_n'") }
        }
      }
    }

    N_historicaldata <- length(historical_data)
    # Check about alpha_0
    # If no chosen value for alpha_0 use default setting
    if (is.null(alpha_0)) { alpha_0 <-1/N_historicaldata
    } else {
      if ( length(unlist(alpha_0))>1 ) {
        message("More than one value for 'alpha_0', the first one will only be used")
        if ( !is.numeric(alpha_0) | alpha_0<0 | alpha_0>1) { stop("Invalid 'alpha_0' value")
        } else { alpha_0 <- alpha_0[1] }
      } else { if ( !is.numeric(alpha_0) | alpha_0<0 | alpha_0>1 ) { stop("Invalid 'alpha_0' value") } }
    }
    # Process historical data
    # Power Prior parameters
    a0_PowerP <- a0 + alpha_0*sum(historical_data)
    b0_PowerP  <- b0 + alpha_0*( sum(historical_n) - sum(historical_data) )
    # Keep similar notation as input
    a0 <- a0_PowerP ; b0 <- b0_PowerP
  }

  ### PCC implementation
  # Sum of observations
  dataSum <- cumsum(data)[seq(1, length(data))]
  # Sum of rates
  sumN <- cumsum(n)[ seq(1, length(data)) ]
  # Posterior parameters
  a0_Post <- a0 + dataSum
  b0_Post <- b0 + sumN - dataSum


  # Control limits
  if (!FIR) { CL <- t( mapply( function(NN, A0, B0, FD = FAR) {
    hd <- betabinom_HM( cover = 1-FD, n = NN, a = A0, b = B0 )
    return( c( hd$lower.bound, hd$upper.bound ) )
  }, NN = n, A0 = a0_Post, B0 = b0_Post ) )
  } else { CL <- t( mapply( function(NN, A0, B0, FD) {
    hd <- betabinom_HM( cover = 1-FD, n = NN, a = A0, b = B0 )
    return( c( hd$lower.bound, hd$upper.bound ) )
  }, NN = n, A0 = a0_Post, B0 = b0_Post, FD = FAR ) )
  }


  CL <- rbind( c(NA, NA), CL )
  CL <- CL[-nrow(CL), ]


  ## Prior Posterior plot
  if ( PriorPosterior_PLOT | pdf_report ) {


    PrPostPLOT <-
      ggplot( data.frame(x = c(0, 1)), aes(x = c(0, 1)) ) +
      stat_function(fun = dbeta, args = list(shape1 = a0_Post[N], shape2 = b0_Post[N]), aes(colour = "Posterior", linetype = "Posterior"), size = 1) +
      {if(a0 == 1/2 & b0 == 1/2) { stat_function(fun = function(x) { ( 1/( sqrt(x * (1-x)) ) ) / (dbeta( .95 * (a0_Post[N]/(a0_Post[N] + b0_Post[N])), shape1 = a0_Post[N], shape2 = b0_Post[N] )) }, aes(colour = "Prior", linetype = "Prior"), size = 1)
      } else { stat_function(fun = dbeta, args = list(shape1 = a0, shape2 = b0), aes(colour = "Prior", linetype = "Prior"), size = 1) } } +
      scale_x_continuous(name = "") +
      scale_y_continuous(name = "Density") +
      scale_linetype_manual(values = c("solid", "dashed"), guide = "none") +
      scale_colour_manual(values = c("#3CB371", "#FF4500"),
                                   labels = c( bquote("Prior: Beta(" ~ theta ~ "|" ~ .(round(a0, digits = 1)) ~ ", " ~ .(round(b0, digits = 1)) ~ ")" ),
                                               bquote("Posterior: Beta(" ~ theta ~ "|" ~ .(round(a0_Post[N], digits = 1)) ~ ", " ~ .(round(b0_Post[N], digits = 1)) ~ ")" )),
                                   guide = guide_legend(override.aes = list( color = c("#FF4500", "#3CB371"),
                                                                             linetype = c("dashed", "solid"),
                                                                             size = c(.5, .5)), title = NULL)) +
      ggtitle(expression(atop("PCC Binomial likelihood - unknown probability"~theta, "Prior/Posterior distribution"))) +
      {if(a0 == 1/2 & b0 == 1/2) {
        geom_point( aes(x = a0_Post[N]/(a0_Post[N] + b0_Post[N]), y = 0), color = "#3CB371",  show.legend = FALSE, shape = 4, size = 3, stroke = 1.5, na.rm = TRUE )
      } else { geom_point(aes(x = c(a0/(a0 + b0), a0_Post[N]/(a0_Post[N] + b0_Post[N])), y = c(0, 0)), color = c("#FF4500", "#3CB371"),  show.legend = FALSE, shape = 4, size = 3, stroke = 1.5, na.rm = TRUE) } } +
      {if(a0 == 1/2 & b0 == 1/2) {
        annotate("text", x = a0_Post[N]/(a0_Post[N] + b0_Post[N]), y = 0, label = paste(expression(mu[post])), color = "#3CB371", size = 6, parse = TRUE, vjust = 1.25)
      } else { annotate("text", x = c(a0/(a0 + b0), a0_Post[N]/(a0_Post[N] + b0_Post[N])), y = c(0, 0), label = paste(expression(mu[prior], mu[post])),
                                 color = c("#FF4500", "#3CB371"), size = 6, parse = TRUE, vjust = 1.25) } } +
      theme( legend.position = "bottom",
             axis.line = element_line(size=1, colour = "black"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             panel.background = element_blank(),
             plot.title = element_text(size = 18, hjust = 0.5),
             text = element_text(size = 15),
             axis.text.x = element_text(colour="black", size = 12),
             axis.text.y = element_text(colour="black", size = 12))

    if ( PriorPosterior_PLOT ) { print(PrPostPLOT) }

  }


  ####################################################################
  ####################################################################
  ## END (1) Only the above bit changes from function to function   ##
  ####################################################################
  ####################################################################


  ## Output
  # Construction of 'In' and 'Out' of control column for return results
  States <- rep("", times = N)
  States[ifelse(data < CL[, 1], TRUE, FALSE)] <- "Alarm (LL)" ; States[ifelse(data > CL[, 2], TRUE, FALSE)] <- "Alarm (UL)"
  # Return results
  PCC_summary <- data.frame( data = data, HPrD_LL = CL[, 1], HPrD_UL= CL[, 2], Alarms = States )


  ## Dynamic recalculation of PCC plot's y axis
  # PCC y axis limits allowance
  Ratio <- ( PCC_summary$HPrD_UL-PCC_summary$HPrD_LL)/min(PCC_summary$HPrD_UL-PCC_summary$HPrD_LL, na.rm = T )

  # Y axis limits
  AdjustedYlim <- c( min(PCC_summary$data, PCC_summary$HPrD_LL[which(Ratio<=2.5)], na.rm = T),
                     max(PCC_summary$data, PCC_summary$HPrD_UL[which(Ratio<=2.5)], na.rm = T) )


  ### Output of function
  ## PCC plot
  if ( PCC_PLOT | pdf_report ) {
    # Creation of PCC plot
    PCC_PlotSummary <- cbind( Indices = 1:N, PCC_summary )
    PCC <- ggplot( PCC_PlotSummary, aes(PCC_PlotSummary[, "Indices"], PCC_PlotSummary[, "data"]) ) +
      geom_line( aes(x = PCC_PlotSummary[, "Indices"], y = PCC_PlotSummary[, "data"]), na.rm = TRUE ) +
      geom_line( aes(x = PCC_PlotSummary[, "Indices"], y = PCC_PlotSummary[, "HPrD_UL"]), color = "red", linetype = "solid", size = 1, na.rm = TRUE ) +
      geom_line( aes(x = PCC_PlotSummary[, "Indices"], y = PCC_PlotSummary[, "HPrD_LL"]), color = "red", linetype = "solid", size = 1, na.rm = TRUE ) +
      geom_ribbon( aes(x = PCC_PlotSummary[, "Indices"], ymin = PCC_PlotSummary[, "HPrD_UL"], ymax = PCC_PlotSummary[, "HPrD_LL"], fill = TRUE), alpha = 0.25, show.legend = FALSE ) +
      scale_fill_manual( values = c("TRUE" = "green") ) +
      geom_point( aes(group = PCC_PlotSummary[, "Indices"], color = as.factor(PCC_PlotSummary[, "Alarms"]), stroke = 1.5), show.legend = FALSE, na.rm = TRUE ) +
      scale_color_manual( values = c("black", "red", "red"), na.value = "black" ) +
      coord_cartesian( ylim = AdjustedYlim ) +
      labs( title = main, x = xlab, y = ylab ) +
      theme( legend.position = "top",
             legend.title = element_blank(),
             axis.line = element_line( colour = "black", size = 0.5, linetype = "solid" ),
             panel.background = element_blank(),
             plot.title = element_text( hjust = 0.5 ) )

    # Creation of PCC plot if historical data are chosen to be on the plot
    if ( !is.null(historical_data) ) {
      PCC_summary_historicaldata <- data.frame(  data = c(historical_data, data), HPrD_LL = c(rep(NA, times = N_historicaldata), CL[, 1]),
                                                 HPrD_UL = c(rep(NA, times = N_historicaldata), CL[, 2]), Alarms = c(rep("", times = N_historicaldata), States)  )
      PCC_PlotSummaryHist <- cbind( Indices = c(-N_historicaldata:(-1), 1:N),
                                    TypeOfdata = c(rep("Historical", times = N_historicaldata), rep("Current", times = N)),
                                    PCC_summary_historicaldata )
      PCC_historical <- ggplot( PCC_PlotSummaryHist, aes(PCC_PlotSummaryHist[, "Indices"], PCC_PlotSummaryHist[, "data"]) ) +
        geom_line( aes(x = PCC_PlotSummaryHist[, "Indices"], y = PCC_PlotSummaryHist[, "data"], linetype = as.factor(PCC_PlotSummaryHist[, "TypeOfdata"])), na.rm = TRUE ) +
        geom_segment( aes(x = 0, y = min(PCC_PlotSummaryHist[, "HPrD_LL"], na.rm = TRUE), xend = 0, yend = max(PCC_PlotSummaryHist[, "HPrD_UL"], na.rm = TRUE)) ) +
        geom_line( aes(x = PCC_PlotSummaryHist[, "Indices"], y = PCC_PlotSummaryHist[, "HPrD_UL"]), color = "red", linetype = "solid", size = 1, na.rm = TRUE ) +
        geom_line( aes(x = PCC_PlotSummaryHist[, "Indices"], y = PCC_PlotSummaryHist[, "HPrD_LL"]), color = "red", linetype = "solid", size = 1, na.rm = TRUE ) +
        geom_ribbon( aes(x = PCC_PlotSummaryHist[, "Indices"], ymin = PCC_PlotSummaryHist[, "HPrD_UL"], ymax = PCC_PlotSummaryHist[, "HPrD_LL"], fill = TRUE), alpha = 0.25, show.legend = FALSE ) +
        scale_fill_manual( values = c("TRUE" = "green") ) +
        geom_point( aes(group = PCC_PlotSummaryHist[, "Indices"], shape = as.factor(PCC_PlotSummaryHist[, "TypeOfdata"]), color = as.factor(PCC_PlotSummaryHist[, "Alarms"]), stroke = 1.5),
                    show.legend = FALSE, na.rm = TRUE ) +
        scale_color_manual( values=c("black", "red", "red"), na.value = "black" ) +
        scale_linetype_manual( values = c("Historical" = "dotted", "Current" = "solid") ) +
        scale_shape_manual( values = c("Historical" = 1, "Current" = 19) ) +
        coord_cartesian( ylim = AdjustedYlim ) +
        labs( title = main, x = xlab, y = ylab ) +
        theme( legend.position = "top",
               legend.title = element_blank(),
               axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
               panel.background = element_blank(),
               plot.title = element_text( hjust = 0.5 ) )
    }

    if ( historical_data_PLOT ) { print(PCC_historical)
    } else { if ( PCC_PLOT) { print(PCC) } }

  }

  # List of results
  if ( summary_list ) { return(PCC_summary) }


  # List of results return in pdf
  if ( pdf_report ) {

    # save pdf
    pdf(
      paste0( path_pdf_report, "\\", "PCC_results_", paste0( unlist(strsplit(date(), " "))[c(1,2,3,5)], collapse = "_" ), "_",
              paste0( unlist(strsplit( unlist(strsplit(date(), " "))[4], ":" )), collapse = "." ),
              ".pdf" ),
      height = 8.264, width = 11.694 )

    # PCC plot on pdf
    print(PCC)

    # Prior Posterior plot on pdf
    print(PrPostPLOT)

    # Results matrix on pdf
    # Chunk of code to split results matrix to different pages - Set a default number based on pdf height/width
    NRowsPerPage <- 25
    if(NRowsPerPage > nrow(PCC_summary)){ FloatingRow <- nrow(PCC_summary) } else { FloatingRow <- NRowsPerPage }
    sapply(1:ceiling(nrow(PCC_summary)/NRowsPerPage), function(index) {
      if (index==1) { StartingRow <- 1 }
      grid.newpage()
      grid.table(PCC_summary[StartingRow:FloatingRow, ])
      StartingRow <<- FloatingRow + 1
      if( sum(NRowsPerPage, FloatingRow) < nrow(PCC_summary)){ FloatingRow <<-  NRowsPerPage + FloatingRow } else { FloatingRow <<- nrow(PCC_summary) }
    })

    dev.off()

  }


}




