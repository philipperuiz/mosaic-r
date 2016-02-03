#' Plotting method for survFitTKTD objects
#' 
#' @param x An object of class \code{survFitTKTD}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Concentrations}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate}.
#' @param main A main title for the plot.
#' @param ci if \code{TRUE}, draws the 95 \% credible limits of the fitted curve.
#' @param one.plot if \code{TRUE}, draw all the estimeted curves in one plot.
#' @param addlegend if \code{TRUE}, adds a default legend to the plot.
#' @param style Graphical method: \code{generic} or \code{ggplot}.
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @keywords plot 
#' @export
#' 
#' @import ggplot2
#' @import grDevices
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#' 
plot.survFitTKTD <- function(x,
                             xlab = "Time",
                             ylab = "Survival rate",
                             main = NULL,
                             ci = FALSE,
                             one.plot = TRUE,
                             addlegend = FALSE,
                             style = "generic", ...) {
  
  data <- survFitPlotDataTKTD(x)
  
  conf.int <- if (ci) { survTKTDConfInt(x) } else NULL
  if (ci) data$dobs <- cbind(data$dobs, qinf95 = conf.int["qinf95",],
                             qsup95 = conf.int["qsup95",])
  dataCI <- if (ci && !one.plot) { survFitPlotCITKTD(x) } else NULL
  dataCIm <- if (ci && !one.plot) { melt(dataCI,
                                         id.vars = c("conc", "time")) } else NULL
  
  if (style == "generic") {
    survFitPlotTKTDGeneric(data, xlab, ylab, main, one.plot, ci, dataCIm,
                           addlegend)
  }
  else if (style == "ggplot") {
    survFitPlotTKTDGG(data, xlab, ylab, main, one.plot, ci, dataCI, dataCIm,
                      addlegend)
  }
  else stop("Unknown style")
}

Surv <- function(Cw, time, ks, ke, NEC, m0)
  # Fonction S ecrite en R pour la validation en simu ensuite
  # Cw est la concentration dans le milieu
{
  S <- exp(-m0*time) # survie de base avec mortalite naturelle seule
  if (Cw > NEC) {
    tNEC <- -(1/ke)*log(1 - NEC/Cw)
    if (time > tNEC) {
      # ajoute de la mortalite due au toxique
      S <- S * exp( ks/ke*Cw*(exp(-ke*tNEC) -exp(-ke*time))
                    - ks*(Cw-NEC)*(time - tNEC) )
    }
  }
  return(S)
}

Survm00 <- function (Cw, time, ks, ke, NEC)
  # Fonction S ecrite en R pour la validation en simu ensuite
  # Cw est la concentration dans le milieu
{
  S <- 1 # survie de base avec mortalite naturelle seule
  if (Cw > NEC) {
    tNEC <- -(1/ke)*log(1 - NEC/Cw)
    if (time > tNEC) {
      # ajoute de la mortalite due au toxique
      S <- S * exp( ks/ke*Cw*(exp(-ke*tNEC) -exp(-ke*time))
                    - ks*(Cw-NEC)*(time - tNEC) )
    }
  }
  return(S)
}

SurvkeInf <- function (Cw, time, ks, NEC, m0)
  # Fonction S ecrite en R pour la validation en simu ensuite
  # Cw est la concentration dans le milieu
{
  S <- exp(-m0*time) # survie de base avec mortalite naturelle seule
  if (Cw > NEC) {
    # ajoute de la mortalite due au toxique
    S <- S * exp(- ks*(Cw-NEC)*(time) )
  }
  return(S)
}

SurvkeInfm00 <- function (Cw, time, ks, NEC)
  # Fonction S ecrite en R pour la validation en simu ensuite
  # Cw est la concentration dans le milieu
{
  S <- 1 # survie de base avec mortalite naturelle seule
  if (Cw > NEC) {
    # ajoute de la mortalite due au toxique
    S <- S * exp(- ks*(Cw-NEC)*(time ) )
  }
  return(S)
}

survFitPlotDataTKTD <- function(x) {
  # INPUT
  # x : An object of class survFitTKTD
  # OUTPUT
  # A list of - dobs : observed values
  #           - dtheo : estimated values
  npoints <- 100
  dtheo <- data.frame(conc = numeric(), t = numeric(), psurv = numeric())
  
  concobs <- unique(x$transformed.data$conc)
  tfin <- seq(0, max(x$jags.data$t), length.out = npoints)
  
  # parameters
  ks <- x$estim.par["ks", "median"]
  if (x$ke) ke <- x$estim.par["ke", "median"]
  nec <- x$estim.par["nec", "median"]
  if (x$m0) m0 <- x$estim.par["m0", "median"]
  
  for (i in 1:length(concobs)) {
    for (j in 1:npoints) {
      if (x$ke && x$m0) {
        psurv <- Surv(Cw = concobs[i], time = tfin[j],
                      ks = ks, ke = ke,
                      NEC = nec,
                      m0 = m0)
      } else if (x$ke && !x$m0) {
        psurv <- Survm00(Cw = concobs[i], time = tfin[j],
                         ks = ks, ke = ke,
                         NEC = nec)
      } else if (!x$ke && x$m0) {
        psurv <- SurvkeInf(Cw = concobs[i], time = tfin[j],
                           ks = ks, NEC = nec,
                           m0 = m0)
      } else {
        psurv <- SurvkeInfm00(Cw = concobs[i], time = tfin[j],
                              ks = ks, NEC = nec)
      }
      
      dtheo <- rbind(dtheo, data.frame(conc = concobs[i],
                                       t = tfin[j],
                                       psurv = psurv))
    }
  }
  
  dobs <- data.frame(conc = x$transformed.data$conc,
                     t = x$transformed.data$time, 
                     psurv = x$transformed.data$N_alive / x$transformed.data$N_init)
  
  return(list(dtheo = dtheo,
              dobs = dobs))
}

#' @importFrom stats aggregate binom.test
survTKTDConfInt <- function(x) {
  # create confidente interval on observed data
  # binomial model by a binomial test
  # INPUT:
  # - x : object of class survFitTT
  # OUTPUT:
  # - ci : confidente interval
  
  ci <- apply(x$transformed.data, 1, function(x) {
    binom.test(x["N_alive"], x["N_init"])$conf.int
  })
  rownames(ci) <- c("qinf95", "qsup95")
  colnames(ci) <- x$transformed.data$conc
  
  return(ci)
}

survFitPlotCITKTD <- function(x) {
  # INPUT
  # x : An object of class survFitTKTD
  # OUTPUT
  # A list of - dobs : observed values
  #           - dtheo : estimated values
  npoints <- 100
  
  concobs <- unique(x$transformed.data$conc)
  tfin <- seq(0, max(x$jags.data$t), length.out = npoints)
  
  # prameters
  mctot <- do.call("rbind", x$mcmc)
  sel <- sample(nrow(mctot))[1:ceiling(nrow(mctot) / 50)]
  ks <- 10^mctot[, "log10ks"][sel]
  if (x$ke) ke <- 10^mctot[, "log10ke"][sel]
  if (x$m0) m0 <- 10^mctot[, "log10m0"][sel]
  nec <- 10^mctot[, "log10NEC"][sel]
  
  # all theorical
  dtheo = list()
  for (k in 1:length(concobs)) {
    dtheo[[k]] <- array(data = NA, dim = c(npoints, length(nec)))
    for (i in 1:length(nec)) {
      for (j in 1:npoints) {
        if (x$ke && x$m0) {
          dtheo[[k]][j, i] <- Surv(Cw = concobs[k], time = tfin[j],
                                   ks = ks[i], ke = ke[i],
                                   NEC = nec[i],
                                   m0 = m0[i])
        } else if (x$ke && !x$m0) {
          dtheo[[k]][j, i] <- Survm00(Cw = concobs[k], time = tfin[j],
                                      ks = ks[i], ke = ke[i],
                                      NEC = nec[i])
        } else if (!x$ke && x$m0) {
          dtheo[[k]][j, i] <- SurvkeInf(Cw = concobs[k], time = tfin[j],
                                        ks = ks[i], NEC = nec[i],
                                        m0 = m0[i])
        } else {
          dtheo[[k]][j, i] <- SurvkeInfm00(Cw = concobs[k], time = tfin[j],
                                           ks = ks[i], NEC = nec[i])
        }
      }
    }
  }
  
  dtheof <- do.call("rbind", dtheo)
  dtheof <- as.data.frame(cbind(rep(concobs, rep(npoints, length(concobs))),
                                rep(tfin, length(concobs)),
                                dtheof))
  names(dtheof) <- c("conc", "time", paste0("X", 1:length(sel)))
  
  # quantile
  qinf95 = NULL
  qsup95 = NULL
  q50 = NULL
  
  for (i in 1:dim(dtheof)[1]) {
    qinf95[i] <- quantile(dtheof[i, 3:length(dtheof)],
                          probs = 0.025, na.rm = TRUE)
    qsup95[i] <- quantile(dtheof[i, 3:length(dtheof)],
                          probs = 0.975, na.rm = TRUE)
    q50[i] <- quantile(dtheof[i, 3:length(dtheof)],
                       probs = 0.5, na.rm = TRUE)
  }
  
  dtheof <- cbind(qinf95, qsup95, q50, dtheof)
  
  return(dtheof)
}

survFitPlotTKTDGeneric <- function(data, xlab, ylab, main, one.plot, ci, dataCIm,
                                   addlegend) {
  # vector color
  data[["dobs"]]$color <- as.numeric(as.factor(data[["dobs"]][["conc"]]))
  data[["dtheo"]]$color <- as.numeric(as.factor(data[["dtheo"]][["conc"]]))
  
  if (one.plot) {
    survFitPlotTKTDGenericOnePlot(data, xlab, ylab, main, addlegend)
  } else {
    par(mfrow = plotMatrixGeometry(length(unique(data[["dobs"]][["conc"]]))))
    
    survFitPlotTKTDGenericNoOnePlot(data, xlab, ylab, ci, dataCIm)
    
    par(mfrow = c(1, 1))
  }
}

survFitPlotTKTDGenericOnePlot <- function(data, xlab, ylab, main, addlegend) {
  plot(data[["dobs"]][["t"]],
       data[["dobs"]][["psurv"]],
       xlab = xlab,
       ylab = ylab,
       pch = 16,
       col = data[["dobs"]]$color,
       main = main)
  # one line by replicate
  by(data[["dtheo"]], list(data[["dtheo"]]$conc),
     function(x) {
       lines(x$t, x$psurv, # lines
             col = x$color)
     })
  if (addlegend) {
    legend("bottomleft",
           legend = unique(data[["dobs"]]$conc),
           pch = 16,
           bty = "n",
           cex = 1,
           ncol = 2,
           col = unique(data[["dobs"]]$color),
           title = "Concentrations")
  }
}

survFitPlotTKTDGenericNoOnePlot <- function(data, xlab, ylab, ci, dataCIm) {
  if (ci) {
    survFitPlotTKTDGenericNoOnePlotCi(data, xlab, ylab, dataCIm)
  } else {
    survFitPlotTKTDGenericNoOnePlotNoCi(data, xlab, ylab)
  }
}

survFitPlotTKTDGenericNoOnePlotCi <- function(data, xlab, ylab, dataCIm) {
  # one line by replicate
  by(data[["dtheo"]], list(data[["dtheo"]]$conc),
     function(x) {
       plot(x[, "t"],
            x[, "psurv"],
            xlab = xlab,
            ylab = ylab,
            type = "n",
            ylim = c(0, 1),
            col = x[, "color"])
       lines(x[, "t"], x[, "psurv"], # lines
             col = x[, "color"])
     })
}

survFitPlotTKTDGenericNoOnePlotNoCi <- function(data, xlab, ylab) {
  # one line by replicate
  by(data[["dtheo"]], list(data[["dtheo"]]$conc),
     function(x) {
       plot(x[, "t"],
            x[, "psurv"],
            xlab = xlab,
            ylab = ylab,
            type = "n",
            ylim = c(0, 1),
            col = x[, "color"],
            main = paste0("Concentration : ", unique(x[, "conc"])))
       lines(x[, "t"], x[, "psurv"],
             col = x[, "color"])
     })
}

survFitPlotTKTDGG <- function(data, xlab, ylab, main, one.plot, ci, dataCI,
                              dataCIm, addlegend) {

  if (one.plot) {
    if (ci) warning("Credible intervals are only evalables in grid plot !")
    survFitPlotTKTDGGOnePlot(data, xlab, ylab, main, addlegend)
  } else {
    survFitPlotTKTDGGNoOnePlot(data, xlab, ylab, main, ci, dataCI, dataCIm, conf.int)
  }
}

survFitPlotTKTDGGOnePlot <- function(data, xlab, ylab, main, addlegend) {
  gf <- ggplot(data$dobs, aes(x = t, y = psurv, colour = factor(conc))) +
    geom_point() +
    geom_line(data = dataCIm)
    geom_line(data = data$dtheo) +
    labs(x = xlab, y = ylab) + ggtitle(main) +
    ylim(c(0, 1)) +
    theme_minimal()
  if (addlegend) {
    gf + scale_color_discrete("Concentrations")
  } else {
    gf + scale_color_discrete(guide = "none")
  }
}

survFitPlotTKTDGGNoOnePlot <- function(data, xlab, ylab, main, ci, dataCI,
                                       dataCIm, conf.int) {
  if (ci) {
    survFitPlotTKTDGGNoOnePlotCi(data, xlab, ylab, main, dataCI, dataCIm, conf.int)
  } else {
    survFitPlotTKTDGGNoOnePlotNoCi(data, xlab, ylab, main)
  }
}

survFitPlotTKTDGGNoOnePlotCi <- function(data, xlab, ylab, main, dataCI,
                                         dataCIm, conf.int) {
  ggplot(data$dobs,
         aes(x = t, y = psurv, colour = factor(conc))) +
    geom_point(color = "black") +
    geom_segment(aes(x = t, xend = t, y = qinf95, yend = qsup95),
                 arrow = arrow(length = unit(0.25, "cm"), angle = 90,
                               ends = "both"), data$obs, color = "gray",
                 size = 0.5) +
    geom_line(data = dataCIm, aes(x = time, y = value, group = variable),
              alpha = 0.05) +
    geom_line(data = dataCI, aes(x = time, y = q50), linetype = 'dashed', color = "black") +
    geom_line(data = data$dtheo, color = "red") +
    geom_line(data = dataCI, aes(x = time, y = qinf95), linetype = 'dashed', color = "black") +
    geom_line(data = dataCI, aes(x = time, y = qsup95), linetype = 'dashed', color = "black") +
    facet_wrap(~conc) +
    labs(x = xlab, y = ylab) + ggtitle(main) +
    ylim(c(0, 1)) +
    theme_minimal() +
    scale_color_discrete(guide = "none")
}

survFitPlotTKTDGGNoOnePlotNoCi <- function(data, xlab, ylab, main) {
  ggplot(data$dobs,
         aes(x = t, y = psurv, colour = factor(conc))) +
    geom_point() +
    geom_line(data = data$dtheo) +
    facet_wrap(~conc) +
    labs(x = xlab, y = ylab) + ggtitle(main) +
    ylim(c(0, 1)) +
    theme_minimal() +
    scale_color_discrete(guide = "none")
}

