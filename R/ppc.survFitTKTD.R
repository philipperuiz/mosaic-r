#' Posterior predictive check plot for survFitTKTD objects
#'
#' The \code{ppc} function plot the observed versus predicted values for the
#' \code{survFitTKTD} objects.
#'
#' @param x An object of class \code{survFitTKTD}
#' @param style Graphical package method: \code{generic} or \code{ggplot}.
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @examples
#'
#' # (1) Load the data
#' data(cadmium1)
#'
#' # (2) Create an object of class "survData"
#' dat <- survData(cadmium1)
#'
#' \dontrun{
#' # (3) Run the survFitTKTD function with the TKTD model
#' out <- survFitTKTD(dat)
#'
#' # (4) Plot observed versus predicted values
#' ppc(out)
#' }
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot
#' 
#' @export
ppc.survFitTKTD <- function(x, style = "generic", ...) {
  if (!is(x, "survFitTKTD"))
    stop("x is not of class 'survFitTKTD'!")
  
  xlab <- "Observed Nbr. of survivor"
  ylab <- "Predicted Nbr. of survivor"
  
  ppc_gen(EvalsurvTKTDPpc(x), style, xlab, ylab)
}

ppc_gen <- function(tab, style, xlab, ylab) {
  
  if (style == "generic") PpcGeneric(tab, xlab, ylab)
  else if (style == "ggplot") PpcGG(tab, xlab, ylab)
  else stop("Unknown style")
}

#' @importFrom stats rbinom quantile
EvalsurvTKTDPpc <- function(x) {

  tot.mcmc <- do.call("rbind", x$mcmc)

  ke <- 10^sample(tot.mcmc[, "log10ke"], 5000)
  ks <- 10^sample(tot.mcmc[, "log10ks"], 5000)
  nec <- 10^sample(tot.mcmc[, "log10NEC"], 5000)
  m0 <- 10^sample(tot.mcmc[, "log10m0"], 5000)
  
  n <- x$jags.data$ndat
  xconc <- x$jags.data$x
  t <- x$jags.data$t
  tprec <- x$jags.data$tprec
  Nprec <- x$jags.data$Nprec
  NsurvObs <- x$jags.data$y
  Nprec <- x$jags.data$Nprec
  bigtime <- x$jags.data$bigtime
  NsurvPred <- matrix(NA, nrow = 5000, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:length(ke)) {
      xcor <- ifelse(xconc[i] > 0, xconc[i], 10)
      R <- ifelse(xconc[i] > nec[j], nec[j]/xcor, 0.1)
      tNEC <- ifelse(xconc[i] > nec[j], -1 / ke[j] * log(1 - R), bigtime)
      tref <- max(tprec[i], tNEC)
      psurv <- exp(-m0 * (t[i] - tprec[i]) +
                     if (t[i] > tNEC) {
                       -ks * ((xconc[i] - nec[j]) * (t[i] - tref) +
                                xconc[i]/ke[j] * (exp(-ke[j] * t[i]) - exp(-ke[j] * tref)))
                     } else {
                       0
                     })
    }
    NsurvPred[, i] <- rbinom(5000, Nprec[i], psurv)
  }

  QNsurvPred <- t(apply(NsurvPred, 2, quantile,
                        probs = c(2.5, 50, 97.5) / 100, na.rm = TRUE))
  tab <- data.frame(QNsurvPred,
                    Nprec, NsurvObs,
                    col = ifelse(QNsurvPred[,"2.5%"] > NsurvObs |
                                   QNsurvPred[,"97.5%"] < NsurvObs,
                                 "red", "green"))
  colnames(tab) <- c("P2.5", "P50", "P97.5", "Nprec", "Obs", "col")
  
  return(tab)
}

stepCalc <- function(obs_val) {
  # calculation of steps coordinate
  sObs <- sort(obs_val)
  stepX <- c(0, sapply(2:length(sObs), function(i) {
    sObs[i-1] + (sObs[i] - sObs[i-1]) / 2}), max(sObs))
  return(list(sObs = sObs, stepX = stepX))
}

jitterObsGenerator <- function(stepX, tab, obs_val) {
  # uniform jittering of observed values
  allSpaceX <- sapply(2:length(stepX),
                      function(i) { stepX[i] - stepX[i-1] })
  spaceX <- min(allSpaceX[which(allSpaceX != 0)])/2
  lengthX <- table(tab[, "Obs"])
  jitterObs <- mapply(function(x, y) {
    if (y == 1) {
      seq(x, x, length.out = y)
    } else {
      seq(x - spaceX + (2 * spaceX / (y + 1)),
          x + spaceX - (2 * spaceX / (y + 1)), length.out = y)
    }
  }, x = sort(obs_val), y = lengthX)
  return(list(spaceX = spaceX,
              jitterObs = unlist(jitterObs)))
}

#' @importFrom graphics abline segments
PpcGeneric <- function(tab, xlab, ylab) {
  obs_val <- unique(tab[, "Obs"])
  sObs <- stepCalc(obs_val)$sObs
  stepX <- stepCalc(obs_val)$stepX
  jittered_obs <- jitterObsGenerator(stepX, tab, obs_val)$jitterObs
  spaceX <- jitterObsGenerator(stepX, tab, obs_val)$spaceX
  
  plot(c(0, max(tab[, "P97.5"])),
       c(0, max(tab[, "P97.5"])),
       type = "n",
       xlab = xlab,
       ylab = ylab,
       xaxt = "n",
       yaxt = "n")
  
  # axis
  axis(side = 2, at = if (max(tab[, "Obs"]) == 1) {
    c(0, 1)
  } else {
    pretty(c(0, max(tab[, "P97.5"])))
  })
  axis(side = 1, at = if (max(tab[, "Obs"]) == 1) {
    c(0, 1)
  } else {
    pretty(c(0, max(tab[, "P97.5"])))
  })
  
  if (max(sObs) < 20) {
    sapply(1:length(sObs), function(i) {
      segments(sObs[i] - (spaceX * 1.25), sObs[i],
               sObs[i] + (spaceX * 1.25), sObs[i])
    })
  } else {
    abline(0, 1)
  }
  
  tab0 <- tab[order(tab$Obs),]
  delta <- 0.01 * (max(obs_val) - min(obs_val))
  segments(jittered_obs, tab0[, "P2.5"],
           jittered_obs, tab0[, "P97.5"],
           col = as.character(tab0[, "col"]))
  segments(jittered_obs - delta, tab0[, "P2.5"],
           jittered_obs + delta, tab0[, "P2.5"],
           col = as.character(tab0[, "col"]))
  segments(jittered_obs - delta, tab0[, "P97.5"],
           jittered_obs + delta, tab0[, "P97.5"],
           col = as.character(tab0[, "col"]))
  
  points(jittered_obs, tab0[, "P50"],
         pch = 16)
}

#' @import ggplot2
#' @importFrom  grid arrow unit
PpcGG <- function(tab, xlab, ylab) {
  obs_val <- unique(tab[, "Obs"])
  sObs <- stepCalc(obs_val)$sObs
  stepX <- stepCalc(obs_val)$stepX
  jittered_obs <- jitterObsGenerator(stepX, tab, obs_val)$jitterObs
  spaceX <- jitterObsGenerator(stepX, tab, obs_val)$spaceX
  
  tab0 <- cbind(tab[order(tab$Obs),], jittered_obs)
  
  df <- data.frame(sObs, spaceX)
  
  if (max(sObs) < 20) {
    gf1 <- ggplot(df) +
      geom_segment(aes(x = sObs - (spaceX * 1.25),
                       xend = sObs + (spaceX * 1.25),
                       y = sObs, yend = sObs)) +
      scale_x_continuous(breaks = c(0, tab0[, "Obs"]),
                         labels = c(0, tab0[, "Obs"])) +
      scale_y_continuous(breaks = c(0, tab0[, "Obs"]),
                         labels = c(0, tab0[, "Obs"]))
  } else {
    gf1 <- ggplot(tab0) +
      geom_abline(intercept = 0, slope = 1)
  }
  
  gf2 <- gf1 +
    geom_segment(aes(x = jittered_obs, xend = jittered_obs,
                     y = P2.5, yend = P97.5), data = tab0,
                 arrow = arrow(length = unit(0.1, "cm"), angle = 90,
                               ends = "both"),
                 color = tab0$col) +
    geom_point(aes(x = jittered_obs, y = P50), tab0) +
    expand_limits(y = 0) +
    expand_limits(x = 0) +
    labs(x = xlab, y = ylab) +
    theme_minimal()
  
  return(gf2)
}
