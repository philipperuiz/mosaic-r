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

