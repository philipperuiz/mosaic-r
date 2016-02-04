#' Cheks the convergence of the MCMC chains
#' 
#' The \code{convergence} function checks the convergence of the MCMC chains
#' from the JAGS estimate with the Gelman and Rubin convergence diagnostic
#' (Gelman and Rubin, 1992). It summarizes the \code{mcmc} or \code{mcmc.list}
#' object with a trace of the sampled output, a density estimate and an
#' autocorrelation plot for each parameter in the chain.
#' 
#' @param out An object of class \code{reproFitTT}, \code{survFitTT} or 
#' \code{survFitTKTD}.
#' @param trace If \code{TRUE}, the function traces the sampled output estimate
#' for each parameter in the chain.
#' @param density If \code{TRUE}, the function plots the density estimate for
#' each parameter in the chain.
#' @param autocorr If \code{TRUE}, the function plots the autocorrelation for
#' each parameter in each chain.
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' 
#' @return The function returns an object of class list with the point estimate
#' of the multivariate potential scale reduction factor and the point estimate
#' of the potential scale reduction factor (Rhat) for each parameter of the
#' Gelman and Rubin test (Gelman and Rubin, 1992). A value close to 1 is
#' expected when convergence is reached. See the \code{\link[coda]{gelman.diag}}
#' help for more details.
#' 
#' @note When \code{style = "ggplot"}, the function calls packages \code{ggmcmc}
#' and \code{gridExtra} and returns a graphical object of class \code{ggplot}.
#' 
#' @author Marie Laure Delignette-Muller
#' <marielaure.delignettemuller@@vetagro-sup.fr>, Philippe Ruiz
#' <philippe.ruiz@@univ-lyon1.fr>
#' 
#' @seealso \code{\link{reproFitTt}} and \code{\link[coda]{gelman.diag}},
#' \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{autocorr.plot}} from the
#' \code{rjags} package and \code{\link[ggmcmc]{ggs_traceplot}},
#' \code{\link[ggmcmc]{ggs_density}} and
#' \code{\link[ggmcmc]{ggs_autocorrelation}} from the \code{ggmcmc} package
#' (\url{http://xavier-fim.net/packages/ggmcmc})
#' 
#' @references Gelman, A. and Rubin, D.B. (1992) \emph{Inference from iterative
#' simulation using multiple sequences}, Statistical Science, 7, 457-511.
#' 
#' @keywords mcmc-analysis
#' 
#' @examples
#' 
#' # (1) Load the data
#' data(zinc)
#' 
#' # (2) Create an object of class "reproData"
#' dat <- reproData(zinc)
#' 
#' \dontrun{
#' # (3) Run the reproFitTT function
#' out <- reproFitTT(dat)
#' 
#' # (4) Check the convergence
#' convergence(out, trace = TRUE, density = FALSE,
#' autocorr = TRUE)
#' 
#' # (5) Check the convergence using the "ggmcmc" package
#' convergence(out, trace = TRUE, density = TRUE, 
#' autocorr = TRUE, style = "ggplot")
#' }
#' 
#' @export
#' 
#' @importFrom coda autocorr.plot gelman.diag
#' @importFrom ggmcmc ggs ggs_traceplot ggs_density ggs_autocorrelation
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' 
convergence <- function(out,
                        trace = TRUE,
                        density = TRUE,
                        autocorr = TRUE,
                        style = "generic") {
  # test class object
  if (class(out) != "reproFitTT" && class(out) != "survFitTT" &&
      class(out) != "survFitTKTD")
    stop("The object passed in argument 'out' is not of class 'reproFitTt', 'survFitTt' or 'survFitTKTD' !\n")
  
  # keep MCMC list
  mcmc <- out$mcmc
  #  Gelamn and Rubin diagnostic
  GelRubmulti <- gelman.diag(mcmc)$mpsrf
  GelRubesti <- gelman.diag(mcmc)$psrf[, "Point est."]
  
  # return the psrf and mprsf value
  cat("Gelman and Rubin:\n")
  cat("Potential scale reduction factor for each parameter:\n")
  print(GelRubesti)
  cat("\nMultivariate potential scale reduction factor:\n", GelRubmulti,"\n")
  
  # generic plot
  if (style == "generic") {
    # trace and density
    if (trace == TRUE || density == TRUE) {
      plot(mcmc, trace = trace, density = density)
    }
    # autocorrelation
    if (autocorr == TRUE) {
      if (trace == TRUE || density == TRUE){
        if (Sys.getenv("RSTUDIO") == "") dev.new() # create a new page plot
        # when not use RStudio
      }
      autocorr.plot(mcmc, ask = TRUE)
    }
  }
  
  # ggplot
  if (style == "ggplot") {
    # creat ggs objects
    D <- ggs(mcmc)
    if (trace == TRUE) trp <- ggs_traceplot(D) # trace plot
    if (density == TRUE) dns <- ggs_density(D) # density plot
    if (autocorr == TRUE) atc <- ggs_autocorrelation(D) #autocorr plot
    if(trace == FALSE || density == FALSE || autocorr == FALSE) {
      blank <- grid.rect(gp = gpar(col = "white"))
    }
    
    
    
    do.call(grid.arrange, list(if (trace == TRUE) trp else blank,
                               if (density == TRUE) dns else blank,
                               if (autocorr == TRUE) atc else blank,
                               ncol = 2))
  }
  
  return(invisible(list(mpsrf = GelRubmulti,
                        psrf = GelRubesti)))
}
