#' @importFrom dplyr filter
survTKTDCreateJagsData <- function(data, distr = "norm", bond = "01",
                                   m0 = FALSE, ke = FALSE) {
  # Creates the parameters to define the prior of the TKTD model
  # INPUTS
  # data : object of class survData
  # distr : normal or uniform priors
  # bond : 0.01 and 0.99 or 0.001 and 0.999 % of the highest tested
  # concentration
  # OUTPUT
  # jags.data : list of data required for the jags.model function
  
  data <- data[data$time != 0, ]
  
  # Parameter calculation of concentration min and max
  concmin <- min(data$conc[data$conc != 0])
  concmax <- max(data$conc)
  
  tmin <- min(data$time)
  tmax <- max(data$time)
  conc <- sort(unique(data$conc))
  
  deltaCmin = NULL
  for (i in 2:length(conc)) {
    deltaCmin[i - 1] <- conc[i] - conc[i - 1]
  }
  deltaCmin <- min(deltaCmin)
  
  # ks parameters
  if (bond == "01") {
    ksmax <- -log(0.01) / (tmin * deltaCmin)
    ksmin <- -log(0.99) / (tmax * (concmax - concmin))
  } else if (bond == "001") {
    ksmax <- -log(0.001) / (tmin * deltaCmin)
    ksmin <- -log(0.999) / (tmax * (concmax - concmin))
  } else {
    stop("Wrong bond")
  }
  meanlog10ks <- (log10(ksmax) + log10(ksmin)) / 2
  sdlog10ks <- (log10(ksmax) - log10(ksmin)) / 4
  taulog10ks <- 1 / sdlog10ks^2
  
  if (ke) {
    # ke parameters
    if (bond == "01") {
      kemax <- -log(0.01) / tmin
      kemin <- -log(0.99) / tmax
    } else if (bond == "001") {
      kemax <- -log(0.001) / tmin
      kemin <- -log(0.999) / tmax
    } else {
      stop("Wrong bond")
    }
    meanlog10ke <- (log10(kemax) + log10(kemin)) / 2
    
    sdlog10ke <- (log10(kemax) - log10(kemin)) / 4
    taulog10ke <- 1 / sdlog10ke^2
  }
  
  if (m0) {
    # m0 parameters
    m0max <- -log(0.5) / tmin
    if (bond == "01") {
      m0min <- -log(0.99) / tmax
    } else if (bond == "001") {
      m0min <- -log(0.999) / tmax
    } else {
      stop("Wrong bond")
    }
    
    meanlog10m0 <- (log10(m0max) + log10(m0min)) / 2
    sdlog10m0 <- (log10(m0max) - log10(m0min)) / 4
    taulog10m0 <- 1/ sdlog10m0^2
  }
  
  # nec parameters
  meanlog10nec <- (log10(concmax) + log10(concmin))/2
  sdlog10nec <- (log10(concmax) - log10(concmin)) / 4 
  taulog10nec <- 1/ sdlog10nec^2
  
  if (distr == "unif") {
    return(list( x = data$conc, y = data$N_alive,
                 t = data$time, tprec = data$tprec,
                 t2prec = if (m0 && !ke) {data$t2prec} else NULL,
                 Nprec = data$Nprec,
                 minlog10conc = log10(concmin), maxlog10conc = log10(concmax),
                 minlog10ks = log10(ksmin), maxlog10ks = log10(ksmax),
                 minlog10ke = if (ke) {log10(kemin)} else NULL,
                 maxlog10ke = if (ke) {log10(kemax)} else NULL,
                 minlog10m0 = if (m0) {log10(m0min)} else NULL,
                 maxlog10m0 = if (m0) {log10(m0max)} else NULL,
                 ndat = length(data$conc),
                 bigtime = if (ke) {max(data$time) + 10} else NULL))
  } else if (distr == "norm") {
    return(list( x = data$conc, y = data$N_alive,
                 t = data$time, tprec = data$tprec,
                 t2prec = if (m0 && !ke) {data$t2prec} else NULL,
                 Nprec = data$Nprec,
                 meanlog10ks = meanlog10ks, taulog10ks = taulog10ks,
                 meanlog10ke = if (ke) {meanlog10ke} else NULL,
                 taulog10ke = if (ke) {taulog10ke} else NULL,
                 meanlog10m0 = if (m0) {meanlog10m0} else NULL,
                 taulog10m0 = if (m0) {taulog10m0} else NULL,
                 meanlog10nec = meanlog10nec, taulog10nec = taulog10nec,
                 ndat = length(data$conc),
                 bigtime = if (ke) {max(data$time) + 10} else NULL))
  }
}
 
modelTKTDUnif <- "model {
#########priors 
log10ks ~ dunif(minlog10ks, maxlog10ks)
log10NEC ~ dunif(minlog10conc, maxlog10conc)
log10ke ~ dunif(minlog10ke, maxlog10ke)
log10m0 ~ dunif(minlog10m0, maxlog10m0)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
ke <- 10**log10ke
m0 <- 10**log10m0

##########Computation of the likelihood
for (i in 1:ndat)
{
  tNEC[i] <- ifelse(x[i] > NEC, -1 / ke * log(1 - R[i]), bigtime)
  R[i] <- ifelse(x[i] > NEC, NEC/xcor[i], 0.1)
  xcor[i] <- ifelse(x[i] > 0, x[i], 10)
  tref[i] <- max(tprec[i], tNEC[i])
  
  psurv[i] <- exp(-m0 * (t[i] - tprec[i]) + ifelse(t[i] > tNEC[i], -ks * ((x[i] - NEC) * (t[i] - tref[i]) + x[i]/ke * ( exp(-ke * t[i]) - exp(-ke * tref[i]))), 0))
  
  y[i] ~ dbin(psurv[i] , Nprec[i]) 
}
}"

modelTKTDUnifm00 <- "model {
#########priors 
log10ks ~ dunif(minlog10ks, maxlog10ks)
log10NEC ~ dunif(minlog10conc, maxlog10conc)
log10ke ~ dunif(minlog10ke, maxlog10ke)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
ke <- 10**log10ke
eps <- 0.00000001

##########Computation of the likelihood
for (i in 1:ndat)
{
  psurv[i] <- exp(-ks * (ifelse((x[i] - NEC) >= 0, x[i] - NEC, 0) + eps) * (t[i] - tprec[i]) + (x[i] / ke) * (exp(-ke * t[i]) - exp(-ke * tprec[i])))
  y[i] ~ dbin(psurv[i], ifelse(Nprec[i] > 0, Nprec[i], 1))
}
}"

modelTKTDUnifm00v2 <- "model {
#########priors 
log10ks ~ dunif(minlog10ks, maxlog10ks)
log10NEC ~ dunif(minlog10conc, maxlog10conc)
log10ke ~ dunif(minlog10ke, maxlog10ke)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
ke <- 10**log10ke
eps <- 0.00000001

##########Computation of the likelihood
for (i in 1:ndat)
{
  tNEC[i] <- ifelse(x[i] > NEC, -1 / ke * log(1 - R[i]), bigtime)
  R[i] <- ifelse(x[i] > NEC, NEC/xcor[i], 0.1)
  xcor[i] <- ifelse(x[i] > 0, x[i], 10)
  tref[i] <- max(tprec[i], tNEC[i])
  
  psurv[i] <- exp(ifelse(t[i] > tNEC[i], -ks * ((x[i] - NEC) * (t[i] - tref[i]) + x[i]/ke * ( exp(-ke * t[i]) - exp(-ke * tref[i]))), eps))
  
  y[i] ~ dbin(psurv[i], ifelse(Nprec[i] > 0, Nprec[i], 1))
}
}"

modelTKTDUnifkeInf <- "model {
#########priors 
log10ks ~ dunif(minlog10ks, maxlog10ks)
log10NEC ~ dunif(minlog10conc, maxlog10conc)
log10m0 ~ dunif(minlog10m0, maxlog10m0)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
m0 <- 10**log10m0
eps <- 0.00000001

##########Computation of the likelihood
for (i in 1:ndat)
{
  psurv[i] <- exp(m0 * (tprec[i] - t[i]) + ks * ((ifelse((x[i] - NEC) >= 0, x[i] - NEC, 0) + eps) * (tprec[i] - t[i])))
  y[i] ~ dbin(psurv[i], ifelse(Nprec[i] > 0, Nprec[i], 1))
}
}"

# model m0 = 0 ke = Inf new
modelTKTDUnifm00keInf <- "model {
#########priors 
log10ks ~ dunif(minlog10ks, maxlog10ks)
log10NEC ~ dunif(minlog10conc, maxlog10conc)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
eps <- 0.00000001

##########Computation of the likelihood
for (i in 1:ndat)
{
  
  psurv[i] <- exp(ks * (ifelse((x[i] - NEC) >= 0, x[i] - NEC, 0) + eps) * (tprec[i] - t[i]))
  
  y[i] ~ dbin(psurv[i], ifelse(Nprec[i] > 0, Nprec[i], 1))
}
}"

modelTKTDNorm <- "model {
#########priors 
log10ks ~ dnorm(meanlog10ks, taulog10ks)
log10NEC ~ dnorm(meanlog10nec, taulog10nec)
log10ke ~ dnorm(meanlog10ke, taulog10ke)
log10m0 ~ dnorm(meanlog10m0, taulog10m0)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
ke <- 10**log10ke
m0 <- 10**log10m0

##########Computation of the likelihood
for (i in 1:ndat)
{
  tNEC[i] <- ifelse(x[i] > NEC, -1/ke * log( 1- R[i]), bigtime)
  R[i] <- ifelse(x[i] > NEC, NEC/xcor[i], 0.1)
  xcor[i] <- ifelse(x[i] > 0, x[i], 10)
  tref[i] <- max(tprec[i], tNEC[i])
  
  psurv[i] <- exp(-m0 * (t[i] - tprec[i]) + ifelse(t[i] > tNEC[i], -ks * ((x[i] - NEC) * (t[i] - tref[i]) + x[i]/ke * ( exp(-ke * t[i]) - exp(-ke * tref[i]))), 0))
  
  y[i] ~ dbin(psurv[i] , Nprec[i]) 
}
}"

modelTKTDNormm00 <- "model {
#########priors 
log10ks ~ dnorm(meanlog10ks, taulog10ks)
log10NEC ~ dnorm(meanlog10nec, taulog10nec)
log10ke ~ dnorm(meanlog10ke, taulog10ke)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
ke <- 10**log10ke
eps <- 0.00000001

##########Computation of the likelihood
for (i in 1:ndat)
{
  psurv[i] <- exp(-ks * (ifelse((x[i] - NEC) >= 0, x[i] - NEC, 0) + eps) * (t[i] - tprec[i]) + (x[i] / ke) * (exp(-ke * t[i]) - exp(-ke * tprec[i])))
  y[i] ~ dbin(psurv[i], ifelse(Nprec[i] > 0, Nprec[i], 1))
}
}"

modelTKTDNormm00v2 <- "model {
#########priors 
log10ks ~ dnorm(meanlog10ks, taulog10ks)
log10NEC ~ dnorm(meanlog10nec, taulog10nec)
log10ke ~ dnorm(meanlog10ke, taulog10ke)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
ke <- 10**log10ke
eps <- 0.00000001

##########Computation of the likelihood
for (i in 1:ndat)
{
  tNEC[i] <- ifelse(x[i] > NEC, -1/ke * log( 1- R[i]), bigtime)
  R[i] <- ifelse(x[i] > NEC, NEC/xcor[i], 0.1)
  xcor[i] <- ifelse(x[i] > 0, x[i], 10)
  tref[i] <- max(tprec[i], tNEC[i])
  
  psurv[i] <- exp(ifelse(t[i] > tNEC[i], -ks * ((x[i] - NEC) * (t[i] - tref[i]) + x[i]/ke * ( exp(-ke * t[i]) - exp(-ke * tref[i]))), eps))
  
  y[i] ~ dbin(psurv[i], ifelse(Nprec[i] > 0, Nprec[i], 1))
}
}"


# model ke = Inf
modelTKTDNormkeInf <- "model {
#########priors 
log10ks ~ dnorm(meanlog10ks, taulog10ks)
log10NEC ~ dnorm(meanlog10nec, taulog10nec)
log10m0 ~ dnorm(meanlog10m0, taulog10m0)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
m0 <- 10**log10m0
eps <- 0.00000001

##########Computation of the likelihood
for (i in 1:ndat)
{
  psurv[i] <- exp(m0 * (tprec[i] - t[i]) + ks * ((ifelse((x[i] - NEC) >= 0, x[i] - NEC, 0) + eps) * (tprec[i] - t[i])))
  y[i] ~ dbin(psurv[i], ifelse(Nprec[i] > 0, Nprec[i], 1))
}
}"

# model m0 = 0 ke = Inf new
modelTKTDNormm00keInf <- "model {
#########priors 
log10ks ~ dnorm(meanlog10ks, taulog10ks)
log10NEC ~ dnorm(meanlog10nec, taulog10nec)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
eps <- 0.00000001

##########Computation of the likelihood
for (i in 1:ndat)
{
  
  psurv[i] <- exp(ks * (ifelse((x[i] - NEC) >= 0, x[i] - NEC, 0) + eps) * (tprec[i] - t[i]))
  
  y[i] ~ dbin(psurv[i], ifelse(Nprec[i] > 0, Nprec[i], 1))
}
}"

survTKTDPARAMS <- function(mcmc, m0, ke) {
  # create the table of posterior estimated parameters
  # for the survival analyses
  # INPUT:
  # - mcmc:  list of estimated parameters for the model with each item representing
  # a chain
  # OUTPUT:
  # - data frame with 3 columns (values, CIinf, CIsup) and 3-4rows (the estimated
  # parameters)
  
  # Retrieving parameters of the model
  res.M <- summary(mcmc)
  
  if (ke) {
    ke <- 10^res.M$quantiles["log10ke", "50%"]
    keinf <- 10^res.M$quantiles["log10ke", "2.5%"]
    kesup <- 10^res.M$quantiles["log10ke", "97.5%"]
  }
  ks <- 10^res.M$quantiles["log10ks", "50%"]
  ksinf <- 10^res.M$quantiles["log10ks", "2.5%"]
  kssup <- 10^res.M$quantiles["log10ks", "97.5%"]
  nec <- 10^res.M$quantiles["log10NEC", "50%"]
  necinf <- 10^res.M$quantiles["log10NEC", "2.5%"]
  necsup <- 10^res.M$quantiles["log10NEC", "97.5%"]
  if (m0) {
    m0 <- 10^res.M$quantiles["log10m0", "50%"]
    m0inf <- 10^res.M$quantiles["log10m0", "2.5%"]
    m0sup <- 10^res.M$quantiles["log10m0", "97.5%"]
  }
  
  # Definition of the parameter storage and storage data
  if (ke && m0) {
    rownames <- c("ke", "ks", "nec", "m0")
    params <- c(ke, ks, nec, m0)
    CIinf <- c(keinf, ksinf, necinf, m0inf)
    CIsup <- c(kesup, kssup, necsup, m0sup)
  } else if (!ke && m0) {
    rownames <- c("ks", "nec", "m0")
    params <- c(ks, nec, m0)
    CIinf <- c(ksinf, necinf, m0inf)
    CIsup <- c(kssup, necsup, m0sup)
  } else if (ke && !m0) {
    rownames <- c("ke", "ks", "nec")
    params <- c(ke, ks, nec)
    CIinf <- c(keinf, ksinf, necinf)
    CIsup <- c(kesup, kssup, necsup)
  } else {
    rownames <- c("ks", "nec")
    params <- c(ks, nec)
    CIinf <- c(ksinf, necinf)
    CIsup <- c(kssup, necsup)
  }
  
  res <- data.frame(median = params, Q2.5 = CIinf, Q97.5 = CIsup,
                    row.names = rownames)
  
  return(res)
}

#' Fit a Bayesian TKTD model for survival analysis among time
#' 
#' The \code{survFitTKTD} function estimates the parameters of a TKTD model
#' for survival analysis using Bayesian inference.
#' 
#' The function returns
#' parameter estimates of the TKTD model.
#' 
#' @param data An object of class \code{survData}.
#' @param distr If \code{"norm"} gives normal distributions for prior distribution
#' parameters, else if \code{"unif"} gives uniform distributions for prior
#' distribution parameters.
#' @param bond If \code{"01"} the extreme case where the survival remains at 99\%
#' at the highest tested concentration, else if \code{"001"} it's 99.9\%.
#' @param m0 If \code{TRUE}, m0 is in the model.
#' @param ke If \code{TRUE}, ke is in the model.
#' @param n.chains Number of MCMC chains. The minimum required number of chains
#' is 2.
#' @param quiet If \code{TRUE}, make silent all prints and progress bars of
#' JAGS compilation.
#' 
#' @return The function returns an object of class \code{survFitTKTD}. A list
#' of 9 objects:
#' \item{estim.par}{A table of the estimated parameters as medians and 95 \%
#' credible intervals.}
#' \item{mcmc}{An object of class \code{mcmc.list} with the posterior
#' distributions.}
#' \item{model}{A JAGS model object.}
#' \item{parameters}{A list of the parameters names used in the model.}
#' \item{n.chains}{An integer value corresponding to the number of chains used
#' for the MCMC computation.}
#' \item{n.iter}{A list of two numerical value corresponding to the beginning
#' and the end of monitored iterations.}
#' \item{n.thin}{A numerical value corresponding to the thinning interval.}
#' \item{jags.data}{A list of data used by the internal \code{\link[rjags]{jags.model}}
#' function. This object is intended for the case when the user wishes to use
#' the \code{\link[rjags]{rjags}} package instead of the automatied estimation
#' function.}
#' \item{transformed.data}{The \code{survData} object.
#' See \code{\link{survData}} for details.}
#' 
#' @author Marie Laure Delignette-Muller
#' <marielaure.delignettemuller@@vetagro-sup.fr>, Philippe Ruiz
#' <philippe.ruiz@@univ-lyon1.fr>
#' 
#' @references Plummer, M. (2013) JAGS Version 3.4.0 user manual.
#' \url{http://sourceforge.net/projects/mcmc-jags/files/Manuals/3.x/jags_user_manual.pdf/download}
#'
#' Spiegelhalter, D., N. Best, B. Carlin, and A. van der Linde (2002) Bayesian
#' measures of model complexity and fit (with discussion).  \emph{Journal of
#' the Royal Statistical Society}, Series B 64, 583-639.
#'
#' @keywords estimation
#
#' @examples
#' 
#' # (1) Load the survival data
#' data(cadmium1)
#' 
#' # (2) Create an object of class "survData"
#' dat <- survData(cadmium1)
#' 
#' \dontrun{
#' # (3) Run the survFitTKTD function
#' out <- survFitTKTD(dat)
#' }
#' 
#' @export
#' @import rjags
#' @importFrom dplyr group_by summarise filter
#' 
survFitTKTD <- function(data,
                        distr = "norm",
                        bond = "01",
                        m0 = TRUE,
                        ke = TRUE,
                        n.chains = 3,
                        quiet = FALSE) {
  # test class object
  if(!is(data, "survData"))
    stop("survFitTKTD: object of class survData expected")
  
#   # Choose model by testing mortality in the control
#   if (!m0) {
#     control <- filter(data, conc == 0)
#     if (any(control$Nsurv < control$Ninit)) {
#       m0 <- TRUE
#       message("m0 is turned TRUE because there is mortality in the control !")
#     }
#   }
#   
  # data transformation
  data <- summarise(group_by(data, conc, time), N_alive = sum(Nsurv))

  n <- nrow(data)
  data$tprec <- NA
  data$t2prec <- NA
  data$Nprec <- NA
  data$N_init <- NA
  for (i in 1:n)
  {
    if (data$time[i] != 0)
    {
      data$tprec[i] <- data$time[i - 1]
      data$Nprec[i] <- data$N_alive[i - 1]
      data$N_init[i] <- data$N_alive[data$conc == data$conc[i] & data$time == 0]
    }
    if (data$time[i] != 0 && data$time[i] != min(data$time[data$time != 0])) {
      data$t2prec[i] <- data$tprec[i - 1]
    }
  }
  
  # control
  datasurv0 <- subset(data, time == min(data$time[data$time != 0]))
  datasurv0$time <- 0
  datasurv0$N_alive <- datasurv0$N_init
  data[is.na(data$tprec),
       c("tprec", "Nprec", "N_init")] <- datasurv0[, c("tprec", "Nprec", "N_init")]

  jags.data <- survTKTDCreateJagsData(data, distr, bond, m0, ke)
  jags.data <- jags.data[!sapply(jags.data, is.null)]

  # Define model
  if (distr == "norm") {
    if (m0 && ke) {
      model <- survLoadModel(model.program = modelTKTDNorm,
                             data = jags.data, n.chains,
                             Nadapt = 3000, quiet)
    } else if (!m0 && ke) {
      model <- survLoadModel(model.program = modelTKTDNormm00v2,
                             data = jags.data, n.chains,
                             Nadapt = 3000, quiet)
    } else if (m0 && !ke) {
      model <- survLoadModel(model.program = modelTKTDNormkeInf,
                             data = jags.data, n.chains,
                             Nadapt = 3000, quiet)
    } else {
      model <- survLoadModel(model.program = modelTKTDNormm00keInf,
                             data = jags.data, n.chains,
                             Nadapt = 3000, quiet)
    }
  } else if (distr == "unif") {
    if (m0 && ke) {
      model <- survLoadModel(model.program = modelTKTDUnif,
                             data = jags.data, n.chains,
                             Nadapt = 3000, quiet)
    } else if (!m0 && ke) {
      model <- survLoadModel(model.program = modelTKTDUnifm00v2,
                             data = jags.data, n.chains,
                             Nadapt = 3000, quiet)
    } else if (m0 && !ke) {
      model <- survLoadModel(model.program = modelTKTDUnifkeInf,
                             data = jags.data, n.chains,
                             Nadapt = 3000, quiet)
    } else {
      model <- survLoadModel(model.program = modelTKTDUnifm00keInf,
                             data = jags.data, n.chains,
                             Nadapt = 3000, quiet)
    }
  }

  # Determine sampling parameters
  parameters <- if (m0 && ke) {
    c("log10ke", "log10NEC","log10ks", "log10m0")
  } else if (!m0 && ke) {
    c("log10ke", "log10NEC","log10ks")
  } else if (m0 && !ke) {
    c("log10NEC","log10ks", "log10m0")
  } else {
    c("log10NEC","log10ks")
  }
  
  sampling.parameters <- modelSamplingParameters(model,
                                                 parameters, n.chains, quiet)
  
  tktd.dic <- calcDIC(model, sampling.parameters, quiet)
  
  # Sampling
  prog.b <- ifelse(quiet == TRUE, "none", "text")
  
  mcmc <- coda.samples(model, parameters,
                       n.iter = sampling.parameters$niter,
                       thin = sampling.parameters$thin,
                       progress.bar = prog.b)
  
  # summarize estime.par et CIs
  # calculate from the estimated parameters
  estim.par <- survTKTDPARAMS(mcmc, m0, ke)
  
  #OUTPUT
  OUT <- list(estim.par = estim.par,
              DIC = tktd.dic,
              mcmc = mcmc,
              model = model,
              ke = ke,
              m0 = m0,
              parameters = parameters,
              n.chains = summary(mcmc)$nchain,
              n.iter = list(start = summary(mcmc)$start,
                            end = summary(mcmc)$end),
              n.thin = summary(mcmc)$thin,
              distr = distr,
              jags.data = jags.data,
              transformed.data = data)
  
  class(OUT) <- "survFitTKTD"
  return(OUT)
}

