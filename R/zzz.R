GamaWelcomeMessage <- function() {

  msg <- c(paste0(
    "GAMA - Genetic Approach to Maximize clustering criterion, version ", packageVersion("gama")),
    "\nType 'citation(\"gama\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg) {
  # startup message
  msg <- GamaWelcomeMessage()
  if(!interactive())
    msg[1] <- paste("Package 'gama', version", packageVersion("gama"))
  packageStartupMessage(msg)
  invisible()
}
