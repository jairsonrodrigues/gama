GAStartupMessage <- function()
{
  # Startup message obtained as
  # > figlet GA
  msg <- c(paste0(
    "GAMA - Genetic Approach to MAximize clustering criterion, version ", packageVersion("gama")),
    "\nType 'citation(\"gama\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # unlock .ga.default variable allowing its modification
  #unlockBinding(".gama.default", asNamespace("gama"))
  # startup message
  msg <- GAStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'gama' version", packageVersion("GA"))
  packageStartupMessage(msg)
  invisible()
}
