get.os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
    else if (os == 'linux')
      os <- "linux"
    else if (os == 'windows')
      os <- "windows"
  } else {

    ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
    if (grepl("windows", R.version$os))
      os <- "windows"
  }

  tolower(os)
}

