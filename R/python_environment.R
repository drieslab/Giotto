


#' @title checkGiottoEnvironment
#' @name checkGiottoEnvironment
#' @param verbose be verbose
#' @details Checks if a miniconda giotto environment can be found.
#' Can be installed with \code{\link{installGiottoEnvironment}}.
#' @export
checkGiottoEnvironment =  function(verbose = TRUE) {

  ## get operating system
  os_specific_system = get_os()

  ## check if giotto environment is already installed
  conda_path = reticulate::miniconda_path()
  if(os_specific_system == 'osx') {
    full_path = paste0(conda_path, "/envs/giotto_env/bin/pythonw")
  } else if(os_specific_system == 'windows') {
    full_path = paste0(conda_path, "\\envs\\giotto_env\\python.exe")
  } else if(os_specific_system == 'linux') {
    full_path = paste0(conda_path, "/envs/giotto_env/bin/python")
  }

  if(file.exists(full_path)) {
    if(verbose) cat('\n giotto environment found at \n',
                    full_path, '\n')
    return(TRUE)

  } else {
    return(FALSE)
  }

}

#' @title return_giotto_environment_path
#' @description returns the path to the detected giotto environment
#' @return string path
#' @keywords internal
return_giotto_environment_path = function() {

  ## get operating system
  os_specific_system = get_os()

  ## check if giotto environment is already installed
  conda_path = reticulate::miniconda_path()
  if(os_specific_system == 'osx') {
    full_path = paste0(conda_path, "/envs/giotto_env")
  } else if(os_specific_system == 'windows') {
    full_path = paste0(conda_path, "\\envs\\giotto_env")
  } else if(os_specific_system == 'linux') {
    full_path = paste0(conda_path, "/envs/giotto_env")
  }

  if(file.exists(full_path)) {
    return(full_path)
  } else {
    return(NA)
  }

}

#' @title return_giotto_environment_path_executable
#' @description returns the path to the detected giotto environment executable
#' @return string path
#' @keywords internal
return_giotto_environment_path_executable = function() {

  ## get operating system
  os_specific_system = get_os()

  ## check if giotto environment is already installed
  conda_path = reticulate::miniconda_path()
  if(os_specific_system == 'osx') {
    full_path = paste0(conda_path, "/envs/giotto_env/bin/pythonw")
  } else if(os_specific_system == 'windows') {
    full_path = paste0(conda_path, "\\envs\\giotto_env\\python.exe")
  } else if(os_specific_system == 'linux') {
    full_path = paste0(conda_path, "/envs/giotto_env/bin/python")
  }

  if(file.exists(full_path)) {
    return(full_path)
  } else {
    return(NA)
  }

}


#' @title install_giotto_environment_specific
#' @description installation of giotto environment
#' @keywords internal
install_giotto_environment_specific = function(packages_to_install = c('pandas', 'networkx', 'python-igraph',
                                                                       'leidenalg', 'python-louvain', 'python.app',
                                                                       'scikit-learn'),
                                               verbose = TRUE) {

  ## install Giotto environment
  if(verbose) cat('\n |---- install giotto environment ----| \n')
  conda_path = reticulate::miniconda_path()

  ## 3. identify operating system and adjust the necessary packages
  os_specific_system = get_os()

  if(os_specific_system != 'osx') {
    packages_to_install = packages_to_install[packages_to_install != 'python.app']
  }


  ## for unix-like systems ##
  if(.Platform[['OS.type']] == 'unix') {

    conda_full_path = paste0(conda_path,'/','bin/conda')
    reticulate::conda_create(envname = 'giotto_env',
                             conda = conda_full_path)


    full_envname = paste0(conda_path,'/envs/giotto_env')

    if(os_specific_system == 'osx') {
      python_full_path = paste0(conda_path, "/envs/giotto_env/bin/pythonw")
    } else if(os_specific_system == 'linux') {
      python_full_path = paste0(conda_path, "/envs/giotto_env/bin/python")
    }


    reticulate::py_install(packages = packages_to_install,
                           envname = 'giotto_env',
                           method = 'conda',
                           conda = conda_full_path,
                           python_version = '3.6')

    reticulate::py_install(packages = 'smfishhmrf',
                           envname = full_envname,
                           method = 'conda',
                           conda = conda_full_path,
                           pip = TRUE,
                           python_version = '3.6')


    ## for windows systems ##
  } else if(.Platform[['OS.type']] == 'windows') {

    conda_full_path = paste0(conda_path,'/','condabin/conda.bat')
    reticulate::conda_create(envname = 'giotto_env',
                             conda = conda_full_path)


    full_envname = paste0(conda_path,'/envs/giotto_env')
    python_full_path = paste0(conda_path, "/envs/giotto_env/python.exe")

    reticulate::py_install(packages = packages_to_install,
                           envname = 'giotto_env',
                           method = 'conda',
                           conda = conda_full_path,
                           python_version = '3.6',
                           channel = c('conda-forge', 'vtraag'))

    reticulate::py_install(packages = 'smfishhmrf',
                           envname = full_envname,
                           method = 'conda',
                           conda = conda_full_path,
                           pip = TRUE,
                           python_version = '3.6')
  }


}


#' @title install_giotto_environment
#' @description installation options of giotto environment
#' @keywords internal
install_giotto_environment = function(force_environment = FALSE,
                                      packages_to_install = c('pandas', 'networkx', 'python-igraph',
                                                              'leidenalg', 'python-louvain', 'python.app',
                                                              'scikit-learn'),
                                      verbose = TRUE) {

  # first see if Giotto is already installed
  giotto_installed = checkGiottoEnvironment(verbose = verbose)


  if(giotto_installed == TRUE & force_environment == FALSE) {
    # do nothing if already installed and no force required

    if(verbose) cat('Giotto environment is already installed, set force_environment = TRUE to reinstall \n')

  } else if(giotto_installed == TRUE & force_environment == TRUE) {
    # reinstall giotto if force required

    # first remove giotto environment, then install
    reticulate::conda_remove(envname = 'giotto_env')

    install_giotto_environment_specific(packages_to_install = packages_to_install,
                                        verbose = verbose)

  } else {

    # install giotto if nothing is found
    install_giotto_environment_specific(packages_to_install = packages_to_install,
                                        verbose = verbose)

  }

}




#' @title installGiottoEnvironment
#' @description Installs a giotto environment
#' @param packages_to_install all python modules (packages) that should be installed for Giotto to work
#' @param force_miniconda force reinstallation of miniconda
#' @param force_environment force reinstallation of the giotto environment
#' @param verbose be verbose
#' @return installs a giotto environment using the reticulate miniconda system
#' @details This function will install a local giotto environment using
#' the miniconda system as implemented by reticulate. Once this giotto environment is
#' installed it will be automatically detected when you run the Giotto toolbox. If you want to use
#' your own python path then you can set the python_path in the \code{\link{createGiottoInstructions}}
#' and provide the instructions to the \code{\link{createGiottoObject}} function.
#' @export
#' @examples
#' \dontrun{
#'
#'   # this command will install r-miniconda
#'   # and a giotto environment with all necessary python modules
#'   installGiottoEnvironment()
#' }
#'
installGiottoEnvironment =  function(packages_to_install = c('pandas', 'networkx', 'python-igraph',
                                                             'leidenalg', 'python-louvain', 'python.app',
                                                             'scikit-learn'),
                                     force_miniconda = FALSE,
                                     force_environment = FALSE,
                                     verbose = TRUE) {


  ## 1. check and install miniconda locally if necessary
  conda_path = reticulate::miniconda_path()

  if(!file.exists(conda_path) | isTRUE(force_miniconda)) {
    if(verbose) cat('\n |---- install local miniconda ----| \n')
    reticulate::install_miniconda(force = force_miniconda)
  }


  ## 2. install giotto environment
  install_giotto_environment(force_environment = force_environment,
                             packages_to_install = packages_to_install,
                             verbose = verbose)


}



#' @title removeGiottoEnvironment
#' @name removeGiottoEnvironment
#' @param verbose be verbose
#' @details Removes a previously installed giotto environment.
#' See \code{\link{installGiottoEnvironment}}.
#' @export
removeGiottoEnvironment = function(verbose = TRUE) {

  # first see if Giotto is already installed
  giotto_installed = checkGiottoEnvironment(verbose = verbose)

  if(giotto_installed == FALSE) {
    cat('Giotto environment is not found and probably never installed')
  }

  # first remove giotto environment, then install
  reticulate::conda_remove(envname = 'giotto_env')

}


#' @title set_giotto_python_path
#' @name set_giotto_python_path
#' @description sets the python path
#' @keywords internal
set_giotto_python_path = function(python_path = NULL) {


  # check if giotto environment exists
  giotto_environment_installed = checkGiottoEnvironment(verbose = FALSE)

  ## if a python path is provided, use that path
  if(!is.null(python_path)) {
    cat('\n external python path provided and will be used \n')
    python_path = as.character(python_path)
    reticulate::use_python(required = T, python = python_path)

  } else if(giotto_environment_installed == TRUE) {

    cat('\n no external python path was provided, but a giotto python environment was found and will be used \n')
    python_path = return_giotto_environment_path_executable()
    reticulate::use_python(required = T, python = python_path)

  } else {

    cat('\n no external python path or giotto environment was specified, will check if a default python path is available \n')

    if(.Platform[['OS.type']] == 'unix') {
      python_path = try(system('which python3', intern = T))
    } else if(.Platform[['OS.type']] == 'windows') {
      python_path = try(system('where python3', intern = T))
    }

    if(class(python_path) == "try-error") {
      cat('\n no default python path found, install python and/or use strategy 1 or 2 \n')
      python_path = '/need/to/set/path/to/python'
    } else {
      python_path = python_path
      reticulate::use_python(required = T, python = python_path)
      cat('\n A default python path was found: ', python_path, ' and will be used\n')
    }

    cat('\n If this is not the correct python path, either\n')
    cat('\n 1. use installGiottoEnvironment() to install a local miniconda python environment along with required modules \n')
    cat('\n 2. provide an existing python path to python_path to use your own python path which has all modules installed \n')

  }

  return(python_path)
}

