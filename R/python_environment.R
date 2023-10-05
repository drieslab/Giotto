


#' @title checkGiottoEnvironment
#' @name checkGiottoEnvironment
#' @param mini_install_path (optional) path to miniconda or conda directory within which the giotto environment lives
#' If not provided, automatically determined by reticulate::miniconda_path()
#' Note the required input format:
#' Correct format --> mini_install_path = "C:/my/conda/lives/here" OR "C:\\my\\conda\\lives\\here"
#' INCORRECT formats --> mini_install_path = "C:/my/conda/lives/here/" AND "C:\\my\\conda\\lives\\here\\"
#' @param verbose be verbose
#' @details Checks if a miniconda giotto environment can be found.
#' Can be installed with \code{\link{installGiottoEnvironment}}.
#' @export
checkGiottoEnvironment =  function(mini_install_path = NULL, verbose = TRUE) {

  ## get operating system
  os_specific_system = get_os()

  ## check if giotto environment is already installed
  if (is.null(mini_install_path)){
    conda_path = reticulate::miniconda_path()
  } else if (!dir.exists(mini_install_path)) {
    if(verbose) wrap_msg(" Unable to find miniconda directory", mini_install_path,
                         "\nPlease ensure the directory exists and is provided as a string.")
    return(FALSE)
  } else {
    conda_path = mini_install_path
  }

  if(os_specific_system == 'osx') {
    full_path = paste0(conda_path, "/envs/giotto_env/bin/pythonw")
  } else if(os_specific_system == 'windows') {
    full_path = paste0(conda_path, "\\envs\\giotto_env\\python.exe")
  } else if(os_specific_system == 'linux') {
    full_path = paste0(conda_path, "/envs/giotto_env/bin/python")
  }

  if(file.exists(full_path)) {
    if(verbose) message('\n giotto environment found at \n',
                        full_path, '\n')
    return(TRUE)

  } else {
    if(verbose) wrap_msg('\n giotto environment was expected, but NOT found at \n', full_path, '\n')
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
                                               python_version = '3.10.2',
                                               mini_install_path = NULL,
                                               verbose = TRUE) {

  ## install Giotto environment
  if(verbose) message('\n |---- install giotto environment ----| \n')
  conda_path = mini_install_path

  ## 3. identify operating system and adjust the necessary packages
  os_specific_system = get_os()


  if(os_specific_system != 'osx') {
    packages_to_install = packages_to_install[!grepl(pattern = 'python.app', x = packages_to_install)]
  }

  py_lou = packages_to_install[grepl(pattern = 'python-louvain',x = packages_to_install)]

  pip_packages = c("smfishhmrf", "session-info", py_lou)

  # python-louvain must be installed with pip, not with conda-forge
  packages_to_install = packages_to_install[!grepl(pattern = 'python-louvain',x = packages_to_install)]

  ## for unix-like systems ##
  if(.Platform[['OS.type']] == 'unix') {

    conda_full_path = paste0(conda_path,'/','bin/conda')
    # If this does not exist, check for alternative conda config
    # i.e. an env created through .condarc config
    if (!file.exists(conda_full_path)) conda_full_path = paste0(conda_path,"/conda.exe")
    reticulate::conda_create(envname = 'giotto_env',
                             conda = conda_full_path,
                             python_version = python_version)

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
                           python_version = python_version,
                           channel = c('conda-forge', 'vtraag'))

    reticulate::py_install(packages = pip_packages,
                           envname = full_envname,
                           method = 'conda',
                           conda = conda_full_path,
                           pip = TRUE,
                           python_version = python_version)


    ## for windows systems ##
  } else if(.Platform[['OS.type']] == 'windows') {

    conda_full_path = paste0(conda_path,'/','condabin/conda.bat')
    reticulate::conda_create(envname = 'giotto_env',
                             conda = conda_full_path,
                             python_version = python_version)


    full_envname = paste0(conda_path,'/envs/giotto_env')
    python_full_path = paste0(conda_path, "/envs/giotto_env/python.exe")

    reticulate::py_install(packages = packages_to_install,
                           envname = 'giotto_env',
                           method = 'conda',
                           conda = conda_full_path,
                           python_version = python_version,
                           channel = c('conda-forge', 'vtraag'))

    reticulate::py_install(packages = pip_packages,
                           envname = full_envname,
                           method = 'conda',
                           conda = conda_full_path,
                           pip = TRUE,
                           python_version = python_version)
  }

}


#' @title install_giotto_environment
#' @description installation options of giotto environment
#' @keywords internal
install_giotto_environment = function(force_environment = FALSE,
                                      packages_to_install = c('pandas', 'networkx', 'python-igraph',
                                                              'leidenalg', 'python-louvain', 'python.app',
                                                              'scikit-learn'),
                                      python_version = '3.10.2',
                                      mini_install_path = NULL,
                                      verbose = TRUE) {

  # first see if Giotto is already installed
  giotto_installed = checkGiottoEnvironment(mini_install_path = mini_install_path,
                                            verbose = verbose)


  if(giotto_installed == TRUE & force_environment == FALSE) {
    # do nothing if already installed and no force required

    if(verbose) wrap_msg('Giotto environment is already installed, set force_environment = TRUE to reinstall \n')

  } else if(giotto_installed == TRUE & force_environment == TRUE) {
    # reinstall giotto if force required
    if(.Platform[['OS.type']] == 'unix') {
      conda_full_path = paste0(mini_install_path,'/','bin/conda')
      # If this does not exist, check for alternative conda config
      # i.e. an env created through .condarc config
      if (!file.exists(conda_full_path)) conda_full_path = paste0(mini_install_path,"/conda.exe")
    } else if(.Platform[['OS.type']] == 'windows') {
      conda_full_path = paste0(mini_install_path,'/','condabin/conda.bat')
    }

    # first remove giotto environment, then install
    reticulate::conda_remove(envname = 'giotto_env', conda = conda_full_path)

    install_giotto_environment_specific(packages_to_install = packages_to_install,
                                        python_version = python_version,
                                        mini_install_path = mini_install_path,
                                        verbose = verbose)

  } else {

    # install giotto if nothing is found
    install_giotto_environment_specific(packages_to_install = packages_to_install,
                                        python_version = python_version,
                                        mini_install_path = mini_install_path,
                                        verbose = verbose)

  }

}




#' @title installGiottoEnvironment
#' @description Installs a giotto environment
#' @param packages_to_install all python modules (packages) that should be installed for Giotto to work
#' @param python_version python version to use within the giotto conda environment
#' @param mini_install_path (optional) desired location miniconda installation.
#' If not provided, it is chosen by reticulate::install_miniconda()
#' Note the required input format:
#' Correct format --> mini_install_path = "C:/my/conda/lives/here" OR "C:\\my\\conda\\lives\\here"
#' INCORRECT formats --> mini_install_path = "C:/my/conda/lives/here/" AND "C:\\my\\conda\\lives\\here\\"
#' @param force_miniconda force reinstallation of miniconda
#' @param force_environment force reinstallation of the giotto environment
#' @param verbose be verbose
#' @return installs a giotto environment using the reticulate miniconda system
#' @details This function will install a local giotto environment using
#' the miniconda system as implemented by reticulate. Once this giotto environment is
#' installed it will be automatically detected when you run the Giotto toolbox. If you want to use
#' your own python path then you can set the python_path in the \code{\link{createGiottoInstructions}}
#' and provide the instructions to the \code{\link{createGiottoObject}} function.
#'
#' By default, Python v3.10.2 will be used with the following python modules
#' for Giotto Suite implementations:
#' \preformatted{
#'    - pandas==1.5.1
#'    - networkx==2.8.8
#'    - python-igraph==0.10.2
#'    - leidenalg==0.9.0
#'    - python-louvain==0.16
#'    - python.app==1.4
#'    - scikit-learn==1.1.3
#' }
#'
#'  The giotto environment can be custom installed by changing the
#'  python_version parameter and module versions in the
#'  packages_to_install parameter.
#'
#'  For example, this giotto environment works as well, and was the
#'  default environment status for past releases of Giotto.
#'  Python  v3.6
#'  \preformatted{
#'   - pandas==1.1.5
#'   - networkx==2.6.3
#'   - python-igraph==0.9.6
#'   - leidenalg==0.8.7
#'   - python-louvain==0.15
#'   - python.app==2 # macOS only
#'   - scikit-learn==0.24.2
#' }
#'
#' @export
#' @examples
#' \dontrun{
#'
#'   # this command will install r-miniconda
#'   # and a giotto environment with all necessary python modules
#'   installGiottoEnvironment()
#'}
installGiottoEnvironment =  function(packages_to_install = c('pandas==1.5.1',
                                                             'networkx==2.8.8',
                                                             'python-igraph==0.10.2',
                                                             'leidenalg==0.9.0',
                                                             'python-louvain==0.16',
                                                             'python.app==1.4',
                                                             'scikit-learn==1.1.3'),
                                     python_version = '3.10.2',
                                     mini_install_path = NULL,
                                     force_miniconda = FALSE,
                                     force_environment = FALSE,
                                     verbose = TRUE) {


  ## 1. check and install miniconda locally if necessary
  conda_path = NULL
  if (is.null(mini_install_path)){
    conda_path = reticulate::miniconda_path()
  } else if (!dir.exists(mini_install_path)) {
    stop(wrap_msg(paste0(" Unable to install miniconda in ", mini_install_path, "\nPlease ensure the directory has been created and provided as a string.")))
  } else {
    conda_path = mini_install_path

    wrap_msg("NOTICE: Attempting to install the Giotto Environment at a custom path.\n",
    "Please note that multiple .yml files are provided in the repository for advanced installation and convenience.",
    "To install the most up-to-date Giotto environment using a .yml file, open a shell",
    " compatible with conda/miniconda and navigate to the directory containing Giotto.",
    "If you are unsure where Giotto lives on your machine, run the R function",
    " `.libPaths()`, which will return the path(s) at which R packages install on your machine.",
    "Once in the directory containing Giotto, run the following to create your environment in one step:\n\n",
    "conda env create -n giotto_env -f ./python/configuration/genv.yml\n\n",
    "Alternatively, Giotto environment configurations are stored in",
    " the directory Giotto/inst/python/configuration/ on the github repository.")

    manual_install = as.character(readline("Would you prefer to install manually? [y/n] "))

    if(!manual_install %in% c("y","Y","n","N")) stop("Invalid input. Please try again.")

    if (manual_install %in% c("y","Y")) stop(wrap_txt("There is no error; this just stops function execution. Please follow the instructions above for manual installation. Thank you!"))
    else wrap_msg("Continuing with automatic installation...\n")


    if(.Platform[['OS.type']] == 'unix') {
      conda_full_path = paste0(conda_path,'/','bin/conda')
      # If this does not exist, check for alternative conda config
      # i.e. an env created through .condarc config
      if (!file.exists(conda_full_path)) conda_full_path = paste0(conda_path,"/conda.exe")
    } else if(.Platform[['OS.type']] == 'windows') {
      conda_full_path = paste0(conda_path,'/','condabin/conda.bat')
    }
    if (!file.exists(conda_full_path)) force_miniconda = TRUE
  }

  if(!file.exists(conda_path) | isTRUE(force_miniconda)) {
    if(verbose) message('\n |---- install local miniconda ----| \n')
    reticulate::install_miniconda(path = conda_path, force = force_miniconda)
  }


  ## 2. install giotto environment
  install_giotto_environment(force_environment = force_environment,
                             packages_to_install = packages_to_install,
                             python_version = python_version,
                             mini_install_path = conda_path,
                             verbose = verbose)


}



#' @title removeGiottoEnvironment
#' @name removeGiottoEnvironment
#' @param mini_path path to anaconda/miniconda.
#' Default: reticulate::miniconda_path()
#' i.e. "C:/my/conda/lives/here" OR "C:\\my\\conda\\lives\\here"
#' @param verbose be verbose
#' @details Removes a previously installed giotto environment.
#' See \code{\link{installGiottoEnvironment}}.
#' @export
removeGiottoEnvironment = function(mini_path = NULL, verbose = TRUE) {

  if (is.null(mini_path)){
    conda_path = reticulate::miniconda_path()
  } else if (!dir.exists(mini_path)) {
     stop(wrap_msg(paste0("Directory", mini_path, "could not be found. Please ensure it was entered properly and as a string.")))
  } else {
     conda_path = mini_path
  }

  if(.Platform[['OS.type']] == 'unix') {
    conda_full_path = paste0(conda_path,'/','bin/conda')

  } else if(.Platform[['OS.type']] == 'windows') {
    conda_full_path = paste0(conda_path,'/','condabin/conda.bat')
  }

  # first see if Giotto is already installed
  giotto_installed = checkGiottoEnvironment(verbose = verbose)

  if(giotto_installed == FALSE) {
    wrap_msg('Giotto environment is not found and probably never installed')
  }

  # first remove giotto environment, then install
  reticulate::conda_remove(envname = 'giotto_env',
                           conda = conda_full_path)

}


#' @title set_giotto_python_path
#' @name set_giotto_python_path
#' @description sets the python path
#' @param python_path character. Full path to python executable
#' @param verbose be verbose
#' @keywords internal
set_giotto_python_path = function(python_path = NULL,
                                  verbose = TRUE) {

  # If a path is provided by the user and it exists,
  # direct reticulate to said executable and exit immediately
  if(!is.null(python_path) && file.exists(python_path)) {
    if(verbose) wrap_msg('\n external python path provided and will be used \n')
    python_path = as.character(python_path)
    reticulate::use_python(required = T, python = python_path)
    return (python_path)
  }

  # Otherwise, check the OS and if a Giotto Environment exists
  # use that executable
  os_specific_system = get_os()
  conda_path = reticulate::miniconda_path()
  if(os_specific_system == 'osx') {
    python_path = paste0(conda_path, "/envs/giotto_env/bin/pythonw")
  } else if(os_specific_system == 'windows') {
    python_path = paste0(conda_path, "\\envs\\giotto_env\\python.exe")
  } else if(os_specific_system == 'linux') {
    python_path = paste0(conda_path, "/envs/giotto_env/bin/python")
  }


  # check if giotto environment exists
  giotto_environment_installed = checkGiottoEnvironment(mini_install_path = conda_path,
                                                        verbose = FALSE)

  if(giotto_environment_installed == TRUE) {

    if(verbose) wrap_msg('\n no external python path was provided, but a giotto python environment was found and will be used \n')
    #python_path = return_giotto_environment_path_executable()
    reticulate::use_python(required = T, python = python_path)

  } else {

    if(verbose) wrap_msg('\n no external python path or giotto environment was specified, will check if a default python path is available \n')

    if(.Platform[['OS.type']] == 'unix') {
      python_path = try(system('which python3', intern = T))
    } else if(.Platform[['OS.type']] == 'windows') {
      python_path = try(system('where python3', intern = T))
    }

    if(inherits(python_path, 'try-error')) {
      wrap_msg('\n no default python path found, install python and/or use strategy 1 or 2 \n')
      python_path = NULL
    } else {
      python_path = python_path
      reticulate::use_python(required = T, python = python_path)
      if(verbose) wrap_msg('\n A default python path was found: ', python_path, ' and will be used\n')
    }

    wrap_msg('\n If this is not the correct python path, either')
    wrap_msg('\n 1. use installGiottoEnvironment() to install a local miniconda python environment along with required modules')
    wrap_msg('\n 2. provide an existing python path to python_path to use your own python path which has all modules installed')

  }

  return(python_path)
}

#' @title Prompt User for Python Install
#' @name py_install_prompt
#' @param package python package/github url
#' @param env environment into which package will be installed
#' @description prompts user to install a package
#' @keywords internal
py_install_prompt <- function (package = NULL,
                               env = NULL){
  if(is.null(package) | is.null(env)) {
    stop(GiottoUtils::wrap_txt("Incorrect Usage.\n", errWidth = TRUE))
  }

  install_py_pkg_msg = paste0("Python package `",
                              package,
                              "` is required and not installed.\n")
  warning(install_py_pkg_msg, immediate. = TRUE)
  install_py_pkg_msg = paste0("Enter 0 to skip installation and quit.\n")
  install_py_pkg_msg = paste0(install_py_pkg_msg, "Enter any other number ",
                              "to install\n`",
                              package,
                              "`\nto python environment: `",
                              env, '`\n\n')
  resp = as.integer(readline(prompt = install_py_pkg_msg))
  return(resp)

}

#' @title Install Package from GitHub Link
#' @name install_github_link_pip
#' @param link link to github repository containing a python package, 
#' e.g. `git+https://github.com/TencentAILabHealthcare/pysodb.git`
#' @param env conda environment to which `link` will be installed via pip
#' @description
#' Installs `link` to python `env`
#' @keywords internal
install_github_link_pip <- function(link = NULL,
                                            env = NULL){
  # Guard 
  if (is.null(link) | is.null(env)) stop(GiottoUtils::wrap_txt("Incorrect Usage.", errWidth = TRUE))

  config <- reticulate::py_discover_config(use_environment=env)
  # system commands return 0 if ran successfully, 1 otherwise
  successful_install = !as.logical(system2(config$python, 
                                          c("-m", 
                                            "pip",
                                            "install",
                                            link)))
  if (successful_install) {
    return(TRUE)
  } else {
    git_url_err_msg = "Provided GitHub URL `"
    git_url_err_msg = paste0(git_url_err_msg, github_package_URL, "`")
    git_url_err_msg = paste0(git_url_err_msg, " Could not be installed.\n")
    git_url_err_msg = paste0(git_url_err_msg, 
                            "Please try again with a different URL.")
    stop(GiottoUtils::wrap_txt(git_url_err_msg, errWidth = TRUE))
  }
}

#' @title Install Python Package with Reticulate
#' @name install_py_pkg_reticulate
#' @param package name of python package
#' @param env name of the environment into which the python 
#' package should be installed.
#' @details 
#' Installs `package` to python `env` after prompting user.
#' Installation is done via `py_install` from the 
#' `reticulate` package.
#' @keywords internal
install_py_pkg_reticulate <- function(package = NULL,
                                      env = NULL){
  resp = py_install_prompt(package = package,
                           env = env)
  if(resp != 0){
    try_install = tryCatch(expr = { 
      reticulate::py_install(package = package, 
                             envname = env)
      return(TRUE)

    }, error = function(e){
      return(FALSE)
    })

    if (!try_install){
      cannot_install_msg = paste0("Could not install `",package)
      cannot_install_msg = paste0(cannot_install_msg, "` using `reticulate::py_install()`\n")

      GiottoUtils::wrap_msg(cannot_install_msg,
                            errWidth = TRUE)
      return(FALSE)
    }
  } else {
    stop(GiottoUtils::wrap_txt("Package not installed.\n", errWidth = TRUE))
  }
}

#' @title Check Python Package Installation
#' @name checkPythonPackage
#' @param package_name name of python package. See details.
#' @param github_package_url URL linking to github repository containing
#' a python package that may be installed with pip,
#' e.g. `git+https://github.com/TencentAILabHealthcare/pysodb.git`;
#' see details.
#' @param env_to_use name of the environment into which the python 
#' package should be installed.
#' @description checks python environment for a 
#' provided package, installs if it is not found.
#' @details
#' Parameter `github_package_url` takes precedent over 
#' `package_name`, i.e. if both are provided, only the github
#' URL will be installed. This function should only be provided
#' one parameter, or the other. 
#' 
#' @keywords export
checkPythonPackage <- function (package_name = NULL,
                              github_package_url = NULL,
                              env_to_use = "giotto_env"){
  # Guard clauses
  if (is.null(package_name) & is.null(github_package_url)){
    null_input_err_msg = "A python package name must be provided, e.g. `scanpy==1.9.0`"
    null_input_err_msg = paste0(null_input_err_msg, 
                                "\nAlternatively, provide a github package URL, ")
    null_input_err_msg = paste0(null_input_err_msg, 
                                "e.g. `git+https://github.com/TencentAILabHealthcare/pysodb.git` ")
    stop(GiottoUtils::wrap_txt(null_input_err_msg, 
                                errWidth = TRUE))
  }
  # Find path to currently initialized python env
  path_to_env = reticulate::py_config()$pythonhome

  if (!grepl(env_to_use, path_to_env)){
    env_err_msg = paste0("Provided python environment `", 
                         env_to_use,"` is not initialized.")
    env_err_msg = paste0(env_err_msg, 
                         "\nThe following python environment is in use: `", 
                         path_to_env,"`")
    env_err_msg = paste0(env_err_msg, 
                         "\nTo initialize `",env_to_use,"`, you must ",
                         "restart your R session.")

    stop(GiottoUtils::wrap_txt(env_err_msg,
                               errWidth = TRUE))
  }

  env_str_location = stringr::str_locate(path_to_env, env_to_use)[2]
  # Change env_to_use from name of environment 
  # to the full environment path
  env_to_use = substr(path_to_env, 1, env_str_location)

  # If a github link is provided, install it and exit
  if (!is.null(github_package_url)){
    resp = py_install_prompt(package = github_package_url,
                             env = env_to_use)
    if (resp != 0){
      install_status = install_github_link_pip(link = github_package_url,
                                               env = env_to_use)
      return(install_status) 
    } else {
       stop(GiottoUtils::wrap_txt("Package not installed.\n", errWidth = TRUE))
    }
  }

  package_config = reticulate::py_list_packages()
  pkgs_in_py_env = package_config$package
  versions = package_config$version
  
  # package installed, right version --> exit
  # package installed, but wrong version --> prompt for install
  # package not installed --> prompt for install
  
  version_number = NULL

  contains_version_number = grepl("==", package_name) 
  if (contains_version_number){
    split_package_version = stringr::str_split(package_name, 
                                               pattern = "==")
    package_name = split_package_version[[1]][1]
    version_number = split_package_version[[1]][2]
  }

  if (package_name %in% pkgs_in_py_env){
    if (!contains_version_number){
      # If a version number is not provided,
      # and the package exists within the 
      # reticulate package list, exit, 
      # since it is already installed.
      return(TRUE)

    } else {
      # Check that the version numbers match, if provided
      idx = which(pkgs_in_py_env == package_name)
      version_match = (version_number == versions[idx])

      if (version_match){
        # if the versions match, the right version
        # is installed
        return(TRUE)
      } else{ 
        # Otherwise, install the provided version
        inst_result = install_py_pkg_reticulate(package = paste0(package_name,
                                                                 "==",
                                                                 version_number),
                                                env = env_to_use)
        return(inst_result)
      }
    }
  } else {
    if (!contains_version_number) {
      # If it is not installed, and has no version
      # number, install it.
      inst_result = install_py_pkg_reticulate(package = package_name,
                                              env = env_to_use)
    } else{
      # If it is not installed, and has a version
      # number, concatenate the package and verion
      # strings, and install
      inst_result = install_py_pkg_reticulate(package = paste0(package_name,
                                                                 "==",
                                                                 version_number),
                                                env = env_to_use)
    }
    return(inst_result)
  }
}
