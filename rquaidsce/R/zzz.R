# Package startup and python backend loading

.pyquaidsce_env <- new.env(parent = emptyenv())

.find_generic_python <- function() {
  # 1. User explicit environment variable
  env_py <- Sys.getenv("PYQUAIDSCE_PYTHON", Sys.getenv("RETICULATE_PYTHON", ""))
  if (nzchar(env_py) && file.exists(env_py)) {
    return(env_py)
  }

  # 2. System PATH
  w_py <- Sys.which("python")
  if (nzchar(w_py) && file.exists(w_py)) {
    return(w_py)
  }
  w_py3 <- Sys.which("python3")
  if (nzchar(w_py3) && file.exists(w_py3)) {
    return(w_py3)
  }

  # 3. Standard OS install directories (generic without user-specific hardcoded paths)
  if (.Platform$OS.type == "windows") {
    local_app <- Sys.getenv("LOCALAPPDATA")
    prog_files <- Sys.getenv("ProgramFiles")
    prog_files86 <- Sys.getenv("ProgramFiles(x86)")
    vers <- c("Python314", "Python313", "Python312", "Python311", "Python310", "Python39")

    cands <- c(
      file.path(local_app, "Programs", "Python", vers, "python.exe"),
      file.path(local_app, "Python", "pythoncore-3.14-64", "python.exe"),
      file.path(prog_files, vers, "python.exe"),
      file.path(prog_files86, vers, "python.exe")
    )
    cands <- cands[nzchar(cands) & file.exists(cands)]
    if (length(cands) > 0) {
      return(cands[1])
    }
  } else {
    cands <- c("/usr/bin/python3", "/usr/local/bin/python3", "/opt/homebrew/bin/python3")
    cands <- cands[file.exists(cands)]
    if (length(cands) > 0) {
      return(cands[1])
    }
  }

  return(NULL)
}

.onLoad <- function(libname, pkgname) {
  # Configure RETICULATE_PYTHON early before reticulate binds
  if (!nzchar(Sys.getenv("RETICULATE_PYTHON", ""))) {
    cand <- .find_generic_python()
    if (!is.null(cand)) {
      Sys.setenv(RETICULATE_PYTHON = cand)
    }
  }

  if (requireNamespace("reticulate", quietly = TRUE)) {
    reticulate::configure_environment(pkgname)
  }
}

.get_pyquaidsce <- function() {
  if (exists("module", envir = .pyquaidsce_env, inherits = FALSE)) {
    return(get("module", envir = .pyquaidsce_env))
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required to use 'rquaidsce'. Please install it with `install.packages('reticulate')`.", call. = FALSE)
  }

  # Import pyquaidsce
  py <- tryCatch({
    reticulate::import("pyquaidsce", delay_load = FALSE)
  }, error = function(e) {
    # If standard import failed, try with generic discovered python
    cand <- .find_generic_python()
    if (!is.null(cand)) {
      tryCatch({
        reticulate::use_python(cand, required = TRUE)
        reticulate::import("pyquaidsce", delay_load = FALSE)
      }, error = function(e2) NULL)
    } else {
      NULL
    }
  })

  if (is.null(py)) {
    stop("The 'pyquaidsce' Python module could not be loaded.\n",
         "Please ensure 'pyquaidsce' is installed in your active Python environment (`pip install pyquaidsce`).\n",
         "If your Python environment is in a custom path, specify it before loading the package:\n",
         "  Sys.setenv(RETICULATE_PYTHON = 'path/to/python')\n",
         "  library(rquaidsce)",
         call. = FALSE)
  }

  # Check minimum version requirement
  ver <- tryCatch(as.character(py$`__version__`), error = function(e) "0.0.0")
  if (utils::compareVersion(ver, "1.3.0") < 0) {
    warning(sprintf("Detected pyquaidsce version %s. Version >= 1.3.0 is recommended.", ver), call. = FALSE)
  }

  assign("module", py, envir = .pyquaidsce_env)
  return(py)
}
