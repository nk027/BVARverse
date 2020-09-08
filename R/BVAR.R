
#
# Internal BVAR functions that we use here. Basically the contents of
# BVAR's 11..., 12... and 15... scripts. At some point we should include
# these directly from the main repository.
#


#' Check numeric scalar
#'
#' Check whether an object is bounded and coercible to a numeric value.
#'
#' @param x Numeric scalar.
#' @param min Numeric scalar. Minimum value of \emph{x}.
#' @param max Numeric scalar. Maximum value of \emph{x}.
#' @param fun Function to apply to \emph{x} before returning.
#' @param msg String fed to \code{\link[base]{stop}} if an error occurs.
#'
#' @return Returns \code{fun(x)}.
#'
#' @noRd
num_check <- function(
  x, min = 0, max = Inf,
  msg = "Please check the numeric parameters.",
  fun = as.numeric) {

  if(!is.numeric(x) || length(x) != 1 || x < min || x > max) {stop(msg)}

  return(fun(x))
}


#' @noRd
int_check <- function(
  x, min = 0L, max = Inf,
  msg = "Please check the integer parameters.") {

  num_check(x, min, max, msg, fun = as.integer)
}

#' Lag a matrix
#'
#' Compute a lagged version of a matrix to be used in vector autoregressions.
#' Higher lags are further to the right.
#'
#' @param x Matrix (\eqn{N * M}) to lag.
#' @param lags Integer scalar. Number of lags to apply.
#'
#' @return Returns an \eqn{N * (M * lags)} matrix with consecutive lags on the
#' right. The elements of the first \emph{lags} rows are 0.
#'
#' @noRd
lag_var <- function(x, lags) {

  x_rows <- nrow(x)
  x_cols <- ncol(x)

  x_lagged <- matrix(0, x_rows, lags * x_cols)
  for(i in 1:lags) {
    x_lagged[(lags + 1):x_rows, (x_cols * (i - 1) + 1):(x_cols * i)] <-
      x[(lags + 1 - i):(x_rows - i), (1:x_cols)]
  }

  return(x_lagged)
}


#' Compute gamma coefficients
#'
#' Compute the shape \emph{k} and scale \emph{theta} of a Gamma
#' distribution via the mode and standard deviation.
#'
#' @param mode Numeric scalar.
#' @param sd Numeric scalar.
#'
#' @return Returns a list with shape \emph{k} and scale parameter \emph{theta}.
#'
#' @noRd
gamma_coef <- function(mode, sd) {

  mode_sq <- mode ^ 2
  sd_sq <- sd ^ 2
  k <- (2 + mode_sq / sd_sq + sqrt((4 + mode_sq / sd_sq) * mode_sq / sd_sq)) / 2
  theta <- sqrt(sd_sq / k)

  return(list("k" = k, "theta" = theta))
}


#' Name hyperparameters
#'
#' Function to help name hyperparameters. Accounts for multiple occurences
#' of \emph{psi} by adding sequential numbers.
#'
#' @param x Character vector. Parameter names.
#' @param M Integer scalar. Number of columns in the data.
#'
#' @return Returns a character vector of adjusted parameter names.
#'
#' @noRd
name_pars <- function(x, M) {

  out <- Reduce(c, sapply(x, function(y) {
    if(y == "psi") {paste0(y, 1:M)} else {y}}))

  return(out)
}


#' Fill credible intervals
#'
#' Helper function to fill data, colours or similar things based on credible
#' intervals. These are used in \code{\link{plot.bvar_irf}} and
#' \code{\link{plot.bvar_fcast}}.
#'
#' Note that transparency may get appended to recycled HEX colours. Also note
#' that no, i.e. a length 0 central element is required when drawing polygons.
#'
#' @param x Scalar or vector. The central element.
#' @param y Scalar or vector. Value(s) to surround the central element with.
#' The first value is closest, values may get recycled.
#' @param P Odd integer scalar. Number of total bands.
#'
#' @return Returns a vector or matrix (if \emph{x} is a vector) of \emph{x},
#' surrounded by \emph{y}.
#'
#' @noRd
fill_ci <- function(x, y, P) {

  n_y <- if(P %% 2 == 0) {
    stop("No central position for x found.")
  } else {P %/% 2}

  fill <- rep(y, length.out = n_y)

  if(length(x) > 1) { # Matrix
    n_row <- length(x)
    return(cbind(t(rev(fill))[rep(1, n_row), ], x, t(fill)[rep(1, n_row), ]))
  } else { # Vector
    return(c(rev(fill), x, fill))
  }
}


#' @noRd
fill_ci_na <- function(x, P) {

  # Corner case when quantiles are missing (t_back or conditional forecasts)
  if(P == 2) {return(if(length(x > 1)) {cbind(x, NA)} else {c(x, NA)})}

  fill_ci(x = x, y = NA, P = P)
}


#' @noRd
fill_ci_col <- function(x, y, P) {

  # Apply transparency to HEX colours
  if(length(y) == 1 && is_hex(y, alpha = FALSE)) {
    y <- paste0(y, alpha_hex(P))
  }

  fill_ci(x = x, y = y, P = P)
}


#' Get a transparency HEX code
#'
#' @param P Integer scalar. Number of total bands.
#'
#' @return Returns a character vector of transparency codes.
#'
#' @importFrom grDevices rgb
#'
#' @noRd
alpha_hex <- function(P) {

  n_trans <- P %/% 2
  out <- switch(n_trans, # Handpicked with love
    "FF", c("FF", "80"), c("FF", "A8", "54"),
    c("FF", "BF", "80", "40"), c("FF", "CC", "99", "66", "33"))

  if(is.null(out)) { # Let rgb() sort it out otherwise
    out <- substr(rgb(1, 1, 1, seq(1, 0, length.out = n_trans)), 8, 10)
  }

  return(out)
}


#' Check valid HEX colour
#'
#' @param x Character scalar or vector. String(s) to check.
#' @param alpha Logical scalar. Whether the string may contain alpha values.
#'
#' @return Returns a logical scalar or vector.
#'
#' @noRd
is_hex <- function(x, alpha = FALSE) {

  if(alpha) return(grepl("^#[0-9a-fA-F]{6,8}$", x))

  return(grepl("^#[0-9a-fA-F]{6,6}$", x))
}


#' Get variable positions
#'
#' Helper functions to aid with variable selection, e.g. in
#' \code{\link{plot.bvar_irf}} and \code{\link{plot.bvar_fcast}}.
#'
#' @param vars Numeric or character vector of variables to subset to.
#' @param variables Character vector of all variable names. Required if
#' \emph{vars} is provided as character vector.
#' @param M Integer scalar. Count of all variables.
#'
#' @return Returns a numeric vector with the positions of desired variables.
#'
#' @noRd
pos_vars <- function(vars, variables, M) {

  if(is.null(vars) || length(vars) == 0L) {
    return(1:M) # Full set
  }
  if(is.numeric(vars)) {
    return(vapply(vars, int_check, # By position
      min = 1, max = M, msg = "Variable(s) not found.", integer(1)))
  }
  if(is.character(vars) && !is.null(variables)) {
    out <- do.call(c, lapply(vars, grep, variables)) # By name
    if(length(out) > 0) {return(out)}
  }

  stop("Variable(s) not found.")
}


#' Name dependent / explanatory variables
#'
#' @param variables Character vector of all variable names.
#' @param M Integer scalar. Count of all variables.
#' @param lags Integer scalar. Number of lags applied.
#'
#' @return Returns a character vector of variable names.
#'
#' @noRd
name_deps <- function(variables, M) {

  if(is.null(variables)) {
    variables <- paste0("var", seq(M))
  } else if(length(variables) != M) {
    stop("Vector with variables is incomplete.")
  }

  return(variables)
}


#' @noRd
name_expl <- function(variables, M, lags) {

  if(is.null(variables)) {
    variables <- name_deps(variables, M)
  }
  explanatories <- c("constant", paste0(rep(variables, lags), "-lag",
    rep(seq(lags), each = length(variables))))

  return(explanatories)
}


#' Compute log distribution function of Inverse Gamma
#'
#' @param x Numeric scalar. Draw of the IG-distributed variable.
#' @param shape Numeric scalar.
#' @param scale Numeric scalar.
#'
#' @return Returns the log Inverse Gamma distribution function.
#'
#' @noRd
p_log_ig <- function(x, shape, scale) {

  return(shape * log(scale) - (shape + 1) * log(x) - scale / x - lgamma(shape))
}


#' Compute companion matrix
#'
#' Compute the companion form of the VAR coefficients.
#'
#' @param beta Numeric (\eqn{K * M}) matrix with VAR coefficients.
#' @param K Integer scalar. Number of columns in the independent data.
#' @param M Integer scalar. Number of columns in the dependent data.
#' @param lags Integer scalar. Number of lags applied.
#'
#' @return Returns a numeric (\eqn{K - 1 * K -1}) matrix with \emph{beta} in
#' companion form.
#'
#' @noRd
get_beta_comp <- function(beta, K, M, lags) {

  beta_comp <- matrix(0, K - 1, K - 1)

  beta_comp[1:M, ] <- t(beta[2:K, ]) # Kick constant
  if(lags > 1) { # Add block-diagonal matrix beneath VAR coefficients
    beta_comp[(M + 1):(K - 1), 1:(K - 1 - M)] <- diag(M * (lags - 1))
  }

  return(beta_comp)
}


#' Check whether a package is installed
#'
#' @param package Character scalar.
#'
#' @noRd
has_package <- function(package) {

  if(!requireNamespace(package, quietly = TRUE)) {
    stop("Package \'", package, "\' required for this method.", call. = FALSE)
  }

  return(NULL)
}


#' Generate quantiles
#'
#' Check a vector of confidence bands and create quantiles from it.
#'
#' @param conf_bands Numeric vector of probabilities (\eqn{(0, 1)}).
#'
#' @return Returns a sorted, symmetric vector of quantiles.
#'
#' @noRd
quantile_check <- function(conf_bands) {

  conf_bands <- sapply(conf_bands, num_check,
    min = 0 + 1e-16, max = 1 - 1e-16, msg = "Confidence bands misspecified.")

  # Allow only returning the median
  if(length(conf_bands) == 1 && conf_bands == 0.5) {return(conf_bands)}

  # Sort and make sure we have no duplicates (thank mr float)
  quants <- sort(c(conf_bands, 0.5, (1 - conf_bands)))
  quants <- quants[!duplicated(round(quants, digits = 12L))]

  return(quants)
}


#' Prepare BVAR data for methods
#'
#' Helper function to retrieve hyperparameters or coefficient values based on
#' name / position. Also supports multiple \code{bvar} objects and may be used
#' to check them for similarity.
#'
#' @param x A \code{bvar} object, obtained from \code{\link{bvar}}.
#' @param vars Character vector used to select variables. Elements are matched
#' to hyperparameters or coefficients. Coefficients may be matched based on
#' the dependent variable (by providing the name or position) or the
#' explanatory variables (by providing the name and the desired lag). See the
#' example section for a demonstration. Defaults to \code{NULL}, i.e. all
#' hyperparameters.
#' @param vars_response,vars_impulse Optional character or integer vectors used
#' to select coefficents. Dependent variables are specified with
#' \emph{vars_response}, explanatory ones with \emph{vars_impulse}. See the
#' example section for a demonstration.
#' @param chains List with additional \code{bvar} objects. Contents are then
#' added to trace and density plots.
#' @param check_chains Logical scalar. Whether to check \emph{x} and
#' \emph{chains} for similarity.
#' @param ... Fed to \code{\link{chains_fit}}.
#'
#' @return Returns a named list with:
#' \itemize{
#'   \item \code{data} - Numeric matrix with desired data.
#'   \item \code{vars} - Character vector with names for the desired data.
#'   \item \code{chains} - List of numeric matrices with desired data.
#'   \item \code{bounds} - Numeric matrix with optional boundaries.
#' }
#'
#' @noRd
prep_data <- function(
  x,
  vars = NULL,
  vars_response = NULL, vars_impulse = NULL,
  chains = list(),
  check_chains = FALSE, ...) {

  if(!inherits(x, "bvar")) {stop("Please provide a `bvar` object.")}

  if(inherits(chains, "bvar")) {chains <- list(chains)}
  lapply(chains, function(x) {if(!inherits(x, "bvar")) {
    stop("Please provide `bvar` objects to the chains parameter.")
  }})

  if(check_chains) {chains_fit(x, chains, ...)}


  # Prepare selection ---

  vars_hyp <- c("ml", colnames(x[["hyper"]]))
  vars_dep <- x[["variables"]]
  vars_ind <- x[["explanatories"]]
  if(is.null(vars_ind)) { # Compatibility to older versions (<= 0.2.2)
    vars_ind <- name_expl(vars_dep,
      M = x[["meta"]][["M"]], lags = x[["meta"]][["lags"]])
  }

  if(is.null(vars) && is.null(vars_impulse) && is.null(vars_response)) {
    vars <- vars_hyp
  }

  choice_hyp <- vars_hyp[unique(do.call(c, lapply(vars, grep, vars_hyp)))]

  choice_dep <- if(is.null(vars_response)) {
    # Interpret numbers as positions, exclude independents
    vars_dep[unique(c(as.integer(vars[grep("^[0-9]+$", vars)]),
      do.call(c, lapply(vars[!grepl("(^const|lag[0-9]+$)", vars)],
        grep, vars_dep))))]
  } else {pos_vars(vars_response, vars_dep, M = x[["meta"]][["M"]])}
  choice_dep <- choice_dep[!is.na(choice_dep)]

  choice_ind <- if(is.null(vars_impulse)) {
    # Limit to ones with "-lag#" or "constant" to separate from dependents
    vars_ind[unique(do.call(c, lapply(vars[grep("(^const|lag[0-9]+$)", vars)],
      grep, vars_ind)))]
  } else {pos_vars(vars_impulse, vars_ind, M = x[["meta"]][["K"]])}

  if(all(c(length(choice_hyp), length(choice_dep), length(choice_ind)) == 0)) {
    stop("No matching data found.")
  }


  # Build up required outputs ---

  out <- out_vars <- out_bounds <- out_chains <- list()
  N <- x[["meta"]][["n_save"]]

  if(length(choice_hyp) > 0) { # Hyperparameters
    out[["hyper"]] <- cbind("ml" = x[["ml"]], x[["hyper"]])[seq(N), choice_hyp]
    out_vars[["hyper"]] <- choice_hyp
    out_bounds[["hyper"]] <- vapply(choice_hyp, function(z) {
      if(z == "ml") {c(NA, NA)} else {
        c(x[["priors"]][[z]][["min"]], x[["priors"]][[z]][["max"]])
      }}, double(2))
    out_chains[["hyper"]] <- lapply(chains, function(z) {
      cbind("ml" = z[["ml"]], z[["hyper"]])[seq(N), choice_hyp]
    })
  } else {
    out_chains[["hyper"]] <- rep(list(NULL), length(chains))
  }

  if(length(choice_dep) > 0 || length(choice_ind) > 0) { # Betas
    pos_dep <- pos_vars(choice_dep,
      variables = vars_dep, M = x[["meta"]][["M"]])
    pos_ind <- pos_vars(choice_ind,
      variables = vars_ind, M = x[["meta"]][["K"]])
    K <- length(pos_dep) * length(pos_ind)

    out[["betas"]] <- grab_betas(x, N, K, pos_dep, pos_ind)
    out_vars[["betas"]] <- paste0(
      rep(vars_dep[pos_dep], length(pos_ind)), "_",
      rep(vars_ind[pos_ind], each = length(pos_dep)))
    out_bounds[["betas"]] <- matrix(NA, ncol = K, nrow = 2)
    out_chains[["betas"]] <- lapply(chains, grab_betas, N, K, pos_dep, pos_ind)
 } else {
   out_chains[["betas"]] <- rep(list(NULL), length(chains))
 }

  # Merge stuff and return ---

  out_data <- cbind(out[["hyper"]], out[["betas"]])
  out_vars <- c(out_vars[["hyper"]], out_vars[["betas"]])
  out_chains <- mapply(cbind,
    out_chains[["hyper"]], out_chains[["betas"]], SIMPLIFY = FALSE)
  out_chains <- lapply(out_chains, `colnames<-`, out_vars)
  colnames(out_data) <- out_vars

  out <- list(
    "data" = out_data, "vars" = out_vars, "chains" = out_chains,
    "bounds" = cbind(out_bounds[["hyper"]], out_bounds[["betas"]]))

  return(out)
}


#' Grab draws of certain betas
#'
#' Helper function for \code{\link{prep_data}}.
#'
#' @param x A \code{bvar} object, obtained from \code{\link{bvar}}.
#' @param N,K Integer scalars. Number of rows and columns to return.
#' @param pos_dep,pos_ind Numeric vectors. Positions of desired variables.
#'
#' @return Returns a matrix with the requested data.
#'
#' @noRd
grab_betas <- function(x, N, K, pos_dep, pos_ind) {
  data <- matrix(NA, nrow = N, ncol = K)
  k <- 1
  for(i in pos_ind) {for(j in pos_dep) {
    data[, k] <- x[["beta"]][seq(N), i, j] # seq() for longer chains
    k <- k + 1
  }}
  return(data)
}


#' Check equalities across chains
#'
#' Function to help check whether \code{bvar} objects are close enough to
#' compare. Accessed via \code{\link{prep_data}}.
#'
#' @param x A \code{bvar} object, obtained from \code{\link{bvar}}.
#' @param chains List with additional \code{bvar} objects.
#' @param Ms Logical scalar. Whether to check equality of
#' \code{x[["meta"]][["M"]]}.
#' @param n_saves Logical scalar. Whether to check equality of
#' \code{x[["meta"]][["n_save"]]}.
#' @param hypers Logical scalar. Whether to check equality of
#' \code{x[["priors"]][["hyper"]]}.
#'
#' @return Returns \code{TRUE} or throws an error.
#'
#' @noRd
chains_fit <- function(
  x, chains,
  Ms = TRUE,
  n_saves = FALSE,
  hypers = FALSE) {

  if(is.null(chains) || length(chains) == 0) {return(TRUE)}

  if(Ms) {
    Ms <- c(x[["meta"]][["M"]],
      vapply(chains, function(x) {x[["meta"]][["M"]]}, integer(1)))
    if(!all(duplicated(Ms)[-1])) {stop("Number of variables does not match.")}
  }
  if(n_saves) {
    n_saves <- c(x[["meta"]][["n_save"]],
      vapply(chains, function(x) {x[["meta"]][["n_save"]]}, integer(1)))
    if(!all(duplicated(n_saves)[-1])) {
      stop("Number of stored iterations does not match.")
    }
  }
  if(hypers) {
    hypers <- vapply(chains, function(z) {
      x[["priors"]][["hyper"]] == z[["priors"]][["hyper"]]}, logical(1))
    if(!all(hypers)) {stop("Hyperparameters do not match.")}
  }

  return(TRUE)
}