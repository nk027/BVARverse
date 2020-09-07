
#' Function for tidying up Bayesian VAR outputs
#'
#' Function for tidy summarizing output obtained from \code{\link{bvar}}
#' in tibble format.
#' All or a subset of the available variables can be summarized.
#'
#' @param x A \code{bvar} object, obtained from \code{\link{bvar}}.
#' @param vars Character vector used to select variables. Elements are matched
#' to hyperparameters or coefficients. Coefficients may be matched based on
#' the dependent variable (by providing the name or position) or the
#' explanatory variables (by providing the name and the desired lag). See the
#' example section for a demonstration. Defaults to \code{NULL}, i.e. all
#' hyperparameters.
#' @param vars_impulse,vars_response Optional character or integer vectors used
#' to select coefficents. Dependent variables are specified with
#' \emph{vars_response}, explanatory ones with \emph{vars_impulse}. Default to
#' \code{NULL}, indicating that no coefficients will be processed.
#' @param variables Optional character vector. Names of all variables in the
#' object. Used to subset and title. Taken from \code{x$variables} if available.
#' @param quants Optional numeric vector for computing quantiles of parameter
#' draws.
#' @param chains List of \code{bvar} objects. Contents are then added to trace
#' and density plots to help assessing covergence.
#'
#' @return Returns a \code{\link[tibble]{tibble}} with relevant information
#' for further processing.
#'
#' @export
tidy.bvar <- function(x, conf_bands = 0.16, ...) {

  coefs <- coef.bvar(x, conf_bands = conf_bands)

  df <- fortify.bvar(coefs, conf_bands = conf_bands)

  out <- pivot_wider(df, names_from = quantile, names_prefix = "q")

  return(out)
}


#' @rdname tidy.bvar
tidy.bvar_fcast <- function(x, t_back = 0L, ...) {

  df <- fortify.bvar_fcast(x, t_back = t_back)

  out <- pivot_wider(df, names_from = quantile, names_prefix = "q")

  return(out)
}


#' @rdname tidy.bvar
tidy.bvar_irf <- function(x, ...) {

  df <- fortify.bvar_irf(x)

  out <- pivot_wider(df, names_from = quantile, names_prefix = "q")

  return(out)
}
