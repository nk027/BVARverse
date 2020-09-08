
#' @importFrom generics augment
#' @export
generics::augment


#' Augment BVAR outputs and convert into a tibble
#'
#' Turn the outputs of a Bayesian VAR (see \code{\link[BVAR]{bvar}}) into a
#' an augmented tibble. Methods are available for \code{bvar} objects (will
#' yield coefficients and their quantiles), \code{bvar_fcast} objects (with
#' predictions, their quantiles and optionally real datapoints), and
#' \code{bvar_irf} objects (with impulse responses).
#'
#' @param x A \code{bvar} or derived object to turn into a tibble.
#' @param conf_bands Numeric vector. Credible intervals of coefficients to
#' include in the tibble.
#' @param t_back Integer scalar. Whether to include actual datapoints in the
#' tidied forecast.
#' @param ... Not used.
#'
#' @return Returns a \code{\link[tibble]{tibble}} with relevant information;
#' quantiles can be found in the columns.
#'
#' @import BVAR
#' @importFrom rlang .data
#' @importFrom tidyr pivot_wider
#' @importFrom stats coef
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Access a subset of the fred_qd dataset
#' data <- fred_qd[, c("CPIAUCSL", "UNRATE", "FEDFUNDS")]
#' # Transform it to be stationary
#' data <- fred_transform(data, codes = c(5, 5, 1), lag = 4)
#'
#' # Estimate a BVAR using one lag, default settings and very few draws
#' x <- bvar(data, lags = 1, n_draw = 1000L, n_burn = 200L, verbose = FALSE)
#'
#' # Create tibbles from the outputs
#' augment(x)
#' augment(irf(x))
#' augment(predict(x))
#' }
augment.bvar <- function(x, conf_bands = 0.16, ...) {

  coefs <- coef(x, conf_bands = conf_bands)

  df <- tidy.bvar_coefs(coefs, conf_bands = conf_bands)

  out <- pivot_wider(df, names_from = .data$quantile, names_prefix = "q")

  return(out)
}


#' @rdname augment.bvar
#' @export
augment.bvar_fcast <- function(x, t_back = 0L, ...) {

  df <- tidy.bvar_fcast(x, t_back = t_back)

  out <- pivot_wider(df, names_from = .data$quantile, names_prefix = "q")

  return(out)
}


#' @rdname augment.bvar
#' @export
augment.bvar_irf <- function(x, ...) {

  df <- tidy.bvar_irf(x)

  out <- pivot_wider(df, names_from = .data$quantile, names_prefix = "q")

  return(out)
}
