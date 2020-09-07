
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
fortify.bvar <- function(x,
  vars = NULL, vars_response = NULL, vars_impulse = NULL,
  chains = list(), ...) {

  if(!inherits(x, "bvar")) {
    if(inherits(x[[1]], "bvar")) { # Allow chains to x
      chains <- x
      x <- x[[1]]
      chains[[1]] <- NULL
    } else {stop("Please provide a `bvar` object.")}
  }

  x <- prep_data(x,
    vars = vars, vars_response = vars_response, vars_impulse = vars_impulse,
    chains = chains, check_chains = FALSE)

  df <- as.data.frame.table(x[["data"]],
    stringsAsFactors = FALSE, responseName = "value")
  names(df)[1:2] <- c("chain", "variable")
  df[["chain"]] <- 1

  dfs <- vector("list", length(x[["chains"]]) + 1)
  dfs[[1]] <- df
  for(i in seq_along(x[["chains"]])) {
    df_chain <- as.data.frame.table(x[["chains"]][[i]],
      stringsAsFactors = FALSE, responseName = "value")
    names(df_chain) [1:2] <- c("chain", "variable")
    df_chain[["chain"]] <- i + 1
    dfs[[i + 1]] <- df_chain
  }

  out <- do.call(rbind, dfs)
  out[["draw"]] <- seq(nrow(x[["data"]]))
  out[["chain"]] <- as.factor(out[["chain"]])
  out <- out[c(2, 4, 1, 3)]

  return(out)
}


#' @rdname fortify.bvar
fortify.bvar_coefs <- function(x, ...) {

  out <- as.data.frame.table(x,
    stringsAsFactors = FALSE, responseName = "value")

  has_quants <- length(dim(x)) == 3
  if(has_quants) {
    out[["quantile"]] <- as.factor(gsub("([0-9]+)%", "\\1", out[[1]]))
    out[[1]] <- NULL
  } else {
    out[["quantile"]] <- 0.5
  }
  out <- out[, c(2, 1, 3, 4)]
  names(out)[1:2] <- c("variable", "term")

  return(out)
}


#' @rdname fortify.bvar
fortify.bvar_fcast <- function(x, t_back = 0L, ...) {

  H <- x[["setup"]][["horizon"]]
  variables <- x[["variables"]]
  M <- length(variables)

  quants <- x[["quants"]]

  has_quants <- length(dim(quants)) == 3L
  if(has_quants) {
    P <- dim(quants)[1]
    bands <- dimnames(quants)[[1]]
  } else {
    P <- 2 # We make quants 3-dimensional so filling with t_back is easier
    quants <- array(NA, c(2, dim(x[["quants"]])))
    quants[1, , ] <- x[["quants"]]
    bands <- c("50%", "NA")
  }

  # Add t_back actual datapoints
  t_back <- int_check(t_back, 0, Inf, msg = "Issue with t_back.")
  if(t_back != 0) {
    data <- tail(x[["data"]], t_back)
    # Extend the quants array with data, quantiles are set to NA
    quants <- vapply(seq(M), function(i) {
      t(rbind(fill_ci_na(data[, i], P), t(quants[, , i])))
    }, matrix(0, P, t_back + H), USE.NAMES = FALSE)
    dimnames(quants)[[1]] <- bands
  }

  # Names as identifier
  dimnames(quants)[[2]] <- seq.int(-t_back + 1, H)
  dimnames(quants)[[3]] <- variables

  # Prep dataframe
  out <- as.data.frame.table(quants,
    stringsAsFactors = FALSE, responseName = "value")

  # Mess with it
  names(out)[1:3] <- c("quantile", "time", "variable")
  out[["time"]] <- as.integer(out[["time"]])
  out <- out[out[["quantile"]] != "NA", ] # Clean dummy dimension from before
  out <- out[, c(3, 2, 4, 1)]
  out[["quantile"]] <- as.factor(gsub("([0-9]+)%", "\\1", out[["quantile"]]))
  out <- out[order(out[["variable"]], out[["time"]]), ]

  return(out)
}


#' @rdname fortify.bvar
fortify.bvar_irf <- function(x, ...) {

  H <- x[["setup"]][["horizon"]]
  variables <- x[["variables"]]
  quants <- x[["quants"]]

  has_quants <- length(dim(quants)) == 4L
  if(has_quants) {
    P <- dim(quants)[1]
    bands <- dimnames(quants)[[1]]
  } else {
    P <- 2 # We make quants 3-dimensional so filling with t_back is easier
    quants <- array(NA, c(2, dim(x[["quants"]])))
    quants[1, , , ] <- x[["quants"]]
    bands <- c("50%", "NA")
  }

  # Names as identifier
  dimnames(quants)[[2]] <- variables # Response
  dimnames(quants)[[3]] <- seq.int(1, H)
  dimnames(quants)[[4]] <- variables # Impulse

  # Prep dataframe
  out <- as.data.frame.table(quants,
    stringsAsFactors = FALSE, responseName = "value")

  # Mess with it
  names(out)[1:4] <- c("quantile", "response", "time", "impulse")
  out[["time"]] <- as.integer(out[["time"]])
  out <- out[out[["quantile"]] != "NA", ] # Clean dummy dimension from before
  out <- out[, c(4, 2, 3, 5, 1)]
  out[["quantile"]] <- as.factor(gsub("([0-9]+)%", "\\1", out[["quantile"]]))
  out <- out[order(out[["impulse"]], out[["response"]], out[["time"]]), ]

  return(out)
}
