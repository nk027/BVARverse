#' Quick ggplot2 plots for Bayesian VARs
#'
#' Function to quickly plot outputs from \code{bvar} and derived objects.
#' Supported plots include traces and densities, forecasts, and impulse
#' response functions. For more flexible plots one may use the outputs of
#' \code{\link{tidy.bvar}} and \code{\link{augment.bvar}}.
#'
#' @inheritParams tidy.bvar
#' @param type A string with the type (trace or density) of plot desired.
#' @param orientation A string indicating the desired orientation of trace or
#' density plots
#' @param col Character vector. Colour(s) of the lines delineating credible
#' intervals. Single values will be recycled if necessary. Recycled HEX color
#' codes are varied in transparency if not provided (e.g. "#737373FF"). Lines
#' can be bypassed by setting this to \code{"transparent"}.
#'
#' @return Returns a \code{ggplot} object with a basic structure.
#'
#' @import BVAR
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes labs theme_bw geom_line geom_density
#' @importFrom ggplot2 geom_hline facet_wrap scale_x_continuous
#' @importFrom ggplot2 scale_color_manual
#'
#' @export
#'
#' @examples
#' # Access a subset of the fred_qd dataset
#' data <- fred_qd[, c("CPIAUCSL", "UNRATE", "FEDFUNDS")]
#' # Transform it to be stationary
#' data <- fred_transform(data, codes = c(5, 5, 1), lag = 4)
#'
#' # Estimate a BVAR using one lag, default settings and very few draws
#' x <- bvar(data, lags = 1, n_draw = 1000L, n_burn = 200L, verbose = FALSE)
#'
#' # Plot the outputs - alternatively use ggplot() with fortify()
#' bv_ggplot(x)
#' bv_ggplot(irf(x))
#' bv_ggplot(predict(x))
bv_ggplot <- function(x, ...) {UseMethod("bv_ggplot", x)}


#' @rdname bv_ggplot
#' @export
bv_ggplot.default <- function(x, ...) {
  stop("No methods for class ", paste0(class(x), collapse = " / "), " found.")
}


#' @rdname bv_ggplot
#' @export
bv_ggplot.bvar_chains <- function(x, ...) {
  bv_ggplot.bvar(x, ...)
}


#' @rdname bv_ggplot
#' @export
bv_ggplot.bvar <- function(x,
  type = c("trace", "density"),
  vars = NULL, vars_response = NULL, vars_impulse = NULL,
  orientation = c("horizontal", "vertical"),
  chains = list(),
  ...) {

  if(!inherits(x, "bvar")) {
    if(inherits(x[[1]], "bvar")) { # Allow chains to x
      chains <- x
      x <- x[[1]]
      chains[[1]] <- NULL
    } else {stop("Please provide a `bvar` object.")}
  }

  type <- match.arg(type)
  orientation <- match.arg(orientation)

  df <- tidy.bvar(x,
    vars = vars, vars_response = vars_response, vars_impulse = vars_impulse,
    chains = chains)

  # Trace or density plot
  if(type == "trace") {
    p <- ggplot(df, aes(x = .data$draw)) +
      labs(x = NULL, y = NULL, color = "MCMC chain") +
      theme_bw()

    p <- p + facet_wrap(. ~ .data$variable, scales = "free",
                        nrow = ifelse(orientation == "horizontal", 1, 
                                      length(unique(df[["variable"]]))))

    if(length(chains) == 0) {
      p <- p + geom_line(aes(y = .data$value))
    } else {
      p <- p +
        geom_line(aes(y = .data$value, color = .data$chain), alpha = 0.5)
    }

  } else if(type == "density") {
    p <- ggplot(df, aes(x = .data$value)) +
      labs(x = NULL, y = NULL, fill = "MCMC chain", color = "MCMC chain") +
      theme_bw() +
      facet_wrap(. ~ .data$variable, scales = "free")

    if(length(chains) == 0) {
      p <- p + geom_density(fill = "grey", alpha = 0.5)
    } else {
      p <- p +
        geom_density(aes(color = .data$chain, fill = .data$chain), alpha = 0.5)
    }
  }
  
  # To-do: Add information about quantiles to the plot

  return(p)
}


#' @rdname bv_ggplot
#' @export
bv_ggplot.bvar_irf <- function(x,
  vars_response = NULL,
  vars_impulse = NULL,
  col = "#737373",
  ...) {

  df <- tidy.bvar_irf(x)

  P <- length(unique(df[["quantile"]]))
  variables <- x[["variables"]]
  M <- length(variables)

  pos_imp <- pos_vars(vars_impulse, variables, M)
  pos_res <- pos_vars(vars_response, variables, M)

  # Only keep the ones specified
  df <- df[intersect(which(df[["impulse"]] %in% variables[pos_imp]),
    which(df[["response"]] %in% variables[pos_res])), ]

  df[["impulse"]] <- factor(paste(df[["impulse"]], "shock"))
  df[["response"]] <- factor(paste(df[["response"]], "response"))

  p <- ggplot(df, aes(x = .data$time - 1)) +
    scale_x_continuous(breaks = function(x) {
      unique(floor(pretty(seq(0, (max(x)) * 1.1))))
    }) +
    labs(x = NULL, y = NULL) +
    theme_bw()

  p <- p + facet_wrap(.data$response ~ .data$impulse, scales = "free_y")

  col <- fill_ci_col(x = "#000000", y = col, P = P)
  p <- p + geom_line(aes(y = .data$value, col = .data$quantile)) +
    scale_color_manual(values = col, 
                        breaks = as.character(unique(df$quantile))) +
    geom_hline(yintercept = 0, colour = "darkgray", lty = 2)

  # To-do: Allow ribbons

  return(p)
}


#' @rdname bv_ggplot
#' @export
bv_ggplot.bvar_fcast <- function(x,
  vars = NULL,
  col = "#737373",
  t_back = 1L,
  ...) {

  df <- tidy.bvar_fcast(x, t_back)

  P <- length(unique(df[["quantile"]]))
  variables <- x[["variables"]]
  M <- length(variables)

  pos <- pos_vars(vars, variables, M)

  # Only keep the ones specified
  df <- df[which(df[["variable"]] %in% variables[pos]), ]

  df[["variable"]] <- factor(paste(df[["variable"]], "prediction"))

  p <- ggplot(df, aes(x = .data$time)) +
    scale_x_continuous(breaks = function(x) {
      unique(floor(pretty(seq((-t_back + 1), (max(x)) * 1.1))))
    }) +
    labs(x = NULL, y = NULL) +
    theme_bw()

  p <- p + facet_wrap(. ~ .data$variable, scales = "free_y")

  col <- fill_ci_col(x = "#000000", y = col, P = P)
  p <- p +
    geom_line(aes(y = .data$value, col = .data$quantile), na.rm = TRUE) +
    scale_color_manual(values = col, 
                        breaks = as.character(unique(df$quantile))) +
    geom_hline(yintercept = 0, colour = "darkgray", lty = 2)

  # To-do: Allow ribbons

  return(p)
}
