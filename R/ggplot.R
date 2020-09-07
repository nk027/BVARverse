
#' Quick BVAR plots using ggplot2
#'
#'
#'
#' @param x
#' @param type
#' @param vars
#' @param vars_response,vars_impulse
#' @param chains
#' @param t_back
#' @param ... Not used.
#'
#' @return Returns a \code{ggplot} object.
#'
#' @export
bv_ggplot <- function(x, ...) {UseMethod("bv_ggplot", x)}


#' @rdname bv_ggplot
bv_ggplot.default <- function(x, ...) {
  stop("No methods for class ", paste0(class(x), collapse = " / "), " found.")
}


#' @rdname bv_ggplot
bv_ggplot.bvar <- function(x,
  type = c("trace", "density"),
  vars = NULL, vars_response = NULL, vars_impulse = NULL,
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

  df <- fortify.bvar(x,
    vars = vars, vars_response = vars_response, vars_impulse = vars_impulse,
    chains = chains)

  # Trace or density plot
  if(type == "trace") {
    p <- ggplot(df, aes(x = draw)) +
      labs(x = NULL, y = NULL, color = "MCMC chain") +
      theme_bw()

    p <- p + facet_wrap(. ~ variable, scales = "free")

    if(length(chains) == 0) {
      p <- p + geom_line(aes(y = value))
    } else {
      p <- p + geom_line(aes(y = value, color = chain), alpha = 0.5)
    }

  } else if(type == "density") {
    p <- ggplot(df, aes(x = value)) +
      labs(x = NULL, y = NULL, fill = "MCMC chain", color = "MCMC chain") +
      theme_bw() +
      facet_wrap(. ~ variable, scales = "free")

    if(length(chains) == 0) {
      p <- p + geom_density(fill = "grey", alpha = 0.5)
    } else {
      p <- p + geom_density(aes(color = chain, fill = chain), alpha = 0.5)
    }
  }

  # To-do: Add information about quantiles to the plot

  return(p)
}


#' @rdname bv_ggplot
bv_ggplot.bvar_irf <- function(x,
  vars_response = NULL,
  vars_impulse = NULL,
  # col = "#737373",
  ...) {

  df <- fortify.bvar_irf(x)

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

  p <- ggplot(df, aes(x = time - 1)) +
    scale_x_continuous(breaks = function(x) {
      unique(floor(pretty(seq(0, (max(x)) * 1.1))))
    }) +
    labs(x = NULL, y = NULL) +
    theme_bw()

  p <- p + facet_wrap(response ~ impulse, scales = "free_y")

  # col <- fill_ci_col(x = "#000000", y = col, P = P)
  p <- p + geom_line(aes(y = value, col = quantile)) +
    # scale_colour_manual(values = col) +
    geom_hline(yintercept = 0, colour = "darkgray", lty = 2)

  # To-do: Allow ribbons

  return(p)
}


#' @rdname bv_ggplot
bv_ggplot.bvar_fcast <- function(x,
  vars = NULL,
  # col = "#737373",
  t_back = 1L,
  ...) {

  df <- fortify.bvar_fcast(x, t_back)

  P <- length(unique(df[["quantile"]]))
  variables <- x[["variables"]]
  M <- length(variables)

  pos <- pos_vars(vars, variables, M)

  # Only keep the ones specified
  df <- df[which(df[["variable"]] %in% variables[pos]), ]

  df[["variable"]] <- factor(paste(df[["variable"]], "prediction"))

  p <- ggplot(df, aes(x = time)) +
    scale_x_continuous(breaks = function(x) {
      unique(floor(pretty(seq((-t_back + 1), (max(x)) * 1.1))))
    }) +
    labs(x = NULL, y = NULL) +
    theme_bw()

  p <- p + facet_wrap(. ~ variable, scales = "free_y")

  # col <- fill_ci_col(x = "#000000", y = col, P = P)
  p <- p + geom_line(aes(y = value, col = quantile), na.rm = TRUE) +
    # scale_colour_manual(values = col) +
    geom_hline(yintercept = 0, colour = "darkgray", lty = 2)

  # To-do: Allow ribbons

  return(p)
}
