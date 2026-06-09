#' Automatic recoding into factors
#'
#' @param formula model expression
#' @param data data frame
#' @return recoded data frame
#' @export
auto_factor <- function(formula, data, max_levels = 10) {

  mf <- model.frame(formula, data = data)
  predictor_names <- attr(terms(mf), "term.labels")
  X <- mf[, predictor_names, drop = FALSE]

  for (nm in predictor_names) {

    x <- X[[nm]]

    if (is.ordered(x)) {
      contrasts(x) <- cum_contr(x)
      X[[nm]] <- x
    } else if (!is.factor(x)) {
      nunique <- length(unique(x[!is.na(x)]))
      if (nunique <= max_levels) {
        X[[nm]] <- factor(x)
      }
    }
  }

  yname <- as.character(formula[[2]])
  result <- data.frame(
    mf[[yname]],
    X,
    check.names = FALSE
  )

  names(result)[1] <- yname

  return(result)
}
