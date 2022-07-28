#' print.malt
#'
#' @param x an output of the malt function
#' @param digits the number of significative digits
#' @param ... other parameters
#'
#' @export
print.malt=function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nLast draw:\n")
  print.default(format(x$draw, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\nAcceptance rate: ", paste(format(x$accept, digits = digits)))
    # cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
  #     "\n\n", sep = "")
  # if (length(coef(x))) {
  #   cat("Coefficients:\n")
  #   print.default(format(coef(x), digits = digits), print.gap = 2L,
  #                 quote = FALSE)
  # }
  # else cat("No coefficients\n")
  # cat("\n")
  invisible(x)
}
