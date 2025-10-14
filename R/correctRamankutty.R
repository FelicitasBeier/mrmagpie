#' @title correctRamankutty
#' @description Read Available Land Si
#'
#' @param x magpie object provided by the read function
#'
#' @return magpie object
#' @author Felicitas Beier
#'
#' @seealso
#'   \code{\link{readRamankutty}}
#' @examples
#'
#' \dontrun{
#'   readSource("Ramankutty", convert="onlycorrect")
#' }
correctRamankutty <- function(x) {

  x <- toolConditionalReplace(x, conditions = c("is.na()", "<0"), replaceby = 0)

  return(x)
}
