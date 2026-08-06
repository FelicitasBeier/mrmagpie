#' @title correctGAMI
#' @description Correct the GAMI v2.1 forest-age dataset: replace missing and negative
#' class-fraction values (ocean / non-forest fill) by 0. The object stays on GAMI's native
#' 0.5-degree grid; the projection onto the 67420 lpj-cell grid happens in
#' \code{calcAgeClassDistribution}.
#' @param x magpie object provided by \code{readGAMI}
#' @return magpie object on GAMI's native 0.5-degree grid, cleaned
#' @author Florian Humpenoeder
#' @seealso \code{\link{readGAMI}}
#' @examples
#' \dontrun{
#' readSource("GAMI", convert = "onlycorrect")
#' }
correctGAMI <- function(x) {
  x[is.na(x)] <- 0
  x[x < 0]    <- 0
  return(x)
}
