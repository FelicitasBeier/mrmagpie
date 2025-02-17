#' @title readSchulteUebbing
#' @description read the critical nitrogen surplus dataset of Schulte-Uebbing et al. 2022
#' @author Michael Crawford
#'
#' @return a magpie object containing critical nitrogen surplus
#'
#' @examples
#' \dontrun{
#' readSource("SchulteUebbing")
#' }
readSchulteUebbing <- function() {

  filePath <- file.path(
    "Global_critical_N_surpluses_and_N_inputs_and_their_exceedances",
    "Output_files",
    "Critical N surpluses",
    "nsur_crit_mi_all_ph.asc"
  )

  out <- suppressWarnings(
    magclass::read.magpie(filePath) |>
      magpiesets::addLocation(fillMissing = NA) |>
      magclass::setYears(nm = "y2010") |>
      magclass::collapseDim(dim = 1.4) # remove cell number that's also appended
  )

  return(out)
}
