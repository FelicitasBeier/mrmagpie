#' @title downloadSchulteUebbing
#' @description Download the critical nitrogen surplus values of Schulte-Uebbing et al. 2022.
#'   This function requires the 'archive' package, which is suggested. If 'archive' is not available,
#'   the function will stop with an error.
#' @return A list containing metadata from Schulte-Uebbing et al. 2022.
#' @examples
#' \dontrun{
#'   downloadSchulteUebbing()
#' }

downloadSchulteUebbing <- function() {

  # Check if the 'archive' package is available.
  if (!requireNamespace("archive", quietly = TRUE)) {
    stop("Package 'archive' is required for this function. Please install it or manually load the data.")
  }

  # Define URL and destination for the 7z archive.
  url7z <- "https://zenodo.org/records/6395016/files/Global_critical_N_surpluses_and_N_inputs_and_their_exceedances.7z?download=1" # nolint
  dest7z <- "Global_critical_N_surpluses_and_N_inputs_and_their_exceedances.7z"

  # Download the file if it does not already exist.
  if (!file.exists(dest7z)) {
    tryCatch({
      utils::download.file(url7z, destfile = dest7z, mode = "wb")
    }, error = function(e) {
      stop("Download failed: ", e$message)
    })
  }

  # Create an extraction directory if it doesn't exist.
  extractDir <- "Global_critical_N_surpluses_and_N_inputs_and_their_exceedances"
  if (!dir.exists(extractDir)) {
    dir.create(extractDir)
  }

  # Extract the 7z archive using the archive package.
  archive::archive_extract(dest7z, dir = extractDir)

  # nolint start
  list(
    url = "https://doi.org/10.1038/s41586-022-05158-2",
    doi = "10.1038/s41586-022-05158-2",
    title = "From Planetary to Regional Boundaries for Agricultural Nitrogen Pollution",
    unit = "Various",
    author = list(
      person("Lena", "Schulte-Uebbing"),
      person("Arthur", "Beusen"),
      person("Lex", "Bouwman"),
      person("Wim", "De Vries")
    ),
    release_date = "2022-10-12",
    description = "This study presents a framework for translating global planetary boundaries for agricultural nitrogen pollution into regional boundaries. It provides a detailed analysis of nitrogen inputs and surpluses, addressing their environmental impacts across different regions.",
    license = "Creative Commons Attribution 4.0 International License",
    reference = "Schulte-Uebbing, L., Beusen, A., Bouwman, L., and De Vries, W. (2022). From Planetary to Regional Boundaries for Agricultural Nitrogen Pollution. Nature, 610(7932), 507-512. doi:10.1038/s41586-022-05158-2"
  )
  # nolint end

}
