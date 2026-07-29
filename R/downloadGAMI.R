#' @title downloadGAMI
#' @description Downloads the Global Age Mapping Integration (GAMI) v2.1 forest-age dataset
#' (Besnard et al. 2024, GFZ Data Services, doi:10.5880/GFZ.1.4.2023.006), the 0.5-degree
#' `class_fraction` product read by \code{\link{readGAMI}}. GAMI provides the within-forest
#' age distribution (12 age classes, 20-member ensemble) for 2010 and 2020.
#' @return metadata list describing the downloaded GAMI source
#' @author Florian Humpenoeder
#' @seealso [downloadSource()], \code{\link{readGAMI}}
#' @examples
#' \dontrun{
#' downloadSource("GAMI")
#' }
downloadGAMI <- function() {
  # GFZ Data Services download folder for doi:10.5880/GFZ.1.4.2023.006 (GAMI v2.1)
  baseURL <- "https://datapub.gfz-potsdam.de/download/10.5880.GFZ.1.4.2023.006-VEnuo"
  file    <- "GAMIv2-1_2010-2020_class_fraction_0deg50.nc"

  utils::download.file(file.path(baseURL, file), destfile = file, mode = "wb")

  meta <- list(
    url          = "https://doi.org/10.5880/GFZ.1.4.2023.006",
    doi          = "10.5880/GFZ.1.4.2023.006",
    title        = "Global Age Mapping Integration (GAMI)",
    author       = paste("Besnard, Simon; Santoro, Maurizio; Herold, Martin; Cartus, Oliver;",
                         "Guetter, Jonas; Heinrich, Viola H.A.; Herault, Bruno; Kassi, Justin;",
                         "N'Guessan, Anny; Koirala, Sujan; Neigh, Christopher; Nelson, Jacob A.;",
                         "Poulter, Benjamin; Weber, Ulrich; Zhang, Tao; Carvalhais, Nuno"),
    version      = "2.1",
    release_date = "2025-01-07",
    description  = paste("Global within-forest age distribution (12 age classes) at 0.5 degree for",
                         "2010 and 2020, 20-member ensemble (class_fraction product)."),
    license      = "CC BY 4.0",
    reference    = "Besnard et al. (2024), GFZ Data Services, doi:10.5880/GFZ.1.4.2023.006",
    unit         = "fraction"
  )
  return(meta)
}
