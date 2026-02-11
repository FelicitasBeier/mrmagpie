#' @title downloadMehta2024
#' @description download Global Area Equipped for Irrigation Dataset 1900-2015 from Mehta et al. (2024).
#'              Gridded dataset is created based on (sub-)national statistics from FAOSTAT, AQUASTAT, EUROSTAT
#'              and country's census data downscaled using two alternative gridded irrigation maps
#'              (GMIA from Siebert et al. 2013 and Meier et al. 2018)
#'
#' @param subtype data subtype to be downloaded
#'                Subtypes available:
#'                Combination of version (v3, v4) and data source:
#'                'GMIA': gridded base map for downscaling from Stefan et al. (2013).
#'                        Global Map of Irrigation Areas version 5.
#'                'Meier2018': gridded base map for downscaling from Meier, et al. (2018).
#'                             Global Irrigated Areas.
#'                Separated by "_" (e.g., subtype = "v4_GMIA")
#'
#' @author  Felicitas Beier
#' @seealso [downloadSource()] [readMehta2024()]
#' @examples
#' \dontrun{
#' a <- downloadSource()
#' }
#'
#' @importFrom utils download.file
#' @importFrom withr with_options

downloadMehta2024 <- function(subtype = "v4_GMIA") {

  # extract version and subtype
  subtype <- unlist(strsplit(subtype, split = "_"))
  version <- subtype[1]
  subtype <- subtype[2]

  if (version == "v4") {
    aeiURL <- "https://zenodo.org/records/14219723"
    doi    <- "10.5281/zenodo.6740334"
    releaseDate <- "2024-11-26"
  } else if (version == "v3") {
    aeiURL <- "https://zenodo.org/records/7809342"
    doi    <- "10.5281/zenodo.6740334"
    releaseDate <- "2023-04-07"
  } else {
    stop("Version does not exist or is not implemented in downloadMehta2024.")
  }

  if (subtype == "GMIA") {
    dataname <- "G_AEI_"
  } else if (subtype == "Meier2018") {
    dataname <- "MEIER_G_AEI_"
  } else {
    stop("The selected subtype is not available for downloadMehta2024. Please select 'GMIA' or 'Meier2018'.")
  }

  years <- c(seq(1900, 1970, by = 10),
             seq(1980, 2015, by = 5))

  files <- c(paste0(dataname, years, ".ASC"))

  # repeat with count and maxcount is necessary for case of 'Couldn't connect to server' error
  maxcount <- 10
  count    <- 0
  for (file in files) {
    repeat {

      withr::with_options(list(timeout = NULL),
                          code = try(download.file(paste0(aeiURL, "/files/", file),
                                                   destfile = file, mode = "wb")))

      count <- count + 1
      if (file.exists(file) || count >= maxcount) {
        break
      }
    }
  }

  return(list(url          = aeiURL,
              doi          = doi,
              title        = "Global Area Equipped for Irrigation Dataset 1900-2015",
              revision     = "2024",
              version      = version,
              release_date = releaseDate,
              publication  = "https://doi.org/10.1038/s44221-024-00206-9",
              author       = "Mehta, Piyush;  Siebert, Stefan;
                              Kummu, Matti; Deng, Qinyu; Ali, Tariq;
                              Marston, Landon; Xie, Wei;  Davis, Kyle",
              description  = paste0("Gridded dataset of Areas Equipped for Irrigation. ",
                                    "Created based on (sub-)national statistics from FAOSTAT, ",
                                    "AQUASTAT, EUROSTAT and country's census data ",
                                    "downscaled using two alternative gridded irrigation maps ",
                                    "(GMIA from Siebert et al. 2013 and Meier et al. 2018)"),
              unit         = "ha",
              license      = "Creative Commons Attribution 4.0 International"))
}
