#' @title calcNonLocalProduction
#' @description Calculates grid-level amount of food that would need to be transported,
#' assuming that food produced in the grid cell is first consumed by local population
#' i.e. amount of food greater than local rural demand, split into that which feeds the local
#' urban population, and that which exceeds total local demand and is available to export
#' @param cells magpiecell or lpjcell (default)
#' @author David M Chen
#' @examples
#' \dontrun{
#' calcOutput("NonLocalTransport")
#' }
calcNonLocalProduction <- function(cells = "lpjcell") {

  productionPri  <- calcOutput("Production",  products = "kcr",
                               cellular = TRUE,
                               aggregate = FALSE)
  productionLi   <- calcOutput("Production", products = "kli",
                               cellular = TRUE,
                               aggregate = FALSE)
  productionPast <- calcOutput("Production", products = "pasture", cellular = TRUE,
                               aggregate = FALSE)
  productionPast <- add_dimension(productionPast,   dim = 3.1, add = "Item", nm = "pasture")
  production <- collapseNames(mbind(productionPri, productionLi)[, , "dm"])
  production <- mbind(production, collapseNames(productionPast[, , "dm"], collapsedim  = 2))

  foodDisaggUrb <- calcOutput("FoodDemandGridded", feed = TRUE,
                              cells = cells, aggregate = FALSE)
  # note this also includes feed demand!

  foodDemPrim <- collapseNames(foodDisaggUrb[, , getItems(production, dim = 3)])

  transpRurToUrb <-  ifelse(production > dimSums(foodDemPrim, dim = 3.2),
                            collapseNames(foodDemPrim[, , "urban"]),
                            collapseNames(production - foodDemPrim[, , "rural"]))
  transpRurToUrb[which(transpRurToUrb < 0)] <- 0
  transpRurToUrb <- add_dimension(transpRurToUrb, dim = 3.2, nm = "RurToUrb")


  transpExp <- ifelse(production > dimSums(foodDemPrim, dim = 3.2),
                      collapseNames(production - dimSums(foodDemPrim, dim = 3.2)),
                      0)
  transpExp <- add_dimension(transpExp, dim = 3.2, nm = "Exported")

  out <- mbind(transpRurToUrb, transpExp)
  # all secondary products need to be transported, but we don't know where they are produced.

  return(list(x = out,
              weight = NULL,
              unit = "Mt dm",
              description = "Food produced that is transported to
                          local city, and amount that is transported and exported",
              isocountries = FALSE))
}
