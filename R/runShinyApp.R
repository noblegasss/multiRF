#' Run Shiny App
#'
#' This function launches the Shiny application.
#'
#' @export runShinyApp
runShinyApp <- function() {
  appDir <- system.file("shiny-app", package = "multiRF")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `multiRF`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
