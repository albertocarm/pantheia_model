#' Launch the Pantheia Prognostic Calculator
#'
#' Starts the Shiny application locally in the user's browser.
#'
#' @param host Host address. Default "127.0.0.1" (localhost).
#' @param port Port number. Default lets Shiny choose an available port.
#' @param launch.browser Whether to open the app in the browser. Default TRUE.
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pantheia()
#' }
pantheia <- function(host = "127.0.0.1", port = NULL, launch.browser = TRUE, ...) {
  app_dir <- system.file("app", package = "pantheia_model")
  if (app_dir == "") {
    stop("Could not find the app directory. Try reinstalling the package.")
  }
  shiny::runApp(
    app_dir,
    host = host,
    port = port,
    launch.browser = launch.browser,
    ...
  )
}
