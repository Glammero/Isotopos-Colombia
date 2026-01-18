# app.R (entrypoint para Connect Cloud)
# Carga el archivo real de la app (isotopos_app.R) y devuelve el objeto Shiny.

res <- source("isotopos_app.R", local = TRUE)
app <- res$value

if (inherits(app, "shiny.appobj")) {
  app
} else if (exists("ui", inherits = TRUE) && exists("server", inherits = TRUE)) {
  shiny::shinyApp(ui = ui, server = server)
} else {
  stop("No se encontró shiny.appobj ni ui/server después de source('isotopos_app.R').")
}
