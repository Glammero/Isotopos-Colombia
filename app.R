# app.R  (en la raíz del repo)

library(shiny)

# 1) Cargar rutas (relativas al repo)
source("config.R")

# 2) Cargar todo tu código en R/
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source))

# 3) Wrapper: aquí conectaremos TU función real de predicción
predict_wrapper <- function(lat, lon, elev = NA_real_) {
  predict_series_point(
    lat = lat,
    lon = lon,
    elev = elev,
    PATHS = PATHS
  )
}

ui <- fluidPage(
  titlePanel("Isótopos de precipitación — Colombia (serie temporal por punto)"),
  sidebarLayout(
    sidebarPanel(
      numericInput("lat", "Latitud", value = 4.65, min = -5, max = 15, step = 0.0001),
      numericInput("lon", "Longitud", value = -74.1, min = -82, max = -66, step = 0.0001),
      checkboxInput("use_dem", "Calcular elevación desde DEM", value = TRUE),
      conditionalPanel(
        condition = "input.use_dem == false",
        numericInput("elev", "Elevación (m)", value = 2600, min = -50, max = 6000, step = 1)
      ),
      actionButton("run", "Consultar"),
      br(), br(),
      downloadButton("download_csv", "Descargar CSV")
    ),
    mainPanel(
      verbatimTextOutput("status"),
      plotOutput("ts_plot", height = "420px"),
      tableOutput("ts_table")
    )
  )
)

server <- function(input, output, session) {

  res <- eventReactive(input$run, {
    lat  <- input$lat
    lon  <- input$lon
    elev <- if (isTRUE(input$use_dem)) NA_real_ else input$elev

    tryCatch({
      df <- predict_wrapper(lat, lon, elev)

      if (!is.data.frame(df)) {
        stop("La predicción no devolvió un data.frame.")
      }

      list(ok = TRUE, df = df, err = NULL)
    }, error = function(e) {
      list(ok = FALSE, df = NULL, err = conditionMessage(e))
    })
  }, ignoreInit = TRUE)

  output$status <- renderText({
    req(res())
    if (isTRUE(res()$ok)) {
      paste("OK. Filas:", nrow(res()$df), "| Columnas:", paste(names(res()$df), collapse = ", "))
    } else {
      paste("ERROR:", res()$err)
    }
  })

  output$ts_plot <- renderPlot({
    req(res())
    req(isTRUE(res()$ok))
    df <- res()$df

    date_col <- names(df)[tolower(names(df)) %in% c("date", "fecha", "time", "t")]
    if (length(date_col) == 0) stop("No encuentro columna de fecha (date/fecha).")

    d <- df[[date_col[1]]]
    if (!inherits(d, "Date")) d <- as.Date(d)

    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No encuentro columnas numéricas para graficar.")

    y <- df[[num_cols[1]]]
    plot(d, y, type = "l", xlab = "Fecha", ylab = num_cols[1])
  })

  output$ts_table <- renderTable({
    req(res())
    req(isTRUE(res()$ok))
    head(res()$df, 20)
  })

  output$download_csv <- downloadHandler(
    filename = function() paste0("isotopos_", Sys.Date(), ".csv"),
    content = function(file) {
      req(res())
      if (!isTRUE(res()$ok)) stop(res()$err)
      write.csv(res()$df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
