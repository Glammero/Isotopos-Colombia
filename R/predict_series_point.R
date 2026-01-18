# R/predict_series_point.R

predict_series_point <- function(
  lat, lon,
  elev = NA_real_,
  start_date = as.Date("1971-01-15"),
  end_date   = as.Date("2024-12-15"),
  PATHS
){
  suppressPackageStartupMessages({
    library(raster)
    library(rts)
    library(lubridate)
    library(xgboost)
    library(ncdf4)
    library(sf)
  })

  # Fechas mensuales (día 15)
  fechas <- seq(as.Date(start_date), as.Date(end_date), by = "month")

  # Rutas desde config (PATHS)
  predictor_dir   <- PATHS$predictors_dir
  telecon_path    <- file.path(predictor_dir, "TELE.csv")
  dem_path        <- PATHS$dem_path
  kop_raster_path <- file.path(predictor_dir, "koepclim/koepclim.grd")
  kop_table_path  <- file.path(predictor_dir, "koepclim/koep.cat.csv")

  # Modelos
  modelo18_path <- file.path(PATHS$models_dir, "d18O", "d18O.model1RDS")
  pp18_path     <- file.path(PATHS$models_dir, "d18O", "d18O.pp")
  modelo2H_path <- file.path(PATHS$models_dir, "d2H",  "d2H.model1RDS")
  pp2H_path     <- file.path(PATHS$models_dir, "d2H",  "d2H.pp")

  # Incertidumbres (RMSE)
  RMSE_d18O <- 2.319
  RMSE_d2H  <- 16.311

  # Chequeos mínimos de existencia
  must_exist <- c(
    predictor_dir, telecon_path, dem_path, kop_raster_path, kop_table_path,
    modelo18_path, pp18_path, modelo2H_path, pp2H_path
  )
  for (p in must_exist){
    if (!file.exists(p) && !dir.exists(p)){
      stop("No existe en el entorno: ", p)
    }
  }

  # Carga modelos
  bst18  <- readRDS(modelo18_path)
  pp18   <- readRDS(pp18_path)
  feat18 <- bst18$feature_names

  bst2H  <- readRDS(modelo2H_path)
  pp2H   <- readRDS(pp2H_path)
  feat2H <- bst2H$feature_names

  # ----------------------------
  # Helpers generales
  # ----------------------------
  `%||%` <- function(a,b) if (!is.null(a)) a else b

  # Reemplazo seguro de hasSlot (NO usar methods::hasSlot)
  has_slot <- function(x, s){
    isS4(x) && (s %in% methods::slotNames(x))
  }

  pick_subdir <- function(vn){
    if (startsWith(vn, "CRU."))  return("CRU")
    if (startsWith(vn, "ERA5.")) return("ERA5")
    if (startsWith(vn, "NCEP.")) return("NCEP")
    ""
  }

  # --- FIX: reescribir rutas absolutas Windows embebidas en objetos Raster (.rds) ---
  fix_raster_file_paths <- function(r, predictors_dir){
    if (!inherits(r, "Raster")) return(r)

    fn <- raster::filename(r)
    if (!is.character(fn) || !nzchar(fn)) return(r)

    fnn <- gsub("\\\\", "/", fn)

    # Raíz vieja (Windows) embebida al crear los .rds
    old_root <- "B:/PisoAI_master_data/data/predictors"

    # Caso 1: reemplazar raíz vieja por la del repo
    if (startsWith(fnn, old_root)) {
      rel  <- sub(paste0("^", old_root, "/?"), "", fnn)
      cand <- file.path(predictors_dir, rel)
      cand <- gsub("\\\\", "/", cand)
      if (file.exists(cand)) {
        r@file@name <- cand
        return(r)
      }
    }

    # Caso 2: fallback: buscar en _cache_grd del repo por basename
    if (!file.exists(fnn)) {
      cand2 <- file.path(predictors_dir, "_cache_grd", basename(fnn))
      cand2 <- gsub("\\\\", "/", cand2)
      if (file.exists(cand2)) {
        r@file@name <- cand2
        return(r)
      }
    }

    r
  }

  .rts_cache <- new.env(parent = emptyenv())

  get_rts_layer <- function(vn, fecha, base_dir){
    static_cru <- c("CRU.DEM","CRU.cont","CRU.trng","CRU.atp.e")

    subd <- pick_subdir(vn)
    base <- if (nzchar(subd)) file.path(base_dir, subd) else base_dir

    if (vn %in% static_cru){
      nc <- file.path(base, paste0(vn, ".nc"))
      if (file.exists(nc)) {
        rl <- raster(nc)
        rl <- fix_raster_file_paths(rl, base_dir)
        return(rl)
      }
    }

    rds <- file.path(base, paste0(vn, ".rds"))
    nc  <- file.path(base, paste0(vn, ".nc"))

    if (file.exists(rds)){

      if (!exists(vn, envir = .rts_cache)){
        obj <- readRDS(rds)

        # FIX rutas embebidas en obj@raster si existe
        if (has_slot(obj, "raster")) {
          rr <- tryCatch(obj@raster, error = function(e) NULL)
          if (!is.null(rr) && inherits(rr, "Raster")) {
            obj@raster <- fix_raster_file_paths(rr, base_dir)
          }
        }

        assign(vn, obj, envir = .rts_cache)
      }

      obj <- get(vn, envir = .rts_cache)
      tt  <- index(obj)
      idx <- which(tt == fecha)
      if (!length(idx)) return(NULL)

      rb <- tryCatch(obj@raster, error = function(e) stack(obj))
      rb <- fix_raster_file_paths(rb, base_dir)

      lyr <- raster::subset(rb, idx[1])
      lyr <- fix_raster_file_paths(lyr, base_dir)

      return(lyr)

    } else if (file.exists(nc)){

      rb <- brick(nc)
      rb <- fix_raster_file_paths(rb, base_dir)

      if (nlayers(rb) == 1) return(rb)
      return(NULL)
    }

    NULL
  }

  fuzzy_point_value <- function(rl, lon, lat, fuzz = 0.5){
    if (is.null(rl)) return(NA_real_)
    pts <- rbind(
      c(lon, lat),
      c(lon, lat + fuzz),
      c(lon, lat - fuzz),
      c(lon + fuzz, lat),
      c(lon - fuzz, lat)
    )
    vals <- raster::extract(rl, pts, method = "simple")
    v0   <- vals[1]
    fill <- mean(vals[2:5], na.rm = TRUE)
    if (is.na(v0)) return(fill) else return(v0)
  }

  # --- Reemplazo de predict(preProcess, ...) SIN caret ---
  apply_preprocess_center_scale <- function(pp, newdata){
    x <- newdata

    if (!is.null(pp$mean)) {
      mu <- pp$mean
      common <- intersect(names(mu), names(x))
      for (nm in common) x[[nm]] <- as.numeric(x[[nm]]) - as.numeric(mu[[nm]])
    }

    if (!is.null(pp$std)) {
      sdv <- pp$std
      common <- intersect(names(sdv), names(x))
      for (nm in common) {
        s <- as.numeric(sdv[[nm]])
        if (is.finite(s) && s != 0) x[[nm]] <- as.numeric(x[[nm]]) / s
      }
    }

    x
  }

  # ----------------------------
  # Insumos estáticos: plantilla, DEM, Köppen, TELE, punto
  # ----------------------------

  plantilla_path <- file.path(predictor_dir, "CRU", "CRU.cld.rds")
  if (!file.exists(plantilla_path)) {
    stop("No encuentro plantilla CRU.cld.rds en: ", plantilla_path)
  }

  plantilla_obj <- readRDS(plantilla_path)
  if (has_slot(plantilla_obj, "raster")) {
    rr <- tryCatch(plantilla_obj@raster, error = function(e) NULL)
    if (!is.null(rr) && inherits(rr, "Raster")) {
      plantilla_obj@raster <- fix_raster_file_paths(rr, predictor_dir)
    }
  }

  plantilla_any <- tryCatch(raster::subset(plantilla_obj@raster, 1), error=function(e) {
    raster::subset(stack(plantilla_obj), 1)
  })
  plantilla_any <- fix_raster_file_paths(plantilla_any, predictor_dir)

  crs_target <- tryCatch(sf::st_crs(raster::crs(plantilla_any)), error=function(e) NA)
  if (is.na(crs_target)) crs_target <- sf::st_crs(4326)

  pt_sf <- sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon","lat"), crs = 4326)
  pt_sf <- sf::st_transform(pt_sf, crs_target)
  pt_xy <- as.data.frame(sf::st_coordinates(pt_sf))[1,]

  DEM <- raster(dem_path)

  kop_r   <- raster(kop_raster_path)
  kop_tbl <- read.csv(kop_table_path, stringsAsFactors = FALSE)

  TELE <- read.csv(telecon_path, stringsAsFactors = FALSE)
  TELE$Date <- as.Date(TELE$Date)

  # ----------------------------
  # build_features_point
  # ----------------------------
  build_features_point <- function(features_m, fecha){
    out <- list()

    vars_r     <- features_m[grepl("^(CRU|ERA5|NCEP)\\.", features_m)]
    static_cru <- c("CRU.DEM","CRU.cont","CRU.trng","CRU.atp.e")

    for (vn in vars_r){
      rl <- get_rts_layer(vn, fecha, predictor_dir)

      if (vn %in% static_cru){
        out[[vn]] <- as.numeric(raster::extract(
          rl, matrix(c(pt_xy$X, pt_xy$Y), ncol = 2), method = "bilinear"
        ))
      } else {
        out[[vn]] <- fuzzy_point_value(rl, pt_xy$X, pt_xy$Y, fuzz = 0.5)
      }
    }

    if ("elevation" %in% features_m){
      if (!is.na(elev)) {
        out[["elevation"]] <- as.numeric(elev)
      } else {
        out[["elevation"]] <- as.numeric(raster::extract(
          DEM, matrix(c(pt_xy$X, pt_xy$Y), ncol = 2), method = "simple"
        ))
      }
    }

    if ("Longitude" %in% features_m) out[["Longitude"]] <- as.numeric(pt_xy$X)
    if ("Latitude"  %in% features_m) out[["Latitude"]]  <- as.numeric(pt_xy$Y)
    if ("lat2"      %in% features_m){
      base_lat <- out[["Latitude"]] %||% as.numeric(pt_xy$Y)
      out[["lat2"]] <- base_lat^2
    }

    fila <- TELE[TELE$Date == fecha, , drop = FALSE]
    if (nrow(fila) == 1){
      if ("NINO34" %in% features_m) out[["NINO34"]] <- as.numeric(fila$NINO34[1])
      if ("NINO12" %in% features_m) out[["NINO12"]] <- as.numeric(fila$NINO12[1])
    } else {
      if ("NINO34" %in% features_m) out[["NINO34"]] <- NA_real_
      if ("NINO12" %in% features_m) out[["NINO12"]] <- NA_real_
    }

    m <- month(fecha)
    for (k in 1:12){
      nm <- paste0("month.", k)
      if (nm %in% features_m) out[[nm]] <- as.integer(m == k)
    }

    est <- if (m %in% c(12,1,2)) "DJF" else if (m %in% 3:5) "MAM" else if (m %in% 6:8) "JJA" else "SON"
    for (ss in c("DJF","MAM","JJA","SON")){
      nm <- paste0("season.", ss)
      if (nm %in% features_m) out[[nm]] <- as.integer(ss == est)
    }

    kop_cols <- grep("^climate\\.", features_m, value = TRUE)
    if (length(kop_cols)){
      idv <- raster::extract(kop_r, matrix(c(pt_xy$X, pt_xy$Y), ncol = 2), method = "simple")
      lab <- kop_tbl$climate[match(idv, kop_tbl$ID)]
      for (cn in kop_cols){
        clase <- sub("^climate\\.", "", cn)
        out[[cn]] <- as.integer(!is.na(lab) && lab == clase)
      }
    }

    miss <- setdiff(features_m, names(out))
    if (length(miss)) for (mm in miss) out[[mm]] <- 0

    as.data.frame(out[features_m], check.names = FALSE)
  }

  # ----------------------------
  # Predicción
  # ----------------------------
  col_d2H  <- sprintf("d2H±%s",  formatC(RMSE_d2H,  format = "f", digits = 3))
  col_d18O <- sprintf("d18O±%s", formatC(RMSE_d18O, format = "f", digits = 3))

  res <- vector("list", length(fechas))

  for (i in seq_along(fechas)){
    f <- fechas[i]

    row18 <- build_features_point(feat18, f)
    if (any(!is.finite(unlist(row18)))) next
    row18_s <- apply_preprocess_center_scale(pp18, row18)
    for (j in seq_along(row18_s)) if (!is.numeric(row18_s[[j]])) row18_s[[j]] <- as.numeric(row18_s[[j]])
    pred18  <- as.numeric(predict(bst18, xgb.DMatrix(data = data.matrix(row18_s))))

    row2H <- build_features_point(feat2H, f)
    if (any(!is.finite(unlist(row2H)))) next
    row2H_s <- apply_preprocess_center_scale(pp2H, row2H)
    for (j in seq_along(row2H_s)) if (!is.numeric(row2H_s[[j]])) row2H_s[[j]] <- as.numeric(row2H_s[[j]])
    pred2H  <- as.numeric(predict(bst2H, xgb.DMatrix(data = data.matrix(row2H_s))))

    elev_used <- if (!is.na(elev)) elev else
      as.numeric(raster::extract(DEM, matrix(c(pt_xy$X, pt_xy$Y), ncol = 2), method = "simple"))

    tmp <- data.frame(
      date = f,
      elevation_used = elev_used,
      stringsAsFactors = FALSE
    )
    tmp[[col_d2H]]  <- pred2H
    tmp[[col_d18O]] <- pred18
    res[[i]] <- tmp
  }

  res_df <- do.call(rbind, res)
  if (is.null(res_df) || nrow(res_df) == 0L){
    stop("No se generaron filas (todas las fechas fueron omitidas). Revisa TELE y disponibilidad de capas.")
  }

  res_df
}
