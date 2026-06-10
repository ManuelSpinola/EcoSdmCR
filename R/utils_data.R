# ============================================================
# utils_data.R
# Carga de covariables CHELSA pre-procesadas por resolución H3
# OVS-CR · ICOMVIS · UNA
# ============================================================

# ── Rutas de datos ───────────────────────────────────────────
DATA_DIR <- "data"

# Mapeo resolución → archivo de covariables
covariables_paths <- list(
  actual = list(
    "6" = file.path(DATA_DIR, "bio_chelsa_presente_no_cor_cr_5367_h6.gpkg"),
    "7" = file.path(DATA_DIR, "bio_chelsa_presente_no_cor_cr_5367_h7.gpkg"),
    "8" = file.path(DATA_DIR, "bio_chelsa_presente_no_cor_cr_5367_h8.gpkg")
  ),
  futuro = list(
    "6" = file.path(DATA_DIR, "bio_chelsa_futuro_cr_5367_h6.gpkg"),
    "7" = file.path(DATA_DIR, "bio_chelsa_futuro_cr_5367_h7.gpkg"),
    "8" = file.path(DATA_DIR, "bio_chelsa_futuro_cr_5367_h8.gpkg")
  )
)

# ── Cache en memoria (se carga una sola vez por sesión) ─────
.cov_cache <- new.env(parent = emptyenv())

#' Cargar covariables para una resolución dada
#'
#' @param resolucion character "6", "7" o "8"
#' @param escenario  character "actual" o "futuro"
#' @return sf con h3_address + variables (o NULL si el archivo no existe)
cargar_covariables <- function(resolucion, escenario = "actual") {
  key <- paste0(escenario, "_res", resolucion)

  if (exists(key, envir = .cov_cache)) {
    return(get(key, envir = .cov_cache))
  }

  path <- covariables_paths[[escenario]][[as.character(resolucion)]]

  if (is.null(path) || !file.exists(path)) {
    warning(sprintf(
      "[utils_data] Archivo no encontrado: %s\n  Prepará el gpkg con h3sdm_extract_num() y colocálo en data/",
      path %||% "(ruta no definida)"
    ))
    return(NULL)
  }

  dat <- sf::st_read(path, quiet = TRUE)
  assign(key, dat, envir = .cov_cache)
  dat
}

#' Verificar disponibilidad de covariables
#' @return data.frame con estado de cada archivo
verificar_covariables <- function() {
  res <- do.call(rbind, lapply(c("6", "7", "8"), function(r) {
    data.frame(
      resolucion = r,
      actual     = file.exists(covariables_paths$actual[[r]]),
      futuro     = file.exists(covariables_paths$futuro[[r]]),
      stringsAsFactors = FALSE
    )
  }))
  res
}

#' Operador %||% (null coalescing)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Filtro de outliers ambientales ───────────────────────────

#' Eliminar outliers ambientales en presencias usando distancia de Mahalanobis
#'
#' Calcula la distancia de Mahalanobis de cada presencia al centroide
#' ambiental del conjunto. Registros con D² > chi2(p = 0.975, df = k)
#' se consideran outliers y se excluyen del dataset PA.
#'
#' @param pa        sf con columnas `presence` ("1"/"0") y covariables CHELSA
#' @param vars_cov  character vector con nombres de las covariables a usar
#' @param umbral_p  percentil de la distribución chi² (default 0.975)
#'
#' @return lista con:
#'   - pa_limpio:    sf sin los outliers de presencia
#'   - outliers_df:  data.frame con los registros eliminados y su D²
#'   - n_eliminados: integer
#'   - umbral:       valor numérico del umbral D²
filtrar_outliers_ambientales <- function(pa, vars_cov, umbral_p = 0.975) {

  # Separar presencias
  idx_pres <- which(pa$presence == "1")
  pres_df  <- sf::st_drop_geometry(pa)[idx_pres, vars_cov, drop = FALSE]

  # Eliminar filas con NA en covariables
  completos <- complete.cases(pres_df)
  if (sum(completos) < (length(vars_cov) + 1)) {
    warning("[filtrar_outliers] Muy pocas presencias completas para calcular Mahalanobis.")
    return(list(
      pa_limpio    = pa,
      outliers_df  = data.frame(),
      n_eliminados = 0L,
      umbral       = NA_real_
    ))
  }

  pres_completa <- pres_df[completos, , drop = FALSE]

  # Centroide y covarianza
  centro <- colMeans(pres_completa)
  cov_m  <- tryCatch(
    cov(pres_completa),
    error = function(e) NULL
  )

  if (is.null(cov_m) || det(cov_m) < .Machine$double.eps) {
    warning("[filtrar_outliers] Matriz de covarianza singular — se omite el filtro.")
    return(list(
      pa_limpio    = pa,
      outliers_df  = data.frame(),
      n_eliminados = 0L,
      umbral       = NA_real_
    ))
  }

  # Distancia de Mahalanobis para todas las presencias con datos completos
  d2     <- mahalanobis(pres_completa, center = centro, cov = cov_m)
  umbral <- qchisq(umbral_p, df = length(vars_cov))

  # Índices originales en pa que son outliers
  idx_completos_en_pa <- idx_pres[completos]
  idx_outliers_en_pa  <- idx_completos_en_pa[d2 > umbral]

  outliers_df <- if (length(idx_outliers_en_pa) > 0) {
    df <- sf::st_drop_geometry(pa)[idx_outliers_en_pa, , drop = FALSE]
    df$mahal_d2 <- d2[d2 > umbral]
    df
  } else {
    data.frame()
  }

  pa_limpio <- if (length(idx_outliers_en_pa) > 0) pa[-idx_outliers_en_pa, ] else pa

  list(
    pa_limpio    = pa_limpio,
    outliers_df  = outliers_df,
    n_eliminados = length(idx_outliers_en_pa),
    umbral       = umbral
  )
}

# ── Info de resoluciones ─────────────────────────────────────
info_resoluciones <- data.frame(
  res   = c("6", "7", "8"),
  label = c(
    "6 — 36.1 km² (paisaje local)",
    "7 — 5.2 km² (SDM fino ★)",
    "8 — 0.74 km² (escala de sitio)"
  ),
  area  = c("36.1 km²", "5.2 km²", "0.74 km²"),
  stringsAsFactors = FALSE
)
