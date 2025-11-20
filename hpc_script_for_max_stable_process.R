# ---- Libraries ----  
library(readxl)
library(openxlsx)
library(SpatialExtremes)
library(furrr)
library(future)
library(tidyverse)

# ---------------- SIGNIFICANCE FUNCTION ----------------
check_significance <- function(est, se) {
  if (is.null(se) || all(is.na(se))) {
    return(tibble(
      Estimate = est,
      SE = NA_real_,
      z = NA_real_,
      p = NA_real_,
      Significance = NA_character_
    ))
  }
  z <- est / se
  p <- 2 * (1 - pnorm(abs(z)))
  signif <- ifelse(p < 0.05, "Significant", "Not Significant")

  tibble(Estimate = est, SE = se, z = z, p = p, Significance = signif)
}

# ---- Paths ----
base_path <- "/home/pposa/MelbourneDocuments/4SpatialExtremes/"
results_dir <- file.path(base_path, "results")
dir.create(results_dir, showWarnings = FALSE)

path_data <- file.path(base_path, "data/Qmax_nrm1.xlsx")
path_vars <- file.path(base_path, "data/NRM_info1.xlsx")
path_temp_vars <- file.path(base_path, "data/temp_cov.xlsx")
pics_dir <- file.path(base_path, "fmadograms")
source(file.path(base_path, "Rcodes/formulae_all.R"))

# ---- Load Temp Covariates ----
temp_cov <- read_excel(path_temp_vars, sheet = "ann_max") %>%
  mutate(across(everything(), as.numeric))

# ============================================================
# =============== HARD-CODED BEST MODEL TABLES ================
# ============================================================

bm_s_win <- tribble(
  ~sheet,    ~model,     ~formula_no,
  "CS",      "powexp",     9,
  "EC",      "brown",     63,
  "MB",      "gwhitmat", 103,
  "MN_L",    "brown",    109,
  "MN_R",    "brown",    103,
  "R_R",     "brown",      1,
  "R_L",     "bessel",    14,
  "SSWF_L",  "brown",     81,
  "SSWF_R",  "cauchy",   109,
  "SS1",     "brown",    119,
  "SS2",     "brown",    123,
  "WT",      "brown",     63
) %>% mutate(type = "s")

bm_st_win <- tribble(
  ~sheet,    ~model,     ~formula_no,
  "CS",      "powexp",     9,
  "EC",      "brown",   1992,
  "MB",      "brown",   2533,
  "MN_L",    "brown",   3368,
  "MN_R",    "brown",   3509,
  "R_R",     "brown",    129,
  "R_L",     "powexp",  2354,
  "SSWF_L",  "brown",   2539,
  "SSWF_R",  "cauchy",   109,
  "SS1",     "brown",    119,
  "SS2",     "brown",   3133,
  "WT",      "brown",   3190
) %>% mutate(type = "st")

# Combine both tables
all_best_models <- bind_rows(bm_s_win, bm_st_win)

# ---------------------------------------------------------
sheets1 <- all_best_models$sheet
windows <- 1:20

plan(multisession, workers = 12)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ============================================================
# =============== PROCESSING FUNCTION =========================
# ============================================================

process_region_window <- function(entry, w) {

  sheet <- entry$sheet
  model <- entry$model
  formula_no <- entry$formula_no
  type <- entry$type

  dat <- read_excel(path_data, sheet = sheet)
  spa_cov <- read_excel(path_vars, sheet = sheet)
  locations <- as.matrix(spa_cov[, c("Longitude", "Latitude")])

  scov <- spa_cov %>%
    select(Longitude, Latitude, Ele_RANGE_zonal, meanAnnP, Cat_Area__, Slope_MEAN_zonal) %>%
    setNames(c("lon", "lat", "ele", "meanP", "cat_a", "slope")) %>%
    as.matrix()

  df_all <- dat[, 3:ncol(dat)] %>% mutate(across(everything(), as.numeric))
  if ((w + 29) > nrow(df_all)) return(NULL)

  df <- df_all[w:(w + 29), ] %>%
    sweep(2, scov[, "cat_a"], FUN = "/") %>%
    as.data.frame()

  tcov <- if (type == "st") temp_cov[w:(w + 29), ] %>% as.matrix() else NULL

  fs <- formula_sets[[formula_no]]

  # ---- Fit MSP Model ----
  m1 <- tryCatch({
    fitmaxstab(
      data = as.matrix(df), coord = locations, cov.mod = model,
      loc.form = fs$loc.form, scale.form = fs$scale.form, shape.form = fs$shape.form,
      temp.form.loc   = fs$temp.form.loc   %||% NULL,
      temp.form.scale = fs$temp.form.scale %||% NULL,
      temp.form.shape = fs$temp.form.shape %||% NULL,
      marg.cov = scov, temp.cov = tcov,
      iso = FALSE, fit.marge = TRUE
    )
  }, error = function(e) return(NULL))

  if (is.null(m1)) return(NULL)

  # Extract and format parameters
  params <- tryCatch(m1$param, error = function(e) NULL)

  if (is.null(params) || length(params) == 0) {
    param_df <- tibble(region = sheet, type = type, window = w)
  } else {
    se <- tryCatch({
      if (is.null(m1$var.cov)) stop("var.cov NULL")
      vdiag <- diag(m1$var.cov)
      if (any(is.na(vdiag))) stop("NA diag")
      sqrt(vdiag)
    }, error = function(e) {
      warning(paste("SE not available for", sheet, "win", w, ":", conditionMessage(e)))
      rep(NA_real_, length(params))
    })

    param_df_long <- tibble(
      Parameter = names(params),
      Estimate  = as.numeric(params),
      StdError  = as.numeric(se),
      Significance = ifelse(is.na(se), NA, ifelse(2*(1-pnorm(abs(params/se)))<0.05,"Significant","Not Significant"))
    )

    param_order <- c(
      "nugget", "range", "smooth",
      paste0("locCoeff", 1:6),
      paste0("scaleCoeff", 1:6),
      "shapeCoeff1",
      paste0("tempCoeffLoc", 1:5),
      paste0("tempCoeffScale", 1:5),
      "tempCoeffShape1"
    )

    keep_params <- intersect(param_order, param_df_long$Parameter)

    param_df <- param_df_long %>%
      filter(Parameter %in% keep_params) %>%
      pivot_wider(
        names_from = Parameter,
        values_from = c(Estimate, StdError, Significance),
        names_glue = "{Parameter}_{.value}"
      ) %>%
      rename_with(~ gsub("_Estimate$", "", .x), ends_with("_Estimate")) %>%
      rename_with(~ gsub("_StdError$", "_SE", .x), ends_with("_StdError")) %>%
      rename_with(~ gsub("_Significance$", "_Sig", .x), ends_with("_Significance"))

    col_order <- unlist(lapply(keep_params, function(p)
      c(p, paste0(p, "_SE"), paste0(p, "_Sig"))))

    param_df <- param_df[, col_order, drop = FALSE] %>%
      mutate(region = sheet, type = type, window = w)
  }

  # ---- F-madogram + TIC ----
  spa_dep_ratio <- NA_real_
  tryCatch({
    spa_dep <- as.data.frame(madogram(df, locations, fitted = m1, which = "ext", col = c(1, 2)))
    filtered_results <- spa_dep %>% filter(ext.coeff < 1.5)
    spa_dep_ratio <- round(nrow(filtered_results) / nrow(spa_dep), 2)

    png(file.path(pics_dir, paste0("fmado_", type, "_", sheet, "_win", w, ".png")),
        width = 8, height = 6, units = "in", res = 300)
    par(mar = c(3.5, 4, 2, 2))
    madogram(df, locations, fitted = m1, which = "ext", col = c(1, 2),
             ylab = "", xlab = "", cex = 1.5, cex.axis = 1.5, lwd = 2)
    mtext("Distance", side = 1, line = 2.5, cex = 1.5)
    mtext("ext.coeff", side = 2, line = 2.5, cex = 1.5)
    title(main = paste("Madogram -", sheet, "(win:", w, ", sdr:", spa_dep_ratio, ")"), cex.main = 1.5)
    dev.off()
  }, error = function(e) message("F-madogram error: ", e$message))

  extcoe <- m1$ext.coeff(seq(0.5, 6, 0.5))

  tic_summary <- cbind(
    data.frame(reg = sheet, tic = TIC(m1), sd_ratio = spa_dep_ratio, type = type, win = w),
    setNames(as.data.frame(t(extcoe)), paste0(seq(0.5, 6, 0.5), "deg"))
  )

  rl <- predict(m1, newdata = scov, ret.per = c(2, 5, 25, 50, 75, 100))
  rl <- rl * scov[, "cat_a"]

  rl_df <- as.data.frame(rl) %>%
    mutate(lon = scov[, "lon"], lat = scov[, "lat"],
           region = sheet, type = type, window = w)

  list(rl = rl_df, tic = tic_summary)
}

# ============================================================
# =============== RUN PARALLEL LOOP ==========================
# ============================================================

all_combinations <- expand_grid(entry_id = 1:nrow(all_best_models), win = windows)

results_all <- future_map(1:nrow(all_combinations), function(j) {
  row <- all_combinations[j, ]
  entry <- all_best_models[row$entry_id, ]
  process_region_window(entry, row$win)
}, .options = furrr_options(seed = TRUE)) %>% compact()

# ============================================================
# =============== SAVE RESULTS ================================
# ============================================================

return_levels_all <- bind_rows(map(results_all, "rl"))
tic_summary_all  <- bind_rows(map(results_all, "tic"))

types <- c("s", "st")

for (tpe in types) {

  rl_tpe <- return_levels_all %>% filter(type == tpe)

  for (reg in unique(rl_tpe$region)) {
    wb <- createWorkbook()
    for (w in windows) {
      addWorksheet(wb, paste0("win_", w))
      writeData(wb, paste0("win_", w),
                rl_tpe %>% filter(region == reg, window == w) %>%
                  select(-region, -type, -window))
    }
    saveWorkbook(wb, file.path(results_dir, paste0("RL_", tpe, "_", reg, ".xlsx")), overwrite = TRUE)
  }

  wb_tic <- createWorkbook()
  tdf <- tic_summary_all %>% filter(type == tpe)
  for (w in windows) {
    addWorksheet(wb_tic, paste0("win_", w))
    writeData(wb_tic, paste0("win_", w),
              tdf %>% filter(win == w) %>%
                select(reg, tic, sd_ratio, all_of(paste0(seq(0.5, 6, 0.5), "deg"))))
  }
  saveWorkbook(wb_tic, file.path(results_dir, paste0("TIC_summary_", tpe, ".xlsx")), overwrite = TRUE)
}
