# deconv_functions.R

#' @importFrom magrittr %>%
NULL
#' Split barcodes into pseudo-coordinates  (e.g., "s_008um_00792_00584-1" -> "792, 584")
#'
#' @param barcode_list A character vector of Visium barcodes (e.g., "s_008um_00792_00584-1")
#' @return A data.frame containing the original barcode and extracted row and col values
#' @importFrom dplyr select
#' @importFrom tidyr separate
#' @export
barcodes_to_pseudo_coordinates <- function(barcode_list) {

  if (missing(barcode_list) || !is.vector(barcode_list) || !is.character(barcode_list)) {
    stop("`barcode_list` must be a non-empty character vector.")
  }
  if (any(is.na(barcode_list))) {
    stop("`barcode_list` contains NA values; please remove or replace them.")
  }
  df <- data.frame(barcodes = barcode_list, stringsAsFactors = FALSE)
  # Attempt to separate into four parts: prefix1, prefix2, rows, cols
  parts <- strsplit(barcode_list, "_", fixed = TRUE)
  invalid <- vapply(parts, function(x) length(x) < 4, logical(1))
  if (any(invalid)) {
    stop("Some entries in `barcode_list` do not conform to expected pattern 'prefix1_prefix2_rows_cols[-1]'.")
  }

  df <- df %>%
    separate(barcodes,
             into = c("prefix1", "prefix2", "rows", "cols"),
             sep = "_",
             remove = FALSE) %>%
    select(-prefix1, -prefix2)
  df$cols <- gsub("-1", "", df$cols)

  if (any(is.na(as.integer(df$rows))) || any(is.na(as.integer(df$cols)))) {
    stop("Failed to parse numeric `rows` or `cols` from some barcodes.")
  }
  df$rows <- as.integer(df$rows)
  df$cols <- as.integer(df$cols)

  return(df)
}





#' Identify high-resolution barcodes to discard based on low-resolution barcodes
#'
#' @param barcode_list A character vector of low-resolution barcodes
#' @param low_res Integer value indicating low-resolution bin size (e.g., 64)
#' @param high_res Integer value indicating high-resolution bin size (e.g., 8)
#' @return A character vector of high-resolution barcodes to discard
#' @importFrom stringr str_split
#' @export
barcode_based_deconvolve <- function(barcode_list, low_res, high_res) {

  if (missing(barcode_list) || !is.vector(barcode_list) || !is.character(barcode_list)) {
    stop("`barcode_list` must be a non-empty character vector.")
  }
  if (length(barcode_list) == 0) {
    return(character(0))
  }
  if (!is.numeric(low_res) || !is.numeric(high_res) ||
      length(low_res) != 1 || length(high_res) != 1 ||
      low_res != floor(low_res) || high_res != floor(high_res) ||
      low_res <= 0 || high_res <= 0) {
    stop("`low_res` and `high_res` must be positive integer values.")
  }


  if (low_res < high_res) {
    stop("low_res must be larger than high_res.")
  }
  if (low_res %% high_res != 0) {
    stop("low_res must be a multiple of high_res.")
  }
  factor <- low_res / high_res
  prefix <- sprintf("%03dum", high_res)
  discard_barcodes <- character()

  for (bc in barcode_list) {

    if (!is.character(bc) || is.na(bc) || !grepl("_", bc, fixed = TRUE)) {
      stop("Each element of `barcode_list` must be a non-NA string containing '_' separators.")
    }

    bc_untail <- gsub("-1", "", bc)

    parts <- unlist(str_split(bc_untail, "_", simplify = TRUE))
    if (length(parts) < 4) {
      stop(paste0("Barcode '", bc, "' does not conform to expected pattern."))
    }

    row_low <- as.integer(parts[3])
    col_low <- as.integer(parts[4])
    if (is.na(row_low) || is.na(col_low)) {
      stop(paste0("Failed to parse numeric row/col from barcode '", bc, "'."))
    }
    rows_high <- seq(row_low * factor,
                     row_low * factor + factor - 1)
    cols_high <- seq(col_low * factor,
                     col_low * factor + factor - 1)
    for (r in rows_high) {
      for (c in cols_high) {
        discard_barcodes <- c(discard_barcodes,
                              paste0("s_",
                                     prefix, "_",
                                     sprintf("%05d", r), "_",
                                     sprintf("%05d", c),
                                     "-1"))
      }
    }
  }
  return(unique(discard_barcodes))
}



#' Identify high-resolution barcodes to discard based on spatial coordinates
#'
#' @param barcode_list A character vector of low-resolution barcodes
#' @param parquet_low_res Path to the low-resolution .parquet file
#' @param parquet_high_res Path to the high-resolution .parquet file
#' @param json_low_res Path to the low-resolution scale factors .json file
#' @param low_res Integer value indicating low-resolution bin size (e.g., 64)
#' @param high_res Integer value indicating high-resolution bin size (e.g., 8)
#' @return A character vector of high-resolution barcodes to discard
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr filter
#' @importFrom arrow read_parquet
#' @export
coordinate_based_deconvolve <- function(barcode_list,
                                        parquet_low_res,
                                        parquet_high_res,
                                        json_low_res,
                                        low_res,
                                        high_res) {

  if (missing(barcode_list) || !is.vector(barcode_list) || !is.character(barcode_list)) {
    stop("`barcode_list` must be a non-empty character vector.")
  }
  if (length(barcode_list) == 0) {
    return(character(0))
  }
  if (!is.numeric(low_res) || !is.numeric(high_res) ||
      length(low_res) != 1 || length(high_res) != 1 ||
      low_res != floor(low_res) || high_res != floor(high_res) ||
      low_res <= 0 || high_res <= 0) {
    stop("`low_res` and `high_res` must be positive integer values.")
  }

  if (low_res < high_res) {
    stop("low_res must be larger than high_res.")
  }
  if (low_res %% high_res != 0) {
    stop("low_res must be a multiple of high_res.")
  }

  if (!is.character(parquet_low_res) || !file.exists(parquet_low_res)) {
    stop("`parquet_low_res` must be a valid path to an existing .parquet file.")
  }
  if (!is.character(parquet_high_res) || !file.exists(parquet_high_res)) {
    stop("`parquet_high_res` must be a valid path to an existing .parquet file.")
  }
  if (!is.character(json_low_res) || !file.exists(json_low_res)) {
    stop("`json_low_res` must be a valid path to an existing .json file.")
  }

  data <- fromJSON(json_low_res)
  if (!"microns_per_pixel" %in% names(data)) {
    stop("`json_low_res` does not contain 'microns_per_pixel' element.")
  }

  microns_per_pixel <- data$microns_per_pixel
  if (!is.numeric(microns_per_pixel) || length(microns_per_pixel) != 1 || microns_per_pixel <= 0) {
    stop("`microns_per_pixel` must be a positive numeric value in JSON.")
  }

  parquet_low_res <- arrow::read_parquet(parquet_low_res) %>% filter(in_tissue == 1)
  if (!all(c("barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres") %in% colnames(parquet_low_res))) {
    stop("`parquet_low_res` must contain columns: 'barcode', 'in_tissue', 'pxl_row_in_fullres', 'pxl_col_in_fullres'.")
  }
  parquet_high_res <- arrow::read_parquet(parquet_high_res) %>% filter(in_tissue == 1)
  if (!all(c("barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres") %in% colnames(parquet_high_res))) {
    stop("`parquet_high_res` must contain columns: 'barcode', 'in_tissue', 'pxl_row_in_fullres', 'pxl_col_in_fullres'.")
  }

  discard_barcodes <- character()

  for (br in barcode_list) {

    if (!is.character(br) || is.na(br)) {
      stop("Each element of `barcode_list` must be a non-NA string.")
    }
    ctr <- parquet_low_res %>% filter(barcode == br)

    if (nrow(ctr) == 0) {
      warning(paste0("Barcode '", br, "' not found in low-resolution data; skipping."))
      next
    }

    half_side_px <- (low_res / 2) / 4.649492204882606  # microns_per_pixel in JSON
    rmin <- ctr$pxl_row_in_fullres - half_side_px
    rmax <- ctr$pxl_row_in_fullres + half_side_px
    cmin <- ctr$pxl_col_in_fullres - half_side_px
    cmax <- ctr$pxl_col_in_fullres + half_side_px
    child <- parquet_high_res %>%
      filter(pxl_row_in_fullres >= rmin,
             pxl_row_in_fullres <= rmax,
             pxl_col_in_fullres >= cmin,
             pxl_col_in_fullres <= cmax) %>%
      pull(barcode)
    discard_barcodes <- c(discard_barcodes, child)
  }
  return(unique(discard_barcodes))
}



#' Plot spatial QC overlaying spots on the background image
#'
#' @param spe SpatialExperiment object or similar containing spatial coordinates and images
#' @param sample_id Character string specifying the column in colData that contains sample IDs
#' @param sample Character string specifying which sample to plot (defaults to the first unique value)
#' @param metric Character string specifying the column in colData to color by (e.g., "detected")
#' @param outliers NULL or character string specifying the column in colData that marks outliers/highlight groups. Can be boolean or factor.
#' @param point_size Numeric value for spot size
#' @param colors A character vector of length 2 (gradient) or more for fill colors
#' @param stroke Numeric value for border size of highlighted spots
#' @return A ggplot object
#' @importFrom escheR make_escheR add_fill
#' @importFrom ggplot2 annotation_raster coord_equal geom_point scale_fill_gradient scale_fill_gradientn scale_color_discrete guides guide_legend
#' @export
plot_spatial_qc <- function(spe,
                            sample_id = "sample_id",
                            sample = unique(colData(spe)[[sample_id]])[1],
                            metric = "detected",
                            outliers = NULL,
                            point_size = 2, colors = c("white", "black"), stroke = 1) {

  if (missing(spe) || !("colData" %in% methods::slotNames(spe))) {
    stop("`spe` must be a SpatialExperiment-like object with colData and spatialCoords methods.")
  }
  if (!sample_id %in% colnames(as.data.frame(colData(spe)))) {
    stop(paste0("`sample_id` '", sample_id, "' not found in colData(spe)."))
  }
  sample_values <- unique(colData(spe)[[sample_id]])
  if (!sample %in% sample_values) {
    stop(paste0("`sample` '", sample, "' not found among colData(spe)[['", sample_id, "']]."))
  }
  if (!metric %in% colnames(as.data.frame(colData(spe)))) {
    stop(paste0("`metric` '", metric, "' not found in colData(spe)."))
  }
  if (!is.null(outliers) && !outliers %in% colnames(as.data.frame(colData(spe)))) {
    stop(paste0("`outliers` '", outliers, "' not found in colData(spe)."))
  }
  if (!is.numeric(point_size) || point_size <= 0) {
    stop("`point_size` must be a positive numeric value.")
  }
  if (!is.numeric(stroke) || stroke < 0) {
    stop("`stroke` must be a non-negative numeric value.")
  }
  if (!is.vector(colors) || !is.character(colors) || length(colors) < 2) {
    stop("`colors` must be a character vector of length >= 2.")
  }

  spe_sub <- spe[, colData(spe)[[sample_id]] == sample]
  p <- escheR::make_escheR(spe_sub) |>
    add_fill(var = metric, point_size = point_size)

  df <- as.data.frame(colData(spe_sub))
  coords <- spatialCoords(spe_sub)

  if (ncol(coords) < 2) {
    stop("`spatialCoords(spe_sub)` must return at least two columns (x and y).")
  }

  df$x <- coords[, 1]
  df$y <- coords[, 2]
  df$metric_val <- df[[metric]]

  r <- imgRaster(spe_sub)
  img_mat <- as.raster(r)

  xmin <- min(df$x)
  xmax <- max(df$x)
  ymin <- min(df$y)
  ymax <- max(df$y)

  p <- p +
    annotation_raster(img_mat, xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax) +
    coord_equal()

  p <- p +
    geom_point(data = df, aes(x = x, y = y, fill = metric_val),
               shape = 21, size = point_size, stroke = 0, color = NA)

  if (!is.null(outliers)) {
    out_vec <- df[[outliers]]
    if (is.factor(out_vec) || is.character(out_vec)) {
      out_factor <- if (is.factor(out_vec)) out_vec else factor(out_vec)
      idx <- which(!is.na(out_factor))
      if (length(idx) > 0) {
        df_out <- df[idx, , drop = FALSE]
        df_out$outlier_group <- droplevels(out_factor[idx])
        p <- p +
          geom_point(data = df_out,
                     aes(x = x, y = y, fill = metric_val, color = outlier_group),
                     shape = 21, size = point_size, stroke = stroke) +
          scale_color_discrete(name = outliers) +
          guides(color = guide_legend(override.aes = list(fill = NA)))
      }
    } else if (is.logical(out_vec) || all(out_vec %in% c(0, 1, NA))) {
      idx <- which(isTRUE(out_vec) | (out_vec %in% TRUE))
      if (length(idx) > 0) {
        df_out <- df[idx, , drop = FALSE]
        p <- p +
          geom_point(data = df_out,
                     aes(x = x, y = y),
                     shape = 21, size = point_size, stroke = stroke, color = "red")
      }
    }
  }

  if (length(colors) == 2) {
    p <- p + scale_fill_gradient(low = colors[1], high = colors[2], name = metric)
  } else {
    p <- p + scale_fill_gradientn(colors = colors, name = metric)
  }

  return(p)
}


#' Compare two deconvolution methods by plotting Venn diagram and spatial QC
#'
#' @param m1 A character vector of barcodes from barcode-based deconvolution (result of barcode_based_deconvolve)
#' @param m2 A character vector of barcodes from coordinate-based deconvolution (result of coordinate_based_deconvolve)
#' @param spe_high_res A SpatialExperiment object for high-resolution data (used for plotting)
#' @return No return value; prints a Venn diagram and a spatial QC plot
#' @importFrom ggvenn ggvenn
#' @importFrom ggplot2 ggtitle
#' @export
compare_deconvolved_methods <- function(m1, m2, spe_high_res) {

  if (missing(m1) || !is.vector(m1) || !is.character(m1)) {
    stop("`m1` must be a character vector of barcodes.")
  }
  if (missing(m2) || !is.vector(m2) || !is.character(m2)) {
    stop("`m2` must be a character vector of barcodes.")
  }
  if (missing(spe_high_res) || !("colData" %in% methods::slotNames(spe_high_res))) {
    stop("`spe_high_res` must be a SpatialExperiment-like object with colData and spatialCoords methods.")
  }
  if (!"barcode" %in% colnames(as.data.frame(colData(spe_high_res)))) {
    stop("`spe_high_res` colData must contain a 'barcode' column.")
  }
  p_venn <- ggvenn(list(barcode_based  = m1,
                        coordinate_based = m2))

  s1 <- setdiff(m2, m1)
  s2 <- setdiff(m1, m2)
  s3 <- intersect(m1, m2)

  spe_high_res$deconv_res <- NA
  spe_high_res$deconv_res[which(colData(spe_high_res)$barcode %in% s1)] <- "barcode-based excluded"
  spe_high_res$deconv_res[which(colData(spe_high_res)$barcode %in% s2)] <- "coordinate-based excluded"
  spe_high_res$deconv_res[which(colData(spe_high_res)$barcode %in% s3)] <- "both"

  p_QC <-
    plot_spatial_qc(spe_high_res,
                    metric = "sum",
                    outliers = "deconv_res",
                    point_size = 1,
                    stroke = 0.8) +
    ggtitle("Compare deconvolved results from both methods")

  print(p_venn)
  print(p_QC)
  return(spe_high_res)
}



#' Estimate agreement scores between barcode-based and coordinate-based methods
#'
#' @param df A data.frame that must contain a column named 'barcodes' (result of barcodes_to_pseudo_coordinates)
#' @param parquet_low_res Path to the low-resolution .parquet file
#' @param parquet_high_res Path to the high-resolution .parquet file
#' @param json_low_res Path to the low-resolution scale factors .json file
#' @param low_res Integer value indicating low-resolution bin size
#' @param high_res Integer value indicating high-resolution bin size
#' @return Returns the input data.frame with an added 'scores' column;
#'   also prints a density plot of scores and scatter plots for lowest scores
#' @importFrom ggplot2 ggplot aes geom_density theme_minimal labs geom_point scale_color_manual
#' @importFrom dplyr select slice_min
#' @export
estimate_deconvolved_scores <- function(df,
                                        parquet_low_res,
                                        parquet_high_res,
                                        json_low_res,
                                        low_res,
                                        high_res) {

  if (missing(df) || !is.data.frame(df)) {
    stop("`df` must be a data.frame containing at least a 'barcodes' column.")
  }
  if (!"barcodes" %in% colnames(df)) {
    stop("`df` must contain a column named 'barcodes'.")
  }
  if (!is.numeric(low_res) || !is.numeric(high_res) ||
      length(low_res) != 1 || length(high_res) != 1 ||
      low_res != floor(low_res) || high_res != floor(high_res) ||
      low_res <= 0 || high_res <= 0) {
    stop("`low_res` and `high_res` must be positive integer values.")
  }
  if (low_res < high_res) {
    stop("`low_res` must be larger than `high_res`.")
  }
  if ((low_res %% high_res) != 0) {
    stop("`low_res` must be a multiple of `high_res`.")
  }
  if (!is.character(parquet_low_res) || !file.exists(parquet_low_res)) {
    stop("`parquet_low_res` must be a valid path to an existing .parquet file.")
  }
  if (!is.character(parquet_high_res) || !file.exists(parquet_high_res)) {
    stop("`parquet_high_res` must be a valid path to an existing .parquet file.")
  }
  if (!is.character(json_low_res) || !file.exists(json_low_res)) {
    stop("`json_low_res` must be a valid path to an existing .json file.")
  }

  df$scores <- NA_real_
  barcodes <- df[["barcodes"]]

  for (br in barcodes) {

    if (!is.character(br) || is.na(br)) {
      warning("One of the barcodes in `df` is not a valid non-NA string; setting its score to NA.")
      next
    }

    m1 <- barcode_based_deconvolve(br, low_res = low_res, high_res = high_res)
    m2 <- coordinate_based_deconvolve(br,
                                      parquet_low_res = parquet_low_res,
                                      parquet_high_res = parquet_high_res,
                                      json_low_res = json_low_res,
                                      low_res = low_res,
                                      high_res = high_res)

    score <- NA_real_
    if (length(m1) > 0 || length(m2) > 0) {
      score <- length(intersect(m1, m2)) / length(union(m1, m2))
    }
    df$scores[df$barcodes == br] <- score
  }

  p_density <- ggplot(df, aes(x = scores)) +
    geom_density(fill = "lightblue", alpha = 0.6, size = 1) +
    theme_minimal(base_size = 14) +
    labs(title = "Distribution",
         x     = "Scores",
         y     = "Density")

  low_score_df <- df %>% dplyr::select(barcodes, scores) %>% dplyr::slice_min(order_by = scores, n = 10)

  ps <- vector("list", length = nrow(low_score_df))
  for (i in seq_len(nrow(low_score_df))) {
    bc <- low_score_df$barcodes[i]
    sc <- low_score_df$scores[i]
    m1 <- barcode_based_deconvolve(bc, low_res = low_res, high_res = high_res)
    m2 <- coordinate_based_deconvolve(bc,
                                      parquet_low_res = parquet_low_res,
                                      parquet_high_res = parquet_high_res,
                                      json_low_res = json_low_res,
                                      low_res = low_res,
                                      high_res = high_res)
    m1_df <- barcodes_to_pseudo_coordinates(m1)
    if (nrow(m1_df) > 0) {
      m1_df$labels <- "m1"
      m1_df$rows   <- as.integer(m1_df$rows)
      m1_df$cols   <- as.integer(m1_df$cols)
    }
    m2_df <- barcodes_to_pseudo_coordinates(m2)
    if (nrow(m2_df) > 0) {
      m2_df$labels <- "m2"
      m2_df$rows   <- as.integer(m2_df$rows)
      m2_df$cols   <- as.integer(m2_df$cols)
      title_text <- paste0(bc, "  (score=", round(sc, 3), ")")
    } else {
      title_text <- paste0(bc, "  (score=", round(sc, 3), ") \n Coordinate-Based method found 0 spots.")
    }
    df_all <- rbind(m1_df, m2_df)
    p_scatter <- ggplot(df_all, aes(x = cols, y = rows, color = labels)) +
      geom_point(alpha = 0.6, size = 2) +
      theme_minimal() +
      labs(title = title_text,
           x     = "Column",
           y     = "Row") +
      scale_color_manual(values = c("m1" = "tomato", "m2" = "steelblue"))
    ps[[i]] <- p_scatter
  }

  print(p_density)
  for (p in ps) {
    print(p)
  }
  return(df)
}
