% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func.R
\name{coordinate_based_deconvolve}
\alias{coordinate_based_deconvolve}
\title{Identify high-resolution barcodes to discard based on spatial coordinates}
\usage{
coordinate_based_deconvolve(
  barcode_list,
  parquet_low_res,
  parquet_high_res,
  json_low_res,
  low_res,
  high_res
)
}
\arguments{
\item{barcode_list}{A character vector of low-resolution barcodes}

\item{parquet_low_res}{Path to the low-resolution .parquet file}

\item{parquet_high_res}{Path to the high-resolution .parquet file}

\item{json_low_res}{Path to the low-resolution scale factors .json file}

\item{low_res}{Integer value indicating low-resolution bin size (e.g., 64)}

\item{high_res}{Integer value indicating high-resolution bin size (e.g., 8)}
}
\value{
A character vector of high-resolution barcodes to discard
}
\description{
Identify high-resolution barcodes to discard based on spatial coordinates
}
