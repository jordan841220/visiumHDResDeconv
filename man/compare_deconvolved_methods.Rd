% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func.R
\name{compare_deconvolved_methods}
\alias{compare_deconvolved_methods}
\title{Compare two deconvolution methods by plotting Venn diagram and spatial QC}
\usage{
compare_deconvolved_methods(m1, m2, spe_high_res)
}
\arguments{
\item{m1}{A character vector of barcodes from barcode-based deconvolution (result of barcode_based_deconvolve)}

\item{m2}{A character vector of barcodes from coordinate-based deconvolution (result of coordinate_based_deconvolve)}

\item{spe_high_res}{A SpatialExperiment object for high-resolution data (used for plotting)}
}
\value{
No return value; prints a Venn diagram and a spatial QC plot
}
\description{
Compare two deconvolution methods by plotting Venn diagram and spatial QC
}
