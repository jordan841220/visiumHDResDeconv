# visiumHDResDeconv
This package focuses on **resolution‐level deconvolution** of Visium HD spatial transcriptomics data. It provides two main approaches for identifying high‐resolution spots to discard based on low‐resolution information. The first method (barcode‐based deconvolution) uses the original barcode IDs to generate a list of high‐resolution barcodes covering each low‐resolution spot. The second method (coordinate‐based deconvolution) relies on spatial coordinates and pixel dimensions to find all high‐resolution spots that fall within the same physical area as each low‐resolution spot. In addition, the package offers visualization helpers to overlay spot QC metrics and to compare results from both deconvolution methods using Venn diagrams and spatial QC plots. The plot_spatial_qc function is a modified version of the implementation in the SpotSweeper R package (special thanks to the SpotSweeper authors, https://github.com/MicTott/SpotSweeper).Finally, there is a convenience function to compute agreement scores between the two methods and to generate summary plots of the lowest‐agreement cases.

# DEPENDENCIES
The package imports the following R packages:
- dplyr
- tidyr
- stringr
- jsonlite
- arrow
- escheR
- ggplot2
- ggvenn
- SpotSweeper

# INSTALLATION
To install the development version of visiumHDResDeconv from GitHub, run the following commands in R:
```
install.packages("devtools")
devtools::install_github("jordan841220/visiumHDResDeconv")
```

# INPUTS
- parquet file output by spaceranger
- JSON file output by spaceranger


# EXAMPLE
(to be continued)


# AUTHOR
Chun‐Ting Lin jordan841220@gmail.com

# CONTACT
For issues or questions, please open an issue on GitHub at https://github.com/jordan841220/visiumHDResDeconv or email jordan841220@gmail.com


# LICENSE
A software license specifies how others may use, modify, and distribute your code. visiumResDeconv is released under the MIT License. This means anyone can freely use, copy, modify, merge, and distribute the software, provided that the original copyright and license notice appear in all copies. The full MIT License text can be found in the LICENSE file in this repository.

