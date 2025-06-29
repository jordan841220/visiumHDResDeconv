# visiumHDResDeconv

## PURPOSE
This package is designed for users who wish to **perform quality control at lower resolution and then “map” QC results back to high resolution**, because applying standard QC directly on high-resolution Visium HD data can be problematic: each high-res spot often has very few UMIs, leading to extreme local variance. 

<br>

## METHODS
It provides two complementary methods for identifying high-resolution spots to discard based on low-resolution QC:
- barcode-based deconvolution: uses each low-resolution barcode to generate all corresponding high-resolution barcodes that cover the same physical area.
- coordinate-based deconvolution: relies on spatial coordinates and pixel dimensions (from JSON and Parquet files) to find every high-resolution spot that falls within each low-resolution spot’s bounding box. 

<br>

In addition, the package offers visualization helpers: a modified version of SpotSweeper’s QC-plot routine (many thanks to the SpotSweeper authors) for overlaying spot‐level QC metrics on the tissue image, and a comparison function that draws Venn diagrams and spatial QC overlays to illustrate differences between the two deconvolution approaches. Finally, there is a convenience function to compute agreement scores between the barcode-based and coordinate-based methods, along with summary plots highlighting the lowest-agreement cases.

<br>

## NOTES
However, in practice each method alone may fail to identify all problematic high-resolution spots. Therefore, we recommend taking the **union** of the barcodes identified by both approaches - these barcodes represent the spots to consider filtering in the high-resolution dataset. Furthermore, since the structure of output and format matters, we recommend using the output of spaceranger v3.1.3.

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

# INSTALLATION
To install the development version of visiumHDResDeconv from GitHub, run the following commands in R:
```
install.packages("devtools")
devtools::install_github("jordan841220/visiumHDResDeconv")
```

# INPUTS
- 10X Genomics Visium HD data (SpatialExperiment object)
- tissue_positions.parquet
- scalefactors_json.json


# EXAMPLE
Following workflow demonstrates how I perform usual QC pipeline on low resolution data, and then deconvolve "bad spots" to higher resolution data.

### Library
```
library(SpotSweeper)
library(escheR)
library(jsonlite)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
options(future.globals.maxSize = 1e9)
```

### Data loading
```
# Load your 64um SpatialExperiment object "spe_64um"
# identifying the mitochondrial transcripts
is.mito <- rownames(spe_64um)[grepl("^mt-|^MT-", rownames(spe_64um))]
# calculating QC metrics for each spot using scuttle
spe_64um <- scuttle::addPerCellQCMetrics(spe_64um, subsets = list(Mito = is.mito))
keep <- colData(spe_64um)$sum > 0
spe_64um <- spe_64um[, keep]

# Load your 8um SpatialExperiment object "spe_8um"
# ...same as above
```
### Easy preprocess
You can apply SpotSweeper QC pipeline here. For example, if I want to exclude "poor spots" in 64um which qc_lib_size are < 50
```
qc_lib_size <- colData(spe_64um)$sum < 50
# add QC results to colData
spe_64um$qc_lib_size <- qc_lib_size
spe_64um$global_outliers <- as.logical(qc_lib_size) 
```

### Transfer barcodes to "pseudo coordination""
```
bad_spots <- barcodes_to_pseudo_coordinates(colData(spe_64um)$barcode[spe_64um$global_outliers])
```
barcodes_to_pseudo_coordinates() function will output a dataframe with 3 columns as below. Notes that the format of barcodes have to be similar to "s_064um_00011_00001-1".

![High‐resolution spot QC example](man/figures/barcodes_to_pseudo_coordinates_example.png)


### Deconvolve
```
m1 <- barcode_based_deconvolve(bad_spots$barcodes, low_res = 64, high_res = 8)
m2 <- coordinate_based_deconvolve(bad_spots$barcodes, 
                                 parquet_low_res = "/path/to/64um/tissue_positions.parquet",
                                 parquet_high_res = "/path/to/8um/tissue_positions.parquet",
                                 json_low_res = "/path/to/64um/scalefactors_json.json",
                                 low_res = 64, 
                                 high_res = 8)
```
Here m1/m2 are 2 list of barcodes.


### Estimate the effect of deconvolve
- The INPUT is a dataframe with "barcodes" column
```
bad_spots <- estimate_deconvolved_scores(bad_spots, 
                                 parquet_low_res = "/path/to/64um/tissue_positions.parquet",
                                 parquet_high_res = "/path/to/8um/tissue_positions.parquet",
                                 json_low_res = "/path/to/64um/scalefactors_json.json",
                                 low_res = 64, 
                                 high_res = 8)
```
The estimate_deconvolved_scores() function calculates a score between 0 and 1 for each low‐resolution barcode—measuring agreement between the two deconvolution methods. It then generates a density plot of all scores and displays the spatial locations of the 10 lowest‐scoring (worst) deconvolution events.

![High‐resolution spot QC example2](man/figures/estimate_deconvolved_scores_example1.png)

![High‐resolution spot QC example3](man/figures/estimate_deconvolved_scores_example2.png)

The resulting distribution plot lets you assess overall agreement between the two methods: a score of 1 indicates perfect concordance, while 0 indicates no overlap at all.



### Compare quality between two deconvolved methods
```
spe_8um <- compare_deconvolved_methods(m1 = m1, m2 = m2, spe_high_res = spe_8um)
```
![High‐resolution spot QC example4](man/figures/compare_deconvolved_methods_example1.png)

The compare_deconvolved_methods() function displays which barcodes each method identified (using a Venn diagram) and overlays those spots on the tissue image to visualize where each set of barcodes actually lies.
<br>

The "poor spots" highlighted by SpotSweeper on 64um resolution:
![High‐resolution spot QC example5new](man/figures/compare_deconvolved_methods_example2new.png)

After we deconvolve them to 8um resolution:
![High‐resolution spot QC example6](man/figures/compare_deconvolved_methods_example3.png)

The returned spe_8um will have a new colData named "deconv_res" which labels whether certain "poor barcode" is deconvolved by method 1, 2 or both (good barcodes will be labeled as NA).


# AUTHOR
Chun‐Ting Lin jordan841220@gmail.com

# CONTACT
For issues or questions, please open an issue on GitHub at https://github.com/jordan841220/visiumHDResDeconv or email jordan841220@gmail.com


# LICENSE
A software license specifies how others may use, modify, and distribute your code. visiumResDeconv is released under the MIT License. This means anyone can freely use, copy, modify, merge, and distribute the software, provided that the original copyright and license notice appear in all copies. The full MIT License text can be found in the LICENSE file in this repository.

