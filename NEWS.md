# RTCGAToolbox 2.2.57

## New Features

* Migration from RTCGAToolbox - moving to a new package

## Bug fixes

* Re-added the arraydata = FALSE flag for RNAseq data type

# RTCGAToolbox 2.2.56

## New features

* `readExonFiles` function now available for TCGA exon array level data. Returns
a `GRangesList`.

## Bug fixes and minor improvements 

* Import `content_type` function from the `httr` package

# RTCGAToolbox 2.2.55

## Bug fixes and minor improvements 

* `extract` now returns either `GRangesList` or `ExpressionSet` without
clinical data

# RTCGAToolbox 2.2.53

## Bug fixes and minor improvements

* Reverted to using dashes "-" instead of periods "." as barcode separators

# RTCGAToolbox 2.2.52

## New features

* News feed!

* `matchClinical` provides functionality to match your experiment and clinical
    datasets by using TCGA identifiers

* `makeGRangesList` allows users to create a `GRangesList` object from raw data
    in either `list` or `data.frame` form.

## Bug fixes and minor improvements

* Invalid name creation now fixed (replaces dashes "-" with periods ".").
* hugo symbol values now showing in the created `GRangesList` objects.
* Warning message from `data.table` import now taken care of.
* `NAMESPACE` now created using Roxygen2