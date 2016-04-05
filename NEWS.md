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