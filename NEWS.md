## CHANGES IN VERSION 2.16.0

### New features

* `RNASeq2GeneNorm` slot in the `FirehoseData` class is a `list` now (from
`matrix`)
* Use `tempdir()` as the default directory for downloading data in
`getFirehoseData`


### Bug fixes and minor improvements

* Save all `RNASeq2GeneNorm` datasets within the output object as a list.
Previously, only the last dataset would get returned (#30)
* Read files from the appropriate download location in `getFirehoseData`
* Move static text file references from `canevolve` to GitHub hosted locations

## CHANGES IN VERSION 2.14.0

### Bug fixes and minor improvements

* GISTIC data for `SKCM` and `LAML` is correctly returned by the
`getFirehoseData` function
* Helpers now correctly assign row names in presence of ranged data when 
using `biocExtract`

## CHANGES IN VERSION 2.12.0

### New features

* `biocExtract` can now return GISTIC peak data as well with the 'GISTICP'
option
* Methylation data is now converted to `DelayedArray` and contained in
`SummarizedExperiment` objects via `biocExtract`
* Added an `isEmpty` method to the `FirehoseGISTIC` data class

### Bug fixes and minor improvements

* Set default 'peak' argument to 'wide'
* Functionality from the `TCGAutils` package is now employed such as
`findGRangesCols` and `uniformBuilds` 

## CHANGES IN VERSION 2.9.41

### Bug fixes and minor improvements

* Fixed bug where file rename was not working due to incorrect path location
in `getFirehoseData`

## CHANGES IN VERSION 2.9.40

### New features

* `getBroadSubtypes` and `getGISTICPeaks` now included in the package.
Thanks to @lgeistlinger!
* Updated import directives
* Added peaks data slot to the `FirehoseGISTIC` class
* Numerous improvements to internal helpers
* Added a `GISTIC` argument to `getFirehoseData`

## CHANGES IN VERSION 2.7.44

### New features

* `selectType` function replaces all previous data type extractor functions
* `getData` improved, removed `CN` argument
* Cleaner documentation for `biocExtract`

### Bug fixes and minor improvements

* Included helper functions for enabling full functionality of `biocExtract`
* Added appropriate import directives
* Re-organized documentation for readability

## CHANGES IN VERSION 2.7.21

### New features

* `biocExtract` function now available - convert data objects to Bioconductor
friendly data classes

## CHANGES IN VERSION 2.7.21

### New features

* New maintainer for package
* Updated API - see `FirehoseData` slots for changes
* New extractor functions for slots of `FirehoseData`
* Standardize argument names and slots
* Export functions in package
* Updated NAMESPACE

### Bug fixes and minor improvements

* Cleaner documentation
* Import directives appropriate
* Example dataset updated to new API
