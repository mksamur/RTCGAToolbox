## Changes in version 2.34.0

### Significant user-visible changes

* `getReport` has been removed from the package.

### Bug fixes and minor improvements

* Remove unnecessary imports e.g., from `RCircos`.

## Changes in version 2.32.0

### New features

* The functions `getCNGECorrelation`, `getDiffExpressedGenes`, and `getSurvival`
have been removed from the package.

## Changes in version 2.30.0

### New features

* Include additional options for `miRNASeqGeneType` and `RNAseq2Norm` inputs,
see `?getFirehoseData`.
* The functions `getCNGECorrelation`, `getDiffExpressedGenes`, and `getSurvival`
have been defunct and will be removed from the package (see
`?'RTCGAToolbox-defunct'`).
* The `RNAseq2Norm` argument in `getFirehoseData` allows additional
options: "RSEM_normalized_log2", "raw_counts", "scaled_estimate" from the
'preprocessed' tarballs in Firehose. The "normalized_counts" default remains.

### Bug fixes and minor improvements

* `getFirehoseData` when used with the `miRNASeqGene` argument was
downloading read counts data rather than RPM. This has been fixed with the
`miRNASeqGeneType` argument. "read_count" and "cross-mapped" data are still
available but must be entered explicitly in `miRNASeqGeneType`.

## Changes in version 2.28.0

### New features

* Resolve disparate columns in mutation files from `GBMLGG` (@andreagrioni, #45)
* Update codebase to download https resources from newer layout in
<gdac.broadinstitute.org> (@biostars-nyc, #44).
* Update `makeSummarizedExperimentFromGISTIC` interface with `rownameCol` input,
type checking, and set row names when the are not duplicated.

### Bug fixes and minor improvements

* Use cache directory obtained from `tools::R_user_dir`
* Clean up class membership checks, e.g., with `is(x, "classname")`
* Set missing rows to `""` for downstream compatibility with
`SummarizedExperiment`.

## Changes in version 2.24.0

### Bug fixes and minor improvements

* The deprecated functionality vignette moved to the `vignettes/analysis/`
folder

## Changes in version 2.22.0

### New features

* `makeSummarizedExperimentFromGISTIC` transferred from `TCGAutils`

### Bug fixes and minor improvements

* `biocExtract` now merges datasets from the same platforms
* Re-worked and simplified `getData` method for `FirehoseData` and
`FirehoseGISTIC`
* Remove missing ranges when creating `GRanges` and `GRangesList` from
`DataFrame`

## Changes in version 2.20.0

### New features

* Added the `RNASeq2Gene` slot to the `FirehoseData` class. This data type
mainly obtains RNASeq v2 `raw_counts` from the pipeline (`scaled_estimates`
also available; @mherberg #39)
* Added an `accmini` example dataset as obtained from `getFirehoseData`
* `getLinks` function shows the user some file provenance based on data
requested
* Newly deprecated functions: `getDiffExpressedGenes`, `getCNGECorrelation`,
`getSurvival` and `getReport`. It is no longer possible for the maintainer to
update these functions in a way that would benefit users. A transfer of
responsibility would be required, i.e. to another package.
* Vignettes are updated to reflect changes in the codebase.

### Bug fixes and minor improvements

* Improvements to internal functions for converting tabular data to
Bioconductor classes
* Missing (NA) `seqnames` are removed when converting to `RaggedExperiment`
* Remove static file dependencies from GitHub and use text inside function
(@DavisWeaver, #34)
* Added default values to helper for making `SummarizedExperiment` datasets
* Coerce sample names to character when in (the rare) case they're numeric
* Added an ellipsis argument to `biocExtract` for specifying the `names.field`
in tabular data that will correspond to the row names of a
`SummarizedExperiment`

## Changes in version 2.18.0

### New features

* Warning for Windows users added when file paths are too long
* `getGISTICPeaks` now requires a `FirehoseGISTIC` data object obtained from
`getFirehoseData`

### Bug fixes and minor improvements

* Consolidate GISTIC data download methods in `getFirehoseData` and
`getGISTICPeaks`
* Increase robustness of internal helper functions that work with strands
* 'TCGA' sample column identification is less strict

## Changes in version 2.16.0

### New features

* `RNASeq2GeneNorm` slot in the `FirehoseData` class is a `list` now (from
`matrix`)
* Use `tempdir()` as the default directory for downloading data in
`getFirehoseData`

### Bug fixes and minor improvements

* Save all `RNASeq2GeneNorm` datasets within the output object as a list.
Previously, only the last dataset would get returned (#30)
* Read files from the appropriate download location in `getFirehoseData`
* Move static text file references from 'canevolve.org' to GitHub hosted
locations
* Check file sizes using `httr` instead of 'canevolve.org' query
(@mksamur, #32)

## Changes in version 2.14.0

### Bug fixes and minor improvements

* GISTIC data for `SKCM` and `LAML` is correctly returned by the
`getFirehoseData` function
* Helpers now correctly assign row names in presence of ranged data when
using `biocExtract`

## Changes in version 2.12.0

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

## Changes in version 2.9.41

### Bug fixes and minor improvements

* Fixed bug where file rename was not working due to incorrect path location
in `getFirehoseData`

## Changes in version 2.9.40

### New features

* `getBroadSubtypes` and `getGISTICPeaks` now included in the package.
Thanks to @lgeistlinger!
* Updated import directives
* Added peaks data slot to the `FirehoseGISTIC` class
* Numerous improvements to internal helpers
* Added a `GISTIC` argument to `getFirehoseData`

## Changes in version 2.7.44

### New features

* `selectType` function replaces all previous data type extractor functions
* `getData` improved, removed `CN` argument
* Cleaner documentation for `biocExtract`

### Bug fixes and minor improvements

* Included helper functions for enabling full functionality of `biocExtract`
* Added appropriate import directives
* Re-organized documentation for readability

## Changes in version 2.7.21

### New features

* `biocExtract` function now available - convert data objects to Bioconductor
friendly data classes

## Changes in version 2.7.21

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
