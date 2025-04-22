## CytoSimplex 0.2.0

- Added support for coloring dots by categorical or continuous variable, with options to use customized or built-in color palettes.
- Changed interactive 3D support of quaternary simplex plot from *rgl* to *plotly*. Now the default is `interactive = TRUE` for `plotQuaternary()` and it returns a `plotly` object.
- Added interactive support for ternary simplex plot, using plotly, triggered with `interactive = TRUE`
- Fixed wilcoxon bug
- Added `readH5ADObsNames()`, `readH5ADObsVar()`, `readH5ADUnsSpMat()`, `readVelocytoLoom()` for loading commonly seen necessary information from H5AD and LOOM files of RNA velocity analysis
- Added example tutorial for HSPC analysis in the manuscript
- Fixed error in `plotQuaternary()`

## CytoSimplex 0.1.1

- Fixed vignette GIF picture linking issue

## CytoSimplex 0.1.0

- Initial release version.
