# CRAN Readiness Report - tsnetworks

This report summarizes the current status of the `tsnetworks` R package regarding CRAN readiness, outlining completed tasks and areas requiring further revision.

## Done:

*   `.Rbuildignore` updated: Added `.DS_Store` and `..Rcheck` to `.Rbuildignore` to prevent their inclusion in the package build.
*   `DESCRIPTION` file consolidated: Duplicate `Imports` fields in the `DESCRIPTION` file have been merged into a single, correct `Imports` field.
*   `dtwclust` dependency addressed: The `dtwclust` package, a suggested dependency, has been installed to resolve a check error.
*   Internal function calls corrected: Removed incorrect `tsn::` prefixes from internal function calls (e.g., `prepare_for_qgraph` in `R/network_builder.R`).
*   Missing helper functions defined:
    *   `check_required_packages` function defined in `R/validation.R`.
    *   `.ensure_same_length` function defined in `R/validation.R`.
    *   `.jitter_spike_train` function defined in `R/distance_functions.R`.
    *   `.discretize_series_pair` function defined in `R/discretization.R`.
*   `importFrom` directives and `globalVariables` added: Added necessary `@importFrom` directives and `utils::globalVariables` calls to the Roxygen blocks of various R files (`R/regime_detection.R`, `R/plotting.R`, `R/calculator.R`, `R/distance_functions.R`, `R/ts_utils.R`, `R/rolling_measures.R`, `R/validation.R`, `R/similarity_utils.R`, `R/network_builder.R`, `R/stna.R`, `R/discretization.R`, `R/output_formatting.R`, `R/ranking_utils.R`, `R/transition_utils.R`, `R/visibility_graph.R`) to resolve "no visible global function definition" and "no visible binding for global variable" notes.
*   `saqrsteps` dataset documented: A basic documentation file (`R/data.R`) has been created for the `saqrsteps` dataset.

## Needs Revision:

*   `pdflatex` error (External Issue): The `R CMD check` consistently reports an error related to `pdflatex` not being available. This is an external system dependency and cannot be resolved directly by modifying the R package code. Action Required: The user needs to ensure `pdflatex` is installed and correctly configured on their system for the R package check to complete successfully.
*   Rd Syntax Errors in `man/detect_regime.Rd`: Despite efforts to correct the Roxygen comments in `R/regime_detection.R`, the generated `man/detect_regime.Rd` file still shows "unknown macro '\item'" and "unexpected section header" warnings. This indicates persistent issues with how the Roxygen comments are interpreted.
    *   Action Required: Carefully review the Roxygen comments for the `detect_regime` function in `R/regime_detection.R`, paying close attention to the formatting of `@param` and `@return` sections, especially the use of `\item` within `\itemize` or `\describe` environments. Ensure all braces are correctly matched and that the structure adheres strictly to Roxygen/Rd syntax.
*   Remaining "no visible global function definition" and "Undefined global functions or variables" notes: Although many have been addressed, the `R CMD check` output still lists several of these.
    *   Action Required: Systematically go through the latest `R CMD check` output. For each remaining "no visible" note, identify the function or variable and the file where it is used. Add the appropriate `@importFrom` directive to the Roxygen block of the calling function's file, or use `utils::globalVariables()` for global variables.
*   Undocumented Arguments in Rd Files: Several functions still have undocumented arguments, leading to warnings.
    *   Action Required: For each function listed in the `R CMD check` output under "Undocumented arguments in Rd file", add a `@param` tag with a clear description for every argument in its Roxygen block.
*   Vignette Issues: Warnings persist regarding vignettes not being built correctly and the `VignetteBuilder` field.
    *   Action Required:
        *   Verify that `VignetteBuilder: knitr` is correctly specified in the `DESCRIPTION` file.
        *   Ensure the `.Rmd` files in the `vignettes/` directory have the correct YAML header for R Markdown vignettes (e.g., `VignetteIndexEntry`, `\VignetteEngine{knitr::rmarkdown}`).
        *   Consider running `devtools::build_vignettes()` locally to debug vignette compilation.
*   "License stub is invalid DCF" in `DESCRIPTION`: This is a minor but important `NOTE` for CRAN submission.
    *   Action Required: Review the `License` field in the `DESCRIPTION` file. Ensure it strictly follows the format required by CRAN (e.g., `MIT + file LICENSE`). If a custom license is used, the `LICENSE` file must be present and correctly formatted.

---
**Next Steps:** Address the "Needs Revision" items one by one, starting with the Rd syntax errors and remaining "no visible" issues. After each set of changes, run `R CMD check .` again to verify progress.
