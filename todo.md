# TODO: Make tsnetworks CRAN-Ready

This document outlines the remaining tasks to prepare the `tsnetworks` R package for CRAN submission, addressing all warnings and errors from `R CMD check`.

## 1. Resolve Roxygen and Documentation Issues

*   **Fix `detect_regime.Rd` formatting:**
    *   Review `R/regime_detection.R` to ensure all `@param`, `@return`, `@details`, and `@examples` sections are correctly formatted, especially `\itemize` and `\describe` environments.
    *   Ensure `detect_regime` has a proper `@description` tag.
*   **Document `saqrsteps` dataset:**
    *   Create a new Roxygen block for `data/saqrsteps.rda` in `R/data.R` (or a new file like `R/saqrsteps.R`) with `@docType data`, `@name saqrsteps`, `@usage data(saqrsteps)`, `@format`, and `@description`.
*   **Address undocumented arguments in Rd `\usage` sections:**
    *   Review the Roxygen documentation for the following functions and ensure all parameters are explicitly documented:
        *   `create_id_network_output`
        *   `create_network_output`
        *   `create_windowed_network_output`
        *   `detect_regime` (re-check all params)
        *   `extract_distance_value`
        *   `extract_ts_from_dataframe`
        *   `extract_ts_vector`
        *   `perform_id_partition_analysis`
        *   `perform_windowed_analysis`

## 2. Manage Dependencies and Imports

*   **Explicitly declare all imports:**
    *   For `dtwclust` and `mclust`: Add these to `Imports` in `DESCRIPTION` if they are hard dependencies, or `Suggests` if optional. Ensure all calls use `::` (e.g., `dtwclust::tsclust()`, `mclust::Mclust()`).
    *   For `dplyr` and `ggplot2`: Replace `library()` or `require()` calls with `::` (e.g., `dplyr::mutate()`, `ggplot2::ggplot()`) and ensure they are listed in `Imports` in `DESCRIPTION`.
    *   For `grDevices`: If `grDevices::rainbow()` is used, ensure `grDevices` is in `Imports` and the call is explicit. If not used, remove from `Imports`.
    *   For `aricode::VI`: Verify if `aricode` is correctly listed in `DESCRIPTION` and if `VI` is an exported function from `aricode`. If `VI` is not exported, find an alternative or adjust usage.
    *   Review all "no visible global function definition/binding" notes and add appropriate `@importFrom` tags to the Roxygen blocks of the functions where these occur. This includes functions like `.apply_min_change_constraint`, `.generate_regime_ids`, `.detect_complexity_column_v2`, `.consolidate_similar_regimes_v2`, `.describe_complexity_patterns_v2`, `.detect_changepoints_v2`, `.detect_threshold_regimes_v2`, `.detect_variance_shifts_v2`, `.detect_gradient_changes_v2`, `.detect_smart_combination_v2`, `.discretize_series_pair`, `.ensure_same_length`, `.generate_windows`, `.jitter_spike_train`, `.validate_input`, `aes`, `arrange`, `check_required_packages`, `cor`, `cutree`, `element_blank`, `element_text`, `facet_wrap`, `first`, `geom_line`, `geom_point`, `geom_rect`, `geom_smooth`, `ggplot`, `group_by`, `group_id`, `hclust`, `hist`, `kmeans`, `labs`, `lag`, `last`, `lead`, `mad`, `median`, `mutate`, `na.omit`, `na.pass`, `qnorm`, `quantile`, `rainbow`, `row_number`, `runif`, `scale_fill_manual`, `segment_lag`, `segment_lead`, `state_change`, `summarise`, `sym`, `theme`, `theme_minimal`, `ungroup`, `vars`, `xmax`, `xmin`, `ymax`, `ymin`.

## 3. Clean Up and Final Checks

*   **Create `.Rbuildignore`:**
    *   Add `.DS_Store` and `..Rcheck` to `.Rbuildignore` to prevent them from being included in the package build.
*   **Re-check non-ASCII characters:**
    *   Carefully inspect `R/network_builder.R`, `R/similarity_utils.R`, and `R/ts_network.R` for any remaining non-ASCII characters and replace them with ASCII equivalents or `\uxxxx` escapes.
*   **Verify `LICENSE` file:**
    *   Ensure a `LICENSE` file exists in the package root and its content matches the `MIT` license. The `DESCRIPTION` file already specifies `License: MIT`.
*   **Vignette configuration:**
    *   Ensure the vignettes are correctly configured to build. The `DESCRIPTION` already has `VignetteBuilder: knitr`. I need to ensure the Rmd files themselves are valid and can be knitted.
*   **Run `roxygen2::roxygenise()`:** After all Roxygen comments are updated, run this command again to regenerate `NAMESPACE` and all `.Rd` files.
*   **Run `R CMD check`:** Perform a final `R CMD check .` to verify all issues are resolved.
