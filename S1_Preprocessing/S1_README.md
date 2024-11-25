# S1_Preprocessing

This folder handles the preprocessing of raw data, including cleaning, filtering, and extracting relevant variables for downstream analysis. It ensures all data is standardized and ready for further steps in the pipeline.

## Folder Structure

### 1. `Functions`

Contains custom auxiliary functions used across the preprocessing scripts.

-   `auxiliary_functions.R`: Core helper functions to simplify and modularize preprocessing tasks.

### 2. `Miscellaneous`

A collection of reference materials, metadata, and environmental data.

-   **`Human_footprint/`**: Contains raster files for the global human footprint dataset (e.g., `HFP2009.tif`) used to assess human impact on ecosystems.
-   **`wc10/`**: Contains raster files for WorldClim bioclimatic variables (`bio1`, `bio2`, etc.) used for spatial and environmental analyses.
    -   `dataset_info_all.xlsx` : Descriptive summaries of datasets.

### 3. `Processed`

Contains cleaned datasets, output as `.rds` files, which are ready for analysis in subsequent steps.

### 4. `Raw_data`

Stores raw input datasets in `.xlsx` format. These files are cleaned and transformed during preprocessing. Examples include:

-   Data for biodiversity studies, e.g., `N21FBD_2.xlsx`, `Predicts_cc1_2013_waite_1.xlsx`.
-   Taxon-specific datasets (birds, insects, plants, etc.).

### 5. `Scripts`

Preprocessing scripts used to clean and organize the raw data:

-   `generate_dataset_map.R`: Creates spatial representations of datasets.
-   `pre_processing.R`: General preprocessing workflow for all raw datasets.

### 6. `run_step1.sh`

A batch file used to run the entire preprocessing code on the university's HPC (High-Performance Computing) system. This step is not mandatory but is included here for reference and reproducibility.

## Notes

-   Rasters in the `Human_footprint/` and `wc10/` folders provide spatial environmental variables and indices; no need to manipulate their subfiles directly.
-   The cleaned data in the `Processed` folder is essential for subsequent analyses in the pipeline.
