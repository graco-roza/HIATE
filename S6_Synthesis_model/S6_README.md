# S6_Synthesis_Model

This folder contains data, functions, and scripts used to synthesize and model the relationships between species and traits replacements. It adheres to the **FAIR principles** (Findable, Accessible, Interoperable, and Reusable) by providing structured metadata for datasets and detailed documentation of functions and scripts.

## Folder Structure

``` plaintext
S6_Synthesis_model
├── data
│   ├── BBGDM_SES.xlsx
│   ├── complete_data.xlsx
│   └── synthesis_data.xlsx
├── functions
│   ├── helper_functions_S6.R
│   └── helper_functions_plot_extended.R
└── scripts
    ├── 1-compile_synthesis_data.R
    ├── 2-model_direction.R
    ├── 3-model_magnitude.R
    ├── 4-model_shape.R
    └── 5-model_convergence.R
```

------------------------------------------------------------------------

## Data Folder (`data/`)

### 1. `BBGDM_SES.xlsx`

# BBGDM_SES.xlsx

**Description**:\
Contains standardized effect size (SES) results from null models for species and trait replacement analyses. This dataset provides SES values and associated p-values for identifying patterns of differentiation, homogenization, or randomness, as well as the shapes of replacement patterns along environmental or human pressure gradients.

## Columns

-   **`dataset`**:\
    Unique identifier for each dataset. Each dataset is represented twice (once for differentiation and once for homogenization).

-   **`direction`**:\
    Indicates the type of ecological process:

    -   **Differentiation**: Processes increasing dissimilarity or diversity among communities.\
    -   **Homogenization**: Processes reducing dissimilarity or diversity among communities.

-   **`direction_SES`**:\
    Standardized effect size (SES) for the direction of replacement, quantifying deviation from null expectations.

-   **`direction_pvalue`**:\
    P-value associated with `direction_SES`, indicating statistical significance.

-   **`magnitude_SES`**:\
    SES for the magnitude of replacement, representing the strength of ecological processes.

-   **`magnitude_pvalue`**:\
    P-value associated with `magnitude_SES`, indicating statistical significance.

-   **`Absent_ses`**:\
    SES for the "Absent" replacement shape (no significant replacement pattern along the gradient).

-   **`Absent_pvalue`**:\
    P-value associated with `Absent_ses`.

-   **`Exponential_ses`**:\
    SES for the "Exponential" replacement shape (replacement increases exponentially along the gradient).

-   **`Exponential_pvalue`**:\
    P-value associated with `Exponential_ses`.

-   **`Saturating_ses`**:\
    SES for the "Saturating" replacement shape (replacement increases rapidly at first and then plateaus).

-   **`Saturating_pvalue`**:\
    P-value associated with `Saturating_ses`.

-   **`Revlog_ses`**:\
    SES for the "Reverse Logistic" replacement shape (replacement is strongest at low and high human pressure but plateaus in intermediate values).

-   **`Revlog_pvalue`**:\
    P-value associated with `Revlog_ses`.

## Usage

This file serves as an input for modeling scripts, especially in analyzing SES-based trends and divergences between traits and species replacement patterns along environmental gradients. It is particularly useful in testing hypotheses about dependency, convergence, and shape dynamics in ecological processes.

------------------------------------------------------------------------

### 2. `synthesis_data.xlsx`

# Metadata for `synthesis_data.xlsx`

## Description

This dataset compiles results and metadata used in the analysis of species and trait replacement across various datasets. The data includes biodiversity facets (Taxonomic and Functional), spatial gradients, and human pressure metrics. It is utilized for synthesizing results and testing hypotheses regarding biodiversity patterns along environmental gradients.

## Columns

1.  **dataset**: Unique identifier for each dataset, indicating the specific study or data source.

2.  **facet**: The biodiversity facet analyzed, either "Taxonomic" or "Functional".

3.  **direction**: Indicates the process direction, either "Differentiation" or "Homogenization".

4.  **direction_z**: Z-score for the direction process, representing standardized values for differentiation or homogenization.

5.  **predictor**: Indicates the gradient or variable driving the analysis, e.g., "hfp" for human footprint.

6.  **buffer**: Spatial buffer size used in the analysis (e.g., "5000" meters).

7.  **magnitude_z**: Z-score for the magnitude of biodiversity replacement.

8.  **Absent**: Number of cases where species or traits were absent along the gradient.

9.  **Saturating**: Number of cases with saturating relationships along the gradient, indicating diminishing replacement effects at higher values.

10. **Exponential**: Number of cases with exponential relationships, where replacement effects increase rapidly.

11. **Revlog**: Number of cases with reverse logistic replacement patterns, suggesting stronger replacement at low and high gradient values with a plateau in the middle.

12. **Hfp.min**: Minimum human pressure value within the gradient.

13. **Hfp.mean**: Mean human pressure value within the gradient.

14. **Hfp.max**: Maximum human pressure value within the gradient.

15. **Hfp.range**: Range of human pressure values, calculated as the difference between maximum and minimum.

16. **species.number**: Number of species recorded in the dataset.

17. **spatial.extent**: Spatial extent of the study area in square kilometers.

18. **spatial.min**: Minimum spatial unit size considered in the analysis.

19. **spatial.mean**: Mean spatial unit size considered in the analysis.

20. **spatial.max**: Maximum spatial unit size considered in the analysis.

21. **latitude.mean**: Mean absolute latitude of the study area.

22. **disturbance**: Main land-use disturbance type, e.g., "agriculture", "forest", or "multiple".

23. **taxa**: Taxonomic group analyzed in the dataset, e.g., "plant", "bird", or "invertebrate".

24. **biotic.group**: Higher-level classification of taxa, e.g., "vertebrate" or "invertebrate".

25. **realm**: Ecosystem type or realm, e.g., "terrestrial" or "aquatic".

## Usage

This dataset is used in: -modelling the cross-dataset variation in direction, magnitude and shape of relationship for species and trait replacement

**Important Note**: Ensure proper referencing of predictor and buffer variables when applying statistical models to avoid confounding spatial and gradient effects.

------------------------------------------------------------------------

## Functions Folder (`functions/`)

### 1. `helper_functions_S6.R`

-   **Description**: Custom functions for data processing and model preparation.
-   **Key Functions**:
    -   `process_model_data`: Cleans and formats data for model input.
    -   `extract_model_ses_data`: Merges SES results with synthesis data for specific model types (e.g., direction, shape).

------------------------------------------------------------------------

### 2. `helper_functions_plot_extended.R`

-   **Description**: Custom plotting functions for creating publication-quality visualizations.
-   **Key Functions**:
    -   `plot_posterior_draws`: Visualizes posterior distributions of model results.
    -   `plot_shape_scheme`: Creates diagrams to represent possible relationship shapes.
    -   `generate_ses_figure_data`: Prepares data for plotting SES results.

------------------------------------------------------------------------

## Scripts Folder (`scripts/`)

### 1. `1-compile_synthesis_data.R`

-   **Description**: Prepares the `synthesis_data.xlsx`

------------------------------------------------------------------------

### 2. `2-model_direction.R`

-   **Description**: Fits Bayesian models to analyze the directional replacement processes (differentiation vs. homogenization).

------------------------------------------------------------------------

### 3. `3-model_magnitude.R`

-   **Description**: Fits Bayesian models to assess magnitude differences in replacement processes between traits and species across datasets

------------------------------------------------------------------------

### 4. `4-model_shape.R`

-   **Description**: Models and visualizes the shapes of replacement relationships (e.g., "Exponential", "Saturating").

------------------------------------------------------------------------

### 5. `5-model_convergence.R`

-   **Description**: Evaluates the dependency between traits and species using null model results to identify convergence/divergence patterns.

------------------------------------------------------------------------

## Outputs

-   **Model Results**: Saved in `S7_Model_outputs_figures_and_tables/model/`.
-   **Figures**: Exported as PDFs to `S7_Model_outputs_figures_and_tables/main_figures/` and `S7_Model_outputs_figures_and_tables/extended_data/`.
