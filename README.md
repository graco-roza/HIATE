# HIATE Project Repository

Welcome to the HIATE project repository, associated with the manuscript **"Human Pressure Homogenises Species and Traits Globally"** by Graco-Roza et al. [2025??].
This repository includes all the data, scripts, and documentation required to replicate the analyses and figures presented in the manuscript. 
The repository is designed to align with FAIR (Findable, Accessible, Interoperable, Reusable) principles, ensuring clear organization, metadata, and usability.
If any improvement can be done, please contact me at: **caio.roza@helsinki.fi**


## Repository Structure

The project is organized into a series of steps, each corresponding to a stage of the analysis pipeline. 
Each step contains its own README file and metadata for relevant datasets and scripts.

```plaintext
HIATE:
  - HIATES.Rproj
  - S1_Preprocessing:
      - Functions
      - Miscellaneous
      - Processed
      - Raw_data
      - Scripts
      - S1_README.md
      - run_step1.sh
  - S2_get_beta_diversity:
      - S2_README.md
      - Scripts
      - betadiv_input
      - betadiv_output
      - pre_processed
  - S3_get_best_predictors:
      - Functions
      - S3_README.md
      - Scripts
      - best_predictors
      - relationship_direction
  - S4_run_BBGDM:
      - S4_README.md
      - bbgdm_output
      - functions
      - scripts
  - S5_Null_models:
      - S5_README.md
      - functions
      - null_output
      - scripts
  - S6_Synthesis_model:
      - S6_README.md
      - data
      - functions
      - scripts
  - S7_Model_outputs_figures_and_tables:
      - extended_data
      - main_figures
      - model
```

### Step 1: Preprocessing (`S1_Preprocessing`)
- **Description**: Raw datasets are preprocessed to extract key metrics, clean data, and standardize formats for downstream analyses.
- **Key Contents**:
  - **Scripts**: `pre_processing.R`, `Figure_S1_and_S4.R`
  - **Data**: Raw ecological datasets, auxiliary files (e.g., Human Footprint data).
  - **Processed Output**: Cleaned datasets (`Processed/*.rds`).
- **README**: Details preprocessing steps, tools, and key decisions.

---

### Step 2: Beta Diversity Analysis (`S2_get_beta_diversity`)
- **Description**: Calculates species and trait beta diversity using functional dissimilarity matrices.
- **Key Contents**:
  - **Scripts**: `1-get_functional_dissimilarities.R`, `2-get_beta_diversity.R`
  - **Data**: Preprocessed inputs, dissimilarity matrices, and beta diversity results.
  - **Output**: Functional beta diversity metrics (`betadiv_output/*.rds`).
- **README**: Includes instructions for running beta diversity analysis and an overview of input/output data.

---

### Step 3: Best Predictor Selection (`S3_get_best_predictors`)
- **Description**: Identifies the best predictors for beta diversity trends using **Generalised Dissimilarity Models**.
- **Key Contents**:
  - **Scripts**: `get_best_predictors.R`
  - **Functions**: Custom functions for model fitting (`functions_get_best_predictors.R`).
  - **Output**: Predictor coefficients and R² values (`best_predictors/*.rds`).
- **README**: Explains the statistical approach and provides an overview of outputs.

---

### Step 4: BBGDM Analysis (`S4_run_BBGDM`)
- **Description**: Runs the `BBGDM` (Bayesian Bootstrap Generalized Dissimilarity Modeling) framework to assess community turnover along human pressure gradients.
- **Key Contents**:
  - **Scripts**: `main_script_S4.R`, `plot_figure_1.R`
  - **Functions**: `helper_functions_S4.R`
  - **Output**: BBGDM model results and visualizations.
- **README**: Describes the modeling workflow and instructions for reproducing results.

---

### Step 5: Null Models (`S5_Null_models`)
- **Description**: Generates null models to evaluate the significance of BBGDM results and trait/species replacement patterns.
- **Key Contents**:
  - **Scripts**: `1-run_null_analysis.R`, `2-null_model_results.R`
  - **Functions**: `functions_run_null_bbgdm.R`
  - **Output**: Null model results (`null_output/*.rds`).
- **README**: Details null model assumptions and methods for significance testing.

---

### Step 6: Synthesis and Modeling (`S6_Synthesis_model`)
- **Description**: Combines results from previous steps into a synthesis model to assess the variation in responses across datasets.
- **Key Contents**:
  - **Scripts**: Includes scripts for direction, magnitude, shape, and convergence analyses.
  - **Data**: Synthesis datasets (`data/synthesis_data.xlsx`).
  - **Output**: model results and figures 2-6.
- **README**: Includes metadata for synthesis data and modeling steps.

---

### Step 7: Outputs and Figures (`S7_Model_outputs_figures_and_tables`)
- **Description**: Contains all final models and figures for the manuscript.
- **Key Contents**:
  - **Figures**: Extended data figures, main manuscript figures (`*.pdf`).
  - **Models**: RDS files for final models.
- **README**: Summarizes the structure of outputs and figure generation.

## Usage Instructions

### Prerequisites
- **Software**: R (≥4.4.1), RStudio, and required packages (`tidyverse`, `brms`, `gawdis`, etc.).
- **System Requirements**: Sufficient RAM and storage for large datasets and computationally intensive models.

### Running the Pipeline
1. Clone the repository and set up the working directory.
2. Follow the README in each step to execute scripts in sequence.
3. Use the synthesis model (`S6_Synthesis_model`) to generate results and figures for publication.

## Additional Notes
- **Metadata**: Each dataset includes a corresponding metadata file to ensure reproducibility and alignment with FAIR principles.
- **Collaboration**: Contributions to improve this pipeline are welcome. Please submit a pull request or open an issue.
- **Citation**: If using this repository, cite the manuscript: *Graco-Roza et al., Human Pressure Homogenises Species and Traits Globally (XXXX)*.

For further details, see individual READMEs in each step.

