# Hyperspectral Model Training Framework

A comprehensive machine learning pipeline for predictive modeling using hyperspectral imaging data. While originally developed for predicting volatile compounds and sugar concentrations in papaya fruit, this framework is designed to be adaptable for any hyperspectral imaging regression problem.

## Overview

This framework processes raw hyperspectral absorbance data through various spectral transformations and applies multiple machine learning models to predict target variables of interest. The pipeline is designed for high-performance computing (HPC) environments and can process multiple variables and transformations in parallel, making it suitable for large-scale hyperspectral analysis projects.

## Key Features

- **Flexible spectral preprocessing**: Multiple transformation methods that can be combined
- **Scalable architecture**: Designed for HPC environments with parallel processing
- **Model variety**: Integration of multiple ML algorithms for comprehensive comparison
- **Automated workflow**: Batch processing of all transformation-variable-model combinations
- **Reproducible results**: Seed-based reproducibility and comprehensive logging

## Project Structure

```
Hyperspectral_model_training_framework/
├── Scripts/
│   ├── 00-spectral-preprocessing.R   # Spectral data preprocessing functions
│   ├── 01-data-preparation.R         # Feature selection methods
│   ├── 02-feature-selection.R        # Machine learning models
│   ├── 03-model-functions.R          # Model evaluation functions
│   └── 04-main-pipeline.R            # Main pipeline orchestrator
├── data/
│   ├── absorbance_data.RDS           # Raw hyperspectral data
│   ├── all_response_data.csv         # Target variables
│   └── response_names.txt            # List of target variables
├── input_data/                       # Generated preprocessed spectral data
│   ├── smooth_data/
│   ├── derivative1_data/
│   ├── derivative2_data/
│   └── [other transformation combinations\]/
├── results/                          # Model outputs
└── style/                           # R Markdown styling files
```

## Workflow

### 1. Spectral Data Preprocessing (`00-spectral-preprocessing.R`)

The preprocessing pipeline applies various spectral transformations commonly used in chemometrics:

- **Smoothing**: Noise reduction using moving average filters
- **Standard Normal Variate (SNV)**: Scatter correction and normalization
- **Baseline correction**: Polynomial detrending (1st or 2nd order)
- **Derivatives**: Enhanced feature extraction (1st and 2nd derivatives)

The framework automatically generates all valid combinations of these transformations, providing comprehensive data preparation for
downstream modeling.

### 2. Machine Learning Pipeline (`04-main-pipeline.R`)

The main pipeline orchestrates:
- **Feature selection**: Multiple methods for dimensionality reduction
- **Model training**: Various regression algorithms optimized for spectral data
- **Model evaluation**: Cross-validation and performance metrics
- **Results aggregation**: Systematic collection of model outputs

### 3. High-Performance Computing Integration

Designed for SLURM-based HPC systems, the framework enables:
- Parallel processing of multiple parameter combinations
- Efficient resource utilization across compute nodes
- Scalable execution for large experimental designs
- Automated job management and error handling

## Usage

### Prerequisites

```bash
# Create conda environment with required packages
CONDA_NAME="hs_modelling"
mamba create -n $CONDA_NAME r-keras r-cli r-pacman r-tidyverse r-fs r-glue r-tidymodels r-furrr r-bglr r-randomforest r-caret r-glmnet      
r-coda r-xgboost r-kernlab r-ranger r-nnet r-pls r-getoptlong r-liblinear r::r-mdatools r::r-prospectr r-ga r-tictoc
conda activate $CONDA_NAME
```

### Adapting for Your Data

1. **Prepare your hyperspectral data**: Format as RDS file with samples as rows and wavelengths as columns
2. **Prepare response variables**: CSV file containing target variables for prediction
3. **Update variable names**: Modify `response_names.txt` with your target variables
4. **Adjust preprocessing**: Customize transformations in `00-spectral-preprocessing.R` if needed

### Step 1: Preprocess Spectral Data

```r
# Run the preprocessing script
source("Scripts/00-spectral-preprocessing.R")
```

### Step 2: Single Model Run (for testing)

```bash
Rscript Scripts/04-main-pipeline.R \
  -i input_data/derivative1_data/derivative1.RDS \
  -o results/derivative1_target_variable.RDS \
  -s 5 \
  -v target_variable \
  -r data/all_response_data.csv \
  -c 1 \
  --verbose
```

### Step 3: HPC Batch Processing

```bash
# Generate command file for all combinations
parallel --dry-run "Rscript Scripts/04-main-pipeline.R -i {1} -o results/{1/.}_{2}.RDS -s 5 -v {2} -r data/all_response_data.csv -c
\$SLURM_CPUS_PER_TASK --verbose" \
  ::: $(find input_data/ -name "*.RDS") \
  ::: $(cat data/response_names.txt) > commands.txt

# Submit job array
sbatch -a 1-$(cat commands.txt | wc -l) \
  --job-name=hyperspectral_modeling \
  --cpus-per-task=2 \
  --mem=8G \
  --time=6:00:00 \
  --export=ALL,CMDS_FILE=commands.txt,CONDA_NAME=hs_modelling \
  array.slurm
```

## Command Line Arguments

- `-i, --input`: Input RDS file with preprocessed spectral data
- `-o, --output`: Output RDS file for results
- `-s, --seed`: Random seed for reproducibility
- `-v, --variable`: Target variable name to predict
- `-r, --response`: CSV file containing all response variables
- `-c, --cores`: Number of CPU cores to use
- `-m, --models`: Models to run (default: "all")
- `-f, --features`: Feature selection methods (default: "all")
- `--verbose`: Enable verbose output

## Applications

This framework can be applied to various hyperspectral imaging problems including:

- **Agriculture**: Crop quality assessment, disease detection, nutrient analysis
- **Food Science**: Quality control, composition analysis, authenticity testing
- **Environmental Monitoring**: Vegetation analysis, pollution detection, soil characterization
- **Materials Science**: Composition analysis, quality control, defect detection
- **Medical Imaging**: Tissue analysis, diagnostic applications

## Output

Results include:
- Model performance metrics (RMSE, R², MAE, etc.)
- Feature importance rankings
- Cross-validation statistics
- Prediction accuracy assessments
- Model comparison summaries

## Contributing

This is an early version of the framework with several areas for improvement:

### Current Limitations
- **Model Training Efficiency**: Some models are prone to overfitting and would benefit from better regularization strategies
- **Classification Support**: Currently only supports regression models; classification algorithms are not yet implemented
- **Hyperparameter Optimization**: Limited automated tuning capabilities
- **Memory Management**: Large datasets may require optimization for memory usage
- **Visualization Tools**: Basic plotting functions could be enhanced

### Areas for Development
- Implementation of classification models for categorical predictions
- Advanced cross-validation strategies and nested CV for hyperparameter tuning
- Integration of deep learning approaches for spectral analysis
- Enhanced feature selection methods specific to hyperspectral data
- Improved computational efficiency and memory management
- Comprehensive visualization and reporting tools
- Docker containerization for improved portability

We welcome contributions in any of these areas! Please feel free to submit issues, feature requests, or pull requests.

## Example Use Case

The framework was originally developed for predicting volatile compounds and sugar concentrations in papaya fruit, demonstrating its application in agricultural quality assessment. The same methodology can be readily adapted for other spectroscopic prediction problems.    

## Authors

- Josh Lomax
- Ido Bar

## Citation

If you use this framework in your research, please cite TBA.