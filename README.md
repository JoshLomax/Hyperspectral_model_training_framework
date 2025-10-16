# Hyperspectral Model Training Framework

A comprehensive machine learning pipeline for predictive modeling using hyperspectral imaging data. While originally developed for predicting volatile compounds and sugar concentrations in papaya fruit, this framework is designed to be adaptable for any hyperspectral imaging regression problem.

## Overview

This framework processes raw hyperspectral absorbance data through various spectral transformations and applies multiple machine learning models to predict target variables of interest. The pipeline is designed for high-performance computing (HPC) environments and can process multiple variables and transformations in parallel, making it suitable for large-scale hyperspectral analysis projects.

## Key Features

- **Flexible spectral preprocessing**: Multiple transformation methods that can be combined
- **Multiple feature selection methods**: Including PLS regression coefficients, SPA, genetic algorithms, UVE, CARS, and PCA
- **Comprehensive model suite**: 9 different regression algorithms with hyperparameter optimization
- **Scalable architecture**: Designed for HPC environments with parallel processing
- **Automated workflow**: Batch processing of all transformation-variable-model combinations
- **Reproducible results**: Seed-based reproducibility and comprehensive logging
- **Multiple data types**: Supports flesh, skin, and combined tissue analysis

## Project Structure

```
Hyperspectral_model_training_framework/
├── Scripts/
│   ├── 00-spectral-preprocessing.R    # Spectral data preprocessing functions
│   ├── 01-data-preparation.R          # Data preparation and train-test splitting
│   ├── 02-feature-selection.R         # Feature selection algorithms (6 methods)
│   ├── 03-model-functions.R           # Machine learning models (9 algorithms)
│   └── 04-main-pipeline.R             # Main pipeline orchestrator
├── data/
│   ├── absorbance_data.RDS            # Raw hyperspectral data
│   ├── all_response_data.csv          # Target variables
│   └── response_names.txt             # List of target variables
├── input_data/                        # Generated preprocessed spectral data
│   ├── smooth_data/
│   ├── derivative1_data/
│   ├── derivative2_data/
│   └── [other transformation combinations]/
├── results/                           # Model outputs
└── style/                            # R Markdown styling files
```

## Workflow

### 1. Spectral Data Preprocessing (`00-spectral-preprocessing.R`)

The preprocessing pipeline applies various spectral transformations commonly used in chemometrics:

- **Smoothing**: Noise reduction using moving average filters
- **Standard Normal Variate (SNV)**: Scatter correction and normalization
- **Baseline correction**: Polynomial detrending (1st or 2nd order)
- **Derivatives**: Enhanced feature extraction (1st and 2nd derivatives)

The framework automatically generates all valid combinations of these transformations, providing comprehensive data preparation for downstream modeling.

### 2. Feature Selection Methods (`02-feature-selection.R`)

Six feature selection algorithms are implemented:

- **PLSR B-coefficient filtering**: Selects wavelengths based on PLS regression coefficients
- **Successive Projections Algorithm (SPA)**: Minimizes collinearity between selected variables
- **Genetic Algorithm (GA)**: Evolutionary optimization for feature subset selection
- **Uninformative Variable Elimination (UVE)**: Removes variables with low reliability scores
- **Competitive Adaptive Re-weighted Sampling (CARS)**: Iterative wavelength selection
- **Principal Component Analysis (PCA)**: Dimension reduction and variable importance
- **None**: Uses all available wavelengths

### 3. Machine Learning Models (`03-model-functions.R`)

Nine regression algorithms with hyperparameter optimization:

- **Bayesian Regression** (BayesA): With convergence diagnostics and adaptive parameter tuning
- **Random Forest**: Grid search over mtry and ntree parameters
- **XGBoost**: Gradient boosting with early stopping
- **Support Vector Machine**: RBF kernel with C and sigma optimization
- **Neural Networks**: Feed-forward networks with PCA preprocessing for high-dimensional data
- **Ranger**: Fast random forest implementation
- **Partial Least Squares Regression (PLSR)**: With cross-validation for component selection
- **Linear Support Vector Regression**: L2-regularized linear models
- **Elastic Net/Ridge Regression**: Regularized linear models with alpha optimization

### 4. Main Pipeline (`04-main-pipeline.R`)

The main pipeline orchestrates the entire workflow with command-line interface support.

## Usage

### Prerequisites

```bash
# Create conda environment with required packages
CONDA_NAME="hs_modelling" 
mamba create -n $CONDA_NAME r-keras r-cli r-pacman r-tidyverse r-fs r-glue r-tidymodels r-furrr r-bglr r-randomforest r-caret r-glmnet r-coda r-xgboost r-kernlab r-ranger r-nnet r-pls r-getoptlong r-liblinear r::r-mdatools r::r-prospectr r-ga r-tictoc
conda activate $CONDA_NAME
```

### Adapting for Your Data

1. **Prepare your hyperspectral data**: Format as RDS file with:
   - Column 1: Response variable
   - Column 2: Sample ID
   - Column 3: Type (F=flesh, S=skin, or other tissue types)
   - Columns 4+: Wavelength measurements

2. **Prepare response variables**: CSV file with ID column and target variables
3. **Update variable names**: Create `response_names.txt` with your target variables
4. **Adjust preprocessing**: Customize transformations in `00-spectral-preprocessing.R` if needed

### Step 1: Preprocess Spectral Data

```r
# Run the preprocessing script
source("Scripts/00-spectral-preprocessing.R")
```

This generates multiple preprocessed datasets in `input_data/` directories.

### Step 2: Single Model Run (for testing)

```bash
Rscript Scripts/04-main-pipeline.R \
  --input input_data/derivative1_data/derivative1.RDS \
  --output results/derivative1_target_variable.RDS \
  --seeds 5 \
  --var target_variable \
  --response data/all_response_data.csv \
  --cpus 1 \
  --models "rf,xgboost" \
  --features "spa,ga" \
  --verbose
```

### Step 3: HPC Batch Processing

```bash
# Generate command file for all combinations
parallel --dry-run "Rscript Scripts/04-main-pipeline.R --input {1} --output results/{1/.}_{2}.RDS --seeds 5 --var {2} --response data/all_response_data.csv --cpus \$SLURM_CPUS_PER_TASK --verbose" \
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

The main pipeline (`04-main-pipeline.R`) accepts the following arguments:

- `--input`: Input RDS file with preprocessed spectral data
- `--output`: Output file path for results (RDS/CSV/TSV/TXT)
- `--seeds`: Number of random seeds for reproducibility (default: 5)
- `--var`: Target variable name to predict
- `--response`: CSV file containing all response variables
- `--cpus`: Number of CPU cores to use (default: 1)
- `--models`: Models to run (comma-separated or 'all')
  - Options: `bayes,rf,xgboost,svm,nnet,ranger,plsr,linear_svr,elastic_net`
- `--features`: Feature selection methods (comma-separated, 'all', or 'none')
  - Options: `plsr_bcoef,spa,ga,uve,cars,pca,none`
- `--verbose`: Enable verbose output

## Applications

This framework can be applied to various hyperspectral imaging problems including:

- **Agriculture**: Crop quality assessment, disease detection, nutrient analysis
- **Food Science**: Quality control, composition analysis, authenticity testing
- **Environmental Monitoring**: Vegetation analysis, pollution detection, soil characterization
- **Materials Science**: Composition analysis, quality control, defect detection
- **Medical Imaging**: Tissue analysis, diagnostic applications

## Output

Results include comprehensive metrics for each model:

- **Performance metrics**: RMSE, R², MAE, NRMSE, correlation coefficients
- **Training and test metrics**: Separate evaluation for both datasets
- **Model-specific outputs**: Variable importance, hyperparameters, convergence diagnostics
- **Timing information**: Feature selection and model training times
- **Reproducibility info**: Seeds, preprocessing methods, timestamps

## Model Performance Evaluation

The framework calculates extensive metrics including:
- Mean Absolute Error (MAE)
- Root Mean Square Error (RMSE)
- R-squared and Adjusted R-squared
- Normalized RMSE (NRMSE)
- Ratio of Performance to Deviation (RPD)
- Pearson correlation coefficients
- Residual statistics

## Contributing

This is an early version of the framework with several areas for improvement:

### Current Limitations
- **Model Training Efficiency**: Some models (particularly Bayesian methods) are prone to overfitting and convergence issues. Better regularization strategies and more robust convergence diagnostics are needed.
- **Classification Support**: Currently only supports regression models; classification algorithms are not yet implemented for categorical predictions.
- **Hyperparameter Optimization**: Limited automated tuning capabilities - currently uses grid search which may not find optimal parameters efficiently.
- **Memory Management**: Large datasets with many wavelengths may require optimization for memory usage, especially with neural networks.
- **Cross-validation Strategy**: Currently uses simple train-test splits; nested cross-validation for hyperparameter tuning would be more robust.

### Areas for Development
- Implementation of classification models for categorical predictions (SVM classification, Random Forest classification, etc.)
- Advanced cross-validation strategies and nested CV for hyperparameter tuning
- Integration of deep learning approaches specifically designed for spectral analysis
- Enhanced feature selection methods specific to hyperspectral data (e.g., interval-based methods)
- Improved computational efficiency and memory management for large datasets
- Comprehensive visualization and reporting tools for model comparison
- Docker containerization for improved portability and reproducibility
- Automated model selection based on problem characteristics
- Integration with cloud computing platforms

We welcome contributions in any of these areas! Please feel free to submit issues, feature requests, or pull requests.

## Example Use Case

The framework was originally developed for predicting volatile compounds and sugar concentrations in papaya fruit, demonstrating its application in agricultural quality assessment. The pipeline processes hyperspectral data from both flesh and skin tissues, applies multiple preprocessing transformations, and evaluates numerous modeling approaches to identify the best prediction strategy for each target compound.

## Performance Considerations

- **Computational complexity**: With 6 feature selection methods, 9 models, multiple seeds, and 3 dataset types, a single preprocessing method can generate hundreds of model runs
- **Memory requirements**: Neural networks and some feature selection methods may require substantial RAM for high-dimensional data
- **Runtime**: Bayesian methods and genetic algorithms can be computationally intensive; consider limiting these for initial exploration

## Authors

- Josh Lomax
- Ido Bar

## Citation

If you use this framework in your research, please cite TBA.

## Troubleshooting

### Common Issues

1. **Memory errors with neural networks**: Reduce the number of wavelengths using PCA or feature selection
2. **Bayesian model convergence**: The framework includes adaptive parameter tuning, but very noisy data may still cause issues
3. **Long runtimes**: Start with a subset of models (`--models "rf,plsr"`) and feature selection methods for initial testing
4. **Missing packages**: Ensure all required packages are installed in your conda environment

### Performance Tips

- Use `--features "none,pca"` for initial exploration to reduce computation time
- Start with faster models like `rf,plsr,elastic_net` before running computationally intensive methods
- Monitor memory usage when processing high-dimensional data
- Use appropriate number of CPU cores based on your system capabilities