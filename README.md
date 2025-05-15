# Finding Optimal Parameters for Type-II Censoring Experiments: Kolmogorov-Smirnov Statistic

This repository contains a research project conducted at IIT Kharagpur in 2020 as part of my graduation in Mathematics. The project focused on designing Type-II censoring experiments to understand the underlying data distribution by computing percentile points for the Kolmogorov-Smirnov (KS) statistic.

## Project Overview

Type-II censoring is widely used in survival analysis and reliability testing. This project investigates the optimal parameters for Type-II censoring experiments to accurately estimate the underlying data distribution. The Kolmogorov-Smirnov (KS) statistic is used as a tool to assess the goodness of fit, and its percentile points under Type-II censoring were computed using numerical methods.

## Repository Structure

- `R-codes/`: Contains all R scripts used for analysis and simulations
  - Type-II censoring analysis scripts
  - KS statistic distribution analysis
  - Quantile calculations for various distributions
  - Simulation codes
- `Simulation-outputs/`: Results from various simulation experiments
- `Research_papers/`: Relevant research papers and literature
- `Reference_Material/`: Additional reference materials and documentation
- `Graphics_pdf_charts/`: PDF versions of generated charts and visualizations
- `Charts/`: Original chart files and visualizations

## Key Components

### R Scripts
The `R-codes` directory contains several important scripts:
- `Type2_censoring_final_code.R`: Main implementation for Type-II censoring analysis
- `KS_Statistic_Distribution_Analysis_Final_Code.R`: Analysis of KS statistic distribution
- `Simulation.R`: Core simulation framework
- Various distribution-specific analysis scripts (Log-normal, Uniform, Laplacian)

### Documentation
- Master's thesis and presentation materials are included in the root directory
- Research papers and reference materials are available in their respective directories

## Research Objectives

1. Design Type-II censoring experiments for efficient parameter estimation
2. Calculate percentile points for the Kolmogorov-Smirnov statistic under Type-II censoring
3. Analyze the impact of censoring on the accuracy of distributional inference
4. Develop methods to optimize experiment parameters for better inference of data distributions

## Key Findings

- Statistical properties of Type-II censored data
- Methods for optimizing experiment parameters
- Computed percentile points for the KS statistic under Type-II censoring
- Analysis of distributional inference accuracy

## Applications

This research has potential applications in:
- Survival analysis and reliability testing
- Statistical quality control
- Distributional inference in censored datasets
- Experimental design optimization

## Getting Started

To reproduce the analysis:
1. Navigate to the `R-codes` directory
2. Install required R packages (list to be added)
3. Run the simulation scripts in the following order:
   - `Simulation.R`
   - `Type2_censoring_final_code.R`
   - Distribution-specific analysis scripts as needed

## Contact

For questions or collaboration opportunities, please refer to the contact information in the thesis document.

## License

This project is part of academic research conducted at IIT Kharagpur. Please refer to the thesis document for usage rights and permissions.
