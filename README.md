# Finding Optimal Parameters for Type-II Censoring Experiments: Kolmogorov-Smirnov Statistic

This repository contains a comprehensive research project conducted at IIT Kharagpur in 2020 as part of my Master's degree in Mathematics. The project focused on designing Type-II censoring experiments to understand the underlying data distribution by computing percentile points for the Kolmogorov-Smirnov (KS) statistic under various censoring scenarios.

## Project Overview

Type-II censoring is widely used in survival analysis and reliability testing, where experiments are terminated after a predetermined number of failures occur. This project investigates the optimal parameters for Type-II censoring experiments to accurately estimate the underlying data distribution. The Kolmogorov-Smirnov (KS) statistic is used as a tool to assess the goodness of fit, and its percentile points under Type-II censoring were computed using numerical methods and Brownian bridge simulations.

## Research Objectives

1. **Design Type-II censoring experiments** for efficient parameter estimation
2. **Calculate percentile points** for the Kolmogorov-Smirnov statistic under Type-II censoring
3. **Analyze the impact of censoring** on the accuracy of distributional inference
4. **Develop methods to optimize experiment parameters** for better inference of data distributions
5. **Compare Type-I and Type-II censoring** methodologies and their effects on KS statistic distributions
6. **Investigate multiple probability distributions** including Normal, Log-normal, Uniform, Laplacian, and Logistic distributions

## Repository Structure

### Core Analysis Scripts (`R-codes/`)
- **`Simulation.R`**: Core simulation framework for KS statistic computation and Brownian bridge analysis
- **`Type2_censoring_final_code.R`**: Main implementation for Type-II censoring analysis with quantile calculations
- **`Type1_censoring_final_code.R`**: Type-I censoring analysis for comparison studies
- **`KS_Statistic_Distribution_Analysis_Final_Code.R`**: Comprehensive KS statistic distribution analysis
- **`KS_Statistic_Distribution_Analysis_n_15.R`**: Specific analysis for sample size n=15
- **`KS_Statistic_Location_scale_Unknown.R`**: Analysis when location and scale parameters are unknown

### Distribution-Specific Analysis
- **`KS_Statistic_Quantile_values_uniform.R`**: Uniform distribution analysis
- **`KS_Statistic_Quantile_Differences_Log-normal.R`**: Log-normal distribution quantile differences
- **`KS_Statistic_Quantile_Differences_Laplacian.R`**: Laplacian distribution analysis
- **`ECDF_Quantile_Differences_TypeII_Censoring.R`**: Empirical CDF quantile differences

### Version Control
- Multiple versions of distribution analysis codes (v1-v4) showing the evolution of the methodology

### Output and Results (`Simulation-outputs/`)
- **`KS_Statistic_values_n_15.csv`**: KS statistic values for sample size 15
- **`KS_Statistic_Quantile_values_1.csv`** and **`KS_Statistic_Quantile_values_1_2.csv`**: Quantile values for different scenarios
- **`KS_Statistic_Quantile_Differences_Normal.csv`**: Normal distribution quantile differences
- **`KS_Statistic_Quantile_Differences_Wilcox_n_9.csv`** and **`KS_Statistic_Quantile_Differences_Wilcox_n_15.csv`**: Wilcox-based analysis results
- **`Quartile_values.csv`**: Comprehensive quartile analysis results

### Visualizations (`Charts/`)
- **Quantile difference plots** for various sample sizes (n=8, 12, 15, 25, 30)
- **Distribution-specific visualizations** for Normal, Log-normal, Uniform, Laplacian, and Logistic distributions
- **Type-I vs Type-II censoring comparisons**
- **KS statistic distribution plots** and p-value analyses

### Research Documentation
- **`Research_papers/`**: Relevant research papers and literature review materials
- **`Reference_Material/`**: Additional reference materials and documentation
- **`Graphics_pdf_charts/`**: PDF versions of generated charts and visualizations

### Thesis and Presentation
- **`18MA40019_YASH_SETHI_Masters_project_thesis_2020_final_v1.pdf`**: Complete master's thesis
- **`18MA40019_YASH_SETHI_Masters_project_presentation_27052020_v1.pdf`**: Project presentation
- **`18MA40019_YASH_SETHI_Masters_thesis_latex_code.tex`**: LaTeX source code for the thesis

## Methodology

### Type-II Censoring Implementation
The project implements Type-II censoring by:
1. Generating random samples from various probability distributions
2. Applying censoring by retaining only the first `censor_obj` observations
3. Computing empirical CDF and KS statistic components (D+ and D-)
4. Scaling the KS statistic by √n for asymptotic analysis

### Brownian Bridge Simulation
- Uses Brownian bridge processes to simulate the Kolmogorov distribution
- Compares empirical KS statistics with theoretical Kolmogorov distribution
- Analyzes quantile differences and distributional properties

### Multi-Distribution Analysis
The research covers multiple probability distributions:
- **Normal Distribution**: Standard normal with known/unknown parameters
- **Log-normal Distribution**: Log-normal with various shape parameters
- **Uniform Distribution**: Uniform(0,1) distribution
- **Laplacian Distribution**: Double exponential distribution
- **Logistic Distribution**: Logistic distribution for comparison

## Key Findings

### Statistical Properties
- **Type-II censoring effects** on KS statistic distributions
- **Quantile point calculations** for various censoring levels
- **Comparison between Type-I and Type-II censoring** methodologies
- **Impact of sample size** on distributional inference accuracy

### Optimization Results
- **Optimal experiment parameters** for different distributions
- **Censoring level recommendations** based on distribution type
- **Sample size effects** on KS statistic accuracy
- **Parameter estimation efficiency** under censoring

### Distribution-Specific Insights
- **Normal distribution**: Optimal censoring levels and parameter sensitivity
- **Log-normal distribution**: Shape parameter effects on censoring efficiency
- **Uniform distribution**: Baseline comparison for censoring effects
- **Laplacian distribution**: Heavy-tailed distribution behavior under censoring

## Applications

This research has potential applications in:
- **Survival analysis** and reliability testing
- **Statistical quality control** and process monitoring
- **Distributional inference** in censored datasets
- **Experimental design optimization** for limited resources
- **Clinical trials** and medical research
- **Industrial reliability testing** and failure analysis

## Getting Started

### Prerequisites
- R programming language (version 3.6 or higher recommended)
- Required R packages:
  - `sde` (for Brownian bridge simulation)
  - `stats` (for statistical functions)
  - `graphics` (for plotting)

### Installation and Setup
1. Clone this repository
2. Install required R packages:
   ```r
   install.packages(c("sde", "stats", "graphics"))
   ```
3. Navigate to the `R-codes` directory

### Running the Analysis
1. **Start with core simulation**:
   ```r
   source("Simulation.R")
   ```

2. **Run Type-II censoring analysis**:
   ```r
   source("Type2_censoring_final_code.R")
   ```

3. **Perform distribution-specific analysis**:
   ```r
   source("KS_Statistic_Distribution_Analysis_Final_Code.R")
   ```

4. **Compare with Type-I censoring**:
   ```r
   source("Type1_censoring_final_code.R")
   ```

### Output Interpretation
- Check `Simulation-outputs/` directory for CSV results
- View generated plots in `Charts/` directory
- Refer to thesis document for detailed interpretation

## Research Contributions

1. **Novel methodology** for computing KS statistic percentiles under Type-II censoring
2. **Comprehensive comparison** of Type-I vs Type-II censoring effects
3. **Multi-distribution analysis** providing insights across different probability models
4. **Optimization guidelines** for experimental design under resource constraints
5. **Brownian bridge implementation** for theoretical distribution comparison

## Contact and Collaboration

For questions, collaboration opportunities, or academic inquiries:
- **Author**: Yash Sethi (18MA40019)
- **Institution**: IIT Kharagpur
- **Year**: 2020
- **Department**: Mathematics

## License and Usage

This project is part of academic research conducted at IIT Kharagpur. The code and methodology are available for educational and research purposes. Please refer to the thesis document for detailed usage rights and permissions.

## Citation

If you use this research in your work, please cite:
```
Sethi, Y. (2020). Finding Optimal Parameters for Type-II Censoring Experiments: 
Kolmogorov-Smirnov Statistic. Master's Thesis, IIT Kharagpur.
```

## Acknowledgments

- **Supervisor**: [Supervisor Name] - IIT Kharagpur
- **Department**: Department of Mathematics, IIT Kharagpur
- **Research Group**: [Research Group Name]
- **Funding**: [Funding Source if applicable]
