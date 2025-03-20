# Clinical Trials Statistical Analysis in R

## Overview
This repository showcases advanced statistical analysis techniques applied to **clinical trials data** using R. It contains **R scripts** that focus on **non-parametric statistical tests**, **permutation tests**, and **data visualization**. The analysis includes **clinical trial treatment comparisons**, **blood protein concentration studies**, and **rank-based hypothesis testing**.

## Files Included

### **1. ClinicalTrials_nonpar.R**
- **Purpose:** Performs statistical analysis on clinical trial data using non-parametric methods.
- **Key Features:**
  - Implements **Wilcoxon-Mann-Whitney (WMW) test** for rank-sum comparisons.
  - Uses **permutation tests** to approximate the null distribution.
  - Conducts **Kruskal-Wallis tests** for multi-group comparisons.
  - Analyzes the effects of treatments on **rat body weight and protein concentration**.
  - Generates **rank-based relative effect estimates**.
- **Functions Used:**
  - `wilcox.test()`, `kruskal.test()`
  - `rank()`, `table()`, `mean()`, `sd()`
  - `ggplot2` for visualization

### **2. ClinicalTrials_mctp.R**
- **Purpose:** Provides statistical solutions for comparing treatment groups in clinical trials.
- **Key Features:**
  - Computes **permutation-based probability mass functions (PMF)**.
  - Simulates the **null distribution** of rank-based test statistics.
  - Generates **confidence intervals for relative effects**.
  - Implements **multiple contrast test procedures (MCTP)** to compare treatment conditions.
- **Functions Used:**
  - `sample()`, `sum()`, `mean()`, `quantile()`
  - `pander()`, `ggplot2`, `boxplot()`, `stripchart()`
  - `npar.t.test()` from the `nparcomp` package for relative effects

## Features
- **Statistical Analysis**
  - Wilcoxon-Mann-Whitney test, Kruskal-Wallis test
  - Non-parametric multiple contrast test procedures (MCTP)
  - Rank-based confidence intervals
- **Visualization**
  - Boxplots, density plots, and histograms (`ggplot2`)
  - Permutation distributions and rank-based PMF plots
- **Reproducible Reports**
  - Well-structured code for easy replication of results
  - Uses `pander` and `knitr` for formatted output

## Installation & Requirements
Ensure you have the required R packages installed:
```r
install.packages(c("tidyverse", "ggplot2", "dplyr", "nparcomp", "multcomp", "pander"))
```

## Usage
1. Clone the repository:
   ```sh
   git clone https://github.com/yourusername/Clinical-Trials-Stats-R.git
   cd Clinical-Trials-Stats-R
   ```

2. Run the scripts in RStudio or R console:
   ```r
   source("ClinicalTrials_nonpar.R")
   source("ClinicalTrials_mctp.R")
   ```

## Author
**Abhinav Mishra**  
Master's in Bioinformatics, FU Berlin

## License
MIT License
