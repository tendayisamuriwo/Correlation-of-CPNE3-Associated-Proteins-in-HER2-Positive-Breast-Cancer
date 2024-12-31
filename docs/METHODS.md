# Detailed Methods

## Data Selection and Processing

### Initial Protein Selection
- Based on 52 proteins identified in ERBB2 overexpressing HMLEC 3.6
- Proteins identified as differentially expressed in response to CPNE3 knockdown

### Discovery Cohort Processing (n=40)
- Selected 36 proteins consistently quantified across all samples
- Excluded 16 proteins due to incomplete quantification:
  - AHSG, ALPP, HIST1H1A, HIST1H1C, HIST1H1D, HIST1H3A
  - MGST1, PNPLA7, PPP6R2, PYGM, RPL32, S100A10
  - SLC3A2, TIMELESS, UBC, ZFP28

### Validation Cohort Processing (n=75)
- Selected 45 proteins consistently quantified across all samples
- Excluded 7 proteins due to missing values or detection issues:
  - ALPP, HIST1H3A, HMGA1, PNPLA7
  - PTGES, UBC, ZFP28

## Statistical Analysis Framework

### Data Quality Assessment

1. **Normality Testing**
   - Method: Shapiro-Wilk test (Shapiro and Wilk, 1965)
   - Implementation: stats::shapiro.test()
   - Significance threshold: α = 0.05
   - Results reported as percentage of passed tests

2. **Homoscedasticity Testing**
   - Method: Breusch-Pagan test (Breusch and Pagan, 1979)
   - Implementation: car::ncvTest()
   - Tests for constant variance in regression residuals
   - Critical for validating correlation analysis assumptions

3. **Correlation Analysis**
   - Pearson Correlation (Pearson, 1895):
     - Measures linear relationships
     - Assumes normally distributed variables
     - Range: -1 to +1
   
   - Spearman Rank Correlation (Spearman, 1904):
     - Non-parametric rank-based approach
     - Robust to non-linear relationships
     - Less sensitive to outliers

   Implementation:
   ```R
   # Correlation calculations
   cor_matrix_pearson <- cor(expression_data, method = "pearson")
   cor_matrix_spearman <- cor(expression_data, method = "spearman")
   
   # Visualization
   corrplot(cor_matrix,
           method = "circle",
           order = "hclust",
           addrect = 3,
           tl.col = "black",
           tl.cex = 0.6)
   ```

## Data Structure Requirements

### Input Data Format
- CSV files with consistent headers
- Required columns:
  - Sample identifiers
  - Protein expression values
  - Clinical annotations

### Output Specifications
1. **Correlation Matrices**
   - Format: PDF files
   - Size: 12 x 12 inches
   - Separate files for Pearson and Spearman correlations

2. **Statistical Results**
   - Format: Text files
   - Content: Correlation coefficients and p-values
   - Separate files for discovery and validation cohorts

## Reproducibility Information

### System Requirements
- R version 4.3.2 or higher
- Required R packages and versions:
  ```R
  library(stats)     # Core statistical functions
  library(corrplot)  # Correlation visualization
  library(car)       # Additional statistical tests
  ```

## References

### Statistical Methods
1. Pearson, K. (1895). *Proceedings of the Royal Society of London, 58*, 240-242.
2. Spearman C. (1904). *American Journal of Psychology, 15*(1), 72-101.
3. Shapiro, S.S. and Wilk, M.B. (1965). *Biometrika, 52*, 591-611.
4. Breusch, T.S., and Pagan, A.R. (1979). *Econometrica, 47*(5), 1287-94.
5. Bewick V., Cheek L., Ball J. (2003). *Critical Care, 7*(6):451-9.