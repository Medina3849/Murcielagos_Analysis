Archaeological Sediment Analysis: Cueva de los Murciélagos
Overview
This repository contains the complete R code for the physico-chemical analysis of archaeological sediments from the Cueva de los Murciélagos archaeological site. The analysis examines sediment properties across different archaeological periods (Palaeolithic, Neolithic, Chalcolithic, Bronze Age, and Classical Period) to understand site formation processes and human occupation patterns.

Scientific Context
This code supports the research presented in our scientific publication analyzing sediment characteristics from the Cueva de los Murciélagos site. The analysis provides insights into:

Temporal variations in sediment composition

Human impact on sedimentary records

Site formation processes across archaeological periods

Relationships between different physico-chemical parameters

Data Description
The dataset includes 13 sediment samples with the following variables:

Sample: Sample identifiers (MZ-1 to MZ-13)

Period: Archaeological period (Palaeolithic to Classical Period)

Hum: Humidity content (%)

LOI: Loss on Ignition (%)

EC: Electrical Conductivity (dS/m)

% CO3: Carbonate content (%)

M.S.: Magnetic Susceptibility (10⁻⁸ m³/kg)

pH: Soil pH

OrgC: Organic Carbon (%)

Sand: Sand content (%)

Silt: Silt content (%)

Clay: Clay content (%)

Analysis Methods
The code implements a comprehensive analytical workflow:

1. Statistical Analysis
Descriptive statistics by archaeological period

Kruskal-Wallis tests with Benjamini-Hochberg correction

Multiple comparison analysis

2. Multivariate Analysis
Principal Component Analysis (PCA)

Hierarchical clustering (Ward's method)

Silhouette analysis for optimal cluster determination

3. Correlation Analysis
Pearson correlation matrix

Significance testing (p < 0.01)

Visualization of significant relationships

Requirements
R Packages
tidyverse, ggpubr, FactoMineR, factoextra, viridis, 
patchwork, rstatix, corrplot, scales, ggrepel, 
dendextend, cluster
All required packages are automatically installed and loaded by the script.

Usage
1. Run the complete analysis:
source("murcielagos_analysis.R")
2. The script will:

Install missing packages automatically

Perform all statistical analyses

Generate publication-quality figures

Export results to CSV files

Create a dedicated output directory: Murcielagos_Analysis_Results/

Output
Generated Figures
Figure 1: Boxplots of key parameters across archaeological periods

Figure 2: PCA analysis (sample distribution and variable contributions)

Figure 3: Hierarchical clustering dendrogram

Figure 4: Correlation matrix with significance indicators

Data Outputs
Descriptive statistics by period

Statistical test results

Correlation matrices

Cluster assignments

Publication-ready tables

File Formats
Figures: TIFF format (600 DPI, LZW compression)

Data: CSV format

All outputs organized in the results directory

Interpretation
Key Analytical Insights
PCA: Reveals major gradients in sediment composition

Clustering: Identifies natural groupings of samples

Statistical tests: Highlights significant period-based differences

Correlations: Shows relationships between physico-chemical parameters
License
This project is licensed for academic use. Please contact the authors for commercial applications.

Contact
For questions regarding this analysis or the underlying data, please contact the corresponding author of the associated publication.

Repository Structure
├── murcielagos_analysis.R    # Main analysis script
├── README.md                 # This file
└── Murcielagos_Analysis_Results/  # Generated output directory
    ├── Figure1_Parameters.tiff
    ├── Figure2_PCA.tiff
    ├── Figure3_Clustering.tiff
    ├── Figure4_Correlations.tiff
    ├── descriptive_statistics.csv
    ├── statistical_results.csv
    ├── correlation_matrix.csv
    ├── significant_correlations_99.csv
    ├── cluster_assignments.csv
    └── publication_table.csv

Computational Notes
R Version: Compatible with R 4.0.0+

Processing Time: < 2 minutes on standard hardware

Memory Requirements: Minimal (< 100 MB RAM)

Dependencies: All packages available from CRAN