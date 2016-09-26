# dPCA
Simple interactive PCA application for gene expression analysis

## Input

1. Data matrix with measurements
2. Sample information, such as experimental conditions (optional)
3. Results from differential expression tests, or any other external data that can be used to filter genes (optional)

## Data options

* The following normalization alternatives, and combinations thereof, are avaliable:
  - RPM
  - log2
  - zscore
* The data (genes) can be filtered using the following variables (multiple combinations of filters are possible):
  - Mean
  - Coefficient of variation
  - Any variable supplied using input option (3)
  
## PCA

* Perform and visualize PCA in either 2D or 3D
  - Color data points using sample information supplied using input option (2)
