# dPCA
Simple interactive PCA application for gene expression analysis

## Requirements
* An operating system that can run MATLAB R2016b (see https://se.mathworks.com/support/sysreq/current_release/)
* Optional: A actual MATLAB installation (minimum R2015b in this case)

## Installation
* As a standalone application (recommended): 
  *  Windows: Download and run "Installer.exe"

* As an app within MATLAB:
  * Download and run the installer package "dPCA.mlappinstall". The app can the be selected from the "APPS" menu in MATLAB.

* Run with MATLAB, but without installation:
  * Download and run "dPCA.mlapp"

## Usage
The main application window is divided into 4 sections:

### Input
The following data files are supported as input:

1. Data matrix with measurements. A text file where (by default) rows are assumed to represent genes and columns represent samples. Row and column names are optional. Example:
  ```
      sample1 sample2 ... sampleN
  gene1   1   23  ... 2
  gene2   42  890 ... 0
  .
  .
  .
  geneN   4   26  ... 9
  ```

2. *Optional:* Sample information, such as experimental conditions. For instance structure as follows:
  ```
  sample_id   group_id
  sample1     tumor_type1
  sample2     normal
  .
  .
  .
  sampleN     tumor_type2
  ```

3. *Optional:* Results from differential expression tests, or any other external data that can be used to filter genes. If more than two groups have been compared, results for the additional comparisons can be added as extra columns. Example:
  ```
  gene_names	fold_change_comparison_1    p_val_comparison_1  fold_change_comparison_2    p_val_comparison_2
  gene1	2.9743670001268 0.0666802826837871  0.292958630875216   0.0454353028875967
  gene2	1.16014613348014    0.725913115168573   -0.337631319683019  0.86134019484511
  .
  .
  .
  geneN	0.326978314374701   0.845047901160025   0.0205013616738917  0.420237180190287
  ```

The minimal input needed is a text file containing a data matrix with measurements. Additional sample and gene information can be provided for filtering. For each file, one should specify whether the format is tab or comma delimited and whether the file has row and/or column names. The data matrix is by default assumed to have genes as rows and samples as columns. If this is not the case, the file can be transposed by marking the corresponding checkbox.

(Bonus feature: clicking "Load" in the data matrix input section without specifying a file will load an example dataset: Fisher's iris data, which contains measures of sepal length, sepal width, petal length, and petal width for 150 iris specimens.)

### Data options
**Normalization:** The following normalization alternatives, and combinations thereof, are avaliable:
  - RPM
  - log2
  - zscore
 
**Filtering:** The data (genes) can be filtered using the following variables (multiple combinations of filters are possible):
  - Mean
  - Coefficient of variation (CV)
  - Any variable supplied using input option (3)
 
The "Mean" and "CV" options are not available when data is zscore-transformed (since all genes then have the same mean and standard deviation).

To apply a filter, first select a variable from the "Cutoff variable" dropdown list. Then select a value from the slider. Finally, click "Add filter". By default, all genes with a value lower than the cutoff will be filtered out and ignored in the PCA computation. To instead include the genes below the cutoff and ingore those above, check the box "Reverse" below the slider. For some variables (such as p-value), it is more convenient to display the slider on a log scale. To do this, mark the "Log10 scale" checkbox. Multiple filters can be added this way. To remove one, select the filter from the "Active filters" list and click the "Remove filter" button.

The currenly included genes will be displayed next to the "Active filters" list (by names if this information was included in either the data matrix or the test file input, otherwise by row number). These genes can be exported to a text file by clicking the "Export genes" button below the list.
 
### PCA
- The results from the principal component analysis can be visualized in either 2D or 3D
- The plot is updated each time a filter is changed
- Data points can be colored using sample information supplied with input option (2)
- The figure window is interactive and supports image export to a number of formats (among other things)

### Messages
This text area will print any relevant information about the activities performed by the app, including warnings and error messages.
