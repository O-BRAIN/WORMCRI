# README

## Folder Structure
MATLAB analysis scripts

  * `sa_eLoreta.mat` - eLORETA inverse operator which is computed in the `scripts_minimal/prepare_eLORETA_inverse_operator.m` script
  * `S*_DFA_power_slope.mat` - values of alpha power, DFA and 1/f slope for each subject
  * `results_DFA_power_slope.mat` - combined results (power, DFA, and slope) for all subjects
  * `results_ga_DFA_power_slope.mat` - grand-average results (power, DFA, and slope)
  * `source_dfa_slope.csv` - DFA and slope values exported for analysis in R
  * `results_matrix_DFA.csv` - results of statistical analysis for DFA
  * `results_matrix_oneF.csv` - results of statistical analysis for 1/f slope
  
* `results_minimal` - folder with final results

  * `significant_sources_*.png` - brain plots with significant sources highlighted in red (before/after FDR as specified in the file name)
  * `tvalue_update_*.png` - brain plots of update-ignore t-values 
  * `significant_sources_*.txt` - list of significant sources for each ROI of the Harvard-Oxford atlas


* `toolboxes_minimal/haufe` - functions that are required for source space analysis (leadfield, eLORETA) and visualizations

  * These scripts are also available [in my repository](https://github.com/ctrltz/bci-brain-connectivity/tree/master/toolboxes/haufe)
  
## How to Use

1. [Download](https://www.parralab.org/nyhead/sa_nyhead.mat) the New York 
Head model (~680 MB) and save it as `sa_nyhead.mat` file in the `toolboxes_minimal/haufe` 
folder.

2. Configure all paths to the data and toolboxes in the following scripts:

* `source_space_config.m`: lines 4, 7, 10, 13, 16, 19, 22, 25, 28, 32
* `scripts/stats_in_R/behav_on_source.R`: line 7

3. Run `source_space_main.m` to calculate DFA and 1/f values.

4. Run the R script `scripts/stats_in_R/behav_on_source.R` to perform the statistical analysis.

5. Run `source_space_main.m` to plot the results.
