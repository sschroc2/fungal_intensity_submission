# fungal_intensity_submission
Repository of code, files, and plots for the manuscript "Do fungi look like macroparasites? Quantifying the patterns and mechanisms of aggregation for host-fungal parasite relationships".

## Repository Structure
- `code/`: Contains the following R scripts to reproduce analyses and figures
	- `01_Aggregation_Groups.R`: Contains the code to aggregate raw data into  groups and perform group-level analyses (calculating fungal intensity means and variances, fitting distributions, and calculating individual-group aggregation metrics (Poulin's D and variations of the Coefficient of Variation).
		-- Output: `fungal_intensity/data/formatted/summary_aggregation_data.csv`
		-- Note: Cannot provide raw data files, therefore will not be able to execute code
	- `02_Assign_Epi_Phase.R`: Contains code to classify Sierra sites as "Invasion", "Post-Invasion", or "Enzootic" phase
		-- Output: `fungal_intensity/data/formatted/abund_and_prev.csv`
			   `fungal_intensity/data/formatted/summary_sierra_with_epi_phase.csv`
			   `fungal_intensity/data/formatted/analysis_aggregation_data.csv`
		-- Note: Cannot provide raw data files, therefore will not be able to execute first half of code
	- `03_Fig.2_TPL_Plot.R`: Contains the code used to create Figure 2 and calculate the slope of Taylor's Power Law (TPL) for each dataset. 
	- `04_Fig.3_Species_TPL_Histogram.R`: Contains generalized linear mixed models, modeling the relationship between the log mean fungal load on the log variance of fungal load, incorporating region and species as possible predictor variables. Includes code to create Figure 3.
	- `05_TPL_Life_Stage.R`: Contains generalized linear mixed models, modeling the relationship between the log mean fungal load on the log variance of fungal load, incorporating region and life stage as possible predictor variables. 
	- `06_Fig.4_Distribution_Fit_Boxplot.R`: Contains code to create Figure 4 and compare which distributions fit the data best.
	- `07_Fig.5_Trajectory_Plot.R`: Contains code to produce Figure 5.
	- `08_Mean_Agg_Quadratic_Model.R`: Contains code that fits linear and quadratic models to the Sierra dataset and compares the performance. Also compares if this relationship changes based on life stage or epizoological phase.
	- `09_IPM_Fig.S3-S6.R`: Contains the code for the Integral Projection Model (IPM), as well as different parameters used. Contains code to produce Supplemental Materials Figures S3-S6.
		-- Output: `fungal_intensity/data/model/eco_df.csv`
		--  	   `fungal_intensity/data/model/evo_df.csv`
	- `10_Fig.6_Empirical_Model_Agg_Plot.R`: Contains the code to produce Figure 6.
	- `11_Supp_Table.S1.R`: Contains code to extract values for Supplemental Materials Table S1.
	- `12_Supp_Fig.S1.R`: Contains code to create Supplemental Materials Figure S1. Also contains code to extract information on goodness of fit (GOF) of a normal distribution on log-transformed data and find proportion of log-transformed data that conforms to a normal distribution / proportion of natural data that conforms to a log-normal distribution.
	- `13_Supp_Fig.S2.R`: Contains code to create Supplemental Materials Figure S2.


	- `09_PA_GreenFrog_Plots_FIP.R`: Contains code to create Figure S2.
	- `11_LifeStage_GLM_FIP.R`:
	- `12_Season_GLM_FIP.R`:
	- `13_GOF_Boxplots_FIP.R`: Contains code to create Figure S1.
	- `distribution_boxplot_code.R`: Contains code to create Figure 4.
	- `table_1_code_updated.R`: Contains code used to populate Table 1.
- `data/formatted/`: contains the following csv files:
	- `abund_and_prev.csv`: Contains records of Sierra sites over multiple years and includes information on amphibian abundance, Bd prevalence, mean fungal load, and Poulin's D.
	- `analysis_aggregation_dataset.csv`: Combination of `summary_aggregation_data.csv` (without Sierra records) and `summary_sierra_with_epi_phase.csv`. This is the dataset that is used for all calculations and figures, except the IPM.
	- `summary_aggregation_data.csv`: Contains groups aggregated by dataset, site, species, life stage, and season. This dataset includes aggregated values of mean, variance, Poulin's D, coefficient of variation metrics, and AIC values on fitted distributions.
	- `summary_sierra_with_epi_phase.csv`: Contains just the Sierra data, with the same information in `summary_aggregation_data.csv` plus classification into epizoological phase categories. Used in the IPM.
- `data/model/`: contains the following csv files:
	- `eco_df.csv`: Output of the IPM without host evolution
	- `evo_df.csv`: Output of the IPM with host evolution of resistance
- `results/plots`: Repository for saved figures and includes:
	- `agg_emp_mod.png`: Fig. 6
	- `distr_plot.png`: Fig. 4
	- `Focal_supplement.png`: Fig. S3
	- `gof_plot.png`: Fig. S1
	- `hiLD50_supplement.png`: Fig. S6
	- `map_study_sites.png`: Fig. 1
	- `nearfocal_supplement.png`: Fig. S4
	- `penn_green_full_plot.png`: Fig S2
	- `spec_slope_histo.png`: Fig. 3
	- `Tolerance_supplement.png`: Fig. S5
	- `tpl_lin_reg.png`: Fig. 2
	- `trajectory_plot.png`: Fig. 5
