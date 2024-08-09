# fungal_intensity_submission
Repository of code, files, and plots for the manuscript "Do fungi look like macroparasites? Quantifying the patterns and mechanisms of aggregation for host-fungal parasite relationships".

## Repository Structure
- `code/`: Contains the following R scripts to reproduce analyses and figures
	- `01_Aggregation_Groups.R`: Contains the code to aggregate raw data into  groups and perform group-level analyses (calculating fungal intensity means and variances, fitting distributions, and calculating individual-group aggregation metrics (Poulin's D and variations of the Coefficient of Variation).
		-- Output: `data/formatted/summary_aggregation_data.csv`
		-- Note: Cannot provide raw data files, therefore will not be able to execute code
	- `02_Assign_Epi_Phase.R`: Contains code to classify Sierra sites as "Invasion", "Post-Invasion", or "Enzootic" phase
		-- Output: `data/formatted/abund_and_prev.csv`
			   `data/formatted/summary_sierra_with_epi_phase.csv`
			   `data/formatted/analysis_aggregation_data.csv`
		-- Note: Cannot provide raw data files, therefore will not be able to execute first half of code
	- `03_Fig.1_TPL_Plot.R`: Contains the code used to create Figure 1 and calculate the slope of Taylor's Power Law (TPL) for each dataset. 
	- `04_TPL_Life_Stage.R`: Contains generalized linear mixed models, modeling the relationship between the log mean fungal load on the log variance of fungal load, incorporating region and life stage as possible predictor variables. 
	- `05_Fig.2_Distribution_Fit_Boxplot.R`: Contains code to create Figure 2 and compare which distributions fit the data best.
	- `06_Fig.3_Trajectory_Plot.R`: Contains code to produce Figure 3.
	- `07_Mean_Agg_Quadratic_Model.R`: Contains code that fits linear and quadratic models to the Sierra dataset and compares the performance. Also compares if this relationship changes based on life stage or epizoological phase.
	- `08_IPM_Fig.S5-S8.R`: Contains the code for the Integral Projection Model (IPM), as well as different parameters used. Contains code to produce Supplemental Materials Figures S5-S8.
		-- Output: `data/model/eco_df.csv`
		--  	   `data/model/evo_df.csv`
	- `09_Fig.4_Empirical_Model_Agg_Plot.R`: Contains the code to produce Figure 4.
	- `10_Supp_Table.S1.R`: Contains code to extract values for Supplemental Materials Table S1.
	- `11_Supp_Fig.S2_Spec_TPL.R`: Contains generalized linear mixed models, modeling the relationship between the log mean fungal load on the log variance of fungal load, incorporating region and species as possible predictor variables. Includes code to create Figure 3.
	- `12_Supp_Fig.S3_GOF.R`: Contains code to create Supplemental Materials Figure S3. Also contains code to extract information on goodness of fit (GOF) of a normal distribution on log-transformed data and find proportion of log-transformed data that conforms to a normal distribution / proportion of natural data that conforms to a log-normal distribution.
	- `13_Supp_Fig.S4_PA_GreenFrogs.R`: Contains code to create Supplemental Materials Figure S4.
- `data/formatted/`: contains the following csv files:
	- `abund_and_prev.csv`: Contains records of Sierra sites over multiple years and includes information on amphibian abundance, Bd prevalence, mean fungal load, and Poulin's D.
	- `analysis_aggregation_dataset.csv`: Contains groups aggregated by dataset, site, species, life stage, and season. This dataset includes aggregated values of mean, variance, Poulin's D, coefficient of variation metrics, and AIC values on fitted distributions.
- `data/model/`: contains the following csv files:
	- `eco_df.csv`: Output of the IPM without host evolution
	- `evo_df.csv`: Output of the IPM with host evolution of resistance

