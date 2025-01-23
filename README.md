# Downing-et-al-2018
Code and data for: Downing PA, Griffin AS, &amp; Cornwallis, CK. 2018. Sex differences in helping effort reveal the effect of future reproduction on cooperative behaviour in birds. Proceedings of the Royal Society B, 285: 20181164.


Number of supplementary items: two
1. SD_R_Code.R
2. SD_Table_S1-S4.xlsx


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: SD_R_Code.R

This R script contains all the code needed to replicate the analyses (including packages and functions).
- Data manipulation (lines 68 to 107)
- Mean sex differences in helping effort (lines 115 to 176)
- Publication bias (lines 184 to 194)
- Helping effort and breeding in the natal group (lines 198 to 413)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: SD_Tables_S1-S4.xlsx

This Excel document contains the following sheets:

- Table S1
	+ cooperative bird species excluded from the dataset and the reasons why
	+ 150 data rows = 150 species
	+ column descriptions:\
		A. species = latin binomial of each species (matches the Jetz et al. nomenclature)\
		B. reason for exclusion = why the species was excluded\
		C. reference = study with potentially relevant study\

- Table S2
	+ data on sex differences in provisioning for 20 species
	+ column descriptions:\
		A. English name = English name of each species\
		B. binomial = scientific name of each species (matches the Jetz et al. nomenclature)\
		C. what was measured = the metric used in each study to quantify helping effort\
		D. units = the unit associated with each metric\
		E. femaleX = mean female helper helping effort\
		F. var = variance in mean female helper helping effort\
		G. n.female = number of female helpers associated with femaleX\
		H. maleX = mean male helper helping effort\
		I. var = variance in mean male helper helping effort\
		J. n.male = number of male helpers associated with femaleX\
		K. n(nests) = the number of nests associated with the helper effort measure\
		L. n(groups) = the number of cooperative groups associated with the helper effort measure\
		M. random.effects = did the analysis include random effects?\
		N. covariates = did the analysis include covariates?\
		O. study.type = was the study observational or experimental?\
		P. duration (breeding seasons) = the length of the study\
		Q. nestWatches(mins) = total duration of nest watches\
		R. watches = further details on nest watches\
		S. method = were the nests watched by a human or recorded?\
		T. distance(meters) = the distance from the nest that observations took place\
		U. source = where in the study the data were extracted from\
		V. reference = study from which the data to calculate effect sizes were extracted

- Table S3
	+ details on sex biased philopatry and potentially confounding factors

- Table S4
	+ parameter estimates from statistical models (see SD_R_Code.R for further information)
	+ 36 rows = output from 8 statistical models
	+ variable names:\
	  beta = parameter estimate from the model\
		lwr CI = lower 95% credible interval from the posterior distribution of the model\
		upr CI = upper 95% credible interval from the posterior distribution of the model\


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
