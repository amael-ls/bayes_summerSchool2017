
#### Aim of prog:
## Prepare the data for Matlab:
#		- Clim in 2010, already transformed with the good scaling Growth/Mortality
#		- 2 data tables, one for above, one for below
#		- A data table for the allometries (C0_C1, etc...)
#		- A data table for dbh_star and dbh45m
#		- A data table for the scalings (1 for G, x for mortality), x = nb species
#
#### Units are very messy:
## - Adams2007 (from where the formulae are from):
#		* Growth in cm/yr
# 		* dbh in cm
#		* Height allometry (named H in his paper) in m
## - Purves2008 (from where the values are from):
#		* Growth in cm/yr
# 		* dbh in cm
# 		* Height allometry (named α in his paper) in m
## - On my own
#		* Growth in mm
#		* dbh in mm
# It implies I have to convert growth or dbh in some analytical formulae

library(data.table)
library(stringi)
library(rstan)

rm(list = ls())
graphics.off()

options(max.print = 500)

# growth_dt[, min(mean_temp_period_3)] = 8.226
# growth_dt[, max(mean_temp_period_3)] = 23.194
# growth_dt[, min(tot_pp_for_period_3)] = 29.48
# growth_dt[, max(tot_pp_for_period_3)] = 1631.1

#### Include files
source("../toolFunctions.R")
source("../createData/parametersAllometries.R")

#### Create folder for matlab
if (!dir.exists("./matlab_data"))
	dir.create("./matlab_data")

#### Load climate and species list
## Climate of 2010, averaged of 2006-2010 period (i.e., 5 years)
climate = readRDS("./data_randomForest/clim_2010.rds")

## Rstan results
output = readRDS("./results21species.rds")

## Dividing min_temp_coldest_period by 10 (Steve Vissault, it was stored as x10 to get integer)
climate[, min_temp_coldest_period := min_temp_coldest_period/10]
names(climate)[4:7] = c("K", "k", "PP", "pp")

## Species list, /!\ in the same order than glmm_rstan /!\
ls_14species = c("19481-BET-ALL", "19489-BET-PAP", "28731-ACE-SAC", "28728-ACE-RUB",
	"19462-FAG-GRA", "183295-PIC-GLA", "183397-TSU-CAN", "195773-POP-TRE", "183385-PIN-STR",
	"18032-ABI-BAL", "505490-THU-OCC", "183302-PIC-MAR", "18034-PIC-RUB", "183319-PIN-BAN")

ls_21species = readRDS("../data/ls_21species.rds")

sp_index = which(ls_21species %in% ls_14species)
ls_species = ls_21species[sp_index]
nbSpecies = length(ls_species)

## Subset te output
output = output[sp_index]
(n = length(output))

#### Tool functions
## Extract parameters from LMM Growth
extractParamsGrowth = function(bayes)
{
	meanParams = summary(bayes)$summary[, "mean"]
	meanParams = meanParams[names(meanParams) != "lp__"]
	return (meanParams)
}

## Extract parameters from GLMM results of mortality
extractParamsMortality = function(id, species = ls_species[id], quadratic = FALSE)
{
	if (!quadratic)
		model = readRDS(paste0("../glmm_stanFull/array_", id, "/", species,".rds"))

	if (quadratic)
		model = readRDS(paste0("../glmm_stanFull/array_", id, "/", species,"_quadratic.rds"))

	beta_reg = summary(object = model, pars = "beta")
	intercept = summary(object = model, pars = "alpha")
	fixef = c(intercept = intercept[, "mean"], beta_reg[, "mean"])

	return (fixef)
}

#### Functions to calculate parameters value for a given climate and canopy status
## From R to Matlab (growth params), order: intercept, dbh, dbh^2, K, K^2, PP, PP^2
# Nothing to do

## From R to Matlab (mortality), order: intercept, dbh, (dbh^2), K, (K^2), PP, PP^2
rToMatlabM = function(temp, precip, params, species, ppa_status, quadratic = FALSE)
{
	n = length(temp)

	# Intercept, note that ppa1 is sp-specific according to extractParams fct
	beta_0 = params["intercept"] + ppa_status * params["canopy_statusTRUE"]

	# Phi
	beta_1 = params["dbh"]
	if (quadratic)
		beta_2 = params["I(dbh^2)"]

	# Temp
	beta_3 = params["k"] + ppa_status * params["canopy_statusTRUE:k"]
	if (quadratic)
		beta_4 = params["I(k^2)"] + ppa_status * params["canopy_statusTRUE:I(k^2)"]

	# Precip
	beta_5 = params["pp"] + ppa_status * params["canopy_statusTRUE:pp"]
	beta_6 = params["I(pp^2)"] + ppa_status * params["canopy_statusTRUE:I(pp^2)"]

	if (quadratic)
		return (list(rep(beta_0, n), rep(beta_1, n), rep(beta_2, n), beta_3,
			beta_4, beta_5, beta_6))

	# 1/(1 + exp(-logit_p)) is the sigmoid, the reciprocal function of logit
	return (list(rep(beta_0, n), rep(beta_1, n), beta_3, beta_5, beta_6))
}

#### Run for each species
## Create the data tables
dbh_params = data.table(species_id = ls_species, dbh_star10 = numeric(nbSpecies),
	dbh_infinity = numeric(nbSpecies))

## Non species-specific variables
scalingGrowth = c("min_K" = 8.226, "max_K" = 23.194, "min_pp" = 29.48, "max_pp" = 1631.1)
h_star10 = 10

pars_m = paste0("m", 0:4)

## Species-specific variables
for (i in 1:nbSpecies)
{
	clim_params = copy(climate)
	currentSpecies = ls_species[i]

	print(paste0("species: ", currentSpecies))

	# Mortality parameters
	fixef_mortality = extractParamsMortality(i, currentSpecies)
	scalingMortality = readRDS(paste0("../glmm_stanFull/array_", i,
		"/normalisation_mortality_data.rds"))

	write.csv(scalingMortality, paste0("matlab_data/", currentSpecies, "_sc_mort.csv"),
		row.names = FALSE)

	# Canopy height h* = 10m, associated species-specific dbh* (in mm)
	dbh_star = heightToDbh(height = h_star10,
		a = purves2007_allometries[species == currentSpecies, a],
		b = purves2007_allometries[species == currentSpecies, b])*10 # Times 10 to get it in mm

	# Infinty dbh (here 45 meters)
	dbh45m = heightToDbh(height = 45,
		a = purves2007_allometries[species == currentSpecies, a],
		b = purves2007_allometries[species == currentSpecies, b])*10 # Times 10 to get it in mm

	# Fill dbh_params
	dbh_params[i, c("dbh_star10", "dbh_infinity") := .(dbh_star, dbh45m)]

	# Growth climate scaling
	m = scalingGrowth["min_K"]
	delta = scalingGrowth["max_K"] - scalingGrowth["min_K"]
	clim_params[, K := (K - m)/delta]

	m = scalingGrowth["min_pp"]
	delta = scalingGrowth["max_pp"] - scalingGrowth["min_pp"]
	clim_params[, pp_growth := (pp - m)/delta]

	# Mortality climate scaling
	mu = scalingMortality[var == "k", mu]
	sig = scalingMortality[var == "k", sd]
	clim_params[, k := (k - mu)/sig]

	mu = scalingMortality[var == "pp", mu]
	sig = scalingMortality[var == "pp", sd]
	clim_params[, pp := (pp - mu)/sig]

	# Calculate parameters oversorey
	clim_params[, (pars_m) := rToMatlabM(k, pp, fixef_mortality, currentSpecies, TRUE)]

	write.csv(clim_params, paste0("matlab_data/", currentSpecies, "_over.csv"), row.names = FALSE)
	# writeMat(paste0("matlab_data/test"), clim_params_over = clim_params) # Compatibility?

	# Calculate parameters oversorey
	clim_params[, (pars_m) := rToMatlabM(k, pp, fixef_mortality, currentSpecies, FALSE)]

	write.csv(clim_params, paste0("matlab_data/", currentSpecies, "_under.csv"), row.names = FALSE)
}

#### Species-specific integral bounds
write.csv(dbh_params, "matlab_data/dbh_params.csv", row.names = FALSE)

#### Species-specific allometries
write.csv(C0_C1, "matlab_data/C0_C1.csv", row.names = FALSE)
write.csv(purves2007_allometries[species %in% ls_species],
	"matlab_data/purves2007_allometries.csv", row.names = FALSE)

#### Species list
write.csv(ls_species, "matlab_data/ls_species.csv", row.names = FALSE)

## sd and mean of each column is not necessary 1! The scaling was done on the
#	climatic data related to the tree database, not on the climate data of 2010
#	of the whole north America

#### Growth scaling
write.csv(scalingGrowth, "matlab_data/sc_growth.csv", row.names = FALSE)

#### Checking matlab and R results are in accordance
source("../createData/parametersAllometries.R")

i = 1
currentSpecies = ls_species[i]

dbhToCrownArea( 100, 50, purves2007_allometries[species == currentSpecies, a],
	purves2007_allometries[species == currentSpecies, b],
	purves2007_allometries[species == currentSpecies, T],
	C0_C1, TRUE )

# Function integrang is in R_s_starFixed.R prog, as climate data table
integrand(150, climate[1, K], climate[1, PP], climate[1, k], climate[1, pp],
	fixef_growth, fixef_mortality, currentSpecies, 0.0071,
	100, scalingGrowth, scalingMortality, C0_C1, purves2007_allometries) - 2.119540475517675

survivorship(climate[1, K], climate[1, PP], climate[1, k], climate[1, pp],
	fixef_growth, fixef_mortality, currentSpecies,
	100, 0, scalingGrowth, scalingMortality) - 0.151134147956041
