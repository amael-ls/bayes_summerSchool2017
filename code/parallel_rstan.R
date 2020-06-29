
#### Load libraries
# library(rapportools)
library(data.table)
library(magrittr)
library(parallel)
library(ggmcmc)
library(rstan)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

source("./toolFunctions.R")

#### Load data and some work on them
growth_dt = readRDS("../data/growth_dt.rds")

## Rename
names(growth_dt) = rename(growth_dt,
	c("mean_temp_period_3_lag", "tot_annual_pp_lag", "tot_pp_period3_lag",
	"min_temp_coldest_period_lag"),
	c("mean_temp_period_3", "tot_annual_pp", "tot_pp_for_period_3",
	"min_temp_coldest_period"))

## Remove NA
# The first measure is always NA for deltaYear
growth_dt = growth_dt[!is.na(deltaYear),]

# Some climates variables are undefined at some location
growth_dt = growth_dt[!is.na(mean_temp_period_3) & !is.na(tot_pp_for_period_3), ]

## Calculate canopy status (i.e., ppa), and scale between 0 and 1 temp and precip
growth_dt[, ppa := 0]
growth_dt[height > s_star, ppa := 1]

growth_dt[, mean_temp_period_3 := (mean_temp_period_3 - min(mean_temp_period_3))/(max(mean_temp_period_3) - min(mean_temp_period_3))]
growth_dt[, tot_pp_for_period_3 := (tot_pp_for_period_3 - min(tot_pp_for_period_3))/(max(tot_pp_for_period_3) - min(tot_pp_for_period_3))]

# growth_dt = growth_dt[sample(1:nrow(growth_dt), 2e4)]

## Keep only species that have more than 1e4 individuals
table_species = sort(unique(growth_dt[, .(species_id, tree_id)]) %>%
    .[, table(species_id)])
limInf = 1e4
limSup = 3e4
ls_species = sort(names(table_species[(table_species > limInf) & (table_species < limSup)]))

saveRDS(ls_species, "../data/ls_21species.rds")

growth_dt = growth_dt[species_id %in% ls_species]

print(growth_dt[, table(species_id)])

#### Tool function
## Create stan data for each species
makeDataList = function(growth_dt, ls_species)
{
	n = length(ls_species)
	dataStan_ls = vector("list", n)
	for (i in 1:n)
	{
		species = ls_species[i]
		dataStan_ls[[i]] = list(N = growth_dt[species_id == species, .N],
			T_data = growth_dt[species_id == species, mean_temp_period_3],
			P_data = growth_dt[species_id == species, tot_pp_for_period_3],
			C_data = growth_dt[species_id == species, ppa],
			D_data = growth_dt[species_id == species, dbh],
			Y = growth_dt[species_id == species, growth])
	}
	return (dataStan_ls)
}

## Initialisation function for stan prog
initFct = function(chain_id = 1)
{
	xi = 0.8*(1 - 0.8)/0.1 - 1;

	# Priors, hp1 contains their means, hp2 contains their variance
	pdg = rgamma(n = 1, shape = 5^2/100.0, rate = 5/100.0);

	A = rgamma(n = 1, shape = 1^2/1.0, rate = 1/1.0); # A is the location, I guess it will be positive
	B = rgamma(n = 1, shape = 2^2/1.0, rate = 2/1.0); # Scale, kind of "variance"
	C = rgamma(n = 1, shape = 2^2/1.0, rate = 2/1.0);

	P_opt = rgamma(n = 1, shape = 0.5^2/0.5, rate = 0.5/0.5);
	sigmaP_opt = rgamma(n = 1, shape = 400^2/5000.0, rate = 400/5000.0);

	beta = rbeta(n = 1, shape1 = 0.8*xi, shape2 = xi*(1 - 0.8));

	Phi_opt = rgamma(n = 1, shape = 200^2/10000.0, rate = 200/10000.0);
	sigmaPhi_opt = rgamma(n = 1, shape = 10000^2/100000.0, rate = 10000/100000.0);

	sigma_base = rgamma(n = 1, shape = 5^2/10, rate = .5) + 1e-2;
	sigmaT_opt_sig = rgamma(n = 1, shape = 10, rate = 1) + 1e-2; # 10^2/10
	sigmaP_opt_sig = rgamma(n = 1, shape = 250^2/10, rate = 25) + 1e-2;

	init = list(pdg = pdg, A = A, B = B, C = C,
		P_opt = P_opt, sigmaP_opt = sigmaP_opt, beta = beta, Phi_opt = Phi_opt,
		sigmaPhi_opt = sigmaPhi_opt, sigma_base = sigma_base,
		sigmaT_opt_sig = sigmaT_opt_sig, sigmaP_opt_sig = sigmaP_opt_sig)
	return (init)
}

#### Stan model
## Compile
model = stan_model(file = "growth.stan")

## Define stan variables
maxIter = 2e4

nCores = 4
nChains = 4
print(nCores)
print(nChains)

## List for each species
dataStan = makeDataList(growth_dt, ls_species)

fits = mclapply(1:length(ls_species),
	function(i)
		sampling(object = model, data = dataStan[[i]],
			chains = nChains, iter = maxIter,
			cores = nCores,
			control = list(adapt_delta = 0.85, max_treedepth = 11),
			include = FALSE, pars = c("mu_d", "xi", "sigma_duplicate")),
			mc.cores = 4)

saveRDS(fits, "results21species.rds")
