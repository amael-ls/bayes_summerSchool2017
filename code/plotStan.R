
#### Aims of prog:

#### Load libraries and clear memory
# library(rapportools)
# library(data.table)
library(bayesplot)
# library(parallel)
# library(stringi)
# library(ggmcmc)
library(rstan)
library(coda)
# library(sp)

rm(list = ls())
graphics.off()
options(max.print = 500)

#### Load data, create folders, and define variables
## Load
output = readRDS("./results21species.rds")
ls_species = readRDS("../data/ls_species.rds")
(n = length(output))

## Plot results folder
ifelse(!dir.exists("../plot_results"), dir.create("../plot_results"), FALSE)

# color_scheme_set("darkgray")
source("toSource.R")

for (i in 1:5)
{
	species = ls_species[i]
	ifelse(!dir.exists(paste0("../plot_results/", species)), dir.create(paste0("../plot_results/", species)), FALSE)

	sp_fit = output[[i]]

	# Extract from sp_fit quantities needed for plotting mcmc_nuts_*
	log_p = log_posterior(sp_fit)
	params = nuts_params(sp_fit)

	# combinedchains = mcmc.list(chain, chain2)
	# param_samples = as.data.frame(sp_fit)[,c('beta', 'pdg')]
	# plot(as.mcmc(param_samples))
	#
	# x = as.mcmc(param_samples)
	# y = as.mcmc.list(x)
	# coda::gelman.diag(x)

	# Acceptance rate
	jpeg(paste0("../plot_results/", species, "/acceptance.jpg"))
	mcmc_nuts_acceptance(x = params, lp = log_p)
	dev.off()

	# Divergence
	jpeg(paste0("../plot_results/", species, "/divergence.jpg"))
	mcmc_nuts_divergence(x = params, lp = log_p)
	dev.off()

	# Energy
	jpeg(paste0("../plot_results/", species, "/energy.jpg"))
	mcmc_nuts_energy(x = params)
	dev.off()

	# Stepsize
	jpeg(paste0("../plot_results/", species, "/stepsize.jpg"))
	mcmc_nuts_stepsize(x = params, lp = log_p)
	dev.off()

	# Treedepth
	jpeg(paste0("../plot_results/", species, "/treedepth.jpg"))
	mcmc_nuts_treedepth(x = params, lp = log_p)
	dev.off()

	posterior = as.array(sp_fit)
	dim(posterior) # iterations, chains, 13 parameters + lp__

	# jpeg(paste0("../plot_results/", species, "/pairs.jpg"))
	# color_scheme_set("brightblue")
	# mcmc_pairs(posterior, # pars = c(),
	# 	off_diag_args = list(size = 1, alpha = 1/3),
	# 	condition = pairs_condition(nuts = "accept_stat__"),
	# 	np = params)
	# dev.off()

	jpeg(paste0("../plot_results/", species, "/rhat.jpg"))
	mcmc_rhat(rhat = rhat(sp_fit)) + yaxis_text(hjust = 0)
	dev.off()

	ratios = neff_ratio(sp_fit)

	jpeg(paste0("../plot_results/", species, "/ratio.jpg"))
	mcmc_neff(ratios, size = 2)
	dev.off()

	group = list(others = c("pdg", "beta", "sigma_base"),
		dbh = c("Phi_opt", "sigmaPhi_opt"),
		temperature = c("A", "B", "C", "T_opt", "sigmaT_opt_sig"),
		precipitation = c("P_opt", "sigmaP_opt", "sigmaP_opt_sig"))

	jpeg(paste0("../plot_results/", species, "/autocorel_others.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_acf(posterior, pars = c("pdg", "beta", "sigma_base"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/trace_others.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_trace(posterior, facet_args = list(ncol = 1, strip.position = "left"),
		pars = c("pdg", "beta", "sigma_base"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/dens_overlay_others.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_dens_overlay(posterior, pars = c("pdg", "beta", "sigma_base"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/autocorel_dbh.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_acf(posterior, pars = c("Phi_opt", "sigmaPhi_opt"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/trace_dbh.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_trace(posterior, facet_args = list(ncol = 1, strip.position = "left"),
		pars = c("Phi_opt", "sigmaPhi_opt"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/dens_overlay_dbh.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_dens_overlay(posterior, pars = c("Phi_opt", "sigmaPhi_opt"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/autocorel_temperature.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_acf(posterior, pars = c("A", "B", "C", "T_opt", "sigmaT_opt_sig"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/trace_temperature.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_trace(posterior, facet_args = list(ncol = 1, strip.position = "left"),
		pars = c("A", "B", "C", "T_opt", "sigmaT_opt_sig"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/dens_overlay_temperature.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_dens_overlay(posterior, pars = c("A", "B", "C", "T_opt", "sigmaT_opt_sig"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/autocorel_precipitation.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_acf(posterior, pars = c("P_opt", "sigmaP_opt", "sigmaP_opt_sig"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/trace_precipitation.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_trace(posterior, facet_args = list(ncol = 1, strip.position = "left"),
		pars = c("P_opt", "sigmaP_opt", "sigmaP_opt_sig"))
	dev.off()

	jpeg(paste0("../plot_results/", species, "/dens_overlay_precipitation.jpg"), quality = 100,
		width = 1000, height = 1000)
	mcmc_dens_overlay(posterior, pars = c("P_opt", "sigmaP_opt", "sigmaP_opt_sig"))
	dev.off()

	# for (g in 1:length(group))
	# {
	# 	print(names(group)[g])
	# 	jpeg(paste0("../plot_results/", species, "/autocorel_", names(group)[g], ".jpg"), quality = 100,
	# 		width = 1000, height = 1000)
	# 	mcmc_acf(posterior, pars = group[[g]])
	# 	dev.off()
	#
	# 	jpeg(paste0("../plot_results/", species, "/trace_", names(group)[g], ".jpg"), quality = 100,
	# 		width = 1000, height = 1000)
	# 	mcmc_trace(posterior, facet_args = list(ncol = 1, strip.position = "left"), pars = group[[g]])
	# 	dev.off()
	#
	# 	jpeg(paste0("../plot_results/", species, "/dens_overlay_", names(group)[g], ".jpg"), quality = 100,
	# 		width = 1000, height = 1000)
	# 	mcmc_dens_overlay(posterior, pars = group[[g]])
	# 	dev.off()
	# }
}

summerWine = summary(out)
saveRDS(summerWine, paste0("./summary_", sp_load, "_", outputDate, ".rds"))

monitor(extract(out, permuted = FALSE, inc_warmup = TRUE))
monitor(extract(out, permuted = FALSE, inc_warmup = FALSE))

#### Plots
parsToPlot = c("hp1", "hp2", "pdg", "A", "B", "C", "P_opt", "T_opt", "beta", "Phi_opt", "sigma_base")
dir.create(path = paste0("./", outputDate, "_", sp_load))

for (param in parsToPlot)
{
	out_ggs = ggs(S = out, family = param)
	ggmcmc(out_ggs, file = paste0("./", outputDate, "_", sp_load, "/ggmcmc_", param, ".pdf"))
}


lp_cp <- log_posterior(out)
np_cp <- nuts_params(out)
posterior_cp <- as.array(out)

color_scheme_set("darkgray")
mcmc_parcoord(posterior_cp, np = np_cp)
mcmc_pairs(posterior_cp, np = np_cp) # , pars = c("mu","tau","theta[1]")

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_cp, np = np_cp) +
  xlab("Post-warmup iteration")
