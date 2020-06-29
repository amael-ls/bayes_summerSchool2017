## Required libraries
library(data.table)
library(stringi)

## Tool functions
# Print range of dataframe columns
range_df = function(df, cols = 1:ncol(df), na.rm = TRUE)
{
	n = length(cols)
	range = data.table(nameCol = character(n), min = numeric(n), max = numeric(n),
		mean = numeric(n), var = numeric(n))
	if (sum(class(df) == "data.frame") == 0) # because a dt is a df
	{
		for(i in 1:n)
		{
			if ( is.factor(df[, i]) | is.character(df[, i]) )
			{
				print(paste0(names(df)[i], " is a factor/string"))
				range$nameCol[i] = names(df)[i]
				range$min[i] = NA
				range$max[i] = NA
				range$mean[i] = NA
				range$var[i] = NA
				next
			}
			range[i, nameCol := names(df)[i]]
			range[i, min := min(df[, i], na.rm = na.rm)]
			range[i, max := max(df[, i], na.rm = na.rm)]
			range[i, mean := mean(df[, i], na.rm = na.rm)]
			range[i, var := sd(df[, i], na.rm = na.rm)^2]
		}
	}
	if (is.data.table(df))
	{
		for(i in 1:n)
		{
			if ( is.factor(df[[i]]) | is.character(df[[i]]) )
			{
				print(paste0(names(df)[i], " is a factor/string"))
				range$nameCol[i] = names(df)[i]
				range$min[i] = NA
				range$max[i] = NA
				range$mean[i] = NA
				range$var[i] = NA
				next
			}
			range[i, nameCol := names(df)[i]]
			range[i, min := min(df[[i]], na.rm = na.rm)]
			range[i, max := max(df[[i]], na.rm = na.rm)]
			range[i, mean := mean(df[[i]], na.rm = na.rm)]
			range[i, var := sd(df[[i]], na.rm = na.rm)^2]
		}
		return (range[])
	}

	return (range)
}

# Normalised (rescale) columns of a dataframe and save the coefficients in a file
normalisation = function(df, colnames = names(df), filename = "normalisation.rds", rm_na = TRUE)
{
	if (rm_na)
		print("Warning: rm_na is activated, normalisation won't take NA into account")

	mu_sd = data.table(var = character(length(colnames)), mu = numeric(length(colnames)), sd = numeric(length(colnames)))
	i = 1
	if (is.data.frame(df) & !is.data.table(df))
	{
		for (col in colnames)
		{
			print(paste0("calculus for column ", col))
			mu_sd[i, c("var", "mu", "sd") := .(col, mean(df[[col]], na.rm = rm_na), sd(df[[col]], na.rm = rm_na))]
			df[col] = (df[[col]] - mu)/sd
			i = i + 1
		}
	}
	if (is.data.table(df))
	{
		mu_sd[, c("var", "mu", "sd") := .(colnames, as.matrix(df[, lapply(.SD, mean, na.rm = rm_na), .SDcols = colnames])[1,], as.matrix(df[, lapply(.SD, sd, na.rm = rm_na), .SDcols = colnames])[1,])]
		
		df[, c(colnames) := lapply(.SD, function(x) (x - mean(x, na.rm = rm_na))/sd(x, na.rm = rm_na)), .SDcols = colnames]
	}

	print("done")
	ext_file = stri_sub(filename, from =  stri_locate_last(filename, fixed = ".")[1])
	if (ext_file == ".rds")
		saveRDS(mu_sd, file = filename)

	if (ext_file == ".csv")
		write.table(mu_sd, file = filename, col.names = TRUE, row.names = FALSE, quote = TRUE)
	print(paste0("files containing coefficients saved at: ", filename))
	return (df)
}

# To rescale a dataframe
rescaling = function(df, scaling)
{
	cols = scaling$var

	if ( !(FALSE %in% (cols %in% names(df))) )
	{
		for (col in cols)
			df[[col]] = (scaling[scaling$var == col,]$sd)*df[[col]] + scaling[scaling$var == col,]$mu
	}
	else
	{
		print("could not rescale, problem of cols names. Check out the following variable(s) in scaling file:")
		print(cols[!(cols %in% names(df))])
	}

	return (df)
}

# 3-D hist plot
# x and y are 2 vectors (intervals [x_min, x_max], [y_min, y_max])
# data on which we calculate the frequences
# nx, number of breaking points in x = [x_min, x_max]
# ny, number of breaking points in y = [y_min, y_max]
# var_x and var_y std::string type, denote the variables from the data
myHist = function(x, y, data, nx, ny, var_x, var_y)
{
	if ( (x[2] < x[1]) | (y[2] < y[1]) )
		stop("*** Error, check x and y args")

	if ( (nx < 2) | (ny < 2) )
		stop("*** Error, check nx and ny args (must be bigger than two)")

	int_x = seq(x[1], x[2], length.out = nx)
	int_y = seq(y[1], y[2], length.out = ny)

	freq = matrix(nrow = nx - 1, ncol = ny - 1)

	# 'inside' (i.e., no left and top boundaries)
	for (i in 1:(nx - 2))
	{
		index_x = which( (int_x[i] <= data[[var_x]]) & (data[[var_x]] < int_x[i + 1]))
		subset_int_x = data[index_x,]
		for (j in 1:(ny - 2))
		{
			index_y = which( (int_y[j] <= subset_int_x[[var_y]]) & (subset_int_x[[var_y]] < int_y[j + 1]))
			subset_int_x_y = subset_int_x[index_y,]
			freq[j, i] = dim(subset_int_x_y)[1]
		}
	}

	# x boundary condition
	i = nx - 1
	index_x = which( (int_x[i] <= data[[var_x]]) & (data[[var_x]] <= int_x[i + 1]) )
	subset_int_x = data[index_x,]
	for (j in 1:(ny - 2))
	{
		index_y = which( (int_y[j] <= subset_int_x[[var_y]]) & (subset_int_x[[var_y]] < int_y[j + 1]) )
		subset_int_x_y = subset_int_x[index_y,]
		freq[j, i] = dim(subset_int_x_y)[1]
	}

	# y boundary condition
	j = ny - 1
	index_y = which( (int_y[j] <= data[[var_y]]) & (data[[var_y]] <= int_y[j + 1]) )
	subset_int_y = data[index_y,]
	for (i in 1:(nx - 2))
	{
		index_x = which( (int_x[i] <= subset_int_y[[var_x]]) & (subset_int_y[[var_x]] < int_x[i + 1]))
		subset_int_x_y = subset_int_y[index_x,]
		freq[j, i] = dim(subset_int_x_y)[1]
	}

	# x-y top left corner
	i = nx - 1
	j = ny - 1

	freq[j, i] = dim(data)[1] - sum(freq[!is.na(freq)])

	if (TRUE %in% (freq < 0))
		stop("*** error function, negative frequences... This means I am stupid!")

	output = list(int_x = int_x[1:(nx - 1)], int_y = int_y[1:(ny - 1)], freq = freq)
	return (output)
}

## To shape averageClim, to be readible for growth.R program (needs same format as mostFreqClimCombination())
reshapeAverageClim = function(temp_vars, pp_vars, file = "./averageClim_growth.rds") # temp/pp_vars from run_para_model for all combinations
{
	selectedClim = data.frame(temp = temp_vars, temp_values = numeric(4),
							  pp = pp_vars, pp_values = numeric(4))
	averageClim = readRDS(file)
	count = 1
	for (temp in temp_vars)
	{
		for (pp in pp_vars)
		{
			selectedClim[count, "temp_values"] = averageClim[temp]
			selectedClim[count, "pp_values"] = averageClim[pp]
			row.names(selectedClim)[count] = paste0(temp, "---", pp)
			count = count + 1
		}
	}
	return (selectedClim)
}

## To create species specific climatic gradient
gradClimSpeciesSpecific = function(spanSpecies, species, climVar, length = 150)
{
	minRow = paste(species, "min", sep = "_")
	maxRow = paste(species, "max", sep = "_")

	minVal = spanSpecies[minRow, climVar]
	maxVal = spanSpecies[maxRow, climVar]

	climGradient = seq(from = minVal, maxVal, length.out = length)

	return(climGradient);
}

otherClimValSpeciesSpecific = function(spanSpecies, species, climVar)
{
	meanRow = paste(species, "mean", sep = "_")

	return(spanSpecies[meanRow, climVar]);
}

relativeAbundance = function(df)
{
	if ("id_spe" %in% names(df))
	{
		df = unique(subset(df, select = c("tree_id", "id_spe")))
		ls_species = unique(df$id_spe)
	}
	if ("species_id" %in% names(df))
	{
		df = unique(subset(df, select = c("tree_id", "species_id")))
		ls_species= unique(df$species_id)
	}
	nbIndiv = dim(df)[1]
	relAb = data.frame(species = ls_species, relativeAbundance = numeric(length(ls_species)), stringsAsFactors = FALSE)
	rownames(relAb)
	for (species in ls_species)
	{
		subDf = subset(df)
		relAb[]
	}
}

# To change some names automatically. Needs to be in the same order
rename = function(df, colNamesToChange, newColNames)
{
	if (length(colNamesToChange) != length(newColNames))
		stop("*** Error (from rename) ***: colNamesToChange and newColNames should be the same size")
	if (FALSE %in% (colNamesToChange %in% names(df)))
		warning("*** Warning (from rename) ***: there are some unknown colNames")

	names(df)[names(df) %in% colNamesToChange] = newColNames

	return (names(df))
}
