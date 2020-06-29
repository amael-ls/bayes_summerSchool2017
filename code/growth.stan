
data
{
	// Size integer
	int<lower = 1> N; // size of response var

	// Vector data
	vector[N] T_data; // temperature, E data
	vector<lower = 0>[N] P_data; // Precipitation, E data
	vector<lower = 0, upper = 1>[N] C_data; // canopy, E data
	vector<lower = 0>[N] D_data; // diameter, I data
	vector<lower = 0>[N] Y; // response var, not 'logarithmised'
}

parameters // IMPORTANT: it worth adding constraints, at least to respect the priors, otherwise, a lot of divergence!
{
	real<lower = 0> pdg; // Potential Diameter Growth

	real A;
	real B;
	real<lower = 0> C;
	real<lower = 0> sigmaT_opt_sig;

	real<lower = 0> P_opt; // Optimum precipitation of each species
	real<lower = 0> sigmaP_opt; // Variance among individuals of optimal P within a species
	real<lower = 0> sigmaP_opt_sig;

	real<lower = 0, upper = 1> beta;

	real<lower = 0> Phi_opt;
	real<lower = 0> sigmaPhi_opt;

	// real<lower = 0> sigma; // Variance of individuals around there species specific mean
	real<lower = 0> sigma_base;
}

transformed parameters
{
	vector[N] mu_d =
			pdg * (C_data + (1 - C_data)*beta) .*
			exp(-(T_data - A)/B) .* exp( (-1 - C) * log( 1 + exp(-(T_data - A)/B) ) ) * C*(1 + 1/C)^(1 + C) .* // Trick: a^b = Exp[b*Log[a]]
			exp(-(P_data - P_opt) .* (P_data - P_opt)/sigmaP_opt^2) .*
			exp(-log(D_data/Phi_opt).*log(D_data/Phi_opt)/sigmaPhi_opt^2);

	real xi = 0.8*(1 - 0.8)/0.1 - 1;

	real T_opt = A + B .* log(C); // calc_Topt(A, B, C);

	// Duplicating sigma to use elementwise division ./ in model
	vector[N] sigma_duplicate = sigma_base * exp(-(T_data - T_opt) .* (T_data - T_opt)/sigmaT_opt_sig^2) .*
	exp(-(P_data - P_opt) .* (P_data - P_opt)/sigmaP_opt_sig^2);
}

model
{
	// Priors, hp1 contains their means, hp2 contains their variance
	pdg ~ gamma(5^2/100.0, 5/100.0);

	A ~ gamma(0.8^2/5, 0.8/5); // A is the location, I guess it will be positive
	B ~ gamma(2^2/1000.0, 2/1000.0); // Scale, kind of "variance"
	C ~ gamma(2^2/10000.0, 2/10000.0); // C is the shape, must be positive

	P_opt ~ lognormal(log(0.4^2/sqrt(10 + 0.4^2)), sqrt(log(1 + 10/0.4^2))); // Log[μ²/Sqrt[σ² + μ²]], Sqrt[Log[1 + σ²/μ²]]
	sigmaP_opt ~ gamma(100^2/5000.0, 100/5000.0);

	beta ~ beta(0.8*xi, xi*(1 - 0.8));

	Phi_opt ~ gamma(200^2/10000.0, 200/10000.0);
	sigmaPhi_opt ~ gamma(300^2/100000.0, 300/100000.0);

	sigma_base ~ pareto_type_2(0.01, 10, 2.4); // growth_dt[, var(growth)] = 2.6
	sigmaT_opt_sig ~ pareto_type_2(0.01, 1000, 2.4);
	sigmaP_opt_sig ~ pareto_type_2(0.01, 1000, 2.4);

	// Growth model
	Y ~ gamma(mu_d .* mu_d ./ sigma_duplicate, mu_d ./ sigma_duplicate);
}
