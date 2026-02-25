This folder contains R code used for evolutionary simulations.

# gen_effect.R
Code to generate the C-matrix, the data matrix containing the effect of the regulators on the phenotypic traits. The output file is being used for all simulations.

# sim_evo.R
Code for evolutionary simulations where all species start from the same expression program and the C-matrix remains constant.

# sim_evo_prior_div.R
Code for evolutionary simulations with different species start from diverse expression programs and the C-matrix remains constant.

# gen_prior_div.R
Code to generate a set of diverged regulator expression programs. The output file will be an input for "sim_evo_prior_div.R".

# sim_evo_evolving_effect.R
Code for evolutionary simulations with all species start from the same expression programs and the C-matrix is evolvable.
