#Number of simulations for the run
n = 1000000

#Size aimed for during simulations (= number of monomers in chromosome studied)
max_size = 5892

#Substitution parameter fixed value
substitution_rate = 1

#Max size of amplification events (1/4 of the size aimed for because we suppose a simulation with under 4 amplification events is not right for sure)
max_amp_size = as.integer(max_size / 4)

#Fixed value for the Transition/Transversion ratio
Ti_Tv_ratio = 1.4

#The maximum order for HOR generated and search for in the simulation and data description scripts 
nb_orders = 40

#Number of similarities used to calculate an average similarity value in the data description script
nb_sims_per_order = 5

#Dataframe used to generate and store the options and parameter values for each simulation we want to launch
#For each parameter one column with the option and one column with values drawn randomly using a uniform for each simulation 
#(runif(n, min = 0, max = x) -> n values from the [0,x[ range , x - runif(n, min = 0, max = x) -> n values from the ]0,x] range)
df_params = data.frame(

option_max_size = rep("-s", n),
sample_max_size = rep(max_size, n),

option_amp_rate = rep("--amplification_rate", n),
sample_amp_rate = runif(n, min= 0.0001, max = 1),

option_monomeric_prob = rep("--monomeric_prob", n),
sample_monomeric_prob = runif(n, min= 0, max = 1),

option_sub_rate = rep("--substitution_rate", n),
sample_sub_rate = rep(substitution_rate, n),

option_alpha_amp = rep("--alpha_amplification_size", n),
sample_alpha_amp = 50 - runif(n, min= 0, max = 50),

option_beta_amp = rep("--beta_amplification_size", n),
sample_beta_amp = runif(n, min = 1, max = 50),

option_min_amp_size = rep("--min_amplification_size", n),
sample_min_amp_size = rep(1, n),

option_max_amp_size = rep("--max_amplification_size", n),
sample_max_amp_size = rep(max_amp_size, n),

option_alpha_HOR = rep("--alpha_HOR_order", n),
sample_alpha_HOR = 50 - runif(n, min= 0, max = 50),

option_beta_HOR = rep("--beta_HOR_order", n),
sample_beta_HOR = runif(n, min = 1, max = 50),

option_min_HOR_order = rep("--min_HOR_order", n),
sample_min_HOR_order = rep(1, n),

option_max_HOR_order = rep("--max_HOR_order", n),
sample_max_HOR_order = rep(nb_orders, n),

option_Ti_Tv_ratio = rep("--transition_transversion_ratio", n),
sample_Ti_Tv_ratio = rep(Ti_Tv_ratio, n),

option_orders_tested = rep("--nb_orders_tested", n),
sample_orders_tested = rep(nb_orders, n)
)

#Writes in a file the line sections containing the options and parameter values for each simulation that will be added to the line launching the simulation
write.table(df_params, "sim_params.tab", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
