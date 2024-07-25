n = 1000000
max_size = 5892

amplification_rate = 1

max_amp_size = as.integer(max_size / 4)

Ti_Tv_ratio = 1.4

nb_orders = 40
nb_sims_per_order = 5


df_params = data.frame(

option_max_size = rep("-s", n),
sample_max_size = rep(max_size, n),

option_amp_rate = rep("--amplification_rate", n),
sample_amp_rate = rep(amplification_rate, n),

option_monomeric_prob = rep("--monomeric_prob", n),
sample_monomeric_prob = runif(n, min= 0, max = 1),

option_sub_rate = rep("--substitution_rate", n),
sample_sub_rate = 200 - runif(n, min=0, max=200),

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

write.table(df_params, "sim_params.tab", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
