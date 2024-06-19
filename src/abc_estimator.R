

setwd("/mnt/data/Stage_M2_Adelin_2024/Data/runs_abc/run_24_06_12_chr_21/sum_stats")

observed_data_path = '/mnt/data/Stage_M2_Adelin_2024/Data/Human_monomers_reoriented/chromosome_21_CP068257_all.fst.aligned_mafft.fst.sum_stats.tab'
simulated_data_path = 'chromosome_21_CP068257.2_all.fst.aligned_mafft.fst.no_gap.fst.hmm.consensus.fst.simulation_parameters_and_data_description_run_24_06_12.tab'

observed_data = read.table(observed_data_path)

simulated_data = read.table(simulated_data_path)


simulated_data[["sum_stat_1"]] = rep("monomeric",nrow(simulated_data))
simulated_data[["monomeric"]] = simulated_data$V34
observed_data[["monomeric"]] = observed_data$V2

simulated_data[["sum_stat_2"]] = rep("small_HORs",nrow(simulated_data))
simulated_data[["small_HORs"]] = rowSums(simulated_data[,35:38])
observed_data[["small_HORs"]] = rowSums(observed_data[,3:6])

simulated_data[["sum_stat_3"]] = rep("medium_HORs",nrow(simulated_data))
simulated_data[["medium_HORs"]] = rowSums(simulated_data[,39:48])
observed_data[["medium_HORs"]] = rowSums(observed_data[,7:16])

simulated_data[["sum_stat_4"]] = rep("big_HORs",nrow(simulated_data))
simulated_data[["big_HORs"]] = rowSums(simulated_data[,49:63])
observed_data[["big_HORs"]] = rowSums(observed_data[,17:31])

simulated_data[["sum_stat_5"]] = rep("very_big_HORs",nrow(simulated_data))
simulated_data[["very_big_HORs"]] = rowSums(simulated_data[,64:73])
observed_data[["very_big_HORs"]] = rowSums(observed_data[,32:41])



simulated_data[["sum_stat_6"]] = rep("small_diff",nrow(simulated_data))
simulated_data[["small_diff"]] = rowSums(simulated_data[,75:80])
observed_data[["small_diff"]] = rowSums(observed_data[,43:48])

simulated_data[["sum_stat_7"]] = rep("medium_diff",nrow(simulated_data))
simulated_data[["medium_diff"]] = rowSums(simulated_data[,81:85])
observed_data[["medium_diff"]] = rowSums(observed_data[,49:53])

simulated_data[["sum_stat_8"]] = rep("big_diff",nrow(simulated_data))
simulated_data[["big_diff"]] = rowSums(simulated_data[,86:105])
observed_data[["big_diff"]] = rowSums(observed_data[,54:73])

simulated_data[["sum_stat_9"]] = rep("very_big_diff",nrow(simulated_data))
simulated_data[["very_big_diff"]] = rowSums(simulated_data[,106:175])
observed_data[["very_big_diff"]] = rowSums(observed_data[,74:143])


simulated_data[["sum_stat_10"]] = rep("small_similarity",nrow(simulated_data))
simulated_data[["small_similarity"]] = rowSums(simulated_data[,177:236])
observed_data[["small_similarity"]] = rowSums(observed_data[,145:204])

simulated_data[["sum_stat_11"]] = rep("medium_similarity",nrow(simulated_data))
simulated_data[["medium_similarity"]] = rowSums(simulated_data[,237:257])
observed_data[["medium_similarity"]] = rowSums(observed_data[,205:225])

simulated_data[["sum_stat_12"]] = rep("big_similarity",nrow(simulated_data))
simulated_data[["big_similarity"]] = rowSums(simulated_data[,258:271])
observed_data[["big_similarity"]] = rowSums(observed_data[,226:239])

simulated_data[["sum_stat_13"]] = rep("very_big_similarity",nrow(simulated_data))
simulated_data[["very_big_similarity"]] = rowSums(simulated_data[,272:277])
observed_data[["very_big_similarity"]] = rowSums(observed_data[,240:245])


distances_data = data.frame(

monomeric = (observed_data[1,"monomeric"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$monomeric * 100 / rowSums(simulated_data[1,34:73])),
small_HORs = (observed_data[1,"small_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$small_HORs * 100 / rowSums(simulated_data[1,34:73])),
medium_HORs = (observed_data[1,"medium_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$medium_HORs * 100 / rowSums(simulated_data[1,34:73])),
big_HORs = (observed_data[1,"big_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$big_HORs * 100 / rowSums(simulated_data[1,34:73])),
very_big_HORs = (observed_data[1,"very_big_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$very_big_HORs * 100 / rowSums(simulated_data[1,34:73])),

small_diff = (observed_data[1,"small_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$small_diff * 100 / rowSums(simulated_data[1,75:175])),
medium_diff = (observed_data[1,"medium_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$medium_diff * 100 / rowSums(simulated_data[1,75:175])),
big_diff = (observed_data[1,"big_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$big_diff * 100 / rowSums(simulated_data[1,75:175])),
very_big_diff = (observed_data[1,"very_big_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$very_big_diff * 100 / rowSums(simulated_data[1,75:175])),

small_similarity = (observed_data[1,"small_similarity"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$small_similarity * 100 / rowSums(simulated_data[1,177:277])),
medium_similarity = (observed_data[1,"medium_similarity"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$medium_similarity * 100 / rowSums(simulated_data[1,177:277])),
big_similarity = (observed_data[1,"big_similarity"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$big_similarity * 100 / rowSums(simulated_data[1,177:277])),
very_big_similarity = (observed_data[1,"very_big_similarity"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$very_big_similarity * 100 / rowSums(simulated_data[1,177:277]))

)


distances_data$square_sum = rowSums(distances_data ^ 2)

distances_data$euclidean_distance = sqrt(distances_data$square_sum)

simulated_data[["distance_name"]] = rep("euclidean_distance",nrow(simulated_data))
simulated_data$euclidean_distance = distances_data$euclidean_distance


write.table(simulated_data, file = paste(simulated_data_path, ".euclidean_distances.tab", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)


simulated_data_sorted = simulated_data[order(simulated_data$euclidean_distance),]



acceptance_proportion = nrow(simulated_data) * (1 / 100)

accepted_simulated_data = simulated_data_sorted[1:acceptance_proportion,]

write.table(accepted_simulated_data, file = paste(simulated_data_path, ".euclidean_distances.tab.accepted.tab", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)

png(paste(simulated_data_path,".HOR_order_distribution.png", sep = ""), width = 900, height = 600)
barplot(as.numeric(accepted_simulated_data[1,34:73]), main= "Estimation of HOR order distribution for the simulation with lowest euclidean distance \nfor the observed data", xlab = "Estimated order", ylab = "Monomers count", col = "#95B9C4", names.arg = 1:40)
dev.off()

png(paste(simulated_data_path,".max_diff_similarity_distribution.png", sep = ""), width = 900, height = 600)
barplot(as.numeric(accepted_simulated_data[1,75:175]), main= "Estimation of proportion between monomers in a monomeric organisation and monomers in a HOR \nfor the simulation with lowest euclidean distance for the observed data", xlab = "Similarity percentage", ylab = "Monomers count", col = "#95B9C4", names.arg = 0:100)
dev.off()

png(paste(simulated_data_path,".random_similarity_distribution.png", sep = ""), width = 900, height = 600)
barplot(as.numeric(accepted_simulated_data[1,177:277]), main= "Estimation of average similarity distribution for the simulation with lowest euclidean distance \nfor the observed data", xlab = "Similarity percentage", ylab = "Monomers count", col = "#95B9C4", names.arg = 0:100)
dev.off()







png(paste(simulated_data_path,".alpha_amplification_size_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V12, breaks = 300, main= "Posterior distribution for the alpha parameter of the beta law used to draw amplification sizes", xlab = "Alpha value", ylab = "Simulations count", col = "#95B9C4")
dev.off()

png(paste(simulated_data_path,".beta_amplification_size_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V14, breaks = 300, main= "Posterior distribution for the beta parameter of the beta law used to draw amplification sizes", xlab = "Beta value", ylab = "Simulations count", col = "#95B9C4")
dev.off()

png(paste(simulated_data_path,".alpha_HOR_order_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V20, breaks = 300, main= "Posterior distribution for the alpha parameter of the beta law used to draw HOR orders", xlab = "Alpha value", ylab = "Simulations count", col = "#95B9C4")
dev.off()

png(paste(simulated_data_path,".beta_HOR_order_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V22, breaks = 300, main= "Posterior distribution for the beta parameter of the beta law used to draw HOR orders", xlab = "Beta value", ylab = "Simulations count", col = "#95B9C4")
dev.off()

png(paste(simulated_data_path,".substitution_rate_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V24, breaks = 250, main= "Posterior distribution for the substitution rate parameter used as \naverage number of substitution per nucleotides in the dataset in the simulations", xlab = "Substitution rate", ylab = "Simulations count", col = "#95B9C4")
dev.off()

