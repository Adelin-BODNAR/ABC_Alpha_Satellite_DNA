# Directory containing the simulated data
setwd("")

#Path to observed data file
observed_data_path = ''

#Simulated data file name
simulated_data_path = ''


observed_data = read.table(observed_data_path)
simulated_data = read.table(simulated_data_path)


simulated_data[["sum_stat_1"]] = rep("monomeric",nrow(simulated_data))
simulated_data[["monomeric"]] = simulated_data$V36
observed_data[["monomeric"]] = observed_data$V2

simulated_data[["sum_stat_2"]] = rep("small_HORs",nrow(simulated_data))
simulated_data[["small_HORs"]] = rowSums(simulated_data[,37:40])
observed_data[["small_HORs"]] = rowSums(observed_data[,3:6])

simulated_data[["sum_stat_3"]] = rep("medium_HORs",nrow(simulated_data))
simulated_data[["medium_HORs"]] = rowSums(simulated_data[,41:50])
observed_data[["medium_HORs"]] = rowSums(observed_data[,7:16])

simulated_data[["sum_stat_4"]] = rep("big_HORs",nrow(simulated_data))
simulated_data[["big_HORs"]] = rowSums(simulated_data[,51:65])
observed_data[["big_HORs"]] = rowSums(observed_data[,17:31])

simulated_data[["sum_stat_5"]] = rep("very_big_HORs",nrow(simulated_data))
simulated_data[["very_big_HORs"]] = rowSums(simulated_data[,66:75])
observed_data[["very_big_HORs"]] = rowSums(observed_data[,32:41])



simulated_data[["sum_stat_6"]] = rep("very_small_diff",nrow(simulated_data))
simulated_data[["very_small_diff"]] = rowSums(simulated_data[,77:82])
observed_data[["very_small_diff"]] = rowSums(observed_data[,43:48])

simulated_data[["sum_stat_7"]] = rep("small_diff",nrow(simulated_data))
simulated_data[["small_diff"]] = rowSums(simulated_data[,83:87])
observed_data[["small_diff"]] = rowSums(observed_data[,49:53])

simulated_data[["sum_stat_8"]] = rep("medium_diff",nrow(simulated_data))
simulated_data[["medium_diff"]] = rowSums(simulated_data[,88:97])
observed_data[["medium_diff"]] = rowSums(observed_data[,54:63])

simulated_data[["sum_stat_9"]] = rep("big_diff",nrow(simulated_data))
simulated_data[["big_diff"]] = rowSums(simulated_data[,98:107])
observed_data[["big_diff"]] = rowSums(observed_data[,64:73])

simulated_data[["sum_stat_10"]] = rep("very_big_diff",nrow(simulated_data))
simulated_data[["very_big_diff"]] = rowSums(simulated_data[,108:177])
observed_data[["very_big_diff"]] = rowSums(observed_data[,74:143])


simulated_data[["sum_stat_11"]] = rep("small_similarity",nrow(simulated_data))
simulated_data[["small_similarity"]] = rowSums(simulated_data[,179:239])
observed_data[["small_similarity"]] = rowSums(observed_data[,145:205])

simulated_data[["sum_stat_12"]] = rep("medium_similarity_1",nrow(simulated_data))
simulated_data[["medium_similarity_1"]] = rowSums(simulated_data[,240:249])
observed_data[["medium_similarity_1"]] = rowSums(observed_data[,206:215])

simulated_data[["sum_stat_13"]] = rep("medium_similarity_2",nrow(simulated_data))
simulated_data[["medium_similarity_2"]] = rowSums(simulated_data[,250:259])
observed_data[["medium_similarity_2"]] = rowSums(observed_data[,216:225])

simulated_data[["sum_stat_14"]] = rep("medium_similarity_3",nrow(simulated_data))
simulated_data[["medium_similarity_3"]] = rowSums(simulated_data[,260:269])
observed_data[["medium_similarity_3"]] = rowSums(observed_data[,226:235])

simulated_data[["sum_stat_15"]] = rep("big_similarity",nrow(simulated_data))
simulated_data[["big_similarity"]] = rowSums(simulated_data[,270:274])
observed_data[["big_similarity"]] = rowSums(observed_data[,236:240])

simulated_data[["sum_stat_16"]] = rep("very_big_similarity",nrow(simulated_data))
simulated_data[["very_big_similarity"]] = rowSums(simulated_data[,275:279])
observed_data[["very_big_similarity"]] = rowSums(observed_data[,241:245])


index_mode_observed_HOR_distrib = which.max(observed_data[1,3:41]) + 2
index_mode_simulated_HOR_distrib = index_mode_observed_HOR_distrib + 34

simulated_data[["sum_stat_17"]] = rep("most_present_HOR",nrow(simulated_data))
simulated_data[["most_present_HOR"]] = simulated_data[,index_mode_simulated_HOR_distrib]
observed_data[["most_present_HOR"]] = observed_data[,index_mode_observed_HOR_distrib]


distances_data = data.frame(

monomeric = (observed_data[1,"monomeric"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$monomeric * 100 / rowSums(simulated_data[1,36:75])),
small_HORs = (observed_data[1,"small_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$small_HORs * 100 / rowSums(simulated_data[1,36:75])),
medium_HORs = (observed_data[1,"medium_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$medium_HORs * 100 / rowSums(simulated_data[1,36:75])),
big_HORs = (observed_data[1,"big_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$big_HORs * 100 / rowSums(simulated_data[1,36:75])),
very_big_HORs = (observed_data[1,"very_big_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$very_big_HORs * 100 / rowSums(simulated_data[1,36:75])),

very_small_diff = (observed_data[1,"very_small_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$small_diff * 100 / rowSums(simulated_data[1,77:177])),
small_diff = (observed_data[1,"small_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$small_diff * 100 / rowSums(simulated_data[1,77:177])),
medium_diff = (observed_data[1,"medium_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$medium_diff * 100 / rowSums(simulated_data[1,77:177])),
big_diff = (observed_data[1,"big_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$big_diff * 100 / rowSums(simulated_data[1,77:177])),
very_big_diff = (observed_data[1,"very_big_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$very_big_diff * 100 / rowSums(simulated_data[1,77:177])),

small_similarity = (observed_data[1,"small_similarity"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$small_similarity * 100 / rowSums(simulated_data[1,179:279])),
medium_similarity_1 = (observed_data[1,"medium_similarity_1"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$medium_similarity_1 * 100 / rowSums(simulated_data[1,179:279])),
medium_similarity_2 = (observed_data[1,"medium_similarity_2"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$medium_similarity_2 * 100 / rowSums(simulated_data[1,179:279])),
medium_similarity_3 = (observed_data[1,"medium_similarity_3"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$medium_similarity_3 * 100 / rowSums(simulated_data[1,179:279])),
big_similarity = (observed_data[1,"big_similarity"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$big_similarity * 100 / rowSums(simulated_data[1,179:279])),
very_big_similarity = (observed_data[1,"very_big_similarity"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$very_big_similarity * 100 / rowSums(simulated_data[1,179:279])),

most_present_HOR = (observed_data[1,"most_present_HOR"] ) - (simulated_data$most_present_HOR )

)


distances_data = distances_data[,1:17] ^ 2

distances_data[,1:17] = as.data.frame(lapply(distances_data[,1:17], function(col){col * sd(col)}))

distances_data$square_sum = rowSums(distances_data[,1:17])

distances_data$euclidean_distance = sqrt(distances_data$square_sum)

simulated_data[["distance_name"]] = rep("euclidean_distance",nrow(simulated_data))
simulated_data$euclidean_distance = distances_data$euclidean_distance

write.table(simulated_data, file = paste(simulated_data_path, ".euclidean_distances.tab", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)


simulated_data_sorted = simulated_data[order(simulated_data$euclidean_distance),]



acceptance_proportion = nrow(simulated_data) * (1 / 100)

accepted_simulated_data = simulated_data_sorted[1:acceptance_proportion,]

write.table(accepted_simulated_data, file = paste(simulated_data_path, ".euclidean_distances.tab.accepted.tab", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)

png(paste(simulated_data_path,".HOR_order_distribution.png", sep = ""), width = 900, height = 600)
barplot(as.numeric(accepted_simulated_data[1,36:75]), main= "Estimation of HOR order distribution for the simulation with lowest euclidean distance \nfor the observed data", xlab = "Estimated order", ylab = "Monomers count", col = "#95B9C4", names.arg = 1:40)
dev.off()

png(paste(simulated_data_path,".max_diff_similarity_distribution.png", sep = ""), width = 900, height = 600)
barplot(as.numeric(accepted_simulated_data[1,77:177]), main= "Estimation of proportion between monomers in a monomeric organisation and monomers in a HOR \nfor the simulation with lowest euclidean distance for the observed data", xlab = "Similarity percentage", ylab = "Monomers count", col = "#95B9C4", names.arg = 0:100)
dev.off()

png(paste(simulated_data_path,".random_similarity_distribution.png", sep = ""), width = 900, height = 600)
barplot(as.numeric(accepted_simulated_data[1,179:279]), main= "Estimation of average similarity distribution for the simulation with lowest euclidean distance \nfor the observed data", xlab = "Similarity percentage", ylab = "Monomers count", col = "#95B9C4", names.arg = 0:100)
dev.off()



png(paste(simulated_data_path,".amplification_rate_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V6, breaks = 250, main= "Posterior distribution for the amplification rate parameter used as \naverage number of amplifications per time unit per monomers in the dataset in the simulations", xlab = "Amplification rate", ylab = "Simulations count", col = "#95B9C4")
dev.off()

png(paste(simulated_data_path,".monomeric_probability_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V8, breaks = 100, main= "Posterior distribution for the parameter used as probability of an amplification event producing a monomeric organization", xlab = "Monomeric probability", ylab = "Simulations count", col = "#95B9C4")
dev.off()

png(paste(simulated_data_path,".alpha_amplification_size_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V14, breaks = 300, main= "Posterior distribution for the alpha parameter of the beta law used to draw amplification sizes", xlab = "Alpha value", ylab = "Simulations count", col = "#95B9C4")
dev.off()

png(paste(simulated_data_path,".beta_amplification_size_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V16, breaks = 300, main= "Posterior distribution for the beta parameter of the beta law used to draw amplification sizes", xlab = "Beta value", ylab = "Simulations count", col = "#95B9C4")
dev.off()

png(paste(simulated_data_path,".alpha_HOR_order_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V22, breaks = 300, main= "Posterior distribution for the alpha parameter of the beta law used to draw HOR orders", xlab = "Alpha value", ylab = "Simulations count", col = "#95B9C4")
dev.off()

png(paste(simulated_data_path,".beta_HOR_order_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V24, breaks = 300, main= "Posterior distribution for the beta parameter of the beta law used to draw HOR orders", xlab = "Beta value", ylab = "Simulations count", col = "#95B9C4")
dev.off()

png(paste(simulated_data_path,".substitution_rate_parameter_distribution.png", sep = ""), width = 900, height = 600)
hist(accepted_simulated_data$V26, breaks = 250, main= "Posterior distribution for the substitution rate parameter used as \naverage number of substitution per nucleotides in the dataset in the simulations", xlab = "Substitution rate", ylab = "Simulations count", col = "#95B9C4")
dev.off()

