#TODO: Find an easier way (maybe using arguments) to set the paths to the files (and maybe weights for sum stats)

# Directory containing the simulated data
setwd("/mnt/data/Stage_M2_Adelin_2024/Data/runs_abc/run_24_07_18_chr_2/sum_stats")

#Path to observed data file
observed_data_path = '/mnt/data/Stage_M2_Adelin_2024/Data/Human_monomers_reoriented/chromosome_2_CP068276.2_all.fst.aligned_mafft.fst.sum_stats.tab'

#Simulated data file name
simulated_data_path = 'chromosome_2_CP068276.2_all.fst.aligned_mafft.fst.no_gap.fst.hmm.consensus.fst.simulation_parameters_and_data_description_run_24_07_18_chr_2.tab'

#Reads data from observed data file
observed_data = read.table(observed_data_path)

#Reads data from simulated data file
simulated_data = read.table(simulated_data_path)

#Creates classes for the HOR order distributions

##Class for the number of monomers in a monomeric organization
simulated_data[["sum_stat_1"]] = rep("monomeric",nrow(simulated_data))
simulated_data[["monomeric"]] = simulated_data$V36
observed_data[["monomeric"]] = observed_data$V2

##Class for the small HOR orders (2 to 5)
simulated_data[["sum_stat_2"]] = rep("small_HORs",nrow(simulated_data))
simulated_data[["small_HORs"]] = rowSums(simulated_data[,37:40])
observed_data[["small_HORs"]] = rowSums(observed_data[,3:6])

##Class for the medium sized HOR orders (6 to 15)
simulated_data[["sum_stat_3"]] = rep("medium_HORs",nrow(simulated_data))
simulated_data[["medium_HORs"]] = rowSums(simulated_data[,41:50])
observed_data[["medium_HORs"]] = rowSums(observed_data[,7:16])

##Class for the big HOR orders (16 to 30)
simulated_data[["sum_stat_4"]] = rep("big_HORs",nrow(simulated_data))
simulated_data[["big_HORs"]] = rowSums(simulated_data[,51:65])
observed_data[["big_HORs"]] = rowSums(observed_data[,17:31])

##Class for the biggest HOR orders available (31 to 40)
simulated_data[["sum_stat_5"]] = rep("very_big_HORs",nrow(simulated_data))
simulated_data[["very_big_HORs"]] = rowSums(simulated_data[,66:75])
observed_data[["very_big_HORs"]] = rowSums(observed_data[,32:41])


#Creates classes for the similarity difference distributions

##Class for the smallest similarity differences available (0 to 5)
simulated_data[["sum_stat_6"]] = rep("very_small_diff",nrow(simulated_data))
simulated_data[["very_small_diff"]] = rowSums(simulated_data[,77:82])
observed_data[["very_small_diff"]] = rowSums(observed_data[,43:48])

##Class for the small similarity differences (6 to 10)
simulated_data[["sum_stat_7"]] = rep("small_diff",nrow(simulated_data))
simulated_data[["small_diff"]] = rowSums(simulated_data[,83:87])
observed_data[["small_diff"]] = rowSums(observed_data[,49:53])

##Class for the medium sized similarity differences (11 to 20)
simulated_data[["sum_stat_8"]] = rep("medium_diff",nrow(simulated_data))
simulated_data[["medium_diff"]] = rowSums(simulated_data[,88:97])
observed_data[["medium_diff"]] = rowSums(observed_data[,54:63])

##Class for the big similarity differences (21 to 30)
simulated_data[["sum_stat_9"]] = rep("big_diff",nrow(simulated_data))
simulated_data[["big_diff"]] = rowSums(simulated_data[,98:107])
observed_data[["big_diff"]] = rowSums(observed_data[,64:73])

##Class for the biggest similarity differences available (31 to 100)
simulated_data[["sum_stat_10"]] = rep("very_big_diff",nrow(simulated_data))
simulated_data[["very_big_diff"]] = rowSums(simulated_data[,108:177])
observed_data[["very_big_diff"]] = rowSums(observed_data[,74:143])


#Creates classes for the random similarity distributions

##Class for the small similarities (0 to 60)
simulated_data[["sum_stat_11"]] = rep("small_similarity",nrow(simulated_data))
simulated_data[["small_similarity"]] = rowSums(simulated_data[,179:239])
observed_data[["small_similarity"]] = rowSums(observed_data[,145:205])

##Class for the first 1/3 of the medium sized similarities (61 to 70)
simulated_data[["sum_stat_12"]] = rep("medium_similarity_1",nrow(simulated_data))
simulated_data[["medium_similarity_1"]] = rowSums(simulated_data[,240:249])
observed_data[["medium_similarity_1"]] = rowSums(observed_data[,206:215])

##Class for the second 1/3 of the medium sized similarities (71 to 80)
simulated_data[["sum_stat_13"]] = rep("medium_similarity_2",nrow(simulated_data))
simulated_data[["medium_similarity_2"]] = rowSums(simulated_data[,250:259])
observed_data[["medium_similarity_2"]] = rowSums(observed_data[,216:225])

##Class for the last 1/3 of the medium sized similarities (81 to 90)
simulated_data[["sum_stat_14"]] = rep("medium_similarity_3",nrow(simulated_data))
simulated_data[["medium_similarity_3"]] = rowSums(simulated_data[,260:269])
observed_data[["medium_similarity_3"]] = rowSums(observed_data[,226:235])

##Class for the big similarities (91 to 95)
simulated_data[["sum_stat_15"]] = rep("big_similarity",nrow(simulated_data))
simulated_data[["big_similarity"]] = rowSums(simulated_data[,270:274])
observed_data[["big_similarity"]] = rowSums(observed_data[,236:240])

##Class for the biggest similarities available (96 to 100)
simulated_data[["sum_stat_16"]] = rep("very_big_similarity",nrow(simulated_data))
simulated_data[["very_big_similarity"]] = rowSums(simulated_data[,275:279])
observed_data[["very_big_similarity"]] = rowSums(observed_data[,241:245])

#TODO: Can be removed entirely because the summary statistic doesn't seem to give good results (results forced to fit) #########################

#Gets the indexes in the observed and simulated data of the most represented HOR order in observed data
index_mode_observed_HOR_distrib = which.max(observed_data[1,3:41]) + 2
index_mode_simulated_HOR_distrib = index_mode_observed_HOR_distrib + 34

#Gets the values for the order chosen
simulated_data[["sum_stat_17"]] = rep("most_present_HOR",nrow(simulated_data))
simulated_data[["most_present_HOR"]] = simulated_data[,index_mode_simulated_HOR_distrib]
observed_data[["most_present_HOR"]] = observed_data[,index_mode_observed_HOR_distrib]

###############################################################################################################################################

#Calculates the average values for each beta distribution to plot as posterior distribution instead of their alpha and beta parameters
simulated_data$amp_mean = ((simulated_data$V14 / (simulated_data$V14 + simulated_data$V16)) * (simulated_data$V12 - simulated_data$V10)) + simulated_data$V10
simulated_data$HOR_mean = ((simulated_data$V22 / (simulated_data$V22 + simulated_data$V24)) * (simulated_data$V20 - simulated_data$V18)) + simulated_data$V18

#Creates a dataframe containing the distances for each summary statistic
#Difference of percentages for each class of the distributions (classes forming summary statistics) with weight (importance of the sum stat)
distances_data = data.frame(

monomeric = (observed_data[1,"monomeric"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$monomeric * 100 / rowSums(simulated_data[1,36:75])) * 5,
small_HORs = (observed_data[1,"small_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$small_HORs * 100 / rowSums(simulated_data[1,36:75])) * 1,
medium_HORs = (observed_data[1,"medium_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$medium_HORs * 100 / rowSums(simulated_data[1,36:75])) * 1,
big_HORs = (observed_data[1,"big_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$big_HORs * 100 / rowSums(simulated_data[1,36:75])) * 1,
very_big_HORs = (observed_data[1,"very_big_HORs"] * 100 / rowSums(observed_data[1,2:41])  ) - (simulated_data$very_big_HORs * 100 / rowSums(simulated_data[1,36:75])) * 1,

very_small_diff = (observed_data[1,"very_small_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$small_diff * 100 / rowSums(simulated_data[1,77:177])) * 1,
small_diff = (observed_data[1,"small_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$small_diff * 100 / rowSums(simulated_data[1,77:177])) * 1,
medium_diff = (observed_data[1,"medium_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$medium_diff * 100 / rowSums(simulated_data[1,77:177])) * 1,
big_diff = (observed_data[1,"big_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$big_diff * 100 / rowSums(simulated_data[1,77:177])) * 1,
very_big_diff = (observed_data[1,"very_big_diff"] * 100 / rowSums(observed_data[1,43:143]) ) - (simulated_data$very_big_diff * 100 / rowSums(simulated_data[1,77:177])) * 1,

small_similarity = ((observed_data[1,"small_similarity"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$small_similarity * 100 / rowSums(simulated_data[1,179:279]))) * 5,
medium_similarity_1 = ((observed_data[1,"medium_similarity_1"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$medium_similarity_1 * 100 / rowSums(simulated_data[1,179:279]))) * 5,
medium_similarity_2 = ((observed_data[1,"medium_similarity_2"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$medium_similarity_2 * 100 / rowSums(simulated_data[1,179:279]))) * 5,
medium_similarity_3 = ((observed_data[1,"medium_similarity_3"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$medium_similarity_3 * 100 / rowSums(simulated_data[1,179:279]))) * 5,
big_similarity = ((observed_data[1,"big_similarity"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$big_similarity * 100 / rowSums(simulated_data[1,179:279]))) * 5,
very_big_similarity = ((observed_data[1,"very_big_similarity"] * 100 / rowSums(observed_data[1,145:245]) ) - (simulated_data$very_big_similarity * 100 / rowSums(simulated_data[1,179:279]))) * 5,

most_present_HOR = ((observed_data[1,"most_present_HOR"] ) - (simulated_data$most_present_HOR )) * 0

)

#Calculates the euclidean distance for each simulation 
distances_data = distances_data[,1:17] ^ 2
##Multiplying each difference by the standard deviation of the differences of this class for some kind of normalization in addition to the percentages
distances_data[,1:17] = as.data.frame(lapply(distances_data[,1:17], function(col){col * sd(col)}))
distances_data$square_sum = rowSums(distances_data[,1:17])
distances_data$euclidean_distance = sqrt(distances_data$square_sum)

#Stores the euclidean distances calculated alongside the simulated data
simulated_data[["distance_name"]] = rep("euclidean_distance",nrow(simulated_data))
simulated_data$euclidean_distance = distances_data$euclidean_distance

#Writes the simulated data and the euclidean distances in a file
write.table(simulated_data, file = paste(simulated_data_path, ".euclidean_distances.tab", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)

#Sorts the simulated data by euclidean distances
simulated_data_sorted = simulated_data[order(simulated_data$euclidean_distance),]

#Keep a percentage of the best simulations (here 1%)
acceptance_proportion = nrow(simulated_data) * (1 / 100)
accepted_simulated_data = simulated_data_sorted[1:acceptance_proportion,]

#Writes only the best simulation data in a file
write.table(accepted_simulated_data, file = paste(simulated_data_path, ".euclidean_distances.tab.accepted.tab", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)


#Plot data description distributions for the best simulation

##HOR orders distribution
png(paste(simulated_data_path,".HOR_order_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
barplot(as.numeric(accepted_simulated_data[1,36:75]), main= "", xlab = "Estimated order", ylab = "Monomers count", col = "#95B9C4", names.arg = 1:40, cex.lab = 2)
dev.off()

##HOR orders distribution as a density plot (area under the curve = 1)
png(paste(simulated_data_path,".HOR_order_density.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
barplot(as.numeric(accepted_simulated_data[1,36:75]/accepted_simulated_data[1,4]), main= "", xlab = "Estimated order", ylab = "Density", col = "#95B9C4", names.arg = 1:40, cex.lab = 2)
dev.off()

#Similarity difference distribution
png(paste(simulated_data_path,".max_diff_similarity_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
barplot(as.numeric(accepted_simulated_data[1,77:177]), main= "", xlab = "Difference of similarity percentage", ylab = "Monomers count", col = "#95B9C4", names.arg = 0:100, cex.lab = 2)
dev.off()

#Similarity difference distribution as a density
png(paste(simulated_data_path,".max_diff_similarity_density.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
barplot(as.numeric(accepted_simulated_data[1,77:177]/accepted_simulated_data[1,4]), main= "", xlab = "Difference of similarity percentage", ylab = "Density", col = "#95B9C4", names.arg = 0:100, cex.lab = 2)
dev.off()

#Random similarity distribution (estimation of similarity distribution)
png(paste(simulated_data_path,".random_similarity_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
barplot(as.numeric(accepted_simulated_data[1,179:279]), main= "", xlab = "Similarity percentage", ylab = "Monomers count", col = "#95B9C4", names.arg = 0:100, cex.lab = 2)
dev.off()

#Random similarity distribution as a density
png(paste(simulated_data_path,".random_similarity_density.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
barplot(as.numeric(accepted_simulated_data[1,179:279]/accepted_simulated_data[1,4]), main= "", xlab = "Similarity percentage", ylab = "Density", col = "#95B9C4", names.arg = 0:100, cex.lab = 2)
dev.off()


#Plots of posterior distributions

##Amplification rate distribution
png(paste(simulated_data_path,".amplification_rate_parameter_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$V6, breaks = 100, main= "", xlab = "Amplification rate", ylab = "Simulations count", col = "#95B9C4", cex.lab = 2)
dev.off()

##Amplification rate distribution as density
png(paste(simulated_data_path,".amplification_rate_parameter_density.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$V6, prob = TRUE, breaks = 100, main= "", xlab = "Amplification rate", ylab = "Density", col = "#95B9C4", cex.lab = 2)
dev.off()

#Probability distribution of amplification event producing monomeric organizations
png(paste(simulated_data_path,".monomeric_probability_parameter_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$V8, breaks = 100, main= "", xlab = "Monomeric probability", ylab = "Simulations count", col = "#95B9C4", cex.lab = 2)
dev.off()

#Probability distribution as density
png(paste(simulated_data_path,".monomeric_probability_parameter_density.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$V8, prob = TRUE, breaks = 100, main= "", xlab = "Monomeric probability", ylab = "Density", col = "#95B9C4", cex.lab = 2)
dev.off()

####Alpha and beta parameter values distributions for amplification size########################################################################################
#TODO: Consider removing since harder to interpret than average value of the Beta distribution
png(paste(simulated_data_path,".alpha_amplification_size_parameter_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$V14, breaks = 100, main= "", xlab = "Alpha value", ylab = "Simulations count", col = "#95B9C4", cex.lab = 2)
dev.off()

png(paste(simulated_data_path,".beta_amplification_size_parameter_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$V16, breaks = 100, main= "", xlab = "Beta value", ylab = "Simulations count", col = "#95B9C4", cex.lab = 2)
dev.off()

#################################################################################################################################################################

#Average value distribution of the beta distribution for amplification size
png(paste(simulated_data_path,".Beta_distribution_mean_accepted_amplification_size_parameter_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$amp_mean, breaks = 100, main= "", xlab = "Average values of Amplification size Beta distributions", ylab = "Simulations count", col = "#95B9C4", cex.lab = 2)
dev.off()

#Prior distribution for that posterior (checking for biases)
png(paste(simulated_data_path,".Beta_distribution_mean_amplification_size_parameter_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(simulated_data$amp_mean, breaks = 100, main= "", xlab = "Average values of Amplification size Beta distributions", ylab = "Simulations count", col = "#95B9C4", cex.lab = 2)
dev.off()

#Same as density
png(paste(simulated_data_path,".Beta_distribution_mean_accepted_amplification_size_parameter_density.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$amp_mean/accepted_simulated_data$V4, prob = TRUE, breaks = 100, main= "", xlab = "Average values of Amplification size Beta distributions", ylab = "Density", col = "#95B9C4", cex.lab = 2)
dev.off()

#Same for prior as density
png(paste(simulated_data_path,".Beta_distribution_mean_amplification_size_parameter_density.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(simulated_data$amp_mean/simulated_data$V4, prob = TRUE, breaks = 100, main= "", xlab = "Average values of Amplification size Beta distributions", ylab = "Density", col = "#95B9C4", cex.lab = 2)
dev.off()



####Alpha and beta parameter values distributions for HOR orders ########################################################################################
#TODO: Consider removing since harder to interpret than average value of the Beta distribution
png(paste(simulated_data_path,".alpha_HOR_order_parameter_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$V22, breaks = 100, main= "", xlab = "Alpha value", ylab = "Simulations count", col = "#95B9C4", cex.lab = 2)
dev.off()

png(paste(simulated_data_path,".beta_HOR_order_parameter_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$V24, breaks = 100, main= "", xlab = "Beta value", ylab = "Simulations count", col = "#95B9C4", cex.lab = 2)
dev.off()

#########################################################################################################################################################

#Average value distribution of the beta distribution for HOR orders
png(paste(simulated_data_path,".Beta_distribution_mean_accepted_HOR_order_parameter_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$HOR_mean, breaks = 100, main= "", xlab = "Average values of HOR order Beta distributions", ylab = "Simulations count", col = "#95B9C4", cex.lab = 2)
dev.off()

#Prior distribution for that posterior (checking for biases)
png(paste(simulated_data_path,".Beta_distribution_mean_HOR_order_parameter_distribution.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(simulated_data$HOR_mean, breaks = 100, main= "", xlab = "Average values of HOR order Beta distributions", ylab = "Simulations count", col = "#95B9C4", cex.lab = 2)
dev.off()

#Same but as density
png(paste(simulated_data_path,".Beta_distribution_mean_accepted_HOR_order_parameter_density.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(accepted_simulated_data$HOR_mean, prob = TRUE, breaks = 100, main= "", xlab = "Average values of HOR order Beta distributions", ylab = "Density", col = "#95B9C4", cex.lab = 2)
dev.off()

#Same for prior but as density
png(paste(simulated_data_path,".Beta_distribution_mean_HOR_order_parameter_density.png", sep = ""), width = 900, height = 300)
par(mar = c(5, 5, 4, 5))
hist(simulated_data$HOR_mean, prob = TRUE, breaks = 100, main= "", xlab = "Average values of HOR order Beta distributions", ylab = "Density", col = "#95B9C4", cex.lab = 2)
dev.off()


#Posterior for substitution rate 
#TODO: Fixed instead of amplification rate for latest runs so maybe remove
#png(paste(simulated_data_path,".substitution_rate_parameter_distribution.png", sep = ""), width = 900, height = 300)
#par(mar = c(5, 5, 4, 5))
#hist(accepted_simulated_data$V26, breaks = 250, main= "Posterior distribution for the substitution rate parameter used as \naverage number of substitution per nucleotides in the dataset in the simulations", xlab = "Substitution rate", ylab = "Simulations count", col = "#95B9C4")
#dev.off()

