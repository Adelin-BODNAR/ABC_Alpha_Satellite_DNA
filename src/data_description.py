import argparse, csv, os, subprocess

def get_args():
	"""Function parsing the arguments passed in the command executing this script and checking if the input file path is correct

	Parameters
	----------
	None

	Returns
	-------
	args
		Parsed and checked arguments
	"""

	#Instanciates the argument parser
	parser = argparse.ArgumentParser()
	
	#Adds the argument expected for the input
	parser.add_argument("-i","--input", help = "Argument defining the input monomers file on which the sampling step will be performed", required = True)
	
	#Adds the arguments expected for the sampling step
	parser.add_argument("--sampling_output_path", help = "Argument defining the path to the output file for the sampling step")
	parser.add_argument("-s","--max_sample_size", help = "Argument defining the size of the sampling performed to get the subset of monomers", type = int, required = True)
	parser.add_argument("--random_seed", help = "Argument defining the seed used to perform the sampling used to get the subset of monomers", type = int, default = 42)

	#Adds the arguments expected for the alignment step
	parser.add_argument("--align_output_path", help = "Argument defining the path to the output file for the alignment step")
	parser.add_argument("--thread", help = "Argument defining the number of threads used to perform the alignment of the subset of monomers", type = int, default = 8)
	
	#Adds the arguments expected for the distance matrix generation step
	parser.add_argument("--distmat_output_path", help = "Argument defining the path to the output file for the distance matrix generation step")

	#Adds the argument expected for the distance distribution generation step
	parser.add_argument("--distribution_output_path", help = "Argument defining the path to the output file for the distance distribution generation step")
	parser.add_argument("-p","--distance_distribution_precision", help = "Argument defining the size of bins for the distance distribution", type = int, default = 0, choices = [0,1,2])
	
	
	#Stores the parsed arguments
	args = parser.parse_args()
	
	print("\n",args.input)
	
	#Checks if the input file passed in arguments exists and prints usage then exits with an error if not
	if os.path.isfile(args.input):
		print("Input file found.\n")
	else :
		print("No such input file.\n")
		parser.print_help()
		exit(1)
		
	print("\n",args.sampling_output_path)


	#Checking and definition of the path to the output file of the sampling step
	if args.sampling_output_path is None :
		print("\nNo sampling output path. Usage of default path.")
		args.sampling_output_path = os.path.join(os.path.dirname(args.input),f"{os.path.basename(args.input)}.sampled_{args.max_sample_size}_seed_{args.random_seed}.fst")
		print(args.sampling_output_path,"\n")
		
	elif os.path.isdir(os.path.dirname(os.path.abspath(args.sampling_output_path))) :
		print("\nSampling output directory found.\n")
		
		if os.path.isfile(args.sampling_output_path) :
			print("\nPre-existing file found and will be replaced with sampling output file.\n")
	
	else :
		print("No such sampling output directory.\n")
		parser.print_help()
		exit(1)
	

	#Checking and definition of the path to the output file of the alignment step
	if args.align_output_path is None :
		print("\nNo alignment output path. Usage of default path.")
		args.align_output_path = os.path.join(os.path.dirname(args.sampling_output_path),f"{os.path.basename(args.sampling_output_path)}.align.fst")
		print(args.align_output_path,"\n")

	elif os.path.isdir(os.path.dirname(os.path.abspath(args.align_output_path))) :
		print("\nAlignment output directory found.\n")
		if os.path.isfile(args.align_output_path) :
			print("\nPre-existing file found and will be replaced with alignment output file.\n")

	else :
		print("No such alignment output directory.\n")
		parser.print_help()
		exit(1)


	#Checking and definition of the path to the output file of the distance matric generation step
	if args.distmat_output_path is None :
		print("\nNo distance matrix output path. Usage of default path.")
		args.distmat_output_path = os.path.join(os.path.dirname(args.align_output_path),f"{os.path.basename(args.align_output_path)}.distance_matrix.distmat")
		print(args.distmat_output_path,"\n")

	elif os.path.isdir(os.path.dirname(os.path.abspath(args.distmat_output_path))) :
		print("\nDistance matrix output directory found.\n")
		if os.path.isfile(args.distmat_output_path) :
			print("\nPre-existing file found and will be replaced with distance matrix output file.\n")

	else :
		print("No such distance matrix output directory.\n")
		parser.print_help()
		exit(1)

	#Checking and definition of the path to the output file of the distance distribution generation step
	if args.distribution_output_path is None :
		print("\nNo distance distribution output path. Usage of default path.")
		args.distribution_output_path = os.path.join(os.path.dirname(args.distmat_output_path),f"{os.path.basename(args.distmat_output_path)}.distance_distribution.csv")
		print(args.distribution_output_path,"\n")

	elif os.path.isdir(os.path.dirname(os.path.abspath(args.distribution_output_path))) :
		print("\nDistance distribution output directory found.\n")
		if os.path.isfile(args.distribution_output_path) :
			print("\nPre-existing file found and will be replaced with distance distribution output file.\n")

	else :
		print("No such distance distribution output directory.\n")
		parser.print_help()
		exit(1)
		
	#Prints the arguments for debugging purposes
	print(args)

	return args
	
##############################################################################################################################################################

def monomers_sampling_step(max_sample_size : int, input_file_path : str, sampled_monomers_output_path : str, random_seed): #TODO : Add logs if necessary
	"""Function performing a sampling on the monomers input file

	Parameters
	----------
	max_sample_size : int
		Number of monomers that should be kept from the input file using the sampling script 
		(if greater than the number of monomers available, the sample will contain less monomers than expected)
	input_file_path : str
		String containing the path to the input file used for this step
	sampled_monomers_output_path : str
		String containing the path to the file where the results of this step will be stored
	random_seed
		Argument passed to the main script and parsed at the beginning that defines the random seed used for the sampling

	Returns
	-------
	None
	"""

	#Prints for execution checking purposes
	print(f"\nStarting sampling of size {max_sample_size} with seed {random_seed} for file")
	print(f"{input_file_path}\n")

	#Performing the sampling using the seqtk sample command
	subprocess.run(f"seqtk sample -s {random_seed} {input_file_path} {max_sample_size} > {sampled_monomers_output_path} ", shell = True)

	#Print for execution checking purposes
	print("\nDone.")

	return
	
##############################################################################################################################################################

def monomers_align_step(sampled_input_file_path : str, align_output_path : str, nb_threads : int): #TODO : Add logs if necessary
	"""Function launching the sampled monomers alignment script

	Parameters
	----------
	sampled_input_file_path : str
		String containing the path to the output file of the sampling step used as input file for this step
	align_output_path : str
		String containing the path to the output file where the results of this step will be stored
	nb_threads
		Number of threads that will be used for the alignment
	
	Returns
	-------
	None
	"""


	#Prints for execution checking purposes
	print(f"\nStarting alignment for file")
	print(f"{sampled_input_file_path}\n")
	
	subprocess.run(f"mafft --thread {nb_threads} {sampled_input_file_path} > {align_output_path}", shell = True)

	#TODO: Choose alignment tool (muscle: thread option ?, clustalo: precision ?)
	#subprocess.run(f"muscle -super5 {sampled_input_file_path} -output {align_output_path}", shell = True)
	#subprocess.run(f"time clustalo -i {sampled_input_file_path} -o {align_output_path} -v --threads={nb_threads}", shell = True)
	
	#Print for execution checking purposes
	print("\nDone.")

	return
	
#############################################################################################################################################################################################

def distance_matrix_step(aligned_input_file_path : str, distance_matrix_output_path : str): #TODO : Add logs if necessary
	"""Function launching the distance matrices computation script

	Parameters
	----------
	aligned_input_file_path : str
		String variable containing the path to the output file of the alignment step used as input file for this step
	distance_matrices_output_path : str
		String containing the path to the output file where the results of this step will be stored
	
	
	Returns
	-------
	None
	"""

	#Prints for execution checking purposes
	print(f"\nStarting distance matrix generation for file")
	print(f"{aligned_input_file_path}\n")

	subprocess.run(f"distmat -sequence {aligned_input_file_path} -outfile {distance_matrix_output_path} -nucmethod 0", shell = True)

	#Print for execution checking purposes
	print("\nDone.")

	return

##############################################################################################################################################################

def distance_matrix_to_distribution(distance_matrix_input_path : str, distance_distribution_output_path : str, precision : int): #TODO : Add logs if necessary
	"""Function launching the conversion of a distance matrix to a distribution

	Parameters
	----------
	distance_matrix_input_path : str
		String variable containing the path to the output file of the distance matrix generation step used as input file for this step
	distance_distribution_output_path : str
		String containing the path to the output file where the results of this step will be stored
	
	
	Returns
	-------
	None
	"""
	
	#Dataframe supposed to contain the parsed data
	distances_distrib = dict()

	#Parsing the distance matrix input file
	with open (distance_matrix_input_path,"r") as distance_matrix_file :
		
		#Reading from line 9 which is where the distance matrix start in this format
		for line in distance_matrix_file.readlines()[9:] :

			#Getting only the distance values and converting them to float
			for j,distance in enumerate([float(e) for e in line.split("\t")[:-1] if e != ""]) :
				
				#Avoiding to store the self against self distances
				if j != 0 :

					#Formating the values to match the precision defined with the corresponding argument
					formated_distance = round(distance,precision)

					#Counting values to generate a distribution
					if formated_distance not in distances_distrib.keys() :
						distances_distrib[formated_distance] = 0
					distances_distrib[formated_distance] += 1

	#Print for debugging purposes
	print(distances_distrib)

	#Saving the distribution to a csv file
	with open(distance_distribution_output_path, 'w') as distribution_output_file:  
		writer = csv.writer(distribution_output_file)
		for key in sorted(distances_distrib.keys()):
			writer.writerow([key, distances_distrib[key]])


	return

##############################################################################################################################################################

def main():
	"""Main function that is executed when the script is executed.
	The script is used to extract statistics describing the data passed as input (fasta file).

	Parameters
	----------
	None

	Returns
	-------
	None
	"""

	#Stores the parsed and checked arguments of the script returned by the get_args function
	args = get_args()
	
	monomers_sampling_step(args.max_sample_size, args.input, args.sampling_output_path, args.random_seed)
	
	monomers_align_step(args.sampling_output_path, args.align_output_path, args.thread)
	
	distance_matrix_step(args.align_output_path, args.distmat_output_path)
	
	distance_matrix_to_distribution(args.distmat_output_path, args.distribution_output_path, args.distance_distribution_precision)

	return

##############################################################################################################################################################

if __name__ == "__main__" :
	main()
