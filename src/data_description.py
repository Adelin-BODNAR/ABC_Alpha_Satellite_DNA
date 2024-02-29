import argparse

def get_args(): #TODO : Verify and comment
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
	parser.add_argument("-r","--random_seed", help = "Argument defining the seed used to perform the sampling used to get the subset of monomers", type = int, default = 42)
	
	#Adds the argument expected for the distance distribution generation step
	parser.add_argument("-p","--distance_distribution_precision", help = "Argument defining the number of bins for the distance distribution", type = int, default = 100, choices = [100,1000,10000])
	
	
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
		
	sampling_output_path = ""
	
	print("\n",args.sampling_output_path)
		
	if args.sampling_output_path is None :
		print("\nNo sampling output path. Usage of default path.")
		sampling_output_path = os.path.join(os.path.dirname(args.input),f"{os.path.basename(args.input)}.sampled_{args.max_sample_size}.fst")
		print(sampling_output_path,"\n")
		
	elif os.path.isdir(os.path.dirname(args.sampling_output_path)) :
		print("\nSampling output directory found.\n")
		
		if os.path.isfile(args.sampling_output_path) :
		
			print("\nPre-existing file found and will be replaced with sampling output file.\n")
	
	else :
		
		print("No such output directory.\n")
		parser.print_help()
		exit(1)
		
		
	#Prints the arguments for debugging purposes
	print(args)

	return args
	
##############################################################################################################################################################

def monomers_sampling_step(max_sample_size : int, input_file_path : str, sampled_monomers_output_path : str, random_seed): #TODO : Finish, verify and comment
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
	sampled_monomers_output_path
		String containing the path to the output file of the monomers search and extraction step
	"""

	sampled_monomers_output_path = f"{sampled_monomers_dir_path}/{os.path.splitext(os.path.basename(input_file_path))[0]}_sampled_{sample_size}_seed_{random_seed}.fst"

	with open(monomers_sampling_log_path, "a") as monomers_sampling_log_file :

		subprocess.run(f"time seqtk sample -s{random_seed} {input_file_path} {sample_size} > {sampled_monomers_output_path} ", shell = True, stdout= monomers_sampling_log_file, stderr = subprocess.STDOUT)
		monomers_sampling_log_file.write("\n\n")

	return sampled_monomers_output_path


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
	
	#TODO : Sampling
	
	#TODO : Alignment
	
	#TODO : Distmat
	
	#TODO : Distribution

	return

##############################################################################################################################################################

if __name__ == "__main__" :
	main()
