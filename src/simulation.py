import argparse, os, ete3, random
from Bio import SeqIO

def get_args(): #TODO: Add args for all steps (Considered args : option between newick and nexus format for tree output)
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
    parser.add_argument("-i","--input", help = "Argument defining the input fasta file containing the ancestor monomer used to perform the simulation", required = True)

    #Adds the argument expected for defining the size of the dataset for which the simulation must stop
    parser.add_argument("-s","--max_size", help = "Argument defining the max size of the monomers dataset before the simulation stops", type = int, required = True)

    #Adds the argument expected for defining the simulation run identifier
    parser.add_argument("-id","--simulation_id", help = "Argument defining the ID used to identify the current run of the simulation", default = "no_id")

    #Adds the arguments expected for the amplification simulation
    parser.add_argument("--amplification_frequency", help = "Argument defining the frequency at which a monomers are amplified in the simulation", type = float, required = True)
    #TODO : (random during simulation => find the right law) parser.add_argument("--amplification_size", help = "Argument defining the number of new copies each monomer involved gets during an amplification event ")
	#TODO : (random during simulation => find the right law) parser.add_argument("--hor_size", help = "Argument defining the number adjacent monomers involved in the an amplification event ")
	#TODO : (at the end or random during simulation => find the right law) parser.add_argument("--amplification_position", help = "Argument defining the number of new copies each monomer involved gets during an amplification event ")
	


    #Adds the arguments expected for the outputs
    parser.add_argument("--monomers_dataset_output_path", help = "Argument defining the path to the monomers dataset output file")
    parser.add_argument("--tree_output_path", help = "Argument defining the path to the phylogenic tree output file")

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
		
    print("\n",args.monomers_dataset_output_path)


    #Checking and definition of the path to the monomers dataset output file
    if args.monomers_dataset_output_path is None :
        print("\nNo monomers dataset output path. Usage of default path.")
        args.monomers_dataset_output_path = os.path.join(os.path.dirname(args.input),f"{os.path.basename(args.input)}.monomers_dataset_simulation_{args.simulation_id}.fst")
        print(args.monomers_dataset_output_path,"\n")
		
    elif os.path.isdir(os.path.dirname(os.path.abspath(args.monomers_dataset_output_path))) :
        print("\nMonomers dataset output directory found.\n")
		
        if os.path.isfile(args.monomers_dataset_output_path) :
            print("\nPre-existing file found and will be replaced with monomers dataset output file.\n")
	
    else :
        print("No such monomers dataset output directory.\n")
        parser.print_help()
        exit(1)
	

    #Checking and definition of the path to the phylogenic tree output file
    if args.tree_output_path is None :
        print("\nNo phylogenic tree output path. Usage of default path.")
        args.tree_output_path = os.path.join(os.path.dirname(args.input),f"{os.path.basename(args.input)}.phylogenic_tree_simulation_{args.simulation_id}.tree")
        print(args.tree_output_path,"\n")

    elif os.path.isdir(os.path.dirname(os.path.abspath(args.tree_output_path))) :
        print("\nPhylogenic tree output directory found.\n")
        if os.path.isfile(args.tree_output_path) :
            print("\nPre-existing file found and will be replaced with phylogenic tree output file.\n")

    else :
        print("No such phylogenic tree output directory.\n")
        parser.print_help()
        exit(1)
		
    #Prints the arguments for debugging purposes
    print(args)

    return args

##############################################################################################################################################################

def amplification_simulation(amplification_frequency) :
    """Function simulating the amplification of one ancestor monomer during evolution

    Parameters
	----------
    None

    Returns
    -------
    tree
        Tree structure from the ete3 library used to store and see the amplifiations through time
    """

    #Instantiates the tree object used for the amplification simulation
    tree = ete3.Tree()

    #Adds the ancestor monomer (no sequence yet) to the tree and to a variable used to store the start of the double chained list used for to store the spacial organisation of the monomers
    head = tree.add_child(ete3.TreeNode(name = 0))

    #Adds the features used for tracking the previous and next element in the double chained list to the Node of the ancestor monomer
    head.add_features(previous = None, next = None)

    #Calculates the length of the branch depending of the probability to get an amplification event
    c = 0
    while random.random() > amplification_frequency :
        c += 1
    head.dist += c

    #TODO : While loop until tree.get_leaves() reaches the right size

        #TODO : Draw random size of amplification and random size of HOR

        #TODO : Draw random monomer from tree.get_leaves()

        #TODO : Get all the monomers needed to be amplified depending on the random monomer chosen and the size of HOR

        #TODO : Create all new monomers and store them in a list to add them to the chained list (First monomers replace the amplified ones in the chained list and the others are inserted in the list elsewhere)



    return tree

##############################################################################################################################################################

def main():
	"""Main function that is executed when the script is executed.
	The script is used to simulate the amplification and mutation of an ancestor monomer ( insertion and deletion too for phase 2).

	Parameters
	----------
	None

	Returns
	-------
	None
	"""

	#Stores the parsed and checked arguments of the script returned by the get_args function
	args = get_args()



	return

##############################################################################################################################################################

if __name__ == "__main__" :
	main()