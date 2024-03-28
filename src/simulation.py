import argparse, os, ete3, random
from Bio import SeqIO

import time

def get_args(): #TODO: Add args for all steps
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
    parser.add_argument("--amplification_rate", help = "Argument defining the frequency at which a monomers are amplified in the simulation", type = float, required = True)
    parser.add_argument("--alpha_amplification_size", help = "Argument defining the alpha of the law allowing to draw randomly the number of new copies each monomer involved gets during an amplification event ", type = float, default = 1.0)
    parser.add_argument("--beta_amplification_size", help = "Argument defining the beta of the law allowing to draw randomly the number of new copies each monomer involved gets during an amplification event ", type = float, default = 1.0)
    parser.add_argument("--min_amplification_size", help = "Argument defining the minimum of the values from which we have to draw randomly the number of new copies each monomer involved gets during an amplification event ", type = int, default = 1)
    parser.add_argument("--max_amplification_size", help = "Argument defining the maximum of the values from which we have to draw randomly the number of new copies each monomer involved gets during an amplification event ", type = int, default = 15000)
    parser.add_argument("--alpha_HOR_order", help = "Argument defining the alpha for the law allowing to draw the number of adjacent monomers involved in the an amplification event ", type = float, default = 1.0)
    parser.add_argument("--beta_HOR_order", help = "Argument defining the beta for the law allowing to draw the number of adjacent monomers involved in the an amplification event ", type = float, default = 1.0)
    parser.add_argument("--min_HOR_order", help = "Argument defining the minimum of the values from which we have to draw the number of adjacent monomers involved in the an amplification event ", type = int, default = 1)
    parser.add_argument("--max_HOR_order", help = "Argument defining the maximum of the values from which we have to draw the number of adjacent monomers involved in the an amplification event ", type = int, default = 50)
    #TODO : (random during simulation => find the right law) parser.add_argument("--amplification_position", help = "Argument defining the position where the new copies each monomer involved are inserted ")
    
    #Adds the arguments expected for the outputs
    parser.add_argument("--monomers_dataset_output_path", help = "Argument defining the path to the monomers dataset output file")
    parser.add_argument("--tree_output_path", help = "Argument defining the path to the phylogenic tree output file")

    #Adds the argument allowing to define the level of verbosity the user wants for the script execution
    parser.add_argument("-v","--verbose", help = "Argument allowing to print the script execution supervision prints (tree and list after each amplification event)", action=argparse.BooleanOptionalAction)

    #Stores the parsed arguments
    args = parser.parse_args()
    
    print("-----------------------------------------------------------------------------------\n\n") #TODO: Remove once benchmarking done

    #Prints the provided path to the input file for execution supervision purposes
    print("\n",args.input)
    
    #Checks if the input file passed in arguments exists and prints usage then exits with an error if not
    if os.path.isfile(args.input):
        print("Input file found.\n")
    else :
        print("No such input file.\n")
        parser.print_help()
        exit(1)
        
    #Prints the monomers dataset output path for execution supervision purposes
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
        
    #Prints the argument values provided for execution supervision purposes
    print(args)

    return args

##############################################################################################################################################################

def amplification_simulation(dataset_size : int, amplification_rate : float, alpha_amplification_size : float, beta_amplification_size : float, min_amplification_size : int, max_amplification_size : int, alpha_HOR_order : float, beta_HOR_order : float, min_HOR_order : int, max_HOR_order : int, verbose : bool) :
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

    #Adds the ancestor monomer (no sequence yet) to a list used to store the spacial organisation of the monomers
    monomers_dataset = [tree]

    #Initializes the values for the name (id of the monomer) and dist (lenght of the branch) of the ancestor monomer
    tree.name = "0"
    tree.dist = 0

    #Adds the features used for tracking the previous and next element in the double chained list to the Node of the ancestor monomer
    tree.add_features(amplification_id = 0, previous = None, next = None, parent_name = -1)

    #Initializes the variables used to give IDs to the new monomers and each amplification event
    monomers_counter = 1
    amplification_counter = 0

    #TODO Remove once bechmarking over
    time_create, time_replace, time_insert = 0 , 0 , 0 

    #Loops until the monomers dataset reached the aimed size
    while (len(monomers_dataset) < dataset_size):

        #TODO Remove once bechmarking over
        start = time.time() 

        amplification_counter += 1

        #Calculates the time before the next amplification depending of the probability to get an amplification event and adds it to all branches of the leaves (monomers_dataset)
        time_event = random.expovariate( 1/ (amplification_rate * len(monomers_dataset)) )
        for monomer in monomers_dataset :
            monomer.dist += time_event

        #Draws random size of amplification and random size of HOR
        HOR_order = min(min_HOR_order + int( random.betavariate(alpha_HOR_order, beta_HOR_order) * (max_HOR_order - min_HOR_order + 1)),len(monomers_dataset))
        amplification_size = min(min_amplification_size + int( random.betavariate(alpha_amplification_size, beta_amplification_size) * (max_amplification_size - min_amplification_size + 1)), int((dataset_size - len(monomers_dataset))/HOR_order)+1 )


        #Draws random index of monomer from the current monomers dataset
        index_head_monomers_to_amplify = random.randint(0,len(monomers_dataset)-HOR_order)

        #Initializes the list of monomers that will that will be involved in the amplification event
        monomers_to_amplify = monomers_dataset[index_head_monomers_to_amplify:index_head_monomers_to_amplify + HOR_order]

        #Draws randomly the index in the monomers dataset where the new monomers must be inserted
        position_new_monomer = random.randint(0,len(monomers_dataset)) #TODO : Maybe choose with a beta distribution for instance instead of uniform

        #Loops to create all the monomers needed to be amplified depending on the random monomer chosen and the size of HOR
        for i in range(amplification_size + 1) :

            #Loops through the monomers that need to be amplified to amplify them in the order of the HOR
            for j, currently_amplified_monomer in enumerate(monomers_to_amplify) :

                #Creates the new copy of the currently amplified monomer and adds it to the tree (if the monomer is the replacement of the origin monomer it keeps the same name)
                new_monomer = currently_amplified_monomer.add_child(ete3.TreeNode(name = str(monomers_counter) if i > 0 else currently_amplified_monomer.name , dist = 0))
                #Adds the name of the parent as a feature of the node for later use during output files generation
                new_monomer.add_features(parent_name = new_monomer.up.name)

                #TODO Remove once bechmarking over
                end = time.time()
                time_create += end-start
                start = time.time() 
                
                #Replaces the origin monomers of the amplification with their replacements in the monomers dataset
                if i == 0 :

                    new_monomer.add_features(amplification_id = new_monomer.up.amplification_id)
                    new_monomer.parent_name = new_monomer.up.parent_name 
                    monomers_dataset[index_head_monomers_to_amplify + j] = new_monomer

                    #TODO Remove once bechmarking over
                    end = time.time()
                    time_replace += end-start
                    start = time.time() 

                #Inserts the new monomers in the monomers dataset
                else :
                    new_monomer.add_features(amplification_id = amplification_counter)
                    monomers_dataset.insert(position_new_monomer,new_monomer)
                    position_new_monomer += 1
                    monomers_counter += 1

                    #TODO Remove once bechmarking over
                    end = time.time()
                    time_insert += end-start
                    start = time.time() 



        
        if (verbose) :
            #Prints the list to be able to see the spatial organistation of the monomers and the tree after each amplification event
            c=0
            for current_monomer in monomers_dataset :
            
                print(current_monomer.name,"\t",current_monomer.amplification_id,"\t",current_monomer.parent_name)
            
                c += 1
        
            print(tree.get_ascii(show_internal= True, attributes= ["name","amplification_id","dist"]),"\n")
        
            print("\n\n")
            
        
    #Calculates the time before the next supposed amplification depending of the probability to get an amplification event and adds it to all branches of the leaves in order to have a length greater than 0 for the branches created last
    time_event = random.expovariate( 1/ amplification_rate * len(monomers_dataset) )
    for monomer in monomers_dataset :
        monomer.dist += time_event

    #TODO : Remove once benchmarking done
    print(f"\n\n------ Create time : {time_create} ------")
    print(f"\n\n------ Replacement time : {time_replace} ------")
    print(f"\n\n------ Insertion time : {time_insert} ------")

    return tree, monomers_dataset

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

    #TODO Remove once bechmarking over
    start = time.time() 

    #Launches the amplification simulation and gets the tree , the head of the chained list for the monomers dataset and the number of amplifications done used during the generation of the visual representation of the spatial organization
    tree, monomers_dataset = amplification_simulation(args.max_size, args.amplification_rate, args.alpha_amplification_size, args.beta_amplification_size, args.min_amplification_size, args.max_amplification_size, args.alpha_HOR_order, args.beta_HOR_order, args.min_HOR_order, args.max_HOR_order, args.verbose)

    #TODO Remove once bechmarking over
    end = time.time()
    print(f"\n\n------ Simulation time : {end-start} ------")
    start = time.time() 

    #Saves the tree in a NHX format file
    tree.write(outfile= args.tree_output_path,features=["amplification_id"], format= 1) #TODO : Output format (NHX ?) 

    #TODO Remove once bechmarking over
    end = time.time()
    print(f"\n\n------ Write time : {end-start} ------\n\n")
    

    #TODO : Mutation (+ indels)

    #TODO : Write monomers with sequences and valuable infos as commentary in a fasta file and/or tabulated file
    
    if args.verbose :
        #Prints for execution supervision purposes
        c=0
        for current_monomer in monomers_dataset :
        
            print(current_monomer.name,"\t",current_monomer.amplification_id,"\t",current_monomer.parent_name)
        
            c += 1
    
        print(tree.get_ascii(show_internal= True, attributes= ["name","amplification_id","dist"]),"\n")
    
        print("\n\n")

    return

##############################################################################################################################################################

if __name__ == "__main__" :
    main()