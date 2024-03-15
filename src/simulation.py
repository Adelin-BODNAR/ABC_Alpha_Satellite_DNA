import argparse, os, ete3, random
from Bio import SeqIO

from mycolorpy import colorlist as mcp

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
    parser.add_argument("--max_amplification_size", help = "Argument defining the maximum of the values from which we have to draw randomly the number of new copies each monomer involved gets during an amplification event ", type = int, default = 9)
    parser.add_argument("--alpha_HOR_order", help = "Argument defining the alpha for the law allowing to draw the number of adjacent monomers involved in the an amplification event ", type = float, default = 1.0)
    parser.add_argument("--beta_HOR_order", help = "Argument defining the beta for the law allowing to draw the number of adjacent monomers involved in the an amplification event ", type = float, default = 1.0)
    parser.add_argument("--min_HOR_order", help = "Argument defining the minimum of the values from which we have to draw the number of adjacent monomers involved in the an amplification event ", type = int, default = 1)
    parser.add_argument("--max_HOR_order", help = "Argument defining the maximum of the values from which we have to draw the number of adjacent monomers involved in the an amplification event ", type = int, default = 3)
    #TODO : (random during simulation => find the right law) parser.add_argument("--amplification_position", help = "Argument defining the number of new copies each monomer involved gets during an amplification event ")
    
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

def amplification_simulation(dataset_size : int, amplification_rate : float, alpha_amplification_size : float, beta_amplification_size : float, min_amplification_size : int, max_amplification_size : int, alpha_HOR_order : float, beta_HOR_order : float, min_HOR_order : int, max_HOR_order : int) :
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
    head = tree

    #Initializes the values for the name (id of the monomer) and dist (lenght of the branch) of the ancestor monomer
    tree.name = "0"
    tree.dist = 0

    #Adds the features used for tracking the previous and next element in the double chained list to the Node of the ancestor monomer
    tree.add_features(amplification_id = 0, previous = None, next = None, parent_name = 0)

    #Initializes the variables used to give IDs to the new monomers and each amplification event
    node_counter = 1
    amplification_counter = 0

    #Loops until the monomers dataset reached the aimed size
    while (len(tree.get_leaves()) < dataset_size):

        amplification_counter += 1

        #Calculates the time before the next amplification depending of the probability to get an amplification event and adds it to all branches of the leaves
        time_event = random.expovariate( 1/ amplification_rate * len(tree.get_leaves()) )
        for monomer in tree.get_leaves() :
            monomer.dist += time_event

        #Draws random size of amplification and random size of HOR
        amplification_size = range(min_amplification_size, max_amplification_size)[int( random.betavariate(alpha_amplification_size, beta_amplification_size) * (max_amplification_size - min_amplification_size))]
        HOR_order = range(min_HOR_order, max_HOR_order)[int( random.betavariate(alpha_HOR_order, beta_HOR_order) * (max_HOR_order - min_HOR_order))]

        #Draws random monomer from the current dataset (tree.get_leaves())
        head_monomers_to_amplify = random.choice(tree.get_leaves())

        #Initializes the list that will contain the new versions of the origin monomers of the amplification
        old_monomers_replacements = []

        #Initializes the list that will contain the new monomers created during the amplification event
        new_monomers = []

        #Loops to create all the monomers needed to be amplified depending on the random monomer chosen and the size of HOR
        for i in range(amplification_size + 1) :

            #Sets the currently amplified monomer to be the randomly drawn monomer as the first in the list of monomers to amplify at the beginning of each repeat of the HOR
            currently_amplified_monomer = head_monomers_to_amplify

            #Loops through the monomers that need to be amplified to amplify them in the order of the HOR
            for j in range(HOR_order) :

                #Creates the new copy of the currently amplified monomer and adds it to the tree (if the monomer is the replacement of the origin monomer it keeps the same name)
                new_monomer = currently_amplified_monomer.add_child(ete3.TreeNode(name = str(node_counter) if i > 0 else currently_amplified_monomer.name , dist = 0))
                #Adds the name of the parent as a feature of the node for later use during output files generation
                new_monomer.add_features(parent_name = new_monomer.up.name)
                
                #Adds the new monomer to the right list depending of if it's a replacement of a brand new one
                if i == 0 :
                    old_monomers_replacements.append(new_monomer)
                else :
                    new_monomers.append(new_monomer)
                    node_counter += 1

                #Checks and corrects the size of the HOR depending of if there are not enough monomers after the head of the amplified monomers 
                if currently_amplified_monomer.next == None :
                    HOR_order = j+1
                    break
                #Sets the currently amplified mnomer to the next one in the chained list representing the spatial organisation of the monomers dataset
                else :
                    currently_amplified_monomer = currently_amplified_monomer.next

            #Draws randomly the monomer after which the new ones must be inserted right after the creation of the origin monomers replacements
            if i == 0 :
                position_new_monomers = random.choice(tree.get_leaves()) #TODO : Maybe choose with a beta distribution for instance instead of uniform
                
        #Replaces the origin monomers of theamplification with their replacements in the chained list
        for i in range(len(old_monomers_replacements)) :
            old_monomers_replacements[i].add_features(amplification_id = old_monomers_replacements[i].up.amplification_id)
            old_monomers_replacements[i].parent_name = old_monomers_replacements[i].up.parent_name 

            if i == 0 :
                old_monomers_replacements[i].add_features(previous = old_monomers_replacements[i].up.previous)
                if old_monomers_replacements[i].previous is None :
                    head = old_monomers_replacements[i]
                else :
                    old_monomers_replacements[i].previous.next = old_monomers_replacements[i]
            else :
                old_monomers_replacements[i].add_features(previous = old_monomers_replacements[i-1])

            if i == len(old_monomers_replacements) -1 :
                old_monomers_replacements[i].add_features(next = old_monomers_replacements[i].up.next)
                if old_monomers_replacements[i].next is not None :
                    old_monomers_replacements[i].next.previous = old_monomers_replacements[i]
            else :
                old_monomers_replacements[i].add_features(next = old_monomers_replacements[i+1])

        #Inserts the new monomers in the chained list
        for i in range(len(new_monomers)) :
            new_monomers[i].add_features(amplification_id = amplification_counter)

            if i == 0 :
                new_monomers[i].add_features(previous = position_new_monomers)
            else :
                new_monomers[i].add_features(previous = new_monomers[i-1])

            
            if i == len(new_monomers) -1 :
                new_monomers[i].add_features(next = position_new_monomers.next)
                position_new_monomers.next = new_monomers[0]
                if new_monomers[i].next is not None :
                    new_monomers[i].next.previous = new_monomers[i]
            else :
                new_monomers[i].add_features(next = new_monomers[i+1])

        #Prints the chained list to be able to see the spatial organistation of the monomers and the tree after each amplification event TODO : Remove once debugging done#####
        current_monomer = head
        c=0
        while current_monomer is not None :
            print(current_monomer.name,"\t",current_monomer.amplification_id,"\t",current_monomer.parent_name,"\t",current_monomer.up.name,"\t",current_monomer.previous.name if current_monomer.previous is not None else current_monomer.previous, current_monomer.next.name if current_monomer.next is not None else current_monomer.next)
            current_monomer = current_monomer.next
            c += 1

        print(tree.get_ascii(show_internal= True, attributes= ["name","amplification_id","dist"]),"\n")

        print("\n\n")
        #######################################

    #Calculates the time before the next supposed amplification depending of the probability to get an amplification event and adds it to all branches of the leaves in order to have a length greater than 0 for the branches created last
    time_event = random.expovariate( 1/ amplification_rate * len(tree.get_leaves()) )
    for monomer in tree.get_leaves() :
        monomer.dist += time_event

    return tree, head, amplification_counter

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

    #Launches the amplification simulation and gets the tree , the head of the chained list for the monomers dataset and the number of amplifications done used during the generation of the visual representation of the spatial organization
    tree, head, nb_amplifications = amplification_simulation(args.max_size, args.amplification_rate, args.alpha_amplification_size, args.beta_amplification_size, args.min_amplification_size, args.max_amplification_size, args.alpha_HOR_order, args.beta_HOR_order, args.min_HOR_order, args.max_HOR_order)

    #Saves the tree to an image    
    tree.render(args.tree_output_path+".png")

    #Saves the tree in a NHX format file
    tree.write(outfile= args.tree_output_path,features=["amplification_id"], format= 1) #TODO : Output format (NHX ?) 

    #TODO : Mutation (+ indels)

    #Prints for debugging purposes###############
    print(tree.get_ascii(show_internal= True, attributes= ["name","amplification_id","dist"]),"\n")

    #TODO : Write monomers with sequences and valuable infos as commentary in a fasta file and/or tabulated file

    #TODO : Tries to generate a good representation of the spatial organization
    current_monomer = head
    c=0
    centromere = ete3.Tree()
    centromere.dist = 0
    nb_leaves = len(tree.get_leaves())
    color1=mcp.gen_color(cmap="gist_ncar",n=nb_leaves*(nb_amplifications+1))
    random.shuffle(color1)
    colors = [color1[i*nb_leaves:(i+1)*nb_leaves] for i in range(nb_amplifications + 1)]
    random.shuffle(colors)
    style = ete3.NodeStyle()
    style["fgcolor"] = "#ffffff"
    centromere.set_style(style)
    while current_monomer is not None :
        print(current_monomer.name,"\t",current_monomer.amplification_id,"\t",current_monomer.parent_name,"\t",current_monomer.previous.name if current_monomer.previous is not None else current_monomer.previous, current_monomer.next.name if current_monomer.next is not None else current_monomer.next)
        m = centromere.add_child(ete3.TreeNode(dist = 0, name = " "+current_monomer.name))
        m.add_features(amplification_id = current_monomer.amplification_id, parent_name = current_monomer.parent_name)
        style = ete3.NodeStyle()
        style["fgcolor"] = colors[int(current_monomer.amplification_id)][int(current_monomer.parent_name)]
        style["shape"] = "square"
        style["size"] = 30
        m.set_style(style)
        current_monomer = current_monomer.next
        c += 1

    print("\n",c,len(tree.get_leaves()))

    ts = ete3.TreeStyle()
    ts.show_leaf_name = True
    ts.rotation = -90
    centromere.show(tree_style= ts)
    #############################################

    return

##############################################################################################################################################################

if __name__ == "__main__" :
    main()