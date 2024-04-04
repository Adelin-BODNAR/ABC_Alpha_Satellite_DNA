/**
* @file simulation.cc
* @brief Monomers amplification and mutation simulation script.
* @author Adelin BODNAR
* @date 02/04/2024
*
* Script used to simulate the successive amplifications and mutations (substitutions and/or indels) of an ancestor monomer and its children.
*
*/

#include <iostream>
#include <stdlib.h>
#include <string>
#include <filesystem>
#include <getopt.h>


/** @brief Parsed arguments structure
* 
* Structure used to store the arguments provided by the user in the command line after parsing.
*/
struct Args {

    int verbose ; /**< Integer used to determine the level of verbosity for the execution of the script */

    std::string input ; /**< String containing the path provided by the user to the input fasta file with the ancestor monomer */
    int max_size ; /**< Number of monomers aimed for in the dataset at the end of the simulation */
    std::string simulation_id ; /**< String used as an identifier of an execution of the simulation script */
    double amplification_rate ; /**< Rate at which each monomer in the dataset is amplified */

    double alpha_amplification_size ; /**< Alpha parameter for the beta distribution used to draw a random value for the amplification size during each amplification event */
    double beta_amplification_size ; /**< Beta parameter for the beta distribution used to draw a random value for the amplification size during each amplification event */
    int min_amplification_size ; /**< Minimum value of the range from which a random value for the amplification size is drawn during each amplification event */
    int max_amplification_size ; /**< Maximum value of the range from which a random value for the amplification size is drawn during each amplification event */

    double alpha_HOR_order ; /**< Alpha parameter for the beta distribution used to draw a random value for the HOR order during each amplification event */
    double beta_HOR_order ; /**< Beta parameter for the beta distribution used to draw a random value for the HOR order during each amplification event */
    int min_HOR_order ; /**< Minimum value of the range from which a random value for the HOR order is drawn during each amplification event */
    int max_HOR_order ; /**< Maximum value of the range from which a random value for the HOR order is drawn during each amplification event */

    std::string monomers_dataset_output_path ; /**< String containing the path to the output fasta file generated by the script with the monomers from the dataset retaining their spatial organization */
    std::string tree_output_path ; /**< String containing the path to the output NHX file generated by the script with the phylogenic tree retracing the amplification history of the monomers */

};




/** @brief Parsing the arguments.
*
* Function parsing the arguments provided by the user in the command line and storing them in the structure defined above.
*
* @param argc The number of arguments provided by the user in the command line.
* @param argv The list of arguments provided by the user in the command line.
*
* @return parsed_args an instance of the Args structure containing the arguments provided by the user and used as parameters for multiple parts of the simulation script
*/
Args get_args(int argc, char **argv) {

    //Initializing parsed_args and default values for some of its attributes

    Args parsed_args ;  /**< Variable used to store and return the arguments provided by the user after parsing */

    parsed_args.verbose = 0;   /**< Initializes the default verbosity level at 0 */

    parsed_args.max_size = 15000 ;          /**< Initializes the default max size for the monomers dataset at 15 000 monomers (middle ground between ~30 000 monomers for the largest human chromosome and ~5 000 for the smallest) */
    parsed_args.simulation_id = "no_id" ;   /**< Initializes the default simulation id */
    parsed_args.amplification_rate = 1.0 ;  /**< Initializes the default amplification rate at 1 amplification per monomer per time unit */

    parsed_args.alpha_amplification_size = 1.0 ;    /**< Initializes the default alpha parameter for the amplification size beta distribution at 1 to mimic an uniform distribution with the default beta parameter */
    parsed_args.beta_amplification_size = 1.0 ;     /**< Initializes the default beta parameter for the amplification size beta distribution at 1 to mimic an uniform distribution with the default alpha parameter */
    parsed_args.min_amplification_size = 1 ;        /**< Initializes the default minimum amplification size to 1 to perform at the very least a duplication for each amplification event */
    parsed_args.max_amplification_size = 10000 ;    /**< Initializes the default maximum amplification size */ //TODO : Find the right value and a justification for it.

    parsed_args.alpha_HOR_order = 1.0 ;     /**< Initializes the default alpha parameter for the HOR order beta distribution at 1 to mimic an uniform distribution with the default beta parameter */
    parsed_args.beta_HOR_order = 1.0 ;      /**< Initializes the default beta parameter for the HOR order beta distribution at 1 to mimic an uniform distribution with the default alpha parameter */
    parsed_args.min_HOR_order = 1 ;         /**< Initializes the default minimum HOR order to 1 to also mimic the amplifications resulting to a monomeric organization (layer of monomers amplified for one origin monomer) */
    parsed_args.max_HOR_order = 50 ;        /**< Initializes the default maximum HOR order to 50 monomers since it is not much above 34 monomers the maximum observed yet in humans */

    std::filesystem::path dir_dataset ;
    std::filesystem::path dir_tree ;

    /**< String containing the help message printed when the -h , --help or a wrong option are provided ny the user */
    const std::string usage = "Usage: \n"
                            "\t./simulation [-h] -i INPUT [-s MAX_SIZE] [--id SIMULATION_ID] [--amplification_rate AMPLIFICATION_RATE] \n"
                            "\t[--alpha_amplification_size ALPHA_AMPLIFICATION_SIZE] [--beta_amplification_size BETA_AMPLIFICATION_SIZE]\n"
                            "\t[--min_amplification_size MIN_AMPLIFICATION_SIZE] [--max_amplification_size MAX_AMPLIFICATION_SIZE]\n"
                            "\t[--alpha_HOR_order ALPHA_HOR_ORDER] [--beta_HOR_order BETA_HOR_ORDER] [--min_HOR_order MIN_HOR_ORDER]\n"
                            "\t[--max_HOR_order MAX_HOR_ORDER] [--monomers_dataset_output_path MONOMERS_DATASET_OUTPUT_PATH]\n"
                            "\t[--tree_output_path TREE_OUTPUT_PATH] [-v | --verbose | --no-verbose]\n\n"

                            "\toptions:\n\n"

                            "\t-h, --help \n\t\tOption showing this help message and exiting\n\n"

                            "\t-i INPUT, --input INPUT\n"
                            "\t\tArgument defining the input fasta file containing the \n"
                            "\t\tancestor monomer used to perform the simulation \n\n"

                            "\t-s MAX_SIZE, --max_size MAX_SIZE \n"
                            "\t\tArgument defining the max size of the monomers dataset \n"
                            "\t\tbefore the simulation stops \n\n"

                            "\t--id SIMULATION_ID, --simulation_id SIMULATION_ID \n"
                            "\t\tArgument defining the ID used to identify the current \n"
                            "\t\trun of the simulation \n\n"

                            "\t--amplification_rate AMPLIFICATION_RATE \n"
                            "\t\tArgument defining the rate at which monomers \n"
                            "\t\tare amplified in the simulation \n\n"

                            "\t--alpha_amplification_size ALPHA_AMPLIFICATION_SIZE \n"
                            "\t\tArgument defining the alpha of the law allowing to \n"
                            "\t\tdraw randomly the number of new copies each monomer \n"
                            "\t\tinvolved gets during an amplification event \n\n"

                            "\t--beta_amplification_size BETA_AMPLIFICATION_SIZE \n"
                            "\t\tArgument defining the beta of the law allowing to draw \n"
                            "\t\trandomly the number of new copies each monomer \n"
                            "\t\tinvolved gets during an amplification event \n\n"

                            "\t--min_amplification_size MIN_AMPLIFICATION_SIZE \n"
                            "\t\tArgument defining the minimum of the values from which \n"
                            "\t\twe have to draw randomly the number of new copies each \n"
                            "\t\tmonomer involved gets during an amplification event \n\n"

                            "\t--max_amplification_size MAX_AMPLIFICATION_SIZE \n"
                            "\t\tArgument defining the maximum of the values from which \n"
                            "\t\twe have to draw randomly the number of new copies each \n"
                            "\t\tmonomer involved gets during an amplification event \n\n"

                            "\t--alpha_HOR_order ALPHA_HOR_ORDER \n"
                            "\t\tArgument defining the alpha for the law allowing to \n"
                            "\t\tdraw the number of adjacent monomers involved in the \n"
                            "\t\tan amplification event \n\n"

                            "\t--beta_HOR_order BETA_HOR_ORDER \n"
                            "\t\tArgument defining the beta for the law allowing to \n"
                            "\t\tdraw the number of adjacent monomers involved in the \n"
                            "\t\tan amplification event \n\n"

                            "\t--min_HOR_order MIN_HOR_ORDER \n"
                            "\t\tArgument defining the minimum of the values from which \n"
                            "\t\twe have to draw the number of adjacent monomers \n"
                            "\t\tinvolved in the an amplification event \n\n"

                            "\t--max_HOR_order MAX_HOR_ORDER \n"
                            "\t\tArgument defining the maximum of the values from which \n"
                            "\t\twe have to draw the number of adjacent monomers \n"
                            "\t\tinvolved in the an amplification event \n\n"

                            "\t--monomers_dataset_output_path MONOMERS_DATASET_OUTPUT_PATH \n"
                            "\t\tArgument defining the path to the monomers dataset \n"
                            "\t\toutput file \n\n"

                            "\t--tree_output_path TREE_OUTPUT_PATH \n"
                            "\t\tArgument defining the path to the phylogenic tree \n"
                            "\t\toutput file \n\n"

                            "\t-v, --verbose, --no-verbose \n"
                            "\t\tOptions allowing to print or not the script execution \n"
                            "\t\tsupervision prints (tree and list after each amplification event) \n\n" ;


    /** @brief Long options defining structure
    *
    * Structure used by the function getopt_long to define, parse and link to their short versions the long options for the simulation script.
    */
    const struct option long_options[] =
        {
          // These options don’t set a flag. We distinguish them by their indices. 

          {"id", required_argument, 0, 0},
          {"simulation_id", required_argument, 0, 0},

          {"amplification_rate", required_argument, 0, 0},

          {"alpha_amplification_size", required_argument, 0, 0},
          {"beta_amplification_size", required_argument, 0, 0},
          {"min_amplification_size", required_argument, 0, 0},
          {"max_amplification_size", required_argument, 0, 0},

          {"alpha_HOR_order", required_argument, 0, 0},
          {"beta_HOR_order", required_argument, 0, 0},
          {"min_HOR_order", required_argument, 0, 0},
          {"max_HOR_order", required_argument, 0, 0},

          {"monomers_dataset_output_path", required_argument, 0, 0},
          {"tree_output_path", required_argument, 0, 0},

          //  These options set flags. 

          {"help", no_argument, 0, 'h'},

          {"input", required_argument, 0, 'i'},
          {"max_size", required_argument, 0, 's'},

          {"verbose", no_argument, 0, 'v'},

          {0, 0, 0, 0}
        };

    // getopt_long stores the option index here.
    int option_index;

    int c;

    while ((c = getopt_long (argc, argv, "hi:s:v", long_options, &option_index)) != -1) {

        switch (c)
        {
            case 'i':

                //std::cout << "Input: " << optarg << std::endl;
                if (std::filesystem::exists(optarg)) {
                    std::cout << "Input file found.\n" << std::endl;
                    parsed_args.input = optarg;
                }else{
                    std::cout << "No such input file.\n" << std::endl;
                    std::cout << usage.c_str();
                    exit(1);
                }

                break;

            case 's':
                
                //std::cout << "Max size: " << optarg << std::endl << std::endl;
                if (std::stoi(optarg) > 0) {
                    parsed_args.max_size = std::stoi(optarg) ;
                }else{
                    std::cout << "Max size value must be > 0." << std::endl << std::endl;
                    std::cout << usage.c_str();
                    exit(1);
                }

                break;

            case 'v':
                
                parsed_args.verbose ++ ;
               
                break;

            case 0 :
                
                switch (option_index){

                    case 0 :
                    case 1 :

                        //std::cout << "Simulation Id: " << optarg << std::endl << std::endl;
                        parsed_args.simulation_id = optarg;

                        break ;
                    
                    case 2 :

                        //std::cout << "Amplification rate: " << optarg << std::endl << std::endl;
                        parsed_args.amplification_rate = std::stod(optarg);

                        if (parsed_args.amplification_rate <= 0) {
                            std::cout << "Amplification rate value must be > 0." << std::endl << std::endl;
                            std::cout << usage.c_str();
                            exit(1);
                        }

                        break;

                    case 3 :

                        //std::cout << "Alpha amplification size: " << optarg << std::endl << std::endl;
                        parsed_args.alpha_amplification_size = std::stod(optarg);

                        if (parsed_args.alpha_amplification_size <= 0) {
                            std::cout << "Alpha amplification size value must be > 0." << std::endl << std::endl;
                            std::cout << usage.c_str();
                            exit(1);
                        }

                        break;

                    case 4 :

                        //std::cout << "Beta amplification size: " << optarg << std::endl << std::endl;
                        parsed_args.beta_amplification_size = std::stod(optarg);

                        if (parsed_args.beta_amplification_size <= 0) {
                            std::cout << "Beta amplification size value must be > 0." << std::endl << std::endl;
                            std::cout << usage.c_str();
                            exit(1);
                        }

                        break;

                    case 5 :

                        //std::cout << "Min amplification size: " << optarg << std::endl << std::endl;
                        parsed_args.min_amplification_size = std::stoi(optarg);

                        if (parsed_args.min_amplification_size <= 0) {
                            std::cout << "Alpha amplification size value must be > 0." << std::endl << std::endl;
                            std::cout << usage.c_str();
                            exit(1);
                        }

                        break;

                    case 6 :

                        //std::cout << "Max amplification size: " << optarg << std::endl << std::endl;
                        parsed_args.max_amplification_size = std::stoi(optarg);

                        break;

                    case 7 :

                        //std::cout << "Alpha HOR order: " << optarg << std::endl << std::endl;
                        parsed_args.alpha_HOR_order = std::stod(optarg);

                        if (parsed_args.alpha_HOR_order <= 0) {
                            std::cout << "Alpha HOR order value must be > 0." << std::endl << std::endl;
                            std::cout << usage.c_str();
                            exit(1);
                        }

                        break;

                    case 8 :

                        //std::cout << "Beta HOR order: " << optarg << std::endl << std::endl;
                        parsed_args.beta_HOR_order = std::stod(optarg);

                        if (parsed_args.beta_HOR_order <= 0) {
                            std::cout << "Beta HOR order value must be > 0." << std::endl << std::endl;
                            std::cout << usage.c_str();
                            exit(1);
                        }

                        break;

                    case 9 :

                        //std::cout << "Min HOR order: " << optarg << std::endl << std::endl;
                        parsed_args.min_HOR_order = std::stoi(optarg);

                        if (parsed_args.min_HOR_order <= 0) {
                            std::cout << "Min HOR order value must be > 0." << std::endl << std::endl;
                            std::cout << usage.c_str();
                            exit(1);
                        }

                        break;

                    case 10 :

                        //std::cout << "Max HOR order: " << optarg << std::endl << std::endl;
                        parsed_args.max_HOR_order = std::stoi(optarg);

                        break;

                    case 11 :

                        //std::cout << "Monomers dataset output path: " << optarg << std::endl << std::endl;
                        dir_dataset = std::filesystem::path(optarg).parent_path();
                        if (std::filesystem::exists(dir_dataset) && std::filesystem::is_directory(dir_dataset)) {
                            if (std::filesystem::exists(optarg)) {
                                std::cout << "Pre-existing file found. It will be replaced by the newly generated output file." << std::endl << std::endl;
                            }
                            parsed_args.monomers_dataset_output_path = optarg;
                        
                        }else{
                            std::cout << "No such monomers dataset output directory." << std::endl << std::endl;
                            std::cout << usage.c_str();
                            exit(1);
                        }

                        break;

                    case 12 :

                        //std::cout << "Tree output path: " << optarg << std::endl << std::endl;
                        dir_tree = std::filesystem::path(optarg).parent_path();
                        if (std::filesystem::exists(dir_tree) && std::filesystem::is_directory(dir_tree)) {
                            if (std::filesystem::exists(optarg)) {
                                std::cout << "Pre-existing file found. It will be replaced by the newly generated output file." << std::endl << std::endl;
                            }
                            parsed_args.tree_output_path = optarg;
                        
                        }else{
                            std::cout << "No such tree output directory." << std::endl << std::endl;
                            std::cout << usage.c_str();
                            exit(1);
                        }

                        break;

                    default:

                        exit(1);

                }

                break;

            case '?':
            case 'h':

                std::cout << std::endl;
                std::cout << usage.c_str();
                
                exit(1);

            default:

                exit(1);
        }

        if (parsed_args.input.empty()){
            std::cout << "Input argument is required.\n" << std::endl;
            std::cout << usage.c_str();
            exit(1);
        }

        if (parsed_args.min_amplification_size > parsed_args.max_amplification_size) {
            std::cout << "Max amplification size value must be >= min amplification size." << std::endl << std::endl;
            std::cout << usage.c_str();
            exit(1);
        }

        if (parsed_args.min_HOR_order > parsed_args.max_HOR_order) {
            std::cout << "Max HOR order value must be >= min HOR order." << std::endl << std::endl;
            std::cout << usage.c_str();
            exit(1);
        }

        if (parsed_args.monomers_dataset_output_path.empty()){
            parsed_args.monomers_dataset_output_path = parsed_args.input + ".monomers_dataset_simulation_" + parsed_args.simulation_id + ".fst" ;
        }

        if (parsed_args.tree_output_path.empty()){
            parsed_args.tree_output_path = parsed_args.input + ".phylogenic_tree_simulation_" + parsed_args.simulation_id + ".tree" ;
        }

    }

    std::cout << "Input: " << parsed_args.input << std::endl;
    std::cout << "Max size: " << parsed_args.max_size << std::endl;
    std::cout << "Simulation Id: " << parsed_args.simulation_id << std::endl;
    std::cout << "Amplification rate: " << parsed_args.amplification_rate << std::endl;

    std::cout << "Min amplification size: " << parsed_args.min_amplification_size << std::endl;
    std::cout << "Max amplification size: " << parsed_args.max_amplification_size << std::endl;
    std::cout << "Alpha amplification size: " << parsed_args.alpha_amplification_size << std::endl;
    std::cout << "Beta amplification size: " << parsed_args.beta_amplification_size << std::endl;

    std::cout << "Min HOR order: " << parsed_args.min_HOR_order << std::endl;
    std::cout << "Max HOR order: " << parsed_args.max_HOR_order << std::endl;
    std::cout << "Alpha HOR order: " << parsed_args.alpha_HOR_order << std::endl;
    std::cout << "Beta HOR order: " << parsed_args.beta_HOR_order << std::endl;

    std::cout << "Monomers dataset output path: " << parsed_args.monomers_dataset_output_path << std::endl;
    std::cout << "Tree output path: " << parsed_args.tree_output_path << std::endl;

    std::cout << "Verbosity level: " << parsed_args.verbose << std::endl;

    return parsed_args;
}

/** @brief Main function
*
* Main function in charge of the overall execution of the script including the execution of the argument parsing function , the simulation functions and the output generation functions.
*
* @param argc Number of arguments provided by the user in the command line.
* @param argv List of arguments provided by the user in the command line.
*
* @return 0 if normal execution, other number if error.
*/
int main(int argc, char **argv) {

    Args args = get_args(argc, argv);

    

    return 0;
}