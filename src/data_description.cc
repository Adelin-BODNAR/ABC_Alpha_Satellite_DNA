/**
* @file data_description.cc
* @brief Observed and simulated data description.
* @author Adelin BODNAR
* @date 13/05/2024
*
* Script used to generate statistics used to describe and compare the observed and simulated data.
*
*/

#include "../include/data_description.hh"

#include <bits/stdc++.h>

#include <algorithm> //std::min
#include <Bpp/Seq/Alphabet/AlphabetTools.h> //bpp::Alphabet, bpp::AlphabetTools
#include <Bpp/Seq/Io/Fasta.h>   //bpp::Fasta
#include <cctype> //std::tolower
#include <filesystem>   //std::filesystem::path, std::filesystem::exists, std::filesystem::is_directory
#include <getopt.h> //getopt_long 
#include <memory>   //std::shared_ptr<>
//#include <string> //std::string


#ifdef STANDALONE

/** @brief Parsed arguments structure
* 
* Structure used to store the arguments provided by the user in the command line after parsing.
*/
struct Args {

    int verbose ; /**< Integer used to determine the level of verbosity for the execution of the script */

    std::string input ; /**< String containing the path provided by the user to the input fasta file with the aligned monomers dataset */

    int nb_orders_tested; /**< Integer used to determine the number of HOR orders for which an average similarity value will be calculated */

    int nb_similarities_calculated_per_order; /**< Integer used to determine the number of monomers that each monomer will be compared to per HOR order value */

    int seed; /**< Integer variable used to initialize the random number generators. */

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

    parsed_args.nb_orders_tested = 40; /**< Initializes the default number of HOR orders for which an average similarity value will be calculated at 40. */

    parsed_args.nb_similarities_calculated_per_order = 5; /**< Initializes the default number of monomers that each monomer will be compared to per HOR order value at 5. */

    parsed_args.seed = std::time(0) ; /**< Initializes the default seed at the current time in seconds. */

    //TODO: Remove if sure that it will not be used
    //std::filesystem::path summary_statistics_dir_path ; /**< Declares the path variable used to store the path to the directory that will contain the summary statistics output file. */


    /**< String containing the help message printed when the -h , --help or a wrong option are provided ny the user */
    const std::string usage = "Usage: \n"
                            "\t./data_description [-h] -i INPUT [--nb_orders_tested NB_ORDERS_TESTED] [--nb_similarities_calculated_per_order NB_SIMILARITIES_CALCULATED_PER_ORDER] [--seed SEED] [-v | --verbose | --no-verbose]\n\n"

                            "\toptions:\n\n"

                            "\t-h, --help \n\t\tOption showing this help message and exiting\n\n"

                            "\t-i INPUT, --input INPUT\n"
                            "\t\tArgument defining the input fasta file containing the \n"
                            "\t\taligned monomers compared to generate the summary statistics \n\n"

                            "\t--seed SEED \n"
                            "\t\tArgument defining the seed used to initialize the random \n"
                            "\t\tnumber generator \n\n"

                            "\t--nb_orders_tested NB_ORDERS_TESTED \n"
                            "\t\tArgument defining the number of HOR orders for which an \n"
                            "\t\taverage similarity value will be calculated \n\n"

                            "\t--nb_similarities_calculated_per_order NB_SIMILARITIES_CALCULATED_PER_ORDER \n"
                            "\t\tArgument defining the number of monomers that each monomer will \n"
                            "\t\tbe compared to per HOR order value \n\n"

                            //TODO: Remove if sure that it will not be used
                            //"\t--sum_stats_output_path SUM_STATS_OUTPUT_PATH \n"
                            //"\t\tArgument defining the path to the phylogenic tree \n"
                            //"\t\toutput file \n\n"

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

          {"seed", required_argument, 0, 0},

          {"nb_orders_tested", required_argument, 0, 0},

          {"nb_similarities_calculated_per_order", required_argument, 0, 0},

          //TODO: Remove if sure that it will not be used
          //{"sum_stats_output_path", required_argument, 0, 0},

          //  These options set flags. 

          {"help", no_argument, 0, 'h'},

          {"input", required_argument, 0, 'i'},

          {"verbose", no_argument, 0, 'v'},

          {0, 0, 0, 0}
        };

    // getopt_long stores the option index here.
    int option_index;

    int c;

    while ((c = getopt_long (argc, argv, "hi:v", long_options, &option_index)) != -1) {

        switch (c)
        {
            case 'i':

                //std::cout << "Input: " << optarg << std::endl;
                if (std::filesystem::exists(optarg)) {
                    //std::cout << "Input file found.\n" << std::endl;
                    parsed_args.input = optarg;
                }else{
                    std::cout << "No such input file.\n" << std::endl;
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

                        //std::cout << "Seed: " << optarg << std::endl << std::endl;
                        parsed_args.seed = std::stoi(optarg);

                        break;

                    case 1 :

                        //std::cout << "Nb orders tested per monomer: " << optarg << std::endl << std::endl;
                        parsed_args.nb_orders_tested = std::stoi(optarg);

                        break;

                    case 2 :

                        //std::cout << "Nb similarities calculated per order for each monomer: " << optarg << std::endl << std::endl;
                        parsed_args.nb_similarities_calculated_per_order = std::stoi(optarg);

                        break;

                    //TODO: Remove if sure that it will not be used
                    //case 3 :

                    //    //std::cout << "Summary statistics output path: " << optarg << std::endl << std::endl;
                    //    summary_statistics_dir_path = std::filesystem::path(optarg).parent_path();
                    //    if (std::filesystem::exists(summary_statistics_dir_path) && std::filesystem::is_directory(summary_statistics_dir_path)) {
                    //        if (std::filesystem::exists(optarg)) {
                    //            std::cout << "Pre-existing file found. It will be replaced by the newly generated output file." << std::endl << std::endl;
                    //        }
                    //        parsed_args.sum_stats_output_path = optarg;
                        
                    //    }else{
                    //        std::cout << "No such tree output directory." << std::endl << std::endl;
                    //        std::cout << usage.c_str();
                    //        exit(1);
                    //    }

                    //    break;

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

    }

    if (parsed_args.input.empty()){
        std::cout << "Input argument is required.\n" << std::endl;
        std::cout << usage.c_str();
        exit(1);
    }


    /**< Prints the simulation parameters to the terminal if needed. */
    if (parsed_args.verbose > 0) {

        std::cout << "Input: " << parsed_args.input << std::endl;

        std::cout << "Seed: " << parsed_args.seed << std::endl;

        std::cout << "Nb orders tested: " << parsed_args.nb_orders_tested << std::endl;

        std::cout << "Nb similarities calculated per order: " << parsed_args.nb_similarities_calculated_per_order << std::endl;

        //TODO: Remove if sure that it will not be used
        //std::cout << "Summary statistics output path: " << parsed_args.sum_stats_output_path << std::endl;

        std::cout << "Verbosity level: " << parsed_args.verbose << std::endl << std::endl;
    }

    return parsed_args;
}

#endif





/** @brief Parsing the arguments.
*
* Function comparing two strings (supposed to be DNA sequences) to return the corresponding similarity percentage .
*
* @param input_file_path String containing the path to the input file containing the aligned monomers ordered dataset.
* @param list_sequences_ptr Pointer allowing to access the vector in which we want to store the parsed sequences from the input file as strings.
*
* @return void
*/
void parse_input_fasta(std::string input_file_path, std::vector<std::string>* list_sequences_ptr ) {

    /**< Initializes the variables used to read and parse the input fasta file */
    bpp::Fasta fasta = bpp::Fasta();
    std::shared_ptr<const bpp::Alphabet> alpha = bpp::AlphabetTools::DNA_ALPHABET;
    bpp::Sequence monomer_sequence = bpp::Sequence(alpha);
    std::ifstream monomers_input_stream(input_file_path);

    /**< Reads the sequences of the input file one by one from the input stream to add them to the variable the pointer list_sequences_ptr points to. */
    while(monomers_input_stream){
        fasta.nextSequence(monomers_input_stream, monomer_sequence);
        list_sequences_ptr->push_back(monomer_sequence.toString());
    }

}

/** @brief Calculate similarity of two sequences.
*
* Function comparing two strings (supposed to be DNA sequences) to return the corresponding similarity percentage .
*
* @param seq_1_ptr Pointer allowing to access the first sequence of the pair we want to compare.
* @param seq_2_ptr Pointer allowing to access the second sequence of the pair we want to compare.
*
* @return the similarity percentage for this pair of sequences as a decimal number between 0 and 1 both included.
*/
double get_similarity(std::string* seq_1_ptr, std::string* seq_2_ptr) {

    /**< Checks if both the sequences are the same length. Since we want to compare aligned sequences, they should be the same length. */
    if (seq_1_ptr->length() == seq_2_ptr->length()){

        /**< Initializes the double (for division purposes) variable used to store the number of differences between the two sequences at 0. */
        double hamming_distance = 0 ;

        /**< Checks if there are any gaps in either of the compared sequences. */
        if (seq_1_ptr->find('-') != std::string::npos || seq_2_ptr->find('-') != std::string::npos) {

            /**< Initializes the integer variable used to store the number of gaps found in either sequence along the way at 0. */
            int gaps = 0;

            /**< Checks for each position in the compared sequences to see if the nucleotides are identical. */
            for(int i = 0; i < seq_1_ptr->length(); i++){

                /**< Checks if there are any gaps at the current position and doesn't count any difference at that position if so but increments the gaps variable */
                if ((*seq_1_ptr)[i] == '-' || (*seq_2_ptr)[i] == '-') {
                    gaps++;
                
                /**< Checks if the nucleotides at the current position are identical and increments the hamming_distance variable if not. */
                }else if ((*seq_1_ptr)[i] != (*seq_2_ptr)[i]){
                    hamming_distance++;
                }
            }

            /**< Calculates the similarity percentage for the two compared sequences as a decimal number using the hamming_distance and the number of gaps. */
            return (seq_1_ptr->length() - gaps - hamming_distance)/(std::max(int(seq_1_ptr->length() - gaps), 1)) ;

        }else{
        
            /**< Checks for each position in the compared sequences to see if the nucleotides are identical. */
            for(int i = 0; i < seq_1_ptr->length(); i++){
                if ((*seq_1_ptr)[i] != (*seq_2_ptr)[i]) hamming_distance++;
            }

            /**< Calculates the similarity percentage for the two compared sequences as a decimal number using the hamming_distance. */
            return (seq_1_ptr->length()-hamming_distance)/(seq_1_ptr->length()) ;

        }

    /**< Exits the script with an error if the sequences are not the same length since they are supposed to be aligned.*/
    }else{

        std::cout << "Different size of sequences" << std::endl;
        exit(1);
    
    }
}

/** @brief Compare each monomer to N next monomers.
*
* Function comparing string from the list of sequences to the N next string in the list in order to return the position with maximum similarity value for each of them .
*
* @param nb_orders_tested Integer defining the number of the orders we want average similarities for.
* @param list_sequences Pointer allowing to access the list of sequences used for the comparisons.
* @param HOR_order_distribution_ptr Pointer allowing to access the list used to store the results for the HOR order distribution.
* @param max_similarity_difference_distribution_ptr Pointer allowing to access the list used to store the results for the max similarity distance distribution.
*
*/
void compare_monomers_to_next (int nb_orders_tested, int nb_similarities_calculated_per_order, std::vector<std::string>* list_sequences, std::vector<int>* HOR_order_distribution_ptr, std::vector<int>* max_similarity_difference_distribution_ptr) {

    /**< Declares the variables used to store the values allowing to generate the HOR order and max similarity difference distributions. */
    double max_average_similarity;
    int estimated_order;
    double max_similarity_difference ;

    /**< Declares the integer variable used to store the number of similarities that were really calculated for an order. */
    int nb_similarities ;

    //TODO: Make it a real option and find a better default value
    /**< Initializes the integer variable used to store the number of similarity values that should be used to calculate the average similarity for each order. */

    /**< Generates values for the HOR order and max similarity difference distributions for each monomer of list_sequences. */
    for(int i = 0; i < list_sequences->size()-1; i++){

        /**< Initializes the vector variable used to store the similarities calculated at each step and allowing to avoid recalculating similarities. Positions with no similarities calculated contain the values -1. */
        std::vector<double> similarities (nb_orders_tested * nb_similarities_calculated_per_order,-1);

        /**< Initializes the vector variable used to store the average similarities calculated for each HOR order we want. Positions with no average similarities calculated yet contain the values 0. */
        std::vector<double> average_similarities (nb_orders_tested,0);

        /**< Initializes or reinitializes estimated_order at 0 (impossible value) for each monomer. */
        estimated_order = 0;

        /**< Initializes or reinitializes max_average_similarity at -1 (impossible value) for each monomer. */
        max_average_similarity = -1;

        /**< Initializes or reinitializes max_similarity_difference at -1 (impossible value) for each monomer. */
        max_similarity_difference = -1 ;

        /**< Generates the average similarities and the differences for the nb_orders_tested */
        for (int k = 0; k < nb_orders_tested; k++){

            /**< Initializes or reinitializes the integer variable used to store the number of similarities that were really calculated for the order k+1 at 0. */
            nb_similarities = 0;

            /**< Loops through the the positions in the list of sequences that correspond to the nb_similarities_calculated_per_order next monomers after the current monomer at a distance k from the previous one. */
            for(int j = k+1+i; j < std::min(((k+1) * (nb_similarities_calculated_per_order+1))+i, int(list_sequences->size())); j += k+1){

                /**< Increments the variable storing the number of similarities really calculated for this order. */
                nb_similarities ++ ;

                /**< Checks if the similarity for this position has already been calculated and skips this step if so. */
                if (similarities[j-i-1] == -1){

                    /**< Stores the similarity calculated by the function get_similarity at the right position in the similarities list. */
                    similarities[j-i-1] = get_similarity(&((*list_sequences)[i]), &((*list_sequences)[j]));
                
                }

                /**< Adds the similarity value to the value stored in the average_similarities list at the right position to get a sum of similarities for the order k+1 that will be used afterwards to calculate the average similarity. */
                average_similarities[k] += similarities[j-i-1] ;

            }

            /**< Uses the sum stores in the average_similarities list at the right position and the number of similarities really calculated stored in nb_similarities for to calculate the average similarity for the HOR order k+1. */
            average_similarities[k] = average_similarities[k] / std::max(nb_similarities, 1);

            /**< Checks if the newly calculated average similarity for the HOR order k+1 is greater than the max average similarity encountered until now for the current monomer. */
            if (average_similarities[k] > max_average_similarity) {

                /**< Stores the estimated order of the HOR that the current monomer would be part of if the newly calculated average similarity really is the maximum average similarity for the current monomer. */
                estimated_order = k+1;

                /**< Stores the newly calculated average similarity as the maximum average similarity for the current monomer. */
                max_average_similarity = average_similarities[k];

            }

            /**< Checks if the difference between the newly calculated average similarity for the HOR order k+1 and the average similarity caculated for the order 1 is greater than the max average similarity difference encountered until now for the current monomer. */
            if ((average_similarities[k] - average_similarities[0]) > max_similarity_difference) {

                
                /**< Stores the newly calculated average similarity difference as the maximum average similarity difference for the current monomer. */
                max_similarity_difference = average_similarities[k] - average_similarities[0] ;

            }

        }

        /**< Checks if the estimated HOR order for this monomer is 0 and adds 1 to the HOR order distribution at the right position if not. */
        if (estimated_order != 0) {
            (*HOR_order_distribution_ptr)[estimated_order-1] += 1;

        /**< Exits the script with an error if the estimated HOR order for the HOR is 0. */
        }else{
            std::cout << "Estimated order error" << std::endl ;
            exit(1);
        }
        
        /**< Adds 1 to the max similarity difference distribution at the right position */
        (*max_similarity_difference_distribution_ptr)[int(max_similarity_difference * 100)] += 1;

    }
    
}

/** @brief Compare each monomer to a random other monomer.
*
* Function comparing string from the list of sequences to a random other string in the list in order to return a similarity value for each of them .
*
* @param seed Integer defining the seed used to initialize the random number generator.
* @param list_sequences Pointer allowing to access the list of sequences used for the comparisons.
* @param similarity_distribution Pointer allowing to access the list used to store the results.
*
*/
void compare_monomers_to_random (int seed, std::vector<std::string>* list_sequences, std::vector<int>* similarity_distribution ) {

    double similarity;
    int index_random_monomer;

    /**< Initializes srand using the seed provided by the user. */
    srand(seed);
    /**< Initializes a mersenne twister type random number generator seeded by a number from the srand seeded by the user. */
    std::mt19937 random_generator(rand());
    /**< Initializes a uniform distribution used to draw random indexes during the simulation. */
    std::uniform_real_distribution<double> uniform_distrib = std::uniform_real_distribution<double>(0,1);

    /**< Caculates a similarity value every monomer and a random monomer each from list_sequences. */
    for(int i = 0; i < list_sequences->size(); i++){

        /**< Draws the index of a random monomer in list_sequences and draws one again if the random monomer and the current monomer happen to be the same monomer. */
        do { index_random_monomer = int(uniform_distrib(random_generator) * (list_sequences->size())); } while ( i == index_random_monomer);

        /**< Stores the similarity between the current monomer and the random monomer returned by the get_similarity function. */
        similarity = get_similarity(&((*list_sequences)[i]), &((*list_sequences)[index_random_monomer]));

        /**< Adds 1 to the similarity distribution at the position corresponding to the previously calculated similarity. */
        (*similarity_distribution)[int(similarity * 100)] += 1;

    }

}


/** @brief Generates summary statistics for the list of sequences provided.
*
* Function launching the functions generating summary statistics for strings from the list of sequences and returning them in one line .
*
* @param nb_orders_tested Integer defining the number of the next sequences in the list we want to compare each monomer to.
* @param seed Integer defining the seed used to initialize the random number generator.
* @param list_sequences Pointer allowing to access the list of sequences used for the comparisons.
*
* @return the string containing the summary statistics generated.
*/
std::string generate_summary_stats (int nb_orders_tested, int nb_similarities_calculated_per_order, int seed, std::vector<std::string>* list_sequences ) {

    /**< Declares and initializes the string variable used to strore the summary statistics generated in one line.*/
    std::string summary_stats = "";

    /**< Declares and initializes the vector used to store the estimated HOR order distribution generated by the compare_monomeres_to_next function before adding it to the string.*/
    std::vector<int> HOR_order_distribution (nb_orders_tested,0);

    /**< Declares and initializes the vector used to store the  distribution of max difference between average similarities with each monomers and the N monomers after and similarity with the monomer after, generated by the compare_monomeres_to_next function before adding it to the string.*/
    std::vector<int> max_similarity_difference_distribution (101,0);

    /**< Generates the data summary statistics for HOR_order_distribution and max_similarity_difference_distribution */
    compare_monomers_to_next (nb_orders_tested, nb_similarities_calculated_per_order, list_sequences, &HOR_order_distribution, &max_similarity_difference_distribution);

    /**< Adds the results from HOR_order_distribution to the line used as output.*/

    summary_stats.append("order_distrib\t");

    for (int i = 0; i < HOR_order_distribution.size(); i++ ){
        summary_stats += std::to_string(HOR_order_distribution[i]) + "\t";
    }

    /**< Adds the results from max_similarity_difference_distribution to the line used as output.*/

    summary_stats.append("max_similarity_diff_distrib\t");

    for (int i = 0; i < max_similarity_difference_distribution.size(); i++ ){
        summary_stats += std::to_string(max_similarity_difference_distribution[i]) + "\t";
    }

    /**< Declares and initializes the vector used to store the distribution of similarities between each monomer and a random one generated by the compare_monomeres_to_random function before adding it to the string.*/
    std::vector<int> similarity_distribution (101,0);

    /**< Generates the data summary statistics for similarity_distribution */
    compare_monomers_to_random (seed, list_sequences, &similarity_distribution);

    /**< Adds the results from similarity_distrib to the line used as output.*/

    summary_stats.append("similarity_distrib\t");

    for (int i = 0; i < similarity_distribution.size(); i++ ){
        summary_stats += std::to_string(similarity_distribution[i]) + "\t";
    }

    return summary_stats;

}


#ifdef STANDALONE

/** @brief Main function
*
* Main function in charge of the overall execution of the script.
*
* @param argc Number of arguments provided by the user in the command line.
* @param argv List of arguments provided by the user in the command line.
*
* @return 0 if normal execution, other number if error.
*/
int main(int argc, char **argv) {

    /**< Parses the arguments passed in the command line used to launch the script.*/
    Args args = get_args(argc, argv);

    /**< Declares the vector variable used to store the sequences from the fasta file passed as input. */
    std::vector<std::string> list_sequences;

    /**< Parses the input fasta file and stores all the sequences in list_sequences. */
    parse_input_fasta(args.input, &list_sequences);

    /**< Launches the generation of all the summary stats and stores them in a string variable. */
    std::string summary_stats = generate_summary_stats(args.nb_orders_tested, args.nb_similarities_calculated_per_order, args.seed, &list_sequences);

    /**< Prints the string containing the summary statistics. */
    std::cout << summary_stats << std::endl;

    

    return 0;
}

#endif