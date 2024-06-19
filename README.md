# ABC_Alpha_Satellite_DNA
"An Approximate Bayesian Computation approach for modeling the evolution of centromeric DNA in Primat" - M2 Internship - 2024 - MNHN

## Simulation

### Base command

```bash
./simulation -i INPUT
```

### Usage

```
./simulation [-h] -i INPUT [-s MAX_SIZE] [--id SIMULATION_ID] [--amplification_rate AMPLIFICATION_RATE] 
    [--alpha_amplification_size ALPHA_AMPLIFICATION_SIZE] [--beta_amplification_size BETA_AMPLIFICATION_SIZE] 
    [--min_amplification_size MIN_AMPLIFICATION_SIZE] [--max_amplification_size MAX_AMPLIFICATION_SIZE]
    [--alpha_HOR_order ALPHA_HOR_ORDER] [--beta_HOR_order BETA_HOR_ORDER] 
    [--min_HOR_order MIN_HOR_ORDER] [--max_HOR_order MAX_HOR_ORDER] 
    [--substitution_rate SUBSTITUTION_RATE] [--transition_transversion_ratio TRANSITION_TRANSVERSION_RATIO] 
    [--seed SEED] [--output_option OUTPUT_OPTION] 
    [--parameters_log_path PARAMETERS_LOG_PATH] [--mono mers_dataset_output_path MONOMERS_DATASET_OUTPUT_PATH] 
    [--tree_output_path TREE_OUTPUT_PATH] [-v | --verbose | --no-verbose]

    options:

    -h, --help
        Option showing this help message and exiting

    -i INPUT, --input INPUT
        Argument defining the input fasta file containing the
        ancestor monomer used to perform the simulation

    -s MAX_SIZE, --max_size MAX_SIZE
        Argument defining the max size of the monomers dataset
        before the simulation stops

    --id SIMULATION_ID, --simulation_id SIMULATION_ID
        Argument defining the ID used to identify the current
        run of the simulation

    --amplification_rate AMPLIFICATION_RATE
        Argument defining the rate at which monomers
        are amplified in the simulation

    --alpha_amplification_size ALPHA_AMPLIFICATION_SIZE
        Argument defining the alpha of the law allowing to
        draw randomly the number of new copies each monomer
        involved gets during an amplification event

    --beta_amplification_size BETA_AMPLIFICATION_SIZE
        Argument defining the beta of the law allowing to draw
        randomly the number of new copies each monomer
        involved gets during an amplification event

    --min_amplification_size MIN_AMPLIFICATION_SIZE
        Argument defining the minimum of the values from which
        we have to draw randomly the number of new copies each
        monomer involved gets during an amplification event

    --max_amplification_size MAX_AMPLIFICATION_SIZE
        Argument defining the maximum of the values from which
        we have to draw randomly the number of new copies each
        monomer involved gets during an amplification event

    --alpha_HOR_order ALPHA_HOR_ORDER
        Argument defining the alpha for the law allowing to
        draw the number of adjacent monomers involved in the
        an amplification event

    --beta_HOR_order BETA_HOR_ORDER
        Argument defining the beta for the law allowing to
        draw the number of adjacent monomers involved in the
        an amplification event

    --min_HOR_order MIN_HOR_ORDER
        Argument defining the minimum of the values from which
        we have to draw the number of adjacent monomers
        involved in the an amplification event

    --max_HOR_order MAX_HOR_ORDER
        Argument defining the maximum of the values from which
        we have to draw the number of adjacent monomers
        involved in the an amplification event

    --substitution_rate SUBSTITUTION_RATE
        Argument defining the rate at which nucleotides are substituted
        in the sequence of a monomer during the simulation

    --transition_transversion_ratio TRANSITION_TRANSVERSION_RATIO
        Argument defining the transition/transversion ratio used to initialize the matrix
        necessary to perform substitutions in the sequence of a monomer during the simulation
        (== 1 is equivalent to the Jukes-Cantor model, != 1 is equivalent to the Kimura 80 model)

    --seed SEED
         Argument defining the seed used to initialize the random
         number generator

    --output_option OUTPUT_OPTION
         Argument defining the output option given by the user and
         used to determine what outputs must be produced
         (0 : Monomers dataset only, 1 : Monomers dataset and phylogenetic tree,
          2 : Monomers dataset and simulation parameters, 3 : Monomers dataset, phylogenetic tree and simulation parameters)

    --parameters_log_path PARAMETERS_LOG_PATH
          Argument defining the path to the log file containing
          the values chosen for the simulation parameters

    --monomers_dataset_output_path MONOMERS_DATASET_OUTPUT_PATH
          Argument defining the path to the monomers dataset
          output file

    --tree_output_path TREE_OUTPUT_PATH
          Argument defining the path to the phylogenic tree
          output file

    -v, --verbose, --no-verbose
          Options allowing to print or not the script execution
          supervision prints (tree and list after each amplification event)

```



