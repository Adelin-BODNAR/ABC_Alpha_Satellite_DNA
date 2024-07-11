#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=abc_run

### Requirements
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

### Email
#SBATCH --mail-user=adelin.bodnar@mnhn.fr
#SBATCH --mail-type=ALL

### Output
#SBATCH --output=/shared/projects/abc_execution/run_logs/%x_%j.out
#SBATCH --account=abc_execution




################################ Job execution #################################

# Printing job ID for info in log file
echo "Job ID: " $SLURM_JOB_ID

INPUT_PATH=$1

# Printing path to input file used by simulations (consensus) for info in log file
echo "Input path:" $INPUT_PATH

# Checks if parameters file exists
if [[ ! -f "$INPUT_PATH" ]]; then
    echo "Input file not found!\nPlease pass the path to a fasta file containing a sequence to be used as ancestor."
    exit 1
fi

OUTPUT_PATH=$INPUT_PATH".simulation_parameters_and_data_description_"$SLURM_JOB_ID".tab"

# Printing path to file used to store output distributions of each simulation from this job for info in log file
echo "Output path:" $OUTPUT_PATH

PARAM_FILE_PATH=$2

# Printing path to file containing all parameters for each simulation in this job for info in log file
echo "Parameters file path:" $PARAM_FILE_PATH

# Checks if parameters file exists
if [[ ! -f "$PARAM_FILE_PATH" ]]; then
    echo "Parameters file not found!"
    exit 1
fi

# Initializing the counter used to generate an ID for each simulation
counter=0

# Loop reading each line of parameters file and launching a simulaton for parameters on that line
while IFS= read -r line
do

	# Every 100 simulations generates the monomers dataset and the phylogenetic tree output files in addition to the regular output file
	if ((counter % 100 == 0)); then 
		OUTPUT_OPTION=1 

	# Only generates the resgular output file for the rest of the simulations
	else 
		OUTPUT_OPTION=-1 
	fi

	# Prints the parameters used before launching the simulation itself
	echo "--id "$SLURM_JOB_ID"_"$counter" "$line

	# Launches the simulation and times it out after 5 minutes if not over of prints the time spent for the execution otherwise
	if ! ( timeout 5m time /shared/projects/abc_execution/simulation_data_description -i $INPUT_PATH --output_option $OUTPUT_OPTION --final_output_path $OUTPUT_PATH --id $SLURM_JOB_ID"_"$counter --monomers_dataset_output_path $INPUT_PATH".monomers_dataset_"$SLURM_JOB_ID"_"$counter".fst" --tree_output_path $INPUT_PATH".tree_"$SLURM_JOB_ID"_"$counter".tree" $line )
	
	# Prints a message if the simulation has been stopped by timeout
	then echo "FAILED"
	fi

	# Increments the counter used for ID generation
	counter=$((counter + 1))

# Passes all the lines from the parameters file to be used by the loop
done < "$PARAM_FILE_PATH"
