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

echo "Job ID: " $SLURM_JOB_ID

INPUT_PATH=$1

echo "Input path:" $INPUT_PATH

OUTPUT_PATH=$INPUT_PATH".simulation_parameters_and_data_description_"$SLURM_JOB_ID".tab"

echo "Output path:" $OUTPUT_PATH

PARAM_FILE_PATH=$2

echo "Parameters file path:" $PARAM_FILE_PATH

if [[ ! -f "$PARAM_FILE_PATH" ]]; then
    echo "Parameters file not found!"
    exit 1
fi

counter=0

while IFS= read -r line
do

	if ((counter % 100 == 0)); then 
		OUTPUT_OPTION=1 
	else 
		OUTPUT_OPTION=-1 
	fi

	echo $line
	time /shared/projects/abc_execution/simulation_data_description -i $INPUT_PATH --output_option $OUTPUT_OPTION --final_output_path $OUTPUT_PATH --id $counter --monomers_dataset_output_path $INPUT_PATH".monomers_dataset_"$SLURM_JOB_ID"_"$counter".fst" --tree_output_path $INPUT_PATH".tree_"$SLURM_JOB_ID"_"$counter".tree" $line

	counter=$((counter + 1))

done < "$PARAM_FILE_PATH"
