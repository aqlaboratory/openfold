#!/bin/bash
#SBATCH -c 8 # Number of cpu cores requested
#SBATCH -t 10 # timelimit in minutes
#SBATCH -p project-gpu # Partition where the job will run 
#SBATCH --gres=gpu:4 # How many gpus to allocate
#SBATCH -o myoutput_train.out # Standard output

# module load cuda/11.4
srun train_run.sh