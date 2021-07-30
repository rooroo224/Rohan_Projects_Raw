#!/usr/local_rwth/bin/zsh
#SBATCH --account=lect0045
# ask for eight tasks
#SBATCH --ntasks=8
# name the job
#SBATCH --job-name=hw2
# declare the merged STDOUT/STDERR file
#SBATCH --output=output.%J.txt

module load DEV-TOOLS scalasca scorep
scan -t $MPIEXEC $FLAGS_MPI_BATCH ./build/hw2
