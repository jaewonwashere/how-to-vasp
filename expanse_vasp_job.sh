#!/bin/bash -l

# This code submits VASP jobs on the SDSC Expanse system.

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --account=<put in your account code>
#SBATCH -t 00:25:00
#SBATCH -o slurm.%N.%j.o
#SBATCH -e slurm.%N.%j.e
#SBATCH --mail-user=jaewon_lee@ucsb.edu
#SBATCH --mail-type=end
#SBATCH --partition=compute
#SBATCH --export=ALL

cd $SLURM_SUBMIT_DIR

module reset
module load cpu/0.17.3b
module load gcc/10.2.0/npcyll4
module load openmpi/4.1.3/oq3qvsv
module load vasp6/6.2.1/tqvsz4d

exe=$(which vasp_std)
export OMP_NUM_THREADS=1
mpirun --mca btl_openib_if_include "mlx5_2:1" --mca btl self,vader $exe >& vasp.log
