#!/bin/bash
#SBATCH --mem=150GB -n 500 -t 01:00:00
module load cesga/2020 gcc/system openmpi/4.0.5_ft3 mpi4py/3.0.3

srun python launch_E_Simulations.py
python launch_E_CheckSimulations.py
srun python launch_E_SS.py
python launch_E_Results.py

python launch_Fitness_File.py
srun python launch_F_Simulations.py
python launch_F_CheckSimulations.py
srun python launch_F_SS.py
python launch_F_Results.py

python launch_Neutral_File.py
srun python launch_N_Simulations.py
python launch_N_CheckSimulations.py
srun python launch_N_SS.py
python launch_N_Results.py


python order.py

