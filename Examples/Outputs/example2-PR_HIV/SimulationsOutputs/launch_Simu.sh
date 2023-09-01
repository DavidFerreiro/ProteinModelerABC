#!/bin/bash
#SBATCH --mem=250GB -n 500 -t 06:00:00
module load cesga/2020 gcc/system openmpi/4.0.5_ft3 mpi4py/3.0.3

STARTTIME=$(date +%s)

echo ""
echo "_______________Empirical substitution models_______________"
echo "Starting simulations based on empirical amino acid substitution models"
srun python launch_E_Simulations.py
python launch_E_CheckSimulations.py
echo "Simulations based on empirical amino acid substitution models ended"
echo ""
echo "Calculating summary statistics of empirical amino acid substitution models simulations"
srun python launch_E_SS.py
python launch_E_Results.py
echo "Summary statistics calculations ended"
echo ""

echo "________________Fitness substitution models________________"
echo "Starting simulations based on Fitness amino acid substitution model"
python launch_Fitness_File.py
srun python launch_F_Simulations.py
python launch_F_CheckSimulations.py
echo "Simulations based on Fitness amino acid substitution model ended"
echo ""
echo "Calculating summary statistics of Fitness amino acid substitution model simulations"
srun python launch_F_SS.py
python launch_F_Results.py
echo "Summary statistics calculations ended"
echo ""

echo "________________Neutral substitution models________________"
echo "Starting simulations based on Neutral amino acid substitution model"
python launch_Neutral_File.py
srun python launch_N_Simulations.py
python launch_N_CheckSimulations.py
echo "Simulations based on Neutral amino acid substitution model ended"
echo ""
echo "Calculating summary statistics of Neutral amino acid substitution model simulations"
srun python launch_N_SS.py
python launch_N_Results.py
echo "Summary statistics calculations ended"
echo ""

echo "_____________________Executing ABC analysis_____________________"

python order.py
ENDTIME=$(date +%s)
echo "ProteinModelerABC finished after $(($ENDTIME - $STARTTIME)) seconds"

