#!/bin/bash
#SBATCH --mem=250GB -n 500 -t 01:30:00

module load cesga/2020 gcc/system openmpi/4.0.5_ft3 mpi4py/3.0.3

STARTTIME_F=$(date +%s)
echo -e "\n"

echo "_______________HIVb model_______________"
echo "Starting simulations based on 'HIVb' amino acid substitution model"
STARTTIME_0=$(date +%s)
srun python -c "import bin.Scripts.Functions; bin.Scripts.Functions.Cluster_Simu(0, 10000, 'HIVb')"
srun python -c "import bin.Scripts.Functions; bin.Scripts.Functions.CheckSimu(0, 10000, 'HIVb')"
echo "Simulations based on HIVb amino acid substitution model ended"
echo -e "\n"

echo "Calculating summary statistics of HIVb amino acid substitution model simulations"
srun python -c "import bin.Scripts.Functions; bin.Scripts.Functions.Cluster_SS(0, 10000)"
echo "Summary statistics calculations of HIVb amino acid substitution model simulations ended"
ENDTIME_0=$(date +%s)
echo -e "\n\n"

echo "_______________Fitness model_______________"
echo "Starting simulations based on 'Fitness' amino acid substitution model"
STARTTIME_1=$(date +%s)
srun python -c "import bin.Scripts.Functions; bin.Scripts.Functions.Cluster_Simu(10000, 20000, 'Fitness')"
srun python -c "import bin.Scripts.Functions; bin.Scripts.Functions.CheckSimu(10000, 20000, 'Fitness')"
echo "Simulations based on Fitness amino acid substitution model ended"
echo -e "\n"

echo "Calculating summary statistics of Fitness amino acid substitution model simulations"
srun python -c "import bin.Scripts.Functions; bin.Scripts.Functions.Cluster_SS(10000, 20000)"
echo "Summary statistics calculations of Fitness amino acid substitution model simulations ended"
ENDTIME_1=$(date +%s)
echo -e "\n\n"

echo "_______________Neutral model_______________"
echo "Starting simulations based on 'Neutral' amino acid substitution model"
STARTTIME_2=$(date +%s)
srun python -c "import bin.Scripts.Functions; bin.Scripts.Functions.Cluster_Simu(20000, 30000, 'Neutral')"
srun python -c "import bin.Scripts.Functions; bin.Scripts.Functions.CheckSimu(20000, 30000, 'Neutral')"
echo "Simulations based on Neutral amino acid substitution model ended"
echo -e "\n"

echo "Calculating summary statistics of Neutral amino acid substitution model simulations"
srun python -c "import bin.Scripts.Functions; bin.Scripts.Functions.Cluster_SS(20000, 30000)"
echo "Summary statistics calculations of Neutral amino acid substitution model simulations ended"
ENDTIME_2=$(date +%s)
echo -e "\n\n"

echo "_____________________Executing ABC analysis_____________________"
python -c "import bin.Scripts.Functions; bin.Scripts.Functions.Order()"

ENDTIME_F=$(date +%s)
echo "ProteinModelerABC has finished!!"

echo "ProteinModelerABC model 1 finished after $(($ENDTIME_0 - $STARTTIME_0)) seconds"

echo "ProteinModelerABC model 2 finished after $(($ENDTIME_1 - $STARTTIME_1)) seconds"

echo "ProteinModelerABC model 3 finished after $(($ENDTIME_2 - $STARTTIME_2)) seconds"

echo "ProteinModelerABC finished after $(($ENDTIME_F - $STARTTIME_F)) seconds"

