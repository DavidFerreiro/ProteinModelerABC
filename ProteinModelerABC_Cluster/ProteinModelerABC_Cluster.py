import os
import bin.Scripts.Variables
import bin.Scripts.LeerSettings
import bin.Scripts.Errores
import bin.Scripts.ChangeVariablesPE
import bin.Scripts.Functions

bin.Scripts.Variables.init()

bin.Scripts.LeerSettings.leer()

bin.Scripts.Errores.err()

bin.Scripts.ChangeVariablesPE.chavar()

bin.Scripts.Functions.RealSS(bin.Scripts.Variables.NameOfPhylipFile)

bin.Scripts.Functions.Head()

for x in range(0, bin.Scripts.Variables.NumberOfProcessors):
    rest = x
    val = bin.Scripts.Variables.NumberOfSimulations
    val = val - x
    remainder = int(val) % int(bin.Scripts.Variables.NumberOfProcessors)
    is_divisible = remainder == 0
    if is_divisible == True:
        num = int(val) / int(bin.Scripts.Variables.NumberOfProcessors)
        break

if bin.Scripts.Variables.LongEM != 0:
    with open('launch_E_Simulations.py', 'w') as lE:
        lE.write('import bin.Scripts.Functions\n\n')
        lE.write('bin.Scripts.Functions.Cesga_Empirical_Simu(' + str(bin.Scripts.Variables.Total_E_simu) + ')\n')
        lE.close()

    with open('launch_E_CheckSimulations.py', 'w') as cE:
        cE.write('import bin.Scripts.Functions\n\n')
        cE.write('bin.Scripts.Functions.CheckSimuE(0, ' + str(bin.Scripts.Variables.Total_E_simu) + ')\n')
        cE.close()

    with open('launch_E_SS.py', 'w') as sE:
        sE.write('import bin.Scripts.Functions\n\n')
        sE.write('bin.Scripts.Functions.Cesga_Empirical_SS(' + str(bin.Scripts.Variables.Total_E_simu) + ')\n')
        sE.close()

    with open('launch_E_Results.py', 'w') as rE:
        rE.write('import bin.Scripts.Functions\n')
        rE.write('import multiprocessing\n\n')
        rE.write('if __name__ == "__main__":\n')
        rE.write('  a_pool = multiprocessing.Pool(' + str(bin.Scripts.Variables.NumberOfProcessors) + ')\n')
        rE.write('  result = a_pool.map(bin.Scripts.Functions.SummaryResultE, range(0, ' + str(bin.Scripts.Variables.Total_E_simu) + '), chunksize=1)\n')
        rE.write('  a_pool.close()\n')
        rE.write('  a_pool.join()\n')
        rE.close()

for x in bin.Scripts.Variables.SubstitutionModel.split(' '):
    if x == 'Fitness':
        with open('launch_Fitness_File.py', 'w') as fF:
            fF.write('import bin.Scripts.Functions\n\n')
            fF.write('bin.Scripts.Functions.FitnessFile()\n')
            fF.close()

        with open('launch_F_Simulations.py', 'w') as lF:
            lF.write('import bin.Scripts.Functions\n\n')
            lF.write('bin.Scripts.Functions.Cesga_Fitness_Simu(' + str(bin.Scripts.Variables.NumberOfSimulations) + ')\n')
            lF.close()

        with open('launch_F_CheckSimulations.py', 'w') as cF:
            cF.write('import bin.Scripts.Functions\n\n')
            cF.write('bin.Scripts.Functions.CheckSimuF(' + str(bin.Scripts.Variables.Total_E_simu) + ', ' + str(int(bin.Scripts.Variables.Total_E_simu) + int(bin.Scripts.Variables.NumberOfSimulations)) + ')\n')
            cF.close()

        with open('launch_F_SS.py', 'w') as sF:
            sF.write('import bin.Scripts.Functions\n\n')
            sF.write('bin.Scripts.Functions.Cesga_Fitness_SS(' + str(bin.Scripts.Variables.NumberOfSimulations) + ')\n')
            sF.close()

        with open('launch_F_Results.py', 'w') as rF:
            rF.write('import bin.Scripts.Functions\n')
            rF.write('import multiprocessing\n\n')
            rF.write('if __name__ == "__main__":\n')
            rF.write('  a_pool = multiprocessing.Pool(' + str(bin.Scripts.Variables.NumberOfProcessors) + ')\n')
            rF.write('  result = a_pool.map(bin.Scripts.Functions.SummaryResultF, range(' + str(bin.Scripts.Variables.Total_E_simu) + ', ' + str(int(bin.Scripts.Variables.Total_E_simu) + int(bin.Scripts.Variables.NumberOfSimulations)) + '), chunksize=1)\n')
            rF.write('  a_pool.close()\n')
            rF.write('  a_pool.join()\n')
            rF.close()

    elif x == 'Neutral':
        with open('launch_Neutral_File.py', 'w') as fN:
            fN.write('import bin.Scripts.Functions\n\n')
            fN.write('bin.Scripts.Functions.NeutralFile()\n')
            fN.close()

        with open('launch_N_Simulations.py', 'w') as lN:
            lN.write('import bin.Scripts.Functions\n\n')
            lN.write('bin.Scripts.Functions.Cesga_Neutral_Simu(' + str(bin.Scripts.Variables.NumberOfSimulations) + ')\n')
            lN.close()

        with open('launch_N_CheckSimulations.py', 'w') as cN:
            cN.write('import bin.Scripts.Functions\n\n')
            cN.write('bin.Scripts.Functions.CheckSimuN(' + str(bin.Scripts.Variables.Total_E_simu + int(bin.Scripts.Variables.NumberOfSimulations)) + ', ' + str(int(bin.Scripts.Variables.Total_E_simu) + int(bin.Scripts.Variables.NumberOfSimulations) + int(bin.Scripts.Variables.NumberOfSimulations)) + ')\n')
            cN.close()

        with open('launch_N_SS.py', 'w') as sN:
            sN.write('import bin.Scripts.Functions\n\n')
            sN.write('bin.Scripts.Functions.Cesga_Neutral_SS(' + str(bin.Scripts.Variables.NumberOfSimulations) + ')\n')
            sN.close()

        with open('launch_N_Results.py', 'w') as rN:
            rN.write('import bin.Scripts.Functions\n')
            rN.write('import multiprocessing\n\n')
            rN.write('if __name__ == "__main__":\n')
            rN.write('  a_pool = multiprocessing.Pool(' + str(bin.Scripts.Variables.NumberOfProcessors) + ')\n')
            rN.write('  result = a_pool.map(bin.Scripts.Functions.SummaryResultN, range(' + str(int(bin.Scripts.Variables.Total_E_simu) + int(bin.Scripts.Variables.NumberOfSimulations)) + ', ' + str(int(bin.Scripts.Variables.Total_E_simu) + int(bin.Scripts.Variables.NumberOfSimulations) + int(bin.Scripts.Variables.NumberOfSimulations)) + '), chunksize=1)\n')
            rN.write('  a_pool.close()\n')
            rN.write('  a_pool.join()\n')
            rN.close()

if rest != 0:
    with open('launch_Rest_E_Simulations.py', 'w') as lRE:
        lRE.write('import bin.Scripts.Functions\n\n')
        lRE.write('bin.Scripts.Functions.Cesga_Rest_Empirical_Simu()\n')
        lRE.close()

    with open('launch_Rest_E_SS.py', 'w') as sRE:
        sRE.write('import bin.Scripts.Functions\n\n')
        sRE.write('bin.Scripts.Functions.Cesga_Rest_Empirical_SS()\n')
        sRE.close()

    with open('launch_Rest_F_Simulations.py', 'w') as lRF:
        lRF.write('import bin.Scripts.Functions\n\n')
        lRF.write('bin.Scripts.Functions.Cesga_Rest_Fitness_Simu()\n')
        lRF.close()

    with open('launch_Rest_F_SS.py', 'w') as sRF:
        sRF.write('import bin.Scripts.Functions\n\n')
        sRF.write('bin.Scripts.Functions.Cesga_Rest_Fitness_SS()\n')
        sRF.close()

    with open('launch_Rest_N_Simulations.py', 'w') as lRN:
        lRN.write('import bin.Scripts.Functions\n\n')
        lRN.write('bin.Scripts.Functions.Cesga_Rest_Neutral_Simu()\n')
        lRN.close()

    with open('launch_Rest_N_SS.py', 'w') as sRN:
        sRN.write('import bin.Scripts.Functions\n\n')
        sRN.write('bin.Scripts.Functions.Cesga_Rest_Neutral_SS()\n')
        sRN.close()

with open('launch_Simu.sh', 'w') as ls:
    ls.write('#!/bin/bash\n')
    ls.write('#SBATCH --mem=150GB -n ' + str(bin.Scripts.Variables.NumberOfProcessors) + ' -t 06:00:00\n')
    ls.write('module load cesga/2020 gcc/system openmpi/4.0.5_ft3 mpi4py/3.0.3\n\n')
    ls.write('STARTTIME=$(date +%s)\n\n')
    ls.write('echo ""\n')
    if bin.Scripts.Variables.LongEM != 0:
        ls.write('echo "_______________Empirical substitution models_______________"\n')
        ls.write('echo "Starting simulations based on empirical amino acid substitution models"\n')
        ls.write('srun python launch_E_Simulations.py\n')
        if rest != 0:
            ls.write('srun --mem=150GB -n ' + str(rest * bin.Scripts.Variables.LongEM) + ' -t 02:00:00 python launch_Rest_E_Simulations.py\n')
        ls.write('python launch_E_CheckSimulations.py\n')
        ls.write('echo "Simulations based on empirical amino acid substitution models ended"\n')
        ls.write('echo ""\n')
        ls.write('echo "Calculating summary statistics of empirical amino acid substitution models simulations"\n')
        ls.write('srun python launch_E_SS.py\n')
        if rest != 0:
            ls.write('srun --mem=150GB -n ' + str(rest * bin.Scripts.Variables.LongEM) + ' -t 02:00:00 python launch_Rest_E_SS.py\n')
        ls.write('python launch_E_Results.py\n')
        ls.write('echo "Summary statistics calculations ended"\n')
        ls.write('echo ""\n\n')

    for x in bin.Scripts.Variables.SubstitutionModel.split(' '):
        if x == 'Fitness':
            ls.write('echo "________________Fitness substitution models________________"\n')
            ls.write('echo "Starting simulations based on Fitness amino acid substitution model"\n')
            ls.write('python launch_Fitness_File.py\n')
            ls.write('srun python launch_F_Simulations.py\n')
            if rest != 0:
                ls.write('srun --mem=150GB -n ' + str(rest) + ' -t 02:00:00 python launch_Rest_F_Simulations.py\n')
            ls.write('python launch_F_CheckSimulations.py\n')
            ls.write('echo "Simulations based on Fitness amino acid substitution model ended"\n')
            ls.write('echo ""\n')
            ls.write('echo "Calculating summary statistics of Fitness amino acid substitution model simulations"\n')
            ls.write('srun python launch_F_SS.py\n')
            if rest != 0:
                ls.write('srun --mem=150GB -n ' + str(rest) + ' -t 02:00:00 python launch_Rest_F_SS.py\n')
            ls.write('python launch_F_Results.py\n')
            ls.write('echo "Summary statistics calculations ended"\n')
            ls.write('echo ""\n\n')

        if x == 'Neutral':
            ls.write('echo "________________Neutral substitution models________________"\n')
            ls.write('echo "Starting simulations based on Neutral amino acid substitution model"\n')
            ls.write('python launch_Neutral_File.py\n')
            ls.write('srun python launch_N_Simulations.py\n')
            if rest != 0:
                ls.write('srun --mem=150GB -n ' + str(rest) + ' -t 02:00:00 python launch_Rest_N_Simulations.py\n')
            ls.write('python launch_N_CheckSimulations.py\n')
            ls.write('echo "Simulations based on Neutral amino acid substitution model ended"\n')
            ls.write('echo ""\n')
            ls.write('echo "Calculating summary statistics of Neutral amino acid substitution model simulations"\n')
            ls.write('srun python launch_N_SS.py\n')
            if rest != 0:
                ls.write('srun --mem=150GB -n ' + str(rest) + ' -t 02:00:00 python launch_Rest_N_SS.py\n')
            ls.write('python launch_N_Results.py\n')
            ls.write('echo "Summary statistics calculations ended"\n')
            ls.write('echo ""\n\n')

    ls.write('echo "_____________________Executing ABC analysis_____________________"\n')
    ls.write('\n')
    ls.write('python order.py')
    ls.write('\n')
    ls.write('ENDTIME=$(date +%s)\n')
    ls.write('echo "ProteinModelerABC finished after $(($ENDTIME - $STARTTIME)) seconds"\n')
    ls.write('\n')
    ls.close()
    
    with open('order.py', 'w') as O:
        O.write('import os\n')
        O.write('import time\n')
        O.write('import bin.Scripts.Variables\n\n')
        if bin.Scripts.Variables.SaveSimulations == 'Yes':
            O.write('os.system("tar -czf Simulations.tar.gz ./Simulations/ ")\n')
            O.write('os.system("rm -r ./Simulations")\n')
            O.write('os.system("rm -r ./Results")\n')
            O.write('os.system("rm -r ./bin/Scripts/__pycache__")\n')
            O.write('os.system("rm -r newalignment.fasta")\n')
            O.write('os.system("rm -r Local_interactions.dat")\n')
            O.write('os.system("rm -r REM.txt")\n')
            O.write('os.system("rm -r E_loc_0.txt")\n')
            O.write('os.system("rm -r seqGMRCA")\n')
            O.write('os.system("rm -r ProteinModelerABC_arguments.txt")\n')
            O.write('os.system("rm -r ProteinModelerABC_S_arguments0.txt")\n')
            O.write('os.system("rm -r ProteinModelerABC_S_arguments1.txt")\n')
            O.write("if os.path.exists('ProteinModelerABC_E-Repeat_arguments.txt') == True:\n")
            O.write('   os.system("rm -r ProteinModelerABC_E-Repeat_arguments.txt")\n')
            O.write("if os.path.exists('ProteinModelerABC_F-Repeat_arguments.txt') == True:\n")
            O.write('   os.system("rm -r ProteinModelerABC_F-Repeat_arguments.txt")\n')
            O.write("if os.path.exists('ProteinModelerABC_N-Repeat_arguments.txt') == True:\n")
            O.write('   os.system("rm -r ProteinModelerABC_N-Repeat_arguments.txt")\n')
            O.write("if os.path.exists('Rest_E_Simulations.txt') == True:\n")
            O.write('   os.system("rm -r Rest_E_Simulations.txt")\n')
            O.write("if os.path.exists('Rest_F_Simulations.txt') == True:\n")
            O.write('   os.system("rm -r Rest_F_Simulations.txt")\n')
            O.write("if os.path.exists('Rest_N_Simulations.txt') == True:\n")
            O.write('   os.system("rm -r Rest_N_Simulations.txt")\n')
            O.write('os.system("rm -r Pop_evol.in")\n')

        else:
            O.write('os.system("rm -r newalignment.fasta")\n')
            O.write('os.system("rm -r ./Results")\n')
            O.write('os.system("rm -r ./bin/Scripts/__pycache__")\n')
            O.write('os.system("rm -r Local_interactions.dat")\n')
            O.write('os.system("rm -r REM.txt")\n')
            O.write('os.system("rm -r E_loc_0.txt")\n')
            O.write('os.system("rm -r seqGMRCA")\n')
            O.write('os.system("rm -r ProteinModelerABC_arguments.txt")\n')
            O.write('os.system("rm -r ProteinModelerABC_S_arguments0.txt")\n')
            O.write('os.system("rm -r ProteinModelerABC_S_arguments1.txt")\n')
            O.write("if os.path.exists('ProteinModelerABC_E-Repeat_arguments.txt') == True:\n")
            O.write('   os.system("rm -r ProteinModelerABC_E-Repeat_arguments.txt")\n')
            O.write("if os.path.exists('ProteinModelerABC_F-Repeat_arguments.txt') == True:\n")
            O.write('   os.system("rm -r ProteinModelerABC_F-Repeat_arguments.txt")\n')
            O.write("if os.path.exists('ProteinModelerABC_N-Repeat_arguments.txt') == True:\n")
            O.write('   os.system("rm -r ProteinModelerABC_N-Repeat_arguments.txt")\n')
            O.write("if os.path.exists('Rest_E_Simulations.txt') == True:\n")
            O.write('   os.system("rm -r Rest_E_Simulations.txt")\n')
            O.write("if os.path.exists('Rest_F_Simulations.txt') == True:\n")
            O.write('   os.system("rm -r Rest_F_Simulations.txt")\n')
            O.write("if os.path.exists('Rest_N_Simulations.txt') == True:\n")
            O.write('   os.system("rm -r Rest_N_Simulations.txt")\n')
            O.write('os.system("rm -r Pop_evol.in")\n')

        if os.path.exists('ABCOutputs') == True:
            O.write('os.system("rm -r ABCOutputs")\n')
        if os.path.exists('SimulationsOutputs') == True:
            O.write('os.system("rm -r SimulationsOutputs")\n')
        O.write('os.system("mkdir ABCOutputs")\n')
        O.write('os.system("mkdir SimulationsOutputs")\n')
    
        O.write('os.system("mv ABCAnalysis.r ABCOutputs/")\n')
        O.write('os.system("cp PSimulations.txt ABCOutputs/")\n')
        O.write('os.system("cp SSSimulations.csv ABCOutputs/")\n')
        O.write('os.system("cp SSRealData.csv ABCOutputs/")\n')

        if bin.Scripts.Variables.SaveSimulations == 'Yes':
            O.write('os.system("mv Simulations.tar.gz SimulationsOutputs/")\n')
        O.write('os.system("mv PSimulations.txt SimulationsOutputs/")\n')
        O.write('os.system("mv SSSimulations.csv SimulationsOutputs/")\n')
        O.write('os.system("mv SSRealData.csv SimulationsOutputs/")\n\n')
        
        O.write('os.system("mv launch_* SimulationsOutputs/")\n')
        O.write('os.system("mv order.py SimulationsOutputs/")\n')
        #O.write('os.system("cd ABCOutputs")\n')
        O.write('os.system("Rscript --vanilla ./ABCOutputs/ABCAnalysis.r")\n')
        O.write('attempt = 1\n')
        O.write('for x in range(0, 11):\n')
        O.write('    if attempt < 10:\n')
        O.write('        with open("./ABCOutputs/Results_text.txt") as f:\n')
        O.write('             if not "Mean model posterior probabilities" in f.read():\n')
        O.write('                attempt = attempt + 1\n')
        O.write('                print("ABC estimation failed")\n')
        O.write('                print("Starting another try: " + str(x + 1))\n')
        O.write('                os.system("Rscript --vanilla ./ABCOutputs/ABCAnalysis.r")\n')
        O.write('             else:\n')
        O.write('                if x == 0:\n')
        O.write('                    print("ABC analysis has finished\\n")\n')
        O.write('                    break\n')
        O.write('                else:\n')
        O.write('                    print("Try " + str(x + 1) + " succeeded")\n')
        O.write('                    print("ABC analysis has finished\\n")\n')
        O.write('                    break\n')
        O.write('    else:\n')
        O.write('        sys.exit("\\n"\n')
        O.write('                " ______                     _ \\n"\n')
        O.write('                "|  ____|                   | |\\n"\n')
        O.write('                "| |__   _ __ _ __ ___  _ __| |\\n"\n')
        O.write('                "|  __| | |__| |__/ _ \| |__| |\\n"\n')
        O.write('                "| |____| |  | | | (_) | |  |_|\\n"\n')
        O.write('                "|______|_|  |_|  \___/|_|  (_)\\n\\n"\n')
        O.write('                "But dont panic!! You dont have to repeat everything, just go to the ./ABCOutputs and re execute ABCAnalysis.r file increasing the tolerance or changing the ABC method to rejection. Check manual for further information")\n')
        O.write('print("ProteinModelerABC has finished!!\\n")\n')
        O.close()

with open("ABCAnalysis.r", "w") as P:
        P.write('############################\n')
        P.write('suppressPackageStartupMessages(library(abc))\n')
        P.write('############################\n\n')
        P.write('#####################################################\n')
        P.write('################### ABC VARIABLES ###################\n')
        P.write('#####################################################\n')
        P.write('ABC_Method <- "' + str(bin.Scripts.Variables.ABCMethod) + '"\n')
        P.write('ABC_Tolerance <- ' + str(bin.Scripts.Variables.ABCTolerance) + '\n')
        P.write('ABC_N_Iterations <- ' + str(bin.Scripts.Variables.ABCIterations) + '\n')
        P.write('#####################################################\n')
        P.write('#####################################################\n')
        P.write('#####################################################\n\n\n')
        P.write('#Path\n')
        P.write('address<-paste("' + os.getcwd() + '/ABCOutputs", sep="")\n')
        P.write('setwd(address)\n')
        P.write('############################\n\n')
        if bin.Scripts.Variables.Coalescent == 'Coal':
            P.write('#Load Priors\n')
            P.write('matrix <- paste("PSimulations.txt", sep="")\n')
            P.write('ThisName <- paste("PSimulations",sep="")\n')
            P.write('FullPriormatrix <- read.table(matrix, header=TRUE, sep=",")\n')
            P.write('############################\n\n')
        P.write('#Load Simulations SS\n')
        P.write('matrix <- paste("SSSimulations.csv", sep="")\n')
        P.write('ThisName <- paste("SSSimulations",sep="")\n')
        P.write('FullSSmatrix <- read.csv(matrix, header=TRUE, sep=",")\n')
        P.write('col_names_matrixSS <- dimnames(FullSSmatrix)\n')
        P.write('col_names_matrixSS <- col_names_matrixSS[[2]]\n')
        P.write('ncol_matrixSS <- length(FullSSmatrix[1,])\n')
        P.write('nrow_matrixSS <- length(FullSSmatrix[,1])\n')
        P.write('############################\n\n')
        P.write('#Load Real Data SS\n')
        P.write('matrixReal <- paste("SSRealData.csv", sep="")\n')
        P.write('ThisNameReal <- paste("SSReal",sep="")\n')
        P.write('FullRealmatrix <- read.csv(matrixReal, header=TRUE, sep=",")\n')
        P.write('col_names_matrixRD <- dimnames(FullRealmatrix)\n')
        P.write('col_names_matrixRD <- col_names_matrixRD[[2]]\n')
        P.write('ncol_matrixRD <- length(FullSSmatrix[1,])\n')
        P.write('nrow_matrixRD <- length(FullSSmatrix[,1])\n')
        P.write('############################\n\n')
        P.write('# make vector with assignemt models\n')
        P.write('ModelsVector <- vector(,nrow_matrixSS)\n')
        P.write('#Assing models in the vector\n')
        P.write('j <- 0\n')
        Models = []
        for x in bin.Scripts.Variables.SubstitutionModel.split(' '):
            Models.append(x)
        for x in range(1, len(Models) + 1):
            P.write('for (j in ' + str((x - 1) * bin.Scripts.Variables.NumberOfSimulations + 1) + ':' + str(x * bin.Scripts.Variables.NumberOfSimulations) + ')\n')
            P.write('    {\n')
            P.write('    ModelsVector[j] <- "' + Models[x - 1] + '"\n')
            P.write('    }\n')

        SS_list = []
        for x in bin.Scripts.Variables.SummaryStatistics.split(' '):
            if x == '1':
                SS_list.append('DGREM_Mean')
            elif x == '2':
                SS_list.append('DGREM_sd')
            elif x == '3':
                SS_list.append('SegSites')
            elif x == '4':
                SS_list.append('Grantham_mean_Position')
            elif x == '5':
                SS_list.append('Grantham_sd_Position')
            elif x == '6':
                SS_list.append('Grantham_sk_Position')
            elif x == '7':
                SS_list.append('Grantham_ku_Position')

        P.write('############################\n\n')
        if bin.Scripts.Variables.Coalescent == 'Coal':
            P.write('#Load Priors\n')
            P.write('figureName1<-paste("./Histogram_Priors.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if bin.Scripts.Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(1,2))\n')
                P.write('if (nrow_matrixSS < 1000)	{\n')
                P.write('	hist (FullPriormatrix[,1], breaks=nrow_matrixSS, main="Substitution Rate Prior", xlab="Substitution Rate")\n')
                P.write('	hist (FullPriormatrix[,2], breaks=nrow_matrixSS, main="Theta Prior", xlab="Theta")\n')
                P.write('	} else {\n')
                P.write('	hist (FullPriormatrix[,1], breaks=1000, main="Substitution Rate Prior", xlab="Substitution Rate")\n')
                P.write('	hist (FullPriormatrix[,2], breaks=1000, main="Theta Prior", xlab="Theta")\n')
                P.write('	}\n')
            else:
                P.write('if (nrow_matrixSS < 1000)	{\n')
                P.write('	hist (FullPriormatrix[,1], breaks=nrow_matrixSS, main="Substitution Rate Prior", xlab="Substitution Rate")\n')
                P.write('	hist (FullPriormatrix[,2], breaks=nrow_matrixSS, main="Theta Prior", xlab="Theta")\n')
                P.write('	} else {\n')
                P.write('	hist (FullPriormatrix[,1], breaks=1000, main="Substitution Rate Prior", xlab="Substitution Rate")\n')
                P.write('	hist (FullPriormatrix[,2], breaks=1000, main="Theta Prior", xlab="Theta")\n')
                P.write('	}\n')
            P.write('invisible(dev.off())\n')
            P.write('############################\n\n')

        if len(SS_list) == 7:
            #P.write('figureName<-paste("Results_SS_1",ThisName,".jpeg",sep="")\n')
            #P.write('jpeg(figureName, w=1200, h=900)\n')
            P.write('figureName1<-paste("./Results_SS_Energy.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if bin.Scripts.Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(1,2))\n')
            else:
                pass
            for x in range(0,2):
                P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
                P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
            P.write('invisible(dev.off())\n\n')
            #P.write('figureName<-paste("Results_SS_2",ThisName,".jpeg",sep="")\n')
            #P.write('jpeg(figureName, w=1200, h=900)\n')
            P.write('figureName1<-paste("./Results_SS_AAReplacements.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if bin.Scripts.Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(2,3))\n')
            else:
                pass
            for x in range(2,7):
                P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
                P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
            P.write('invisible(dev.off())\n')
        elif len(SS_list) < 7 and len(SS_list) > 4:
            P.write('figureName1<-paste("./Results_SS.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if bin.Scripts.Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(2,3))\n')
            else:
                pass
            for x in range(0,len(SS_list)):
                P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
                P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
            P.write('invisible(dev.off())\n')
        elif len(SS_list) < 5 and len(SS_list) > 2:
            P.write('figureName1<-paste("./Results_SS.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if bin.Scripts.Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(2,2))\n')
            else:
                pass
            for x in range(0,len(SS_list)):
                P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
                P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
            P.write('invisible(dev.off())\n')
        else:
            P.write('figureName1<-paste("./Results_SS.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if bin.Scripts.Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(1,2))\n')
            else:
                pass
            for x in range(0,len(SS_list)):
                P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
                P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
            P.write('invisible(dev.off())\n')
        P.write('############################\n\n')

        P.write('#Create output file\n')
        P.write('unlink("Results_text.txt", recursive = FALSE)\n')
        P.write('############################\n\n')
        #P.write('cv.modsel <- NULL\n')
        #P.write('attempt <- 1\n')
        #P.write('while(is.null(cv.modsel) && attempt <= 10){\n')
        #P.write('   attempt <- attempt + 1\n')
        #P.write('   try(cv.modsel <- suppressWarnings(cv4postpr(ModelsVector, FullSSmatrix, nval=ABC_N_Iterations, tol=ABC_Tolerance, method=ABC_Method)))\n')
        #P.write('}\n')
        P.write('cv.modsel <- suppressWarnings(cv4postpr(ModelsVector, FullSSmatrix, nval=ABC_N_Iterations, tol=ABC_Tolerance, method=ABC_Method))\n')
        P.write('write(paste("Model selection with abc - cross validation - based on ' + str(bin.Scripts.Variables.ABCIterations) + ' samples",sep=""),"Results_text.txt",append=T)\n')
        P.write('try(capture.output(summary(cv.modsel), file = "Results_text.txt", append = TRUE), silent=TRUE)\n')
        P.write('write(paste("\\n\\n",sep=""),"Results_text.txt",append=T)\n')
        P.write('figureName1<-paste("./Results_ConfusionMatrix_' + str(bin.Scripts.Variables.ABCIterations) + 'samp_",ThisName,".pdf",sep="")\n')
        P.write('pdf(figureName1)\n')
        #P.write('figureName<-paste("Results_ConfusionMatrix_' + str(bin.Scripts.Variables.ABCIterations) + 'samp_",ThisName,".jpeg",sep="")\n')
        #P.write('jpeg(figureName, w=1200, h=900)\n')
        Models.sort()
        Str_Models = ''
        for x in Models:
            Str_Models = Str_Models + '"' + str(x) + '"' + ', '
        Str_Models = Str_Models[:-2]
        P.write('try(plot(cv.modsel, names.arg=c(' + Str_Models + ')), silent=TRUE)\n')
        col = []
        for y in range(1, len(Models) + 1):
            col.append(int(100/(len(Models)+1)*y))
        colors = []
        for z in col:
            colors.append("grey" + str(z))
        P.write('try(legend(x="topright", inset = c(-0.05, -0.16), legend = c(' + Str_Models + '), xpd = TRUE, bty = "n", fill = c(' + str(colors)[1:-1] + ')), silent=TRUE)\n')
        P.write('invisible(dev.off())\n')
        P.write('############################\n\n')

        P.write('#Posterior probabilities of Real data\n')
        P.write('modsel.results <- postpr(FullRealmatrix, ModelsVector, FullSSmatrix, tol=ABC_Tolerance, method=ABC_Method)\n')
        P.write('write(paste("Model selection with abc - Real data",sep=""),"Results_text.txt",append=T)\n')
        P.write('capture.output(summary(modsel.results), file = "Results_text.txt", append = TRUE)\n')
        P.write('write(paste("\\n\\n",sep=""),"Results_text.txt",append=T)\n')
        P.write('############################\n\n')

        P.write('#Goodness-of-fit of Real data\n')
        P.write('write(paste("Goodness of fit of Real data\\n",sep=""),"Results_text.txt",append=T)\n')
        P.write('figureName1<-paste("./Histogram_GoodnessOfFit.pdf",sep="")\n')
        P.write('pdf(figureName1)\n')
        if bin.Scripts.Variables.MultiPage == 'Yes':
            P.write('par(mfrow=c(1,' + str(len(Models)) + '))\n')
            for x in Models:
                P.write('res.gfit.' + str(x) + '=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="' + str(x) + '",], statistic=mean, nb.replicate=' + str(bin.Scripts.Variables.ABCIterations) + '))\n')
                P.write('write(paste("   -' + str(x) + ':",sep=""),"Results_text.txt",append=T)\n')
                P.write('capture.output(summary(res.gfit.' + str(x) + '), file = "Results_text.txt", append = TRUE)\n')
                #P.write('write(paste("\\n\\n",sep=""),"Results_text.txt",append=T)\n')
                #P.write('figureName<-paste("Hist_' + str(x) + '",ThisName,".jpeg",sep="")\n')
                #P.write('jpeg(figureName, w=1200, h=900)\n')
                P.write('plot(res.gfit.' + str(x) + ', main="Histogram under ' + str(x) + '")\n')
        else:
            for x in Models:
                P.write('res.gfit.' + str(x) + '=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="' + str(x) + '",], statistic=mean, nb.replicate=' + str(bin.Scripts.Variables.ABCIterations) + '))\n')
                P.write('write(paste("   -' + str(x) + ':",sep=""),"Results_text.txt",append=T)\n')
                P.write('capture.output(summary(res.gfit.' + str(x) + '), file = "Results_text.txt", append = TRUE)\n')
                # P.write('write(paste("\\n\\n",sep=""),"Results_text.txt",append=T)\n')
                # P.write('figureName<-paste("Hist_' + str(x) + '",ThisName,".jpeg",sep="")\n')
                # P.write('jpeg(figureName, w=1200, h=900)\n')
                P.write('plot(res.gfit.' + str(x) + ', main="Histogram under ' + str(x) + '")\n')
        P.write('invisible(dev.off())\n\n')
        P.write('############################\n\n')
        P.write('gfitpca = function(target, sumstat, index, cprob=0.1, xlim=NULL, ylim=NULL, ...){\n')
        P.write('    loc2plot = function(x, y, cprob, ...)\n')
        P.write('    {\n')
        P.write('        fit = locfit(~lp(x, y, nn=.2), maxk=200, mint=100, maxit=100)\n')
        P.write('        lev = sort(fitted(fit))[floor(cprob * length(x))]\n')
        P.write('        plot(fit, lev=lev, m=100, drawlabels=FALSE, ...)\n')
        P.write('        return (list(fit=fit, lev=lev))\n')
        P.write('    }\n\n')
        P.write('    # when target is a vector\n')
        P.write('    if (is.vector(target)){\n')
        P.write('        target = t( as.data.frame(target))\n')
        P.write('        if (is.data.frame(sumstat)){colnames(target)=names(sumstat)}\n')
        P.write('        if (is.matrix(sumstat)){colnames(target)=colnames(sumstat)}\n')
        P.write('    }\n\n')
        P.write('    # acp\n')
        P.write('    res.prcomp = prcomp(sumstat, scale=T, center=T)\n\n')
        P.write('    nmod = length(table(index))\n')
        P.write('    theindex = names(table(index))\n\n')
        P.write('    # plot\n')
        P.write('    figureName <- paste("PCA.pdf", sep="")\n')
        P.write('    pdf(figureName)\n')
        P.write('    par(mfrow=c(1, 1))\n')
        P.write('    if (is.null(xlim)){xlim = ylim}\n')
        P.write('    if (is.null(ylim)){ylim=xlim}\n\n')
        P.write('    if (! is.null(xlim)){\n')
        P.write('        plot(0, type="n", xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2")\n')
        P.write('    }\n')
        P.write('    for (i in 1:nmod){\n')
        P.write('        ind = index == theindex[i]\n')
        P.write('        if ((i == 1) & (is.null(xlim))){add=FALSE} else {add=TRUE}\n')
        P.write('        loc2plot(res.prcomp$x[ind, 1], res.prcomp$x[ind, 2], cprob, col = i, lty = 1, lwd = 2, add = add, ...)\n')
        P.write('    }\n\n')
        P.write('    # observed data\n')
        P.write('    points(predict(res.prcomp, target)[1], predict(res.prcomp, target)[2], col=1, cex=2, pch=3, lwd=2)\n\n')
        P.write('    # legend\n')
        P.write('    legend("topright", legend=theindex, cex=1.5, col=c(1: nmod), lty = 0, pch = 15)\n\n')
        P.write('    # saveplot\n')
        P.write('    invisible(dev.off())\n')
        P.write('}\n\n')

        P.write('############################\n')
        P.write('pca <- gfitpca(target=FullRealmatrix, sumstat=FullSSmatrix, index=ModelsVector, cprob=' + str(bin.Scripts.Variables.ABCTolerance) + ')\n')
        #P.write('figureName1<-paste("./PCA.pdf",sep="")\n')
        #P.write('pdf(figureName1)\n')
        #P.write('#figureName<-paste("PCA.jpeg",sep="")\n')
        #P.write('#jpeg(figureName, w=1200, h=900)\n')
        P.write('#plot(pca, main="PCA")\n')
        P.write('#invisible(dev.off())\n')
        P.write('############################\n\n')

        P.write('############################\n')
        P.write('SS_All <- matrix(nrow=' + str(bin.Scripts.Variables.NumberOfSimulations) + ')\n')
        P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
        P.write('  for(p in 1:' + str(len(Models)) + '){\n')
        P.write('    SS <- FullSSmatrix[,i][(' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * (p-1) + 1):(' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p)]\n')
        P.write('    SS_All <- cbind(SS_All, SS)\n')
        P.write('  }\n')
        P.write('}\n\n')
        P.write('SS_All <- SS_All[,-1]\n\n')
        P.write('pdf(file=paste("Histograms_SStats.pdf",sep=""))\n')
        if bin.Scripts.Variables.MultiPage == 'Yes':
            if len(SS_list) == 7:
                P.write('par(mfrow=c(4,3))\n')
                P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
                P.write('  for(p in 1:' + str(len(Models)) + '){\n')
                P.write('    hist(SS_All[,(i * ' + str(len(Models)) + ' + (-' + str(len(Models)) + ' + p))], main=ModelsVector[' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p], xlab=col_names_matrixSS[i])\n')
                P.write('    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
                P.write('  }\n')
                P.write('}\n')
            elif len(SS_list) < 7 and len(SS_list) > 4:
                P.write('par(mfrow=c(3,3))\n')
                P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
                P.write('  for(p in 1:' + str(len(Models)) + '){\n')
                P.write('    hist(SS_All[,(i * ' + str(len(Models)) + ' + (-' + str(len(Models)) + ' + p))], main=ModelsVector[' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p], xlab=col_names_matrixSS[i])\n')
                P.write('    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
                P.write('  }\n')
                P.write('}\n')
            else:
                P.write('par(mfrow=c(2,3))\n')
                P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
                P.write('  for(p in 1:' + str(len(Models)) + '){\n')
                P.write('    hist(SS_All[,(i * ' + str(len(Models)) + ' + (-' + str(len(Models)) + ' + p))], main=ModelsVector[' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p], xlab=col_names_matrixSS[i])\n')
                P.write('    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
                P.write('  }\n')
                P.write('}\n')
        else:
            P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
            P.write('  for(p in 1:' + str(len(Models)) + '){\n')
            P.write('    hist(SS_All[,(i * ' + str(len(Models)) + ' + (-' + str(len(Models)) + ' + p))], main=ModelsVector[' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p], xlab=col_names_matrixSS[i])\n')
            P.write('    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
            P.write('  }\n')
            P.write('}\n')
        P.write('invisible(dev.off())\n\n')
        #P.write('MinMax <- c()\n')
        #P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
        #P.write('  h <- min(FullSSmatrix[,i])\n')
        #P.write('  MinMax <- c(MinMax, h)\n')
        #P.write('  g <- max(FullSSmatrix[,i])\n')
        #P.write('  MinMax <- c(MinMax, g)\n')
        #P.write('}\n\n')
        #P.write('SS_Abs <- matrix(nrow=' + str(bin.Scripts.Variables.NumberOfSimulations * len(Models)) + ')\n')
        #P.write('Retained_Simu <- ' + str(bin.Scripts.Variables.NumberOfSimulations) + '*ABC_Tolerance\n')
        #P.write('SS_Tol <- matrix(nrow=ABC_N_Iterations)\n\n')
        #P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
        #P.write('  vec <- c()\n')
        #P.write('  Number <- 0\n')
        #P.write('  for(p in 1:' + str(bin.Scripts.Variables.NumberOfSimulations * len(Models)) + '){\n')
        #P.write('    Number <- abs(FullSSmatrix[p,i] - FullRealmatrix[1,i])\n')
        #P.write('    vec = c(vec, Number)\n')
        #P.write('  }\n')
        #P.write('  SS_Abs <- cbind(SS_Abs, FullSSmatrix[,i])\n')
        #P.write('  SS_Abs <- cbind(SS_Abs, Abs=vec)\n')
        #P.write('  SS_Abs <- cbind(SS_Abs, FullPriormatrix[,2])\n')
        #P.write('}\n')
        #P.write('SS_Abs <- SS_Abs[,-1]\n\n')
        #P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
        #P.write('  for(p in 1:' + str(len(Models)) + '){\n')
        #P.write('    SS <- SS_Abs[,(i * ' + str(len(Models)) + ' - ' + str(len(Models)) + '):(i * ' + str(len(Models)) + ')][(' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * (p-1) + 1):(' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p),]\n')
        #P.write('    SS <- SS[order(SS[,2]),]\n')
        #P.write('    SS_Tol <- cbind(SS_Tol, SS[1:Retained_Simu, 1])\n')
        #P.write('    SS_Tol <- cbind(SS_Tol, SS[1:Retained_Simu, ' + str(len(Models)) + '])\n')
        #P.write('  }\n')
        #P.write('}\n')
        #P.write('SS_Tol <- SS_Tol[,-1]\n\n')
        #P.write('pdf(file=paste("SummaryStatistics2.pdf",sep=""))\n')
        #if bin.Scripts.Variables.MultiPage == 'Yes':
            #P.write('par(mfrow=c(1,3))\n')
        #else:
            #pass
        #P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
        #P.write('  for(p in 1:3){\n')
        #P.write('    plot(y=SS_Tol[,(i * 3 * 2 + (-7 + p*2))], x=SS_Tol[,(i * 3 * 2 + (-6 + p*2))], main=ModelsVector[' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p], ylab=col_names_matrixSS[i], xlab="Theta")\n')
        #P.write('    abline(h = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
        #P.write('  }\n')
        #P.write('}\n')
        #P.write('invisible(dev.off())\n\n')
        #P.write('pdf(file=paste("Scaterplots_SStatsVSParams.pdf",sep=""))\n')
        #if bin.Scripts.Variables.MultiPage == 'Yes':
            #if len(SS_list) == 7:
                #P.write('par(mfrow=c(4,3))\n')
                #P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
                #P.write('  for(p in 1:' + str(len(Models)) + '){\n')
                #P.write('    plot(y=SS_Tol[,(i * ' + str(len(Models)) + ' * 2 + (-7 + p*2))], x=SS_Tol[,(i * ' + str(len(Models)) + ' * 2 + (-6 + p*2))], main=ModelsVector[' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p], ylab=col_names_matrixSS[i], xlab="Theta", ylim=c(MinMax[i*2-1], MinMax[i*2]))\n')
               # P.write('    abline(h = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
                #P.write('  }\n')
                #P.write('}\n')
            #elif len(SS_list) < 7 and len(SS_list) > 4:
                #P.write('par(mfrow=c(3,3))\n')
                #P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
                #P.write('  for(p in 1:' + str(len(Models)) + '){\n')
                #P.write('    plot(y=SS_Tol[,(i * ' + str(len(Models)) + ' * 2 + (-7 + p*2))], x=SS_Tol[,(i * ' + str(len(Models)) + ' * 2 + (-6 + p*2))], main=ModelsVector[' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p], ylab=col_names_matrixSS[i], xlab="Theta", ylim=c(MinMax[i*2-1], MinMax[i*2]))\n')
                #P.write('    abline(h = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
                #P.write('  }\n')
                #P.write('}\n')
            #else:
                #P.write('par(mfrow=c(2,3))\n')
                #P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
                #P.write('  for(p in 1:' + str(len(Models)) + '){\n')
                #P.write('    plot(y=SS_Tol[,(i * ' + str(len(Models)) + ' * 2 + (-7 + p*2))], x=SS_Tol[,(i * ' + str(len(Models)) + ' * 2 + (-6 + p*2))], main=ModelsVector[' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p], ylab=col_names_matrixSS[i], xlab="Theta", ylim=c(MinMax[i*2-1], MinMax[i*2]))\n')
                #P.write('    abline(h = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
                #P.write('  }\n')
                #P.write('}\n')
        #else:
            #P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
            #P.write('  for(p in 1:' + str(len(Models)) + '){\n')
            #P.write('    plot(y=SS_Tol[,(i * ' + str(len(Models)) + ' * 2 + (-7 + p*2))], x=SS_Tol[,(i * ' + str(len(Models)) + ' * 2 + (-6 + p*2))], main=ModelsVector[' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * p], ylab=col_names_matrixSS[i], xlab="Theta", ylim=c(MinMax[i*2-1], MinMax[i*2]))\n')
            #P.write('    abline(h = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
            #P.write('  }\n')
            #P.write('}\n')
        #P.write('invisible(dev.off())\n')
        P.close()
  
#elapsed_time = time.time() - start_time
#with open('ProteinModelerABC.out', 'a') as Error:
        #Error.write('ProteinModelerABC finished after: ' + str(elapsed_time) + ' seconds')
        #Error.close()
#print('ProteinModelerABC finished after: ' + str(elapsed_time) + ' seconds')      


