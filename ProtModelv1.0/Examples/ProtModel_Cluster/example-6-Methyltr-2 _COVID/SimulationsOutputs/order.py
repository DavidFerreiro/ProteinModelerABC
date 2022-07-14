import os
import Variables

os.system("rm -r newalignment.fasta")
os.system("rm -r ./Results")
os.system("rm -r ./__pycache__")
os.system("rm -r Local_interactions.dat")
os.system("rm -r REM.txt")
os.system("rm -r E_loc_0.txt")
os.system("rm -r seqGMRCA")
os.system("rm -r ProtModel_arguments.txt")
os.system("rm -r ProtModel_S_arguments0.txt")
os.system("rm -r ProtModel_S_arguments1.txt")
if os.path.exists('ProtModel_E-Repeat_arguments.txt') == True:
   os.system("rm -r ProtModel_E-Repeat_arguments.txt")
if os.path.exists('ProtModel_F-Repeat_arguments.txt') == True:
   os.system("rm -r ProtModel_F-Repeat_arguments.txt")
if os.path.exists('ProtModel_N-Repeat_arguments.txt') == True:
   os.system("rm -r ProtModel_N-Repeat_arguments.txt")
if os.path.exists('Rest_E_Simulations.txt') == True:
   os.system("rm -r Rest_E_Simulations.txt")
if os.path.exists('Rest_F_Simulations.txt') == True:
   os.system("rm -r Rest_F_Simulations.txt")
if os.path.exists('Rest_N_Simulations.txt') == True:
   os.system("rm -r Rest_N_Simulations.txt")
os.system("rm -r Pop_evol.in")
os.system("mkdir ABCOutputs")
os.system("mkdir SimulationsOutputs")
os.system("mv ABCAnalysis.r ABCOutputs/")
os.system("cp PSimulations.txt ABCOutputs/")
os.system("cp SSSimulations.csv ABCOutputs/")
os.system("cp SSRealData.csv ABCOutputs/")
os.system("mv PSimulations.txt SimulationsOutputs/")
os.system("mv SSSimulations.csv SimulationsOutputs/")
os.system("mv SSRealData.csv SimulationsOutputs/")

os.system("mv launch_* SimulationsOutputs/")
os.system("mv order.py SimulationsOutputs/")
