import os
import time
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
os.system("Rscript --vanilla ./ABCOutputs/ABCAnalysis.r")
attempt = 1
for x in range(0, 11):
    if attempt < 10:
        with open("./ABCOutputs/Results_text.txt") as f:
             if not "Mean model posterior probabilities" in f.read():
                attempt = attempt + 1
                print("ABC estimation failed")
                print("Starting another try: " + str(x + 1))
                os.system("Rscript --vanilla ./ABCOutputs/ABCAnalysis.r")
             else:
                if x == 0:
                    print("ABC analysis  has finished\n")
                    break
                else:
                    print("Try " + str(x + 1) + "  succeeded")
                    print("ABC analysis  has finished\n")
                    break
    else:
        sys.exit("\n"
                " ______                     _ \n"
                "|  ____|                   | |\n"
                "| |__   _ __ _ __ ___  _ __| |\n"
                "|  __| | |__| |__/ _ \| |__| |\n"
                "| |____| |  | | | (_) | |  |_|\n"
                "|______|_|  |_|  \___/|_|  (_)\n\n"
                "But dont panic!! You dont have to repeat everything, just go to the ./ABCOutputs and re execute ABCAnalysis.r file increasing the tolerance or changing the ABC method to rejection. Check manual for further information")
print("PROTMODEL has finished!!\n")
