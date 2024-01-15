from mpi4py import MPI
import bin.Scripts.Variables
import bin.Scripts.ReadSettings
import os
import sys
import csv
import pandas as pd
from Bio import SeqIO
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

#Calculate Summary Statistics of real data
def SScalc(zoe, dat):

    # Call variables
    bin.Scripts.ReadSettings.leer()

    if dat == 'real':
        SFormat = zoe[-4:]
        Sequence = zoe[:-4]
    else:
        Sequence = bin.Scripts.Variables.Simulations_Name[zoe]
        SFormat = ''


    with open(str(Sequence) + SFormat, "r") as f:
        # read an store all lines into list
        lines = f.readlines()

    # Create a list with sequences
    # global Sequences
    Sequences = []
    with open(str(Sequence) + SFormat, "r") as fp:
        for linea in fp:
            if len(linea) < 10:
                pass
            else:
                Sequences.append(linea[10:-1])

    # Create a string for every sequence
    # for x in range(0, bin.Scripts.Variables.NumSeq):
    # globals()[f"Seq{x}"] = Sequences[x]

    # Create two lists, the first one without the last sequences and the second one without the first sequence for comparations
    List1 = []
    for x in range(0, len(Sequences) - 1):
        List1.append(Sequences[x])

    List2 = []
    for x in range(1, len(Sequences)):
        List2.append(Sequences[x])

    # Create two names lists, the first one without the last sequences and the second one without the first sequence for comparations
    Names1 = []
    for y in range(0, bin.Scripts.Variables.NumSeq - 1):
        Names1.append("Seq" + str(y))

    Names2 = []
    for y in range(1, bin.Scripts.Variables.NumSeq):
        Names2.append("Seq" + str(y))

    # Check if DeltaGREM will be necessary
    try:
        if dat == 'real':
            with open("Input_DeltaGREM.in", "w") as SequencesStability:
                SequencesStability.write('# A) Input files\n')
                SequencesStability.write('ALI=' + str(Sequence) + '.fasta\n')
                SequencesStability.write('CHAIN=  ' + str(bin.Scripts.Variables.ChainPDB) + '\n')
                SequencesStability.write('PDB=' + bin.Scripts.Variables.Template + '\n\n\n')
                SequencesStability.write('FILE_STR=structures.in\n')
                SequencesStability.write('#MUT=	prot.mut	# list of mutations (optional)\n\n')
                SequencesStability.write('# FASTA file with related proteins (optional)\n')
                SequencesStability.write('#================================================================\n')
                SequencesStability.write('# B) Thermodynamic model\n')
                SequencesStability.write('TEMP=	0.5		# Temperature\n')
                SequencesStability.write('SU1=	0.065		# configurational entropy per res (unfold)\n')
                SequencesStability.write('SC1=  0.065		# configurational entropy per res (misfold)\n')
                SequencesStability.write('SC0=  0.0		# configurational entropy offset (misfold)\n')
                SequencesStability.write('REM=   2		# Use up to 1,2,3r moments of misfolding energy?\n')
        else:
            with open("Input_DeltaGREM" + str(zoe + 1) + ".in", "w") as SequencesStability:
                SequencesStability.write('# A) Input files\n')
                SequencesStability.write('ALI=' + str(Sequence) + '.fasta\n')
                SequencesStability.write('CHAIN=  ' + str(bin.Scripts.Variables.ChainPDB) + '\n')
                SequencesStability.write('PDB=' + bin.Scripts.Variables.Template + '\n\n\n')
                SequencesStability.write('FILE_STR=structures.in\n')
                SequencesStability.write('#MUT=	prot.mut	# list of mutations (optional)\n\n')
                SequencesStability.write('# FASTA file with related proteins (optional)\n')
                SequencesStability.write('#================================================================\n')
                SequencesStability.write('# B) Thermodynamic model\n')
                SequencesStability.write('TEMP=	0.5		# Temperature\n')
                SequencesStability.write('SU1=	0.065		# configurational entropy per res (unfold)\n')
                SequencesStability.write('SC1=  0.065		# configurational entropy per res (misfold)\n')
                SequencesStability.write('SC0=  0.0		# configurational entropy offset (misfold)\n')
                SequencesStability.write('REM=   2		# Use up to 1,2,3r moments of misfolding energy?\n')

    except OSError:
        print('DeltaGREM Failed. Please check DeltaGREM inputs.')
        sys.exit('ERROR!! Please check  DeltaGREM inputs')

    # Read pdb fasta sequence
    warnings.simplefilter('ignore', PDBConstructionWarning)
    PDBFile = './' + bin.Scripts.Variables.Template
    with open(PDBFile, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            if record.annotations["chain"] == str(bin.Scripts.Variables.ChainPDB):
                PDB_Name = ('>' + record.id[:-2])
                PDB_Seq = (record.seq)
                PDB_Seq = str(PDB_Seq)
                PDB_Seq = PDB_Seq.replace('X', '')

    with open("./" + Sequence + ".fasta", "w") as fp:
        # fp.write(PDB_Name)
        # fp.write(str(PDB_Seq))
        for number, line in enumerate(lines):
            if len(line) > 10:
                fp.write(">")
                fp.write(line[0:9])
                fp.write("\n")
                for e in range(10, len(line)):
                    if line[e] == '-':
                        fp.write(str(PDB_Seq)[(e - 10)])
                    else:
                        fp.write(line[e])
        fp.write(PDB_Name + '\n')
        fp.write(str(PDB_Seq) + '\n')

    if dat == 'real':
        DGREM_route = INSTALLATIONROOT + '/bin/DeltaGREM/DeltaGREM '
        os.system(DGREM_route + 'Input_DeltaGREM.in')
        #os.system("./bin/DeltaGREM/DeltaGREM Input_DeltaGREM.in")  # execute program
    else:
        DGREM_route = INSTALLATIONROOT + '/bin/DeltaGREM/DeltaGREM '
        os.system(DGREM_route + 'Input_DeltaGREM'  + str(zoe + 1) + '.in')
        #os.system("./bin/DeltaGREM/DeltaGREM Input_DeltaGREM" + str(zoe + 1) + ".in")  # execute program

    DGREM_RS = []

    try:
        with open(Sequence + ".fasta_DeltaG.dat") as res:
            for line in res:
                line = line.split('	')
                if line[0][0] != "#":  # + str(bin.Scripts.Variables.NumAA):
                    if line[9] != bin.Scripts.Variables.Template[:-4] + "\n":
                        SEQ = line[1]
                        DGREM_RS.append(SEQ)

    except OSError:
        print("DeltaGREM results file can not be readed")
        sys.exit('ERROR!! Please check DeltaGREM inputs')

        # DGREM mean
    DGREM_Mean = 0
    DGREM_Mean_Sum = 0
    for x in range(0, bin.Scripts.Variables.NumSeq):
        DGREM_Mean_Sum = DGREM_Mean_Sum + float(DGREM_RS[x])
    DGREM_Mean = (DGREM_Mean_Sum / bin.Scripts.Variables.NumSeq)
    if bin.Scripts.Variables.ShowInformationScreen == 'Yes':
        if dat == 'real':
            print('  DGREM mean calculated')
        else:
            if zoe == 0 or (int(zoe) + 1) % 50 == 0:
                print('  DGREM mean calculated')

    # DGREM sd
    DGREM_sd = 0
    DGREM_sd_Sum = 0
    for x in range(0, bin.Scripts.Variables.NumSeq):
        DGREM_sd_Sum = DGREM_sd_Sum + ((float(DGREM_RS[x]) - DGREM_Mean) * (float(DGREM_RS[x]) - DGREM_Mean))
    DGREM_sd = (DGREM_sd_Sum / (bin.Scripts.Variables.NumSeq - 1)) ** (0.5)
    if bin.Scripts.Variables.ShowInformationScreen == 'Yes':
        if dat == 'real':
            print('  DGREM sd calculated')
        else:
            if zoe == 0 or (int(zoe) + 1) % 50 == 0:
                print('  DGREM sd calculated')

    # Segregating Sites
    if bin.Scripts.Variables.Indels == 'Ignored':
        SegSites = 0
        for x in range(0, bin.Scripts.Variables.NumAA):
            if List1[0][x] == "-":
                pass
            for y in List2:
                if y[x] == "-":
                    pass
                elif List1[0][x] != y[x]:
                    SegSites = SegSites + 1
                    break
        if bin.Scripts.Variables.ShowInformationScreen == 'Yes':
            if dat == 'real':
                print('  Segregating sites calculated')
            else:
                if zoe == 0 or (int(zoe) + 1) % 50 == 0:
                    print('  Segregating sites calculated')

    if bin.Scripts.Variables.Indels == 'New State':
        SegSites = 0
        for x in range(0, bin.Scripts.Variables.NumAA):
            for y in List2:
                if List1[0][x] != y[x]:
                    SegSites = SegSites + 1
                    break
        if bin.Scripts.Variables.ShowInformationScreen == 'Yes':
            if dat == 'real':
                print('  Segregating sites calculated')
            else:
                if zoe == 0 or (int(zoe) + 1) % 50 == 0:
                    print('  Segregating sites calculated')

    # Grantham Distance per site
    matriz = pd.read_csv(INSTALLATIONROOT + '/bin/Grantham.csv')
    Total_Grantham_Position = []
    Grantham = 0
    for e in range(0, bin.Scripts.Variables.NumAA):

        List2 = []
        for b in range(1, len(Sequences)):
            List2.append(Sequences[b])

        for x in List1:
            for y in List2:
                if x[e] != y[e]:
                    if x[e] == '-' or y[e] == '-':
                        pass
                    else:
                        row = matriz.columns.get_loc(y[e])
                        Grantham = Grantham + matriz._get_value(row, x[e]) ** 2

            List2.remove(List2[0])
        Total_Grantham_Position.append(Grantham)
        Grantham = 0

    for b in range(1, len(Sequences)):
        List2.append(Sequences[b])

    # Common variables for sd, sk and ku
    Grantham_Position_Mean = 0
    for x in range(0, bin.Scripts.Variables.NumAA):
        Grantham_Position_Mean = Grantham_Position_Mean + Total_Grantham_Position[x]
    Grantham_Position_Mean = Grantham_Position_Mean / bin.Scripts.Variables.NumAA
    if bin.Scripts.Variables.ShowInformationScreen == 'Yes':
        if dat == 'real':
            print('  Grantham distance mean calculated')
        else:
            if zoe == 0 or (int(zoe) + 1) % 50 == 0:
                print('  Grantham distance mean calculated')

    # sd
    Grantham_Position_sd_Sum = 0
    for x in range(0, bin.Scripts.Variables.NumAA):
        Grantham_Position_sd_Sum = Grantham_Position_sd_Sum + ((Total_Grantham_Position[x] - Grantham_Position_Mean) * (
                    Total_Grantham_Position[x] - Grantham_Position_Mean))
    Grantham_Position_sd = (Grantham_Position_sd_Sum / bin.Scripts.Variables.NumAA) ** (0.5)
    if bin.Scripts.Variables.ShowInformationScreen == 'Yes':
        if dat == 'real':
            print('  Grantham distance sd calculated')
        else:
            if zoe == 0 or (int(zoe) + 1) % 50 == 0:
                print('  Grantham distance sd calculated')

    # sk
    if Grantham_Position_sd == 0.0:
        Grantham_Position_sk = 0
        if bin.Scripts.Variables.ShowInformationScreen == 'Yes':
            if dat == 'real':
                print('  Grantham distance skewness calculated')
            else:
                if zoe == 0 or (int(zoe) + 1) % 50 == 0:
                    print('  Grantham distance skewness calculated')
    else:
        Grantham_Position_sk_Sum = 0
        for x in range(0, bin.Scripts.Variables.NumAA):
            Grantham_Position_sk_Sum = Grantham_Position_sk_Sum + ((Total_Grantham_Position[x] - Grantham_Position_Mean) * (
                        Total_Grantham_Position[x] - Grantham_Position_Mean) * (Total_Grantham_Position[
                                                                                  x] - Grantham_Position_Mean))
        Grantham_Position_sk = (Grantham_Position_sk_Sum / bin.Scripts.Variables.NumAA) / (
                    Grantham_Position_sd * Grantham_Position_sd * Grantham_Position_sd)
        if bin.Scripts.Variables.ShowInformationScreen == 'Yes':
            if dat == 'real':
                print('  Grantham distance skewness calculated')
            else:
                print('  Grantham distance skewness calculated')

    # ku
    if Grantham_Position_sd == 0.0:
        Grantham_Position_ku = 0
        if bin.Scripts.Variables.ShowInformationScreen == 'Yes':
            if dat == 'real':
                print('  Grantham distance kurtosis calculated')
            else:
                if zoe == 0 or (int(zoe) + 1) % 50 == 0:
                    print('  Grantham distance kurtosis calculated')
    else:
        Grantham_Position_ku_Sum = 0
        for x in range(0, bin.Scripts.Variables.NumAA):
            Grantham_Position_ku_Sum = Grantham_Position_ku_Sum + ((Total_Grantham_Position[x] - Grantham_Position_Mean) * (
                        Total_Grantham_Position[x] - Grantham_Position_Mean) * (Total_Grantham_Position[
                                                                                  x] - Grantham_Position_Mean) * (
                                                                             Total_Grantham_Position[
                                                                                 x] - Grantham_Position_Mean))
        Grantham_Position_ku = (Grantham_Position_ku_Sum / bin.Scripts.Variables.NumAA) / (
                    Grantham_Position_sd * Grantham_Position_sd * Grantham_Position_sd * Grantham_Position_sd)
        if bin.Scripts.Variables.ShowInformationScreen == 'Yes':
            if dat == 'real':
                print('  Grantham distance kurtosis calculated')
            else:
                if zoe == 0 or (int(zoe) + 1) % 50 == 0:
                    print('  Grantham distance kurtosis calculated')

    if dat == 'real':
        SS_HeadersRD = []
        SS_ValRD = []
    else:
        SS_Headers = []
        SS_csv = []

    for x in bin.Scripts.Variables.SummaryStatistics.split(' '):
        if x == '1':
            if dat == 'real':
                SS_HeadersRD.append('DGREM_Mean')
                SS_ValRD.append(str(DGREM_Mean))
            else:
                if zoe == 0:
                    SS_Headers.append('DGREM_Mean')
                SS_csv.append(str(DGREM_Mean))
        elif x == '2':
            if dat == 'real':
                SS_HeadersRD.append('DGREM_sd')
                SS_ValRD.append(str(DGREM_sd))
            else:
                if zoe == 0:
                    SS_Headers.append('DGREM_sd')
                SS_csv.append(str(DGREM_sd))
        elif x == '3':
            if dat == 'real':
                SS_HeadersRD.append('SegSites')
                SS_ValRD.append(str(SegSites))
            else:
                if zoe == 0:
                    SS_Headers.append('SegSites')
                SS_csv.append(str(SegSites))
        elif x == '4':
            if dat == 'real':
                SS_HeadersRD.append('Grantham_mean_Position')
                SS_ValRD.append(str(Grantham_Position_Mean))
            else:
                if zoe == 0:
                    SS_Headers.append('Grantham_mean_Position')
                SS_csv.append(str(Grantham_Position_Mean))
        elif x == '5':
            if dat == 'real':
                SS_HeadersRD.append('Grantham_sd_Position')
                SS_ValRD.append(str(Grantham_Position_sd))
            else:
                if zoe == 0:
                    SS_Headers.append('Grantham_sd_Position')
                SS_csv.append(str(Grantham_Position_sd))
        elif x == '6':
            if dat == 'real':
                SS_HeadersRD.append('Grantham_sk_Position')
                SS_ValRD.append(str(Grantham_Position_sk))
            else:
                if zoe == 0:
                    SS_Headers.append('Grantham_sk_Position')
                SS_csv.append(str(Grantham_Position_sk))
        elif x == '7':
            if dat == 'real':
                SS_HeadersRD.append('Grantham_ku_Position')
                SS_ValRD.append(str(Grantham_Position_ku))
            else:
                if zoe == 0:
                    SS_Headers.append('Grantham_ku_Position')
                SS_csv.append(str(Grantham_Position_ku))

    if dat == 'real':
        Real = open('SSRealData.csv', 'w')
        writer = csv.writer(Real)
        writer.writerow(SS_HeadersRD)
        writer.writerow(SS_ValRD)

        if bin.Scripts.Variables.SaveSimulations == 'Yes':
            os.system('mkdir Simulations')
            os.system('mv Input_DeltaGREM.in ./Simulations')
            os.system('mv ' + str(Sequence) + '.fasta ./Simulations')
            os.system('mv ' + str(Sequence) + '.fasta_DeltaG.dat ./Simulations')
        else:
            os.system('rm -r Input_DeltaGREM.in')
            os.system('rm -r ' + str(Sequence) + '.fasta')
            os.system('rm -r ' + str(Sequence) + '.fasta_DeltaG.dat')
    else:
        if zoe == 0:
            Simu = open('SSSimulations.csv', 'w')
            writer = csv.writer(Simu)
            writer.writerow(SS_Headers)
            #writer.writerow(SS_csv)

        else:
            pass
            #Simu = open('SSSimulations.csv', 'a')
            #writer = csv.writer(Simu)
            #writer.writerow(SS_csv)

        if bin.Scripts.Variables.SaveSimulations == 'Yes':
            os.system('mv ' + str(Sequence) + ' ./Simulations')
            os.system('mv Input_DeltaGREM' + str(zoe + 1) + '.in ./Simulations')
            os.system('mv ' + str(Sequence) + '.fasta ./Simulations')
            os.system('mv ' + str(Sequence) + '.fasta_DeltaG.dat ./Simulations')
        else:
            os.system('rm -r ' + str(Sequence))
            os.system('rm -r Input_DeltaGREM' + str(zoe + 1) + '.in')
            os.system('rm -r ' + str(Sequence) + '.fasta')
            os.system('rm -r ' + str(Sequence) + '.fasta_DeltaG.dat')

        return SS_csv

#Execute simulations
def Cluster_Simu(pv, zoe, model):

    # Open temporary file with ProteinEvolver commands
    with open('ProteinEvolver_Arguments.txt', 'r') as file:
        lines = file.readlines()

    # Initialize MPI communicator
    comm = MPI.COMM_WORLD
    size = comm.Get_size() # new: gives number of ranks in comm
    rank = comm.Get_rank()

    # Calculate the number of simulations per process
    for x in range(0, size):
        val = zoe - pv
        val = val - x
        remainder =int(val) % int(size)
        is_divisible = remainder ==  0
        if is_divisible == True:
            num = int(val) / int(size)
            break

    # Calculate the number of remaining simulations
    rest_S = (zoe - pv) - val

    # Execute the simulations for empirical models
    if model != 'Fitness' and model != 'Neutral':
        for y in range(0, int(num)):
            if rank == 0:
                data = [size * y + z for z in range(size)]
            else:
                data = None

            data = comm.scatter(data, root=0)
            data = data + pv
            
            PE_route = INSTALLATIONROOT + '/bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 '
            os.system(PE_route + lines[data])
            #os.system('./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 ' + lines[data])

        # Execute the remaining simulations
        if rank < rest_S:
            data_r = [size * rank + z for z in range(size)]

            data_r = comm.scatter(data_r, root=0)
            data_r = data_r + pv + val
            
            PE_route = INSTALLATIONROOT + '/bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 '
            os.system(PE_route + lines[data_r])
            #os.system('./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 ' + lines[data_r])

    else:
        if os.path.isfile('Pop_evol.in') == True:
            with open('Pop_evol.in', 'r') as of:
                lineas = of.readlines()
                for l in range(0, len(lineas)):
                    li = lineas[l]
                    linea = lineas[l].split('=')
                    if linea[0] == 'NEUTRAL':
                        if model == 'Neutral':
                            if linea[1].split()[0] == '1':
                                pass
                            else:
                                li = li.replace('NEUTRAL=	 0', 'NEUTRAL=	 1')
                                with open('Pop_evol.in', 'w') as fp:
                                    for u in range(0, len(lineas)):
                                        if u == l:
                                            fp.write(li)
                                        else:
                                            fp.write(lineas[u])
                        elif model == 'Fitness':
                            if linea[1].split()[0] == '0':
                                pass
                            else:
                                li = li.replace('NEUTRAL=	 1', 'NEUTRAL=	 0')
                                with open('Pop_evol.in', 'w') as fp:
                                    for u in range(0, len(lineas)):
                                        if u == l:
                                            fp.write(li)
                                        else:
                                            fp.write(lineas[u])
                of.close()

        for y in range(0, int(num)):
            if rank == 0:
                data = [size * y + z for z in range(size)]
            else:
                data = None
            data = comm.scatter(data, root=0)
            data = data + pv
            
            PE_route = INSTALLATIONROOT + '/bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 '
            os.system(PE_route + lines[data])
            #os.system('./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 ' + lines[data])


        # Execute the remaining simulations
        if rank < rest_S:
            data_r = [size * rank + z for z in range(size)]
            data_r = comm.scatter(data_r, root=0)
            data_r = data_r + pv + val
            
            PE_route = INSTALLATIONROOT + '/bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 '
            os.system(PE_route + lines[data_r])
            #os.system('./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 ' + lines[data_r])

#Check if simulations finished succesfully
def CheckSimu(pv, zoe, model):

    # Call variables
    bin.Scripts.ReadSettings.leer()

    # Open temporary file with ProteinEvolver commands
    with open('ProteinEvolver_Arguments.txt', 'r') as file:
        lines = file.readlines()

    # Initialize MPI communicator
    comm = MPI.COMM_WORLD
    size = comm.Get_size() # new: gives number of ranks in comm
    rank = comm.Get_rank()

    # Calculate the number of simulations per process
    for x in range(0, size):
        val = zoe - pv
        val = val - x
        remainder =int(val) % int(size)
        is_divisible = remainder ==  0
        if is_divisible == True:
            num = int(val) / int(size)
            break

    # Calculate the number of remaining simulations
    rest_S = (zoe - pv) - val

    # Execute the simulations for empirical models
    if model != 'Fitness' and model != 'Neutral':
        for y in range(0, int(num)):
            if rank == 0:
                data = [size * y + z for z in range(size)]
            else:
                data = None

            data = comm.scatter(data, root=0)
            data = data + pv

            num = ''
            file = bin.Scripts.Variables.Simulations_Name[data]
            if os.path.isfile(file) == True:
                if os.path.getsize(file) == 0:
                    num = int(file[9:])

            else:
                num = int(file[9:])

            if num != '':
                print('  Repeting simulation ' + str(num) + ', ' + str(model))
                PE_route = INSTALLATIONROOT + '/bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 '
                os.system(PE_route + lines[data])
                #os.system('./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 ' + lines[data])

        if rank < rest_S:
            data_r = [size * rank + z for z in range(size)]

            data_r = comm.scatter(data_r, root=0)
            data_r = data_r + pv + val

            num = ''
            file = bin.Scripts.Variables.Simulations_Name[data_r]
            if os.path.isfile(file) == True:
                if os.path.getsize(file) == 0:
                    num = int(file[9:])

            else:
                num = int(file[9:])

            if num != '':
                print('  Repeting simulation ' + str(num) + ', ' + str(model))
                PE_route = INSTALLATIONROOT + '/bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 '
                os.system(PE_route + lines[data_r])
                #os.system('./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 ' + lines[data_r])

    else:
        with open('Pop_evol.in', 'r') as of:
            lineas = of.readlines()
            for l in range(0, len(lineas)):
                li = lineas[l]
                linea = lineas[l].split('=')
                if linea[0] == 'NEUTRAL':
                    if model == 'Neutral':
                        if linea[1].split()[0] == '1':
                            pass
                        else:
                            li = li.replace('NEUTRAL=	 0', 'NEUTRAL=	 1')
                            with open('Pop_evol.in', 'w') as fp:
                                for u in range(0, len(lineas)):
                                    if u == l:
                                        fp.write(li)
                                    else:
                                        fp.write(lineas[u])
                    elif model == 'Fitness':
                        if linea[1].split()[0] == '0':
                            pass
                        else:
                            li = li.replace('NEUTRAL=	 1', 'NEUTRAL=	 0')
                            with open('Pop_evol.in', 'w') as fp:
                                for u in range(0, len(lineas)):
                                    if u == l:
                                        fp.write(li)
                                    else:
                                        fp.write(lineas[u])

            of.close()

        for y in range(0, int(num)):
            if rank == 0:
                data = [size * y + z for z in range(size)]
            else:
                data = None

            data = comm.scatter(data, root=0)
            data = data + pv

            num = ''
            file = bin.Scripts.Variables.Simulations_Name[data]
            if os.path.isfile(file) == True:
                if os.path.getsize(file) == 0:
                    num = int(file[9:])

            else:
                num = int(file[9:])

            if num != '':
                print('  Repeting simulation ' + str(num) + ', ' + str(model))
                PE_route = INSTALLATIONROOT + '/bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 '
                os.system(PE_route + lines[data])
                #os.system('./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 ' + lines[data])

        if rank < rest_S:
            data_r = [size * rank + z for z in range(size)]
            data_r = comm.scatter(data_r, root=0)
            data_r = data_r + pv + val

            num = ''
            file = bin.Scripts.Variables.Simulations_Name[data_r]
            if os.path.isfile(file) == True:
                if os.path.getsize(file) == 0:
                    num = int(file[9:])

            else:
                num = int(file[9:])

            if num != '':
                print('  Repeting simulation ' + str(num) + ', ' + str(model))
                PE_route = INSTALLATIONROOT + '/bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 '
                os.system(PE_route + lines[data_r])
                #os.system('./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 ' + lines[data_r])

#Calculate Summary Statistics of simulations (parallel)
def Cluster_SS(pv, zoe):

    results = []
    results_r = []

    # Initialize MPI communicator
    comm = MPI.COMM_WORLD
    size = comm.Get_size()  # new: gives number of ranks in comm
    rank = comm.Get_rank()

    # Calculate the number of SScalc per process
    for x in range(0, size):
        val = zoe - pv
        val = val - x
        remainder = int(val) % int(size)
        is_divisible = remainder == 0
        if is_divisible == True:
            num = int(val) / int(size)
            break

    # Calculate the number of remaining SScalc
    rest_SS = (zoe - pv) - val

    for y in range(0, int(num)):
        if rank == 0:
            data = [size * y + z for z in range(size)]
        else:
            data = None
        data = comm.scatter(data, root=0)
        data = data + pv
        res = SScalc(int(data), 'simu')
        results.append(res)

    all_results = comm.gather(results, root=0)
    if rank == 0:
        with open('SSSimulations.csv', 'a') as Simu:
            writer = csv.writer(Simu)
            for x in all_results:
                for y in x:
                    writer.writerow(y)



    # Execute the remaining simulations
    if rank < rest_SS:
        data_r = [size * rank + z for z in range(size)]
        data_r = comm.scatter(data_r, root=0)
        data_r = data_r + pv + val
        res_r =SScalc(int(data_r), 'simu')
        results_r.append(res_r)

    all_results_r = comm.gather(results_r, root=0)
    if rank == 0:
        with open('SSSimulations.csv', 'a') as Simu:
            writer = csv.writer(Simu)
            for x in all_results_r:
                for y in x:
                    writer.writerow(y)

# Order files
def Order():

    if bin.Scripts.Variables.SaveSimulations == 'Yes':
        os.system('tar -czf Simulations.tar.gz ./Simulations/ ')
        os.system('rm -r ./Simulations')

    os.system('rm -r ./Results')
    os.system('rm newalignment.fasta')
    os.system('rm Local_interactions.dat')
    os.system('rm REM.txt')
    os.system('rm E_loc_0.txt')
    os.system('rm seqGMRCA')
    os.system('rm Pop_evol.in')
    os.system('rm ProteinEvolver_Arguments.txt')


    if os.path.exists('ABCOutputs') == True:
        os.system('rm -r ABCOutputs')
    if os.path.exists('SimulationsOutputs') == True:
        os.system('rm -r SimulationsOutputs')

    os.system('mkdir ABCOutputs')
    os.system('mkdir SimulationsOutputs')


    os.system('mv ABCAnalysis.r ABCOutputs/')
    if bin.Scripts.Variables.Coalescent == 'Coal':
        os.system('mv Histogram_Priors.pdf ABCOutputs/')

    os.system('cp SSSimulations.csv ABCOutputs/')
    os.system('cp SSRealData.csv ABCOutputs/')

    os.system("mv SSSimulations.csv SimulationsOutputs/")
    os.system("mv SSRealData.csv SimulationsOutputs/")
    os.system("mv launch_PMABC.sh SimulationsOutputs/")

    if os.path.exists('Simulations.tar.gz') == True:
        os.system('mv Simulations.tar.gz SimulationsOutputs/')
    if os.path.exists('PSimulations.txt') == True:
        os.system('mv PSimulations.txt SimulationsOutputs/')

    #Execute ABC analysis
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
                        print("ABC analysis has finished\n")
                        break
                    else:
                        print("Try " + str(x + 1) + " succeeded")
                        print("ABC analysis has finished\n")
                        break
        else:
            sys.exit("\n"
                     " ______                              _ \n"
                     "|  ____|                            | |\n"
                     "| |__     _ __   _ __   ___    _ __ | |\n"
                     "|  __|   | '__| | '__| | _ |  | '__|| |\n"
                     "| |____  | |    | |    |(_)|  | |   |_|\n"
                     "|______| |_|    |_|    |___|  |_|   (_)\n\n"
                     "But dont panic!! You dont have to repeat everything, just go to the ./ABCOutputs and re execute ABCAnalysis.r file increasing the tolerance or changing the ABC method to rejection. Check manual for further information")

# Create a plot of the Prior distribution of the substitution rate
def plot():

    SubstitutionRate_LI = []
    for x in bin.Scripts.Variables.SubstitutionRate_L:
        SubstitutionRate_LI.append(float(x))

    # Crear un nuevo archivo PDF
    figure_name = "./Histogram_Priors.pdf"
    pdf_pages = PdfPages(figure_name)

    if bin.Scripts.Variables.Coalescent == 'Coal':
        if bin.Scripts.Variables.MultiPage == 'Yes':
            fig, axs = plt.subplots(1, 2, figsize=(10, 5))

            axs[0].hist(SubstitutionRate_LI, bins=40, color='gray', edgecolor='black')
            axs[0].set_title('Substitution Rate Prior')
            axs[0].set_xlabel('Substitution Rate')

            # Interval
            interval = (max(SubstitutionRate_LI) - min(SubstitutionRate_LI)) / 3
            axs[0].set_xticks(np.arange(min(SubstitutionRate_LI), max(SubstitutionRate_LI) + min(SubstitutionRate_LI), interval))
            axs[0].set_xticklabels([f'{x:.2e}' for x in np.arange(min(SubstitutionRate_LI), max(SubstitutionRate_LI) + min(SubstitutionRate_LI), interval)])


            axs[1].hist(bin.Scripts.Variables.ThetaRate_L, bins=40, color='gray', edgecolor='black')
            axs[1].set_title('Theta Prior')
            axs[1].set_xlabel('Theta')

            # Interval
            interval = (max(bin.Scripts.Variables.ThetaRate_L) - min(bin.Scripts.Variables.ThetaRate_L)) / 3
            axs[1].set_xticks(np.arange(min(bin.Scripts.Variables.ThetaRate_L), max(bin.Scripts.Variables.ThetaRate_L) + min(bin.Scripts.Variables.ThetaRate_L), interval))
            axs[1].set_xticklabels([f'{x:.2f}' for x in np.arange(min(bin.Scripts.Variables.ThetaRate_L),max(bin.Scripts.Variables.ThetaRate_L) + min(bin.Scripts.Variables.ThetaRate_L),interval)])

            # Agregar la figura al PDF
            pdf_pages.savefig(fig)
            plt.close(fig)

        else:
            fig, ax = plt.subplots(figsize=(5, 5))

            ax.hist(SubstitutionRate_LI, bins=40, color='gray', edgecolor='black',)
            ax.set_title('Substitution Rate Prior')
            ax.set_xlabel('Substitution Rate')
            # Interval
            interval = (max(SubstitutionRate_LI) - min(SubstitutionRate_LI)) / 5
            ax.set_xticks(np.arange(min(SubstitutionRate_LI), max(SubstitutionRate_LI) + min(SubstitutionRate_LI), interval))
            ax.set_xticklabels([f'{x:.2e}' for x in np.arange(min(SubstitutionRate_LI), max(SubstitutionRate_LI) + min(SubstitutionRate_LI), interval)])


            # Agregar la figura al PDF
            pdf_pages.savefig(fig)
            plt.close(fig)

            ax.hist(bin.Scripts.Variables.ThetaRate_L, bins=40, color='gray', edgecolor='black')
            ax.set_title('Theta Prior')
            ax.set_xlabel('Theta', fontsize=10)
            # Interval
            interval = (max(bin.Scripts.Variables.ThetaRate_L) - min(bin.Scripts.Variables.ThetaRate_L)) / 5
            ax.set_xticks(np.arange(min(bin.Scripts.Variables.ThetaRate_L), max(bin.Scripts.Variables.ThetaRate_L) + min(bin.Scripts.Variables.ThetaRate_L), interval))
            ax.set_xticklabels([f'{x:.2f}' for x in np.arange(min(bin.Scripts.Variables.ThetaRate_L), max(bin.Scripts.Variables.ThetaRate_L) + min(bin.Scripts.Variables.ThetaRate_L), interval)])


            # Agregar la figura al PDF
            pdf_pages.savefig(fig)
            plt.close(fig)

        # Guardar el histograma como un archivo PDF
        pdf_pages.close()
