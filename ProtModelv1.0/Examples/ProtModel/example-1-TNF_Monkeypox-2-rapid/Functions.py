import os
import sys
from Bio import SeqIO
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import csv
import pandas as pd
import Variables
import LeerSettings

def E_Simu(zoe):
    with open('ProtModel_arguments.txt','r') as f:
        lines = f.readlines()
        os.system("./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 " + lines[zoe])
        #if (zoe + 1) == 1:
            #print(' Simulation # ' + str(zoe + 1))
        #elif (zoe + 1) % 100 == 0:
            #print(' Simulation # ' + str(zoe+1))


def S_Simu_1(zoe):
    Variables.init()
    LeerSettings.leer()
    try:
        with open('Pop_evol.in', 'r') as of:
            for lineas in of:
                lineas = lineas.split('=')
                lineas[0] = lineas[0].strip()
                if lineas[0] == 'NEUTRAL':
                    if Variables.StructuralSubstitutionModel == 'Fitness':
                        if lineas[1][1] == '0':
                            pass
                        else:
                            with open('Pop_evol.in', 'w') as fp:
                                fp.write('PDB= ' + Variables.Template + ' 		# pdb file\n')
                                fp.write('CHAIN=  ' + str(Variables.ChainPDB) + '			# Chain of pdb file\n')
                                fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
                                fp.write('TEMP=	1.8			# Temperature\n')
                                fp.write('S0=	0.05			# s0, configurational entropy per residue (unfolded)\n')
                                fp.write('SC1=    0.05			# configurational entropy per residue (misfolded)\n')
                                fp.write('SC0=    0.0			# configurational entropy offset (misfolded)\n')
                                fp.write('REM3=	0			# Third cumulant in REM calculation\n')
                                fp.write(
                                'NEUTRAL=	 0		# If 1, Neutral landscape, otherwise population size dependent selection \n')
                                fp.write('NPOP=	10			# N_pop, population size\n')
                                fp.write(
                                'TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
                                fp.write(
                                'OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"\n')
                                fp.close()
                    elif Variables.StructuralSubstitutionModel == 'Neutral':
                        if lineas[1][1] == '1':
                            pass
                        else:
                            with open('Pop_evol.in', 'w') as fp:
                                fp.write('PDB= ' + Variables.Template + ' 		# pdb file\n')
                                fp.write('CHAIN=  ' + str(Variables.ChainPDB) + '			# Chain of pdb file\n')
                                fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
                                fp.write('TEMP=	1.8			# Temperature\n')
                                fp.write('S0=	0.05			# s0, configurational entropy per residue (unfolded)\n')
                                fp.write('SC1=    0.05			# configurational entropy per residue (misfolded)\n')
                                fp.write('SC0=    0.0			# configurational entropy offset (misfolded)\n')
                                fp.write('REM3=	0			# Third cumulant in REM calculation\n')
                                fp.write(
                                'NEUTRAL=	 1		# If 1, Neutral landscape, otherwise population size dependent selection \n')
                                fp.write('NPOP=	10			# N_pop, population size\n')
                                fp.write(
                                'TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
                                fp.write(
                                'OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"\n')
                                fp.close()
            of.close()
    except:
        if Variables.StructuralSubstitutionModel == 'Fitness':
            with open('Pop_evol.in', 'w') as fp:
                fp.write('PDB= ' + Variables.Template + ' 		# pdb file\n')
                fp.write('CHAIN=  ' + str(Variables.ChainPDB) + '			# Chain of pdb file\n')
                fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
                fp.write('TEMP=	1.8			# Temperature\n')
                fp.write('S0=	0.05			# s0, configurational entropy per residue (unfolded)\n')
                fp.write('SC1=    0.05			# configurational entropy per residue (misfolded)\n')
                fp.write('SC0=    0.0			# configurational entropy offset (misfolded)\n')
                fp.write('REM3=	0			# Third cumulant in REM calculation\n')
                fp.write(
                'NEUTRAL=	 0		# If 1, Neutral landscape, otherwise population size dependent selection \n')
                fp.write('NPOP=	10			# N_pop, population size\n')
                fp.write(
                'TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
                fp.write(
                'OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"\n')
                fp.close()
        elif Variables.StructuralSubstitutionModel == 'Neutral':
            with open('Pop_evol.in', 'w') as fp:
                fp.write('PDB= ' + Variables.Template + ' 		# pdb file\n')
                fp.write('CHAIN=  ' + str(Variables.ChainPDB) + '			# Chain of pdb file\n')
                fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
                fp.write('TEMP=	1.8			# Temperature\n')
                fp.write('S0=	0.05			# s0, configurational entropy per residue (unfolded)\n')
                fp.write('SC1=    0.05			# configurational entropy per residue (misfolded)\n')
                fp.write('SC0=    0.0			# configurational entropy offset (misfolded)\n')
                fp.write('REM3=	0			# Third cumulant in REM calculation\n')
                fp.write(
                'NEUTRAL=	 1		# If 1, Neutral landscape, otherwise population size dependent selection \n')
                fp.write('NPOP=	10			# N_pop, population size\n')
                fp.write(
                'TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
                fp.write(
                'OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"\n')
                fp.close()


    with open('ProtModel_S_arguments.txt', 'r') as f:
        data = f.readlines()
        os.system("./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 " + data[zoe])
        # print('Simulation # ' + str(data[zoe].split(':')[1][:-1]))
        #if (zoe + 1) == 1:
            #print(' Simulation # ' + str(zoe + 1))
        #elif (zoe + 1) % 100 == 0:
            #print(' Simulation # ' + str(zoe + 1))
        f.close()


def S_Simu_F(zoe):
    Variables.init()
    LeerSettings.leer()
    try:
        with open('Pop_evol.in', 'r') as of:
            for lineas in of:
                lineas = lineas.split('=')
                lineas[0] = lineas[0].strip()
                if lineas[0] == 'NEUTRAL':
                    if lineas[1][1] == '0':
                        pass
                    else:
                        with open('Pop_evol.in', 'w') as fp:
                            fp.write('PDB= ' + Variables.Template + '		# pdb file\n')
                            fp.write('CHAIN=  ' + str(Variables.ChainPDB) + '			# Chain of pdb file\n')
                            fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
                            fp.write('TEMP=	1.8			# Temperature\n')
                            fp.write('S0=	0.05			# s0, configurational entropy per residue (unfolded)\n')
                            fp.write('SC1=    0.05			# configurational entropy per residue (misfolded)\n')
                            fp.write('SC0=    0.0			# configurational entropy offset (misfolded)\n')
                            fp.write('REM3=	0			# Third cumulant in REM calculation\n')
                            fp.write(
                                'NEUTRAL=	0		# If 1, Neutral landscape, otherwise population size dependent selection \n')
                            fp.write('NPOP=	10			# N_pop, population size\n')
                            fp.write(
                                'TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
                            fp.write(
                                'OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"')
                            fp.close()
            of.close()
    except:
        with open('Pop_evol.in', 'w') as fp:
            fp.write('PDB= ' + Variables.Template + '		# pdb file\n')
            fp.write('CHAIN=  ' + str(Variables.ChainPDB) + '			# Chain of pdb file\n')
            fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
            fp.write('TEMP=	1.8			# Temperature\n')
            fp.write('S0=	0.05			# s0, configurational entropy per residue (unfolded)\n')
            fp.write('SC1=    0.05			# configurational entropy per residue (misfolded)\n')
            fp.write('SC0=    0.0			# configurational entropy offset (misfolded)\n')
            fp.write('REM3=	0			# Third cumulant in REM calculation\n')
            fp.write('NEUTRAL=	0		# If 1, Neutral landscape, otherwise population size dependent selection \n')
            fp.write('NPOP=	10			# N_pop, population size\n')
            fp.write('TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
            fp.write(
                'OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"')
            fp.close()

    with open('ProtModel_S_arguments0.txt', 'r') as f:
        data = f.readlines()
        os.system("./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 " + data[zoe])
        # print('Simulation # ' + str(data[zoe].split(':')[1][:-1]))
        #if (zoe + 1) == 1:
            #print(' Simulation # ' + str(zoe + 1))
        #elif (zoe + 1) % 100 == 0:
            #print(' Simulation # ' + str(zoe + 1))
        f.close()


def S_Simu_N(zoe):
    Variables.init()
    LeerSettings.leer()
    try:
        with open('Pop_evol.in', 'r') as of:
            for lineas in of:
                lineas = lineas.split('=')
                lineas[0] = lineas[0].strip()
                if lineas[0] == 'NEUTRAL':
                    if lineas[1][1] == '1':
                        pass
                    else:
                        with open('Pop_evol.in', 'w') as fp:
                            fp.write('PDB= ' + Variables.Template + '		# pdb file\n')
                            fp.write('CHAIN=  ' + str(Variables.ChainPDB) + '			# Chain of pdb file\n')
                            fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
                            fp.write('TEMP=	1.8			# Temperature\n')
                            fp.write('S0=	0.05			# s0, configurational entropy per residue (unfolded)\n')
                            fp.write('SC1=    0.05			# configurational entropy per residue (misfolded)\n')
                            fp.write('SC0=    0.0			# configurational entropy offset (misfolded)\n')
                            fp.write('REM3=	0			# Third cumulant in REM calculation\n')
                            fp.write(
                                'NEUTRAL=	1		# If 1, Neutral landscape, otherwise population size dependent selection \n')
                            fp.write('NPOP=	10			# N_pop, population size\n')
                            fp.write(
                                'TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
                            fp.write(
                                'OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"')
                            fp.close()
            of.close()
    except:
        with open('Pop_evol.in', 'w') as fp:
            fp.write('PDB= ' + Variables.Template + '		# pdb file\n')
            fp.write('CHAIN=  ' + str(Variables.ChainPDB) + '			# Chain of pdb file\n')
            fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
            fp.write('TEMP=	1.8			# Temperature\n')
            fp.write('S0=	0.05			# s0, configurational entropy per residue (unfolded)\n')
            fp.write('SC1=    0.05			# configurational entropy per residue (misfolded)\n')
            fp.write('SC0=    0.0			# configurational entropy offset (misfolded)\n')
            fp.write('REM3=	0			# Third cumulant in REM calculation\n')
            fp.write('NEUTRAL=	1		# If 1, Neutral landscape, otherwise population size dependent selection \n')
            fp.write('NPOP=	10			# N_pop, population size\n')
            fp.write('TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
            fp.write(
                'OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"')
            fp.close()

    with open('ProtModel_S_arguments1.txt', 'r') as f:
        data = f.readlines()
        os.system("./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 " + data[zoe])
        # print('Simulation # ' + str(data[zoe].split(':')[1][:-1]))
        #if (zoe + 1) == 1:
            #print(' Simulation # ' + str(zoe + 1))
        #elif (zoe + 1) % 100 == 0:
            #print(' Simulation # ' + str(zoe + 1))
        f.close()


def CheckSimu():
    Variables.init()
    LeerSettings.leer()
    NumberofRepeats = []
    for x in Variables.Simulations_Name:
        if os.path.isfile(x) == True:
            if os.path.getsize(x) == 0:
                NumberofRepeats.append(x)
                number = x[9:]
                if number[:8] == '00000000':
                    number = number[8:]
                if number[:7] == '0000000':
                    number = number[7:]
                if number[:6] == '000000':
                    number = number[6:]
                if number[:5] == '00000':
                    number = number[5:]
                if number[:4] == '0000':
                    number = number[4:]
                if number[:3] == '000':
                    number = number[3:]
                if number[:2] == '00':
                    number = number[2:]
                if number[:1] == '0':
                    number = number[1:]

                #Eliminar el fichero vacio
                os.system('rm -r ' + x)

                #Crear archivo nuevo para la simulacion. Uno para cada tipo (emp, fit, neu)
                num_b = '-:' + number

                with open("ProtModel_arguments.txt", "r") as openfile:
                    for line in openfile:
                        if num_b in line:
                            with open("ProtModel_Repeat_arguments.txt", "a") as R_file:
                                R_file.write(line)

                with open("ProtModel_S_arguments0.txt", "r") as openfile:
                    for line in openfile:
                        if num_b in line:
                            with open("ProtModel_Repeat_arguments.txt", "a") as R_file:
                                R_file.write(line)

                with open("ProtModel_S_arguments1.txt", "r") as openfile:
                    for line in openfile:
                        if num_b in line:
                            with open("ProtModel_Repeat_arguments.txt", "a") as R_file:
                                R_file.write(line)
        else:
            NumberofRepeats.append(x)
            number = x[9:]
            if number[:8] == '00000000':
                number = number[8:]
            if number[:7] == '0000000':
                number = number[7:]
            if number[:6] == '000000':
                number = number[6:]
            if number[:5] == '00000':
                number = number[5:]
            if number[:4] == '0000':
                number = number[4:]
            if number[:3] == '000':
                number = number[3:]
            if number[:2] == '00':
                number = number[2:]
            if number[:1] == '0':
                number = number[1:]

            # Buscamos en cada archivo y escribimos en uno nuevo las que hay que repetir
            num_b = '-:' + number

            with open("ProtModel_arguments.txt", "r") as openfile:
                for line in openfile:
                    if num_b in line:
                        with open("ProtModel_Repeat_arguments.txt", "a") as R_file:
                            R_file.write(line)

            with open("ProtModel_S_arguments0.txt", "r") as openfile:
                for line in openfile:
                    if num_b in line:
                        with open("ProtModel_Repeat_arguments.txt", "a") as R_file:
                            R_file.write(line)

            with open("ProtModel_S_arguments1.txt", "r") as openfile:
                for line in openfile:
                    if num_b in line:
                        with open("ProtModel_Repeat_arguments.txt", "a") as R_file:
                            R_file.write(line)



    if len(NumberofRepeats) != 0:
        print(len(NumberofRepeats))
        print(NumberofRepeats)

        #Comprobar la clase de simulacion
        with open("ProtModel_Repeat_arguments.txt", "r") as openfile:
            for line in openfile:
                if 'Pop_evol.in' in line:
                    with open("ProtModel_S_arguments0.txt", "r") as openf:
                        linea = openf.readlines()
                        if line in linea:
                            with open('Pop_evol.in', 'r') as f:
                                for lineas in f:
                                    lineas = lineas.split('=')
                                    lineas[0] = lineas[0].strip()
                                    if lineas[0] == 'NEUTRAL':
                                        if lineas[1][2] == '0':
                                            os.system("./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 " + line)
                                        else:
                                            with open ('Pop_evol.in', 'w') as fp:
                                                fp.write('PDB= ' + Variables.Template + ' 		# pdb file\n')
                                                fp.write('CHAIN=  ' + str(Variables.ChainPDB) + '			# Chain of pdb file\n')
                                                fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
                                                fp.write('TEMP=	1.8			# Temperature\n')
                                                fp.write('S0=	0.05			# s0, configurational entropy per residue (unfolded)\n')
                                                fp.write('SC1=    0.05			# configurational entropy per residue (misfolded)\n')
                                                fp.write('SC0=    0.0			# configurational entropy offset (misfolded)\n')
                                                fp.write('REM3=	0			# Third cumulant in REM calculation\n')
                                                fp.write('NEUTRAL=  0		# If 1, Neutral landscape, otherwise population size dependent selection \n')
                                                fp.write('NPOP=	10			# N_pop, population size\n')
                                                fp.write('TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
                                                fp.write('OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"')
                                                fp.close()
                                            os.system("./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 " + line)
                                f.close()
                            openf.close()
                        else:
                            openf.close()

                    with open("ProtModel_S_arguments1.txt", "r") as openf1:
                        linea = openf1.readlines()
                        if line in linea:
                            with open ('Pop_evol.in', 'r') as f:
                                for lineas in f:
                                    lineas = lineas.split('=')
                                    lineas[0] = lineas[0].strip()
                                    if lineas[0] == 'NEUTRAL':
                                        if lineas[1][2] == '1':
                                            os.system("./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 " + line)
                                        else:
                                            with open('Pop_evol.in', 'w') as fp:
                                                fp.write('PDB= ' + Variables.Template + ' 		# pdb file\n')
                                                fp.write('CHAIN=  ' + str(Variables.ChainPDB) + '			# Chain of pdb file\n')
                                                fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
                                                fp.write('TEMP=	1.8			# Temperature\n')
                                                fp.write('S0=	0.05			# s0, configurational entropy per residue (unfolded)\n')
                                                fp.write('SC1=    0.05			# configurational entropy per residue (misfolded)\n')
                                                fp.write('SC0=    0.0			# configurational entropy offset (misfolded)\n')
                                                fp.write('REM3=	0			# Third cumulant in REM calculation\n')
                                                fp.write('NEUTRAL=  1		# If 1, Neutral landscape, otherwise population size dependent selection \n')
                                                fp.write('NPOP=	10			# N_pop, population size\n')
                                                fp.write('TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
                                                fp.write('OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"')
                                                fp.close()
                                            os.system("./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 " + line)
                                f.close()
                            openf.close()
                        else:
                            openf.close()
                else:
                    os.system("./bin/Protein_EvolverABC/ProteinEvolverProtABC1.2.0 " + line)


def RealSummarySta(x):
    # if is the real data get out the extension, if not pass
    SFormat = x[-4:]
    Sequence = x[:-4]
    print('Summary statistics calculation')
    print(' Calculating ' + str(Sequence) + ' summary statistics')

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
    # for x in range(0, Variables.NumSeq):
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
    for y in range(0, Variables.NumSeq - 1):
        Names1.append("Seq" + str(y))

    Names2 = []
    for y in range(1, Variables.NumSeq):
        Names2.append("Seq" + str(y))

    # Check if DeltaGREM will be necessary
    try:
        with open("Input_DeltaGREM.in", "w") as SequencesStability:
            SequencesStability.write('# A) Input files\n')
            SequencesStability.write('ALI=' + str(Sequence) + '.fasta\n')
            SequencesStability.write('CHAIN=  ' + str(Variables.ChainPDB) + '\n')
            SequencesStability.write('PDB=' + Variables.Template + '\n\n\n')
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
    PDBFile = './' + Variables.Template
    with open(PDBFile, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            if record.annotations["chain"] == str(Variables.ChainPDB):
                PDB_Name = ('>' + record.id[:-2])
                PDB_Seq = (record.seq)
                PDB_Seq = str(PDB_Seq)
                PDB_Seq = PDB_Seq.replace('X', '')

    with open("./" + Sequence + ".fasta", "w") as fp:
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

    # Change gaps by AA template

    os.system("./DeltaGREM Input_DeltaGREM.in")  # execute program

    DGREM_RS = []  # DGREM simulations results

    try:
        with open(Sequence + ".fasta_DeltaG.dat") as res:
            for line in res:
                line = line.split('	')
                if line[0][0] != "#":  # + str(Variables.NumAA):
                    if line[7] != Variables.Template + "\n":
                        SEQ = line[1]
                        DGREM_RS.append(SEQ)

    except OSError:
        print("DeltaGREM results file can not be readed")
        sys.exit('ERROR!! Please check DeltaGREM inputs')

        # DGREM mean
    DGREM_Mean = 0
    DGREM_Mean_Sum = 0
    for x in range(0, Variables.NumSeq):
        DGREM_Mean_Sum = DGREM_Mean_Sum + float(DGREM_RS[x])
    DGREM_Mean = (DGREM_Mean_Sum / Variables.NumSeq)
    if Variables.ShowInformationScreen == 'Yes':
        print('  DGREM mean calculated')

    # DGREM sd
    DGREM_sd = 0
    DGREM_sd_Sum = 0
    for x in range(0, Variables.NumSeq):
        DGREM_sd_Sum = DGREM_sd_Sum + ((float(DGREM_RS[x]) - DGREM_Mean) * (float(DGREM_RS[x]) - DGREM_Mean))
    DGREM_sd = (DGREM_sd_Sum / (Variables.NumSeq - 1)) ** (0.5)
    if Variables.ShowInformationScreen == 'Yes':
        print('  DGREM sd calculated')

    # Segregating Sites
    if Variables.Indels == 'Ignored':
        SegSites = 0
        for x in range(0, Variables.NumAA):
            if List1[0][x] == "-":
                pass
            for y in List2:
                if y[x] == "-":
                    pass
                elif List1[0][x] != y[x]:
                    SegSites = SegSites + 1
                    break
        if Variables.ShowInformationScreen == 'Yes':
            print('  Segregating sites calculated')

    if Variables.Indels == 'New State':
        SegSites = 0
        for x in range(0, Variables.NumAA):
            for y in List2:
                if List1[0][x] != y[x]:
                    SegSites = SegSites + 1
                    break
        if Variables.ShowInformationScreen == 'Yes':
            print('  Segregating sites calculated')

    # Grantham Distance per site
    matriz = pd.read_csv("./bin/Grantham.csv")
    Total_Grantham_Position = []
    Grantham = 0
    for e in range(0, Variables.NumAA):

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
    for x in range(0, Variables.NumAA):
        Grantham_Position_Mean = Grantham_Position_Mean + Total_Grantham_Position[x]
    Grantham_Position_Mean = Grantham_Position_Mean / Variables.NumAA
    if Variables.ShowInformationScreen == 'Yes':
        print('  Grantham distance mean calculated')

    # sd
    Grantham_Position_sd_Sum = 0
    for x in range(0, Variables.NumAA):
        Grantham_Position_sd_Sum = Grantham_Position_sd_Sum + ((Total_Grantham_Position[x] - Grantham_Position_Mean) * (
                    Total_Grantham_Position[x] - Grantham_Position_Mean))
    Grantham_Position_sd = (Grantham_Position_sd_Sum / Variables.NumAA) ** (0.5)
    if Variables.ShowInformationScreen == 'Yes':
        print('  Grantham distance sd calculated')

    # sk
    if Grantham_Position_sd == 0.0:
        Grantham_Position_sk = 0
        if Variables.ShowInformationScreen == 'Yes':
            print('  Grantham distance skewness calculated')
    else:
        Grantham_Position_sk_Sum = 0
        for x in range(0, Variables.NumAA):
            Grantham_Position_sk_Sum = Grantham_Position_sk_Sum + ((Total_Grantham_Position[x] - Grantham_Position_Mean) * (
                        Total_Grantham_Position[x] - Grantham_Position_Mean) * (Total_Grantham_Position[
                                                                                  x] - Grantham_Position_Mean))
        Grantham_Position_sk = (Grantham_Position_sk_Sum / Variables.NumAA) / (
                    Grantham_Position_sd * Grantham_Position_sd * Grantham_Position_sd)
        if Variables.ShowInformationScreen == 'Yes':
            print('  Grantham distance skewness calculated')

    # ku
    if Grantham_Position_sd == 0.0:
        Grantham_Position_ku = 0
        if Variables.ShowInformationScreen == 'Yes':
            print('  Grantham distance kurtosis calculated')
    else:
        Grantham_Position_ku_Sum = 0
        for x in range(0, Variables.NumAA):
            Grantham_Position_ku_Sum = Grantham_Position_ku_Sum + ((Total_Grantham_Position[x] - Grantham_Position_Mean) * (
                        Total_Grantham_Position[x] - Grantham_Position_Mean) * (Total_Grantham_Position[
                                                                                  x] - Grantham_Position_Mean) * (
                                                                             Total_Grantham_Position[
                                                                                 x] - Grantham_Position_Mean))
        Grantham_Position_ku = (Grantham_Position_ku_Sum / Variables.NumAA) / (
                    Grantham_Position_sd * Grantham_Position_sd * Grantham_Position_sd * Grantham_Position_sd)
        if Variables.ShowInformationScreen == 'Yes':
            print('  Grantham distance kurtosis calculated')

    SS_HeadersRD = []
    SS_ValRD = []
    for x in Variables.SummaryStatistics.split(' '):
        if x == '1':
            SS_HeadersRD.append('DGREM_Mean')
            SS_ValRD.append(str(DGREM_Mean))
        elif x == '2':
            SS_HeadersRD.append('DGREM_sd')
            SS_ValRD.append(str(DGREM_sd))
        elif x == '3':
            SS_HeadersRD.append('SegSites')
            SS_ValRD.append(str(SegSites))
        elif x == '4':
            SS_HeadersRD.append('Grantham_mean_Position')
            SS_ValRD.append(str(Grantham_Position_Mean))
        elif x == '5':
            SS_HeadersRD.append('Grantham_sd_Position')
            SS_ValRD.append(str(Grantham_Position_sd))
        elif x == '6':
            SS_HeadersRD.append('Grantham_sk_Position')
            SS_ValRD.append(str(Grantham_Position_sk))
        elif x == '7':
            SS_HeadersRD.append('Grantham_ku_Position')
            SS_ValRD.append(str(Grantham_Position_ku))

    Real = open('SSRealData.csv', 'w')
    writer = csv.writer(Real)
    writer.writerow(SS_HeadersRD)
    writer.writerow(SS_ValRD)

    if Variables.SaveSimulations == 'Yes':
        os.system('mkdir Simulations')
        os.system('mv Input_DeltaGREM.in ./Simulations')
        os.system('mv ' + str(Sequence) + '.fasta ./Simulations')
        os.system('mv ' + str(Sequence) + '.fasta_DeltaG.dat ./Simulations')
    else:
        os.system('rm -r Input_DeltaGREM.in')
        os.system('rm -r ' + str(Sequence) + '.fasta')
        os.system('rm -r ' + str(Sequence) + '.fasta_DeltaG.dat')


def SimuSummarySta(zoe):
    Variables.init()
    LeerSettings.leer()
    # time.sleep(1)
    Sequence = Variables.Simulations_Name[zoe]
    SFormat = ''
    # Write a fasta file
    lines = []

    #if zoe == 0:
        #print(' Calculating simulation # ' + str(zoe + 1) + ' summary statistics')
    #else:
        #if (int(zoe) + 1) % 100 == 0:
            #print(' Calculating simulation # ' + str(zoe + 1) + ' summary statistics')

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
    # for x in range(0, Variables.NumSeq):
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
    for y in range(0, Variables.NumSeq - 1):
        Names1.append("Seq" + str(y))

    Names2 = []
    for y in range(1, Variables.NumSeq):
        Names2.append("Seq" + str(y))

    # Check if DeltaGREM will be necessary
    try:
        with open("Input_DeltaGREM" + str(zoe + 1) + ".in", "w") as SequencesStability:
            SequencesStability.write('# A) Input files\n')
            SequencesStability.write('ALI=' + str(Sequence) + '.fasta\n')
            SequencesStability.write('CHAIN=  ' + str(Variables.ChainPDB) + '\n')
            SequencesStability.write('PDB=' + Variables.Template + '\n\n\n')
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

    with open("./" + Sequence + ".fasta", "w") as fp:
        # fp.write(PDB_Name)
        # fp.write(str(PDB_Seq))
        for number, line in enumerate(lines):
            if len(line) > 10:
                fp.write(">")
                fp.write(line[0:9])
                fp.write("\n")
                fp.write(line[10:])
        # Read pdb fasta sequence
        warnings.simplefilter('ignore', PDBConstructionWarning)
        PDBFile = './' + Variables.Template
        with open(PDBFile, 'r') as pdb_file:
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                if record.annotations["chain"] == str(Variables.ChainPDB):
                    PDB_Name = ('>' + record.id[:-2])
                    PDB_Seq = (record.seq)
                    PDB_Seq = str(PDB_Seq)
                    PDB_Seq = PDB_Seq.replace('X', '')
        fp.write(PDB_Name + '\n')
        fp.write(str(PDB_Seq) + '\n')

    os.system("./DeltaGREM Input_DeltaGREM" + str(zoe + 1) + ".in")  # execute program

    DGREM_RS = []  # DGREM simulations results

    try:
        with open(Sequence + ".fasta_DeltaG.dat") as res:
            for line in res:
                line = line.split('	')
                if line[0][0] != "#":  # + str(Variables.NumAA):
                    if line[7] != Variables.Template + "\n":
                        SEQ = line[1]
                        DGREM_RS.append(SEQ)

    except OSError:
        print("DeltaGREM results file can not be readed")
        sys.exit('ERROR!! Please check DeltaGREM inputs')

    # DGREM mean
    DGREM_Mean = 0
    DGREM_Mean_Sum = 0
    for x in range(0, Variables.NumSeq):
        DGREM_Mean_Sum = DGREM_Mean_Sum + float(DGREM_RS[x])
    DGREM_Mean = (DGREM_Mean_Sum / Variables.NumSeq)
    if Variables.ShowInformationScreen == 'Yes':
        print('  DGREM mean calculated')

    # DGREM sd
    DGREM_sd = 0
    DGREM_sd_Sum = 0
    for x in range(0, Variables.NumSeq):
        DGREM_sd_Sum = DGREM_sd_Sum + ((float(DGREM_RS[x]) - DGREM_Mean) * (float(DGREM_RS[x]) - DGREM_Mean))
    DGREM_sd = (DGREM_sd_Sum / (Variables.NumSeq - 1)) ** (0.5)
    if Variables.ShowInformationScreen == 'Yes':
        print('  DGREM sd calculated')

    # Segregating Sites
    if Variables.Indels == 'Ignored':
        SegSites = 0
        for x in range(0, Variables.NumAA):
            if List1[0][x] == "-":
                pass
            for y in List2:
                if y[x] == "-":
                    pass
                elif List1[0][x] != y[x]:
                    SegSites = SegSites + 1
                    break
        if Variables.ShowInformationScreen == 'Yes':
            print('  Segregating sites calculated')

    if Variables.Indels == 'New State':
        SegSites = 0
        for x in range(0, Variables.NumAA):
            for y in List2:
                if List1[0][x] != y[x]:
                    SegSites = SegSites + 1
                    break
        if Variables.ShowInformationScreen == 'Yes':
            print('  Segregating sites calculated')

    # Grantham Distance per site
    matriz = pd.read_csv("./bin/Grantham.csv")
    Total_Grantham_Position = []
    Grantham = 0
    for e in range(0, Variables.NumAA):

        List2 = []
        for b in range(1, len(Sequences)):
            List2.append(Sequences[b])

        for x in List1:
            for y in List2:
                if x[e] != y[e]:
                    row = matriz.columns.get_loc(y[e])
                    Grantham = Grantham + matriz._get_value(row, x[e]) ** 2

            List2.remove(List2[0])
        Total_Grantham_Position.append(Grantham)
        Grantham = 0

    for b in range(1, len(Sequences)):
        List2.append(Sequences[b])

    # Common variables for sd, sk and ku
    Grantham_Position_Mean = 0
    for x in range(0, Variables.NumAA):
        Grantham_Position_Mean = Grantham_Position_Mean + Total_Grantham_Position[x]
    Grantham_Position_Mean = Grantham_Position_Mean / Variables.NumAA
    if Variables.ShowInformationScreen == 'Yes':
        print('  Grantham distance mean calculated')

    # sd
    Grantham_Position_sd_Sum = 0
    for x in range(0, Variables.NumAA):
        Grantham_Position_sd_Sum = Grantham_Position_sd_Sum + ((Total_Grantham_Position[x] - Grantham_Position_Mean) * (
                    Total_Grantham_Position[x] - Grantham_Position_Mean))
    Grantham_Position_sd = (Grantham_Position_sd_Sum / Variables.NumAA) ** (0.5)
    if Variables.ShowInformationScreen == 'Yes':
        print('  Grantham distance sd calculated')

        # sk
    if Grantham_Position_sd == 0.0:
        Grantham_Position_sk = 0
        if Variables.ShowInformationScreen == 'Yes':
            print('  Grantham distance skewness calculated')
    else:
        Grantham_Position_sk_Sum = 0
        for x in range(0, Variables.NumAA):
            Grantham_Position_sk_Sum = Grantham_Position_sk_Sum + ((Total_Grantham_Position[x] - Grantham_Position_Mean) * (
                        Total_Grantham_Position[x] - Grantham_Position_Mean) * (Total_Grantham_Position[
                                                                                  x] - Grantham_Position_Mean))
        Grantham_Position_sk = (Grantham_Position_sk_Sum / Variables.NumAA) / (
                    Grantham_Position_sd * Grantham_Position_sd * Grantham_Position_sd)
        if Variables.ShowInformationScreen == 'Yes':
            print('  Grantham distance skewness calculated')

            # ku
    if Grantham_Position_sd == 0.0:
        Grantham_Position_ku = 0
        if Variables.ShowInformationScreen == 'Yes':
            print('  Grantham distance kurtosis calculated')
    else:
        Grantham_Position_ku_Sum = 0
        for x in range(0, Variables.NumAA):
            Grantham_Position_ku_Sum = Grantham_Position_ku_Sum + ((Total_Grantham_Position[x] - Grantham_Position_Mean) * (
                        Total_Grantham_Position[x] - Grantham_Position_Mean) * (Total_Grantham_Position[
                                                                                  x] - Grantham_Position_Mean) * (
                                                                             Total_Grantham_Position[
                                                                                 x] - Grantham_Position_Mean))
        Grantham_Position_ku = (Grantham_Position_ku_Sum / Variables.NumAA) / (
                    Grantham_Position_sd * Grantham_Position_sd * Grantham_Position_sd * Grantham_Position_sd)
        if Variables.ShowInformationScreen == 'Yes':
            print('  Grantham distance kurtosis calculated')

    SS_csv = []
    for x in Variables.SummaryStatistics.split(' '):
        if x == '1':
            SS_csv.append(str(DGREM_Mean))
        elif x == '2':
            SS_csv.append(str(DGREM_sd))
        elif x == '3':
            SS_csv.append(str(SegSites))
        elif x == '4':
            SS_csv.append(str(Grantham_Position_Mean))
        elif x == '5':
            SS_csv.append(str(Grantham_Position_sd))
        elif x == '6':
            SS_csv.append(str(Grantham_Position_sk))
        elif x == '7':
            SS_csv.append(str(Grantham_Position_ku))

    Simu = open('SSSimulations.csv', 'a')
    writer = csv.writer(Simu)
    writer.writerow(SS_csv)

    if Variables.SaveSimulations == 'Yes':
        os.system('mv ' + str(Sequence) + ' ./Simulations')
        os.system('mv Input_DeltaGREM' + str(zoe + 1) + '.in ./Simulations')
        os.system('mv ' + str(Sequence) + '.fasta ./Simulations')
        os.system('mv ' + str(Sequence) + '.fasta_DeltaG.dat ./Simulations')
    else:
        os.system('rm -r ' + str(Sequence))
        os.system('rm -r Input_DeltaGREM' + str(zoe + 1) + '.in')
        os.system('rm -r ' + str(Sequence) + '.fasta')
        os.system('rm -r ' + str(Sequence) + '.fasta_DeltaG.dat')


def Head():
    SS_Headers = []
    #h = open("SSRealData.ss", "w")
    for x in Variables.SummaryStatistics.split(' '):
        if x == '1':
            SS_Headers.append('DGREM_Mean')
        elif x == '2':
            SS_Headers.append('DGREM_sd')
        elif x == '3':
            SS_Headers.append('SegSites')
        elif x == '4':
            SS_Headers.append('Grantham_mean_Position')
        elif x == '5':
            SS_Headers.append('Grantham_sd_Position')
        elif x == '6':
            SS_Headers.append('Grantham_sk_Position')
        elif x == '7':
            SS_Headers.append('Grantham_ku_Position')
    
    Simu = open('SSSimulations.csv', 'w')
    writer = csv.writer(Simu)
    writer.writerow(SS_Headers)


