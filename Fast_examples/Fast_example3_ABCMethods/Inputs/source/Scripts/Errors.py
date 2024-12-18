import os
import sys
import bin.Scripts.Variables
from Bio import SeqIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import platform

def err():
    # System detection
    with open('ProteinModelerABC.out', 'a') as Error:
        Error.write('________________________System detection________________________\n')
        Error.close()

    sistema = platform.system()
    if sistema == 'Windows':
        with open ('ProteinModelerABC.out', 'a') as Error:
            Error.write('ERROR!!! This program is not executable by Windows\n\n')
            Error.close()
        sys.exit('ERROR!!! This program is not executable by Windows\n')
    else:
        with open ('ProteinModelerABC.out', 'a') as Error:
            Error.write('Working in a ' + str(sistema) + ' machine\n\n')
            Error.close()
        print('Working in a ' + str(sistema) + ' machine')


    ################################
    ## Input Data and Information ##
    ################################
    with open('ProteinModelerABC.out', 'a') as Error:
        Error.write('______________________Checking Settings.txt______________________\n')
        Error.close()

    #MPI file exists
    if bin.Scripts.Variables.MPI == '':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Local version selected\n')
            Error.close()
        print('Local version selected')
    else:
        if os.path.isfile(os.getcwd() + '/' + bin.Scripts.Variables.MPI):
            with open ('ProteinModelerABC.out', 'a') as Error:
                Error.write('Cluster version selected\n')
                Error.close()
            print('Cluster version selected')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('MPI file does not exists\n')
                Error.write('ERROR!!! Please check the MPI file\n')
                Error.close()
            print('MPI file does not exists')
            sys.exit('ERROR!!! Please check the MPI file')

    # NameOfPhylipFile exists
    if os.path.isfile(os.getcwd() + '/' + bin.Scripts.Variables.NameOfPhylipFile):
        with open ('ProteinModelerABC.out', 'a') as Error:
            Error.write('Alignment file exists\n')
            Error.close()
        print('Alignment file exists')
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Alignment file does not exists\n')
            Error.write('ERROR!!! Please check the alignment file\n')
            Error.close()
        print('Alignment file does not exists')
        sys.exit('ERROR!!! Please check the alignment file')

    # NameOfPhylipFile format
    try:
        with open(bin.Scripts.Variables.NameOfPhylipFile) as Filo:
            for line in Filo:
                LongSeq = len(line)
                if LongSeq <= (bin.Scripts.Variables.NumAA + 11):
                    pass
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Sequence name of line must be 10 characters including spaces before sequences starts\n')
                        Error.write('Check if alignment file end with and empty line\n')
                        Error.write('ERROR!!! Sequence name of line must be 10 characters including spaces before sequences starts\n')
                        Error.close()
                    print('Sequence name of line must be 10 characters including spaces before sequences starts')
                    print('Check if alignment file end with and empty line')
                    sys.exit('ERROR!!! Sequence name of line must be 10 characters including spaces before sequences starts')

    except OSError:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Sequence name of line must be 10 characters including spaces before sequences starts\n')
            Error.write('ERROR!!! Sequence name of line must be 10 characters including spaces before sequences starts\n')
            Error.close()
        print('Sequence name of line must be 10 characters including spaces before sequences starts')
        sys.exit('ERROR!!! Sequence name of line must be 10 characters including spaces before sequences starts')


    # NameOfPhylipFile characters
    with open(bin.Scripts.Variables.NameOfPhylipFile, "r") as fp:
    # read an store all lines into list
        lines = fp.readlines()
        #if Indels == 'Default':
        for x in range (0,  bin.Scripts.Variables.NumSeq):
            if len(lines[x]) < 11:
                pass
            else:
                for y in range(10, bin.Scripts.Variables.NumAA + 8):
                    if lines[x][y] == 'A' or lines[x][y] == 'R' or lines[x][y] == 'N' or lines[x][y] == 'D' or lines[x][y] == 'C' or lines[x][y] == 'Q' or lines[x][y] == 'E' or lines[x][y] == 'G' or lines[x][y] == 'H' or lines[x][y] == 'I' or lines[x][y] == 'L' or lines[x][y] == 'K' or lines[x][y] == 'M' or lines[x][y] == 'F' or lines[x][y] == 'P' or lines[x][y] == 'S' or lines[x][y] == 'T' or lines[x][y] == 'W' or lines[x][y] == 'Y' or lines[x][y] == 'V' or lines[x][y] == '-':
                        pass
                    else:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Error in alignment file: ' + lines[x][y] + ' character was found in line ' + str(x+1) + ', position ' + str(y) + '\n')
                            Error.write('Error in alignment file: ' + lines[x][y] + ' character is not allowed\n')
                            Error.close()
                        print('Error in alignment file: ' + lines[x][y] + ' character was found in line ' + str(x+1) + ', position ' + str(y))
                        sys.exit('Error in alignment file: ' + lines[x][y] + ' character is not allowed')


    # Indels
    if bin.Scripts.Variables.Indels == 'Ignored':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Indels are ignored\n')
            Error.close()
        print("Indels are ignored")
    elif bin.Scripts.Variables.Indels == 'New State':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Indels are considered as a new state\n')
            Error.close()
        print("Indels are considered as a new state")
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Indels value are incorrect\n')
            Error.write('ERROR!! Please check the indels value\n')
            Error.close()
        print("Indels value are incorrect")
        sys.exit('ERROR!! Please check the indels value')


    # Template
    if bin.Scripts.Variables.Template != '':
        if bin.Scripts.Variables.Template[-4:] == '.pdb':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write(bin.Scripts.Variables.Template + ' selected as template\n')
                Error.close()
            print(bin.Scripts.Variables.Template + ' selected as template')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('You must select a .pdb file as template\n')
                Error.write('ERROR!! Please check Template value\n')
                Error.close()
            print('You must select a .pdb file as template')
            sys.exit('ERROR!! Please check Template value')
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Template not selected\n')
            Error.write('ERROR!! Please check Template value\n')
            Error.close()
        print('Template not selected')
        sys.exit('ERROR!! Please check Template value')


    #Chain
    if bin.Scripts.Variables.ChainPDB != '':
        if type(bin.Scripts.Variables.ChainPDB)==str:
            if len(bin.Scripts.Variables.ChainPDB) == 1:
                warnings.simplefilter('ignore', PDBConstructionWarning)
                PDBFile = './' + bin.Scripts.Variables.Template
                check = 0
                with open(PDBFile, 'r') as pdb_file:
                    for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                        if record.annotations["chain"] == bin.Scripts.Variables.ChainPDB:
                            check = 1
                        else:
                            pass
                if check == 1:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Chain ' + str(bin.Scripts.Variables.ChainPDB) + ' selected\n')
                        Error.close()
                    print('Chain ' + str(bin.Scripts.Variables.ChainPDB) + ' selected')
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Chain selected not found in the pdb file\n')
                        Error.write('ERROR!! Please check Chain value\n')
                        Error.close()
                    print('Chain selected not found in the pdb file')
                    sys.exit('ERROR!! Please check Chain value')
            else:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Chain must be a unique letter\n')
                    Error.write('ERROR!! Please check Chain value\n')
                    Error.close()
                print('Chain must be a unique letter')
                sys.exit('ERROR!! Please check Chain value')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Chain must be describe with letters from A to Z\n')
                Error.write('ERROR!! Please check Chain value\n')
                Error.close()
            print('Chain must be describe with letters from A to Z')
            sys.exit('ERROR!! Please check Chain value')
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Chain not selected\n')
            Error.write('ERROR!! Please check Chain value\n')
            Error.close()
        print('Chain not selected')
        sys.exit('ERROR!! Please check Chain value')


    #######################################
    ## Settings for the simulation phase ##
    #######################################

    # NumberOfSimulations
    if bin.Scripts.Variables.NumberOfSimulations < 0:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Your number of simulations (' + str(bin.Scripts.Variables.NumberOfSimulations) + ') are not correct\n')
            Error.write('ERROR!!! Number of simulations must be to 1-100000000\n')
            Error.close()
        print("Your number of simulations (" + str(bin.Scripts.Variables.NumberOfSimulations) + ") are not correct")
        sys.exit('ERROR!!! Number of simulations must be to 1-100000000')
    elif bin.Scripts.Variables.NumberOfSimulations > 100000000:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Your number of simulations (' + str(bin.Scripts.Variables.NumberOfSimulations) + ') are not correct\n')
            Error.write('ERROR!!! Number of simulations must be to 1-100000000\n')
            Error.close()
        print("Your number of simulations (" + str(bin.Scripts.Variables.NumberOfSimulations) + ") are not correct")
        sys.exit('ERROR!!! Number of simulations must be to 1-100000000')
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Correct number of simulations\n')
            Error.close()
        print("Correct number of simulations")


    # NumberOfProcessors
    ComputerProcessors = os.cpu_count()

    if bin.Scripts.Variables.NumberOfProcessors == 'Default':
        bin.Scripts.Variables.NumberOfProcessors = ComputerProcessors
    else:
        if int(bin.Scripts.Variables.NumberOfProcessors) < 1:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Number of processor must have a value of 1 at least\n')
                Error.write('ERROR!! Please check NumberOfProcessors value\n')
                Error.close()
            print("Number of processor must have a value of 1 at least")
            sys.exit('ERROR!! Please check NumberOfProcessors value')
        elif int(bin.Scripts.Variables.NumberOfProcessors) <= int(ComputerProcessors):
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write(str(bin.Scripts.Variables.NumberOfProcessors) + ' processors are selected\n')
                Error.close()
            print(str(bin.Scripts.Variables.NumberOfProcessors) + " processors are selected")
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Number of processors introduced is higher than yours computer number of processors\n')
                Error.write('"Please check your computer number of processors\n')
                Error.write('ERROR!! Please check NumberOfProcessors value\n')
                Error.close()
            print("Number of processors introduced is higher than yours computer number of processors")
            print("Please check your computer number of processors")
            sys.exit('ERROR!! Please check NumberOfProcessors value')


    # SaveSimulations
    if bin.Scripts.Variables.SaveSimulations == 'No':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Simulated data is not saved\n')
            Error.close()
        print("Simulated data is not saved")
    elif bin.Scripts.Variables.SaveSimulations == 'Yes':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Simulated data is saved\n')
            Error.close()
        print("Simulated data is saved")
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Simulated data value are incorrect\n')
            Error.write('ERROR!! Please check simulated data value\n')
            Error.close()
        print("Simulated data value are incorrect")
        sys.exit('ERROR!! Please check simulated data value')


    # ShowInformationScreen
    if bin.Scripts.Variables.ShowInformationScreen == 'No':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Running information will not be displayed on the screen\n')
            Error.close()
        print("Running information will not be displayed on the screen")
    elif bin.Scripts.Variables.ShowInformationScreen == 'Yes':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Running information will be displayed on the screen\n')
            Error.close()
        print("Running information will be displayed on the screen")
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Running information value are incorrect\n')
            Error.write('ERROR!! Please check running information value\n')
            Error.close()
        print("Running information value are incorrect")
        sys.exit('ERROR!! Please check running information value')

    #############################    EVOLUTIONARY HISTORY    #############################
    #¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#

    # Coalescent or Tree
    if bin.Scripts.Variables.Coalescent != 'Coal' and bin.Scripts.Variables.Coalescent != 'Phylo':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('User has to perform coalescent or provide a phylogenetic tree input file\n')
            Error.write('ERROR!! Please check CoalescentOrPhylogeny value\n')
            Error.close()
        print("User has to perform coalescent or provide a phylogenetic tree input file")
        sys.exit('ERROR!! Please check CoalescentOrPhylogeny value')


    elif bin.Scripts.Variables.Coalescent == 'Phylo':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('User-specified phylogenetic tree selected\n')
            Error.close()
        print('User-specified phylogenetic tree selected')

        try:
            with open(bin.Scripts.Variables.Tree) as tree:
                tree.close()
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('User-specified phylogenetic tree exists\n')
                Error.write('Note that format and information provided in the user-specified phylogenetic tree is not reviewed so in case of program failure we recommend to check the tree input file\n\n')
                Error.close()

        except OSError:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('User-specified phylogenetic tree does not exist\n')
                Error.write('ERROR!!! Cannot read ' + str(bin.Scripts.Variables.Tree) + '. Check phylogenetic tree input file\n')
                Error.close()
            print('User-specified phylogenetic tree does not exist')
            sys.exit('ERROR!!! Cannot read ' + str(bin.Scripts.Variables.Tree) + '. Check phylogenetic tree input file')


    elif bin.Scripts.Variables.Coalescent == 'Coal':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Coalescent method selected\n')
            Error.close()
        print('Coalescent method selected')

        # HaploidDiploid
        if bin.Scripts.Variables.HaploidDiploid == 1:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Haploid data are selected\n')
                Error.close()
            print("Haploid data are selected")
        elif bin.Scripts.Variables.HaploidDiploid == 2:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Diploid data are selected\n')
                Error.close()
            print("Diploid data are selected")
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Haploid/Diploid information value are incorrect\n')
                Error.write('ERROR!! Please check Haploid/Diploid value\n')
                Error.close()
            print("Haploid/Diploid information value are incorrect")
            sys.exit('ERROR!! Please check Haploid/Diploid value')


        # PopulationSize
        if bin.Scripts.Variables.PopulationSize > 0:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('PopulationSize selected\n')
                Error.close()
            print("PopulationSize selected")
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('PopulationSize information value are incorrect\n')
                Error.write('ERROR!! Please check PopulationSize value\n')
                Error.close()
            print("PopulationSize information value are incorrect")
            sys.exit('ERROR!! Please check PopulationSize value')


        # DatedTips
        if bin.Scripts.Variables.DatedTips == '':
            pass
        else:
            NumOfDatedTips = bin.Scripts.Variables.DatedTips[0]
            DatedTipsList = list(bin.Scripts.Variables.DatedTips.split(" "))
            NumOfTotalDatedTips = int(NumOfDatedTips) * 3 + 1
            if len(DatedTipsList) == int(NumOfTotalDatedTips):
                for i in range(1, int(NumOfDatedTips)):
                    if int(DatedTipsList[(i * 3)]) < bin.Scripts.Variables.NumSeq:
                        pass
                    else:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('The highest value of the ' + str(i) + ' interval seems incorrect\n')
                            Error.write('ERROR!! Please check DatedTips value\n')
                            Error.close()
                        print("The highest value of the " + str(i) + " interval seems incorrect")
                        sys.exit('ERROR!! Please check DatedTips value')
                for e in range(1, (int(NumOfDatedTips) - 1)):
                    if int(DatedTipsList[(e * 3 - 1)]) < int(DatedTipsList[(e * 3 + 2)]):
                        pass
                    else:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('The lowest value of the ' + str(e) + ' interval seems incorrect\n')
                            Error.write('ERROR!! Please check DatedTips value\n')
                            Error.close()
                        print("The lowest value of the " + str(e) + " interval seems incorrect")
                        sys.exit('ERROR!! Please check DatedTips value')
            else:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Number of elements of Dated Tips are incorrect\n')
                    Error.write('ERROR!! Please check DatedTips value\n')
                    Error.close()
                print("Number of elements of Dated Tips are incorrect")
                sys.exit('ERROR!! Please check DatedTips value')


        # GenerationTime
        if bin.Scripts.Variables.GenerationTime == '':
            pass
        else:
            if bin.Scripts.Variables.GenerationTime.split(' ')[0] == 'fix':
                ListGenerationTime = list(bin.Scripts.Variables.GenerationTime.split(" "))
                if len(ListGenerationTime) == 2:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Parameter value is fix with a value of ' + ListGenerationTime[1] + '\n')
                        Error.close()
                    print("Parameter value is fix with a value of " + ListGenerationTime[1])
                elif len(ListGenerationTime) < 2:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Less generation time parameterss than was expected\n')
                        Error.write('ERROR!! Please check GenerationTime value\n')
                        Error.close()
                    print("Less generation time parameters than was expected")
                    sys.exit('ERROR!! Please check GenerationTime value')
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('More generation time parameterss than was expected\n')
                        Error.write('ERROR!! Please check GenerationTime value\n')
                        Error.close()
                    print("More generation time parameterss than was expected")
                    sys.exit('ERROR!! Please check GenerationTime value')
            elif bin.Scripts.Variables.GenerationTime.split(' ')[0] == 'uniform':
                ListGenerationTime = list(bin.Scripts.Variables.GenerationTime.split(" "))
                if len(ListGenerationTime) == 3:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Parameter values sampled from a uniform distribution of ' + ListGenerationTime[1] + ' - ' + ListGenerationTime[2] + '\n')
                        Error.close()
                    print("Parameter values sampled from a uniform distribution of " + ListGenerationTime[1] + "-" + ListGenerationTime[2])
                elif len(ListGenerationTime) < 3:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Less generation time parameters than was expected\n')
                        Error.write('ERROR!! Please check GenerationTime value\n')
                        Error.close()
                    print("Less generation time parameters than was expected")
                    sys.exit('ERROR!! Please check GenerationTime value')
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('More generation time parameters than was expected\n')
                        Error.write('ERROR!! Please check GenerationTime value\n')
                        Error.close()
                    print("More generation time parameters than was expected")
                    sys.exit('ERROR!! Please check GenerationTime value')
            else:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Generation time parameter value error\n')
                    Error.write('ERROR!! Please check GenerationTime value\n')
                    Error.close()
                print("Generation time parameter value error")
                sys.exit('ERROR!! Please check GenerationTime value')


        # GrowthRate
        if bin.Scripts.Variables.GrowthRate.split(' ')[0] == 0:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Exponential growth rate selected\n')
                Error.close()
            print("Exponential growth rate selected")
            if len(bin.Scripts.Variables.GrowthRate.split(' ')) == 2:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Exponential growth rate seems correct\n')
                    Error.close()
                print("Exponential growth rate seems correct")
            elif len(bin.Scripts.Variables.GrowthRate.split(' ')) < 2:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less values than expected in exponential growth rate\n')
                    Error.write('ERROR!! Please check GrowthRate value\n')
                    Error.close()
                print("Less values than expected in exponential growth rate")
                sys.exit('ERROR!! Please check GrowthRate value')
            elif len(bin.Scripts.Variables.GrowthRate.split(' ')) > 2:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More values than expected in exponential growth rate\n')
                    Error.write('ERROR!! Please check GrowthRate value\n')
                    Error.close()
                print("More values than expected in exponential growth rate")
                sys.exit('ERROR!! Please check GrowthRate value')
        elif bin.Scripts.Variables.GrowthRate.split(' ')[0] == 1:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Demographic periods selected\n')
                Error.close()
            print("Demographic periods selected")
            if len(bin.Scripts.Variables.GrowthRate.split(' ')) == (3 * int(bin.Scripts.Variables.GrowthRate.split(' ')[1]) + 2):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Demographic periods seems correct\n')
                    Error.close()
                print("Demographic periods seems correct")
            elif len(bin.Scripts.Variables.GrowthRate.split(' ')) > (3 * int(bin.Scripts.Variables.GrowthRate.split(' ')[1]) + 2):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less values than expected in demographic periods\n')
                    Error.write('ERROR!! Please check GrowthRate value\n')
                    Error.close()
                print("Less values than expected in demographic periods")
                sys.exit('ERROR!! Please check GrowthRate value')
            elif len(bin.Scripts.Variables.GrowthRate.split(' ')) < (3 * int(bin.Scripts.Variables.GrowthRate.split(' ')[1]) + 2):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More values than expected in demographic periods\n')
                    Error.write('ERROR!! Please check GrowthRate value\n')
                    Error.close()
                print("More values than expected in demographic periods")
                sys.exit('ERROR!! Please check GrowthRate value')
        else:
            pass


        # MigrationModel
        if bin.Scripts.Variables.MigrationModel.split(' ')[0] == 1:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Island model selected\n')
                Error.close()
            print("Island model selected")
            if len(bin.Scripts.Variables.MigrationModel.split(' ')) == (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Island model seems correct\n')
                    Error.close()
                print("Island model seems correct")
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) < (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less values than expected in migration model\n')
                    Error.write('ERROR!! Please check MigrationModel value\n')
                    Error.close()
                print("Less values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) > (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More values than expected in migration model\n')
                    Error.write('ERROR!! Please check MigrationModel value\n')
                    Error.close()
                print("More values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
        elif bin.Scripts.Variables.MigrationModel.split(' ')[0] == 2:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Stepping-stone model selected\n')
                Error.close()
            print("Stepping-stone model selected")
            if len(bin.Scripts.Variables.MigrationModel.split(' ')) == (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Stepping-stone seems correct\n')
                    Error.close()
                print("Stepping-stone seems correct")
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) < (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less values than expected in migration model\n')
                    Error.write('ERROR!! Please check MigrationModel value\n')
                    Error.close()
                print("Less values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) > (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More values than expected in migration model\n')
                    Error.write('ERROR!! Please check MigrationModel value\n')
                    Error.close()
                print("More values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
        elif bin.Scripts.Variables.MigrationModel.split(' ')[0] == 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Continent-island model selected\n')
                Error.close()
            print("Continent-island model selected")
            if len(bin.Scripts.Variables.MigrationModel.split(' ')) == (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Continent-island seems correct\n')
                    Error.close()
                print("Continent-island seems correct")
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) < (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less values than expected in migration model\n')
                    Error.write('ERROR!! Please check MigrationModel value\n')
                    Error.close()
                print("Less values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) > (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less values than expected in migration model\n')
                    Error.write('ERROR!! Please check MigrationModel value\n')
                    Error.close()
                print("More values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
        else:
            pass


        # MigrationRate
        if bin.Scripts.Variables.MigrationRate == '':
            pass
        else:
            if len(bin.Scripts.Variables.MigrationRate.split(' ')) == (2 * int(bin.Scripts.Variables.MigrationRate.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Migration rate seems correct\n')
                    Error.close()
                print("Migration rate seems correct")
            elif len(bin.Scripts.Variables.MigrationRate.split(' ')) < (2 * int(bin.Scripts.Variables.MigrationRate.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less values than expected in migration rate\n')
                    Error.write('ERROR!! Please check MigrationRate value\n')
                    Error.close()
                print("Less values than expected in migration rate")
                sys.exit('ERROR!! Please check MigrationRate value')
            elif len(bin.Scripts.Variables.MigrationRate.split(' ')) > (2 * int(bin.Scripts.Variables.MigrationRate.split(' ')[1])):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More values than expected in migration rate\n')
                    Error.write('ERROR!! Please check MigrationRate value\n')
                    Error.close()
                print("More values than expected in migration rate")
                sys.exit('ERROR!! Please check MigrationRate value')


        # ConvergenceDemes
        if bin.Scripts.Variables.ConvergenceDemes == '':
            pass
        else:
            if len(bin.Scripts.Variables.ConvergenceDemes.split(' ')) == (3 * int(bin.Scripts.Variables.ConvergenceDemes(' ')[1]) + 1):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Convergence of demes seems correct\n')
                    Error.close()
                print("Convergence of demes seems correct")
            elif len(bin.Scripts.Variables.ConvergenceDemes.split(' ')) < (3 * int(bin.Scripts.Variables.ConvergenceDemes.split(' ')[1]) + 1):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less values than expected in convergence demes\n')
                    Error.write('ERROR!! Please check MigrationRate value\n')
                    Error.close()
                print("Less values than expected in convergence demes")
                sys.exit('ERROR!! Please check ConvergenceDeme value')
            elif len(bin.Scripts.Variables.ConvergenceDemes.split(' ')) > (3 * int(bin.Scripts.Variables.ConvergenceDemes.split(' ')[1]) + 1):
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More values than expected in convergence demes\n')
                    Error.write('ERROR!! Please check MigrationRate value\n')
                    Error.close()
                print("More values than expected in convergence demes")
                sys.exit('ERROR!! Please check ConvergenceDemes value')


        # SubstitutionRate
        if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'fix':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 2:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Parameter value is fix with a value of ' + ListSubstitutionRate[1] + '\n')
                    Error.close()
                print("Parameter value is fix with a value of " + ListSubstitutionRate[1])
            elif len(ListSubstitutionRate) < 2:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less substitution rate parameters than was expected\n')
                    Error.write('ERROR!! Please check SubstitutionRate value\n')
                    Error.close()
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More substitution rate parameters than was expected\n')
                    Error.write('ERROR!! Please check SubstitutionRate value\n')
                    Error.close()
                print("More substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')

        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'uniform':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 3:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Parameter values sampled from a uniform distribution of ' + ListSubstitutionRate[1] + ' - ' + ListSubstitutionRate[2] + '\n')
                    Error.close()
                print("Parameter values sampled from a uniform distribution of " + ListSubstitutionRate[1] + " - " + ListSubstitutionRate[2])
            elif len(ListSubstitutionRate) < 3:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less substitution rate parameters than was expected\n')
                    Error.write('ERROR!! Please check SubstitutionRate value\n')
                    Error.close()
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More substitution rate parameters than was expected\n')
                    Error.write('ERROR!! Please check SubstitutionRate value\n')
                    Error.close()
                print("More substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            if ListSubstitutionRate[1] > ListSubstitutionRate[2]:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Second value expected to be higher than the previous one\n')
                    Error.write('ERROR!! Please check SubstitutionRate value\n')
                    Error.close()
                print("Second value expected to be higher than the previous one")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            elif ListSubstitutionRate[1] == ListSubstitutionRate[2]:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Second value expected to be higher than the previous one\n')
                    Error.write('ERROR!! Please check SubstitutionRate value\n')
                    Error.close()
                print("Second value expected to be higher than the previous one")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Parameter values sampled from a uniform distribution of ' + ListSubstitutionRate[1] + ' - ' +ListSubstitutionRate[2] + ' seem correct\n')
                    Error.close()
                print("Parameter values sampled from a uniform distribution of " + ListSubstitutionRate[1] + " - " + ListSubstitutionRate[2] + " seem correct")

        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'norm':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 3:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Parameter values sampled from a normal distribution of ' + ListSubstitutionRate[1] + ' - ' + ListSubstitutionRate[2] + '\n')
                    Error.close()
                print("Parameter values sampled from a normal distribution of " + ListSubstitutionRate[1] + " - " + ListSubstitutionRate[2])
            elif len(ListSubstitutionRate) < 3:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less substitution rate parameters than was expected\n')
                    Error.write('ERROR!! Please check SubstitutionRate value\n')
                    Error.close()
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                if len(ListSubstitutionRate) < 6:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Incorrect substitution rate parameters than was expected\n')
                        Error.write('ERROR!! Please check SubstitutionRate value\n')
                        Error.close()
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                if len(ListSubstitutionRate) > 6:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Incorrect substitution rate parameters than was expected\n')
                        Error.write('ERROR!! Please check SubstitutionRate value\n')
                        Error.close()
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                else:
                    if ListSubstitutionRate[3] == 't':
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Truncated distribution detected\n')
                            Error.close()
                        print("Truncated distribution detected")
                        if ListSubstitutionRate[4] > ListSubstitutionRate[5]:
                            with open('ProteinModelerABC.out', 'a') as Error:
                                Error.write('Fifth value expected to be higher than the previous one\n')
                                Error.write('ERROR!! Please check SubstitutionRate value\n')
                                Error.close()
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        elif ListSubstitutionRate[4] == ListSubstitutionRate[5]:
                            with open('ProteinModelerABC.out', 'a') as Error:
                                Error.write('Fifth value expected to be higher than the previous one\n')
                                Error.write('ERROR!! Please check SubstitutionRate value\n')
                                Error.close()
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        else:
                            pass
                    else:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('t value or less parameters were expected\n')
                            Error.write('ERROR!! Please check SubstitutionRate value\n')
                            Error.close()
                        print("t value or less parameters were expected")
                        sys.exit('ERROR!! Please check SubstitutionRate value')

        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'exp':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 2:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Parameter values sampled from a normal distribution of ' + ListSubstitutionRate[1] + '\n')
                    Error.close()
                print("Parameter values sampled from a normal distribution of " + ListSubstitutionRate[1])
            elif len(ListSubstitutionRate) < 2:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less substitution rate parameters than was expected\n')
                    Error.write('ERROR!! Please check SubstitutionRate value\n')
                    Error.close()
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                if len(ListSubstitutionRate) < 5:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Incorrect substitution rate parameters than was expected\n')
                        Error.write('ERROR!! Please check SubstitutionRate value\n')
                        Error.close()
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                if len(ListSubstitutionRate) > 5:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Incorrect substitution rate parameters than was expected\n')
                        Error.write('ERROR!! Please check SubstitutionRate value\n')
                        Error.close()
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                else:
                    if ListSubstitutionRate[2] == 't':
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Truncated distribution detected\n')
                            Error.close()
                        print("Truncated distribution detected")
                        if ListSubstitutionRate[3] > ListSubstitutionRate[4]:
                            with open('ProteinModelerABC.out', 'a') as Error:
                                Error.write('Fourth value expected to be higher than the previous one\n')
                                Error.write('ERROR!! Please check SubstitutionRate value\n')
                                Error.close()
                            print("Fourth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        elif ListSubstitutionRate[3] == ListSubstitutionRate[4]:
                            with open('ProteinModelerABC.out', 'a') as Error:
                                Error.write('Fourth value expected to be higher than the previous one\n')
                                Error.write('ERROR!! Please check SubstitutionRate value\n')
                                Error.close()
                            print("Fourth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        else:
                            pass
                    else:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('t value or less parameters were expected\n')
                            Error.write('ERROR!! Please check SubstitutionRate value\n')
                            Error.close()
                        print("t value or less parameters were expected")
                        sys.exit('ERROR!! Please check SubstitutionRate value')

        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 3:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Parameter values sampled from a normal distribution of ' + ListSubstitutionRate[1] + ' - ' + ListSubstitutionRate[2] + '\n')
                    Error.close()
                print("Parameter values sampled from a normal distribution of " + ListSubstitutionRate[1] + "-" + ListSubstitutionRate[2])
            elif len(ListSubstitutionRate) < 3:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less substitution rate parameters than was expected\n')
                    Error.write('ERROR!! Please check SubstitutionRate value\n')
                    Error.close()
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                if len(ListSubstitutionRate) < 6:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Incorrect substitution rate parameters than was expected\n')
                        Error.write('ERROR!! Please check SubstitutionRate value\n')
                        Error.close()
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                if len(ListSubstitutionRate) > 6:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Incorrect substitution rate parameters than was expected\n')
                        Error.write('ERROR!! Please check SubstitutionRate value\n')
                        Error.close()
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                else:
                    if ListSubstitutionRate[3] == 't':
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Truncated distribution detected\n')
                            Error.close()
                        print("Truncated distribution detected")
                        if ListSubstitutionRate[4] > ListSubstitutionRate[5]:
                            with open('ProteinModelerABC.out', 'a') as Error:
                                Error.write('Fifth value expected to be higher than the previous one\n')
                                Error.write('ERROR!! Please check SubstitutionRate value\n')
                                Error.close()
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        elif ListSubstitutionRate[4] == ListSubstitutionRate[5]:
                            with open('ProteinModelerABC.out', 'a') as Error:
                                Error.write('Fifth value expected to be higher than the previous one\n')
                                Error.write('ERROR!! Please check SubstitutionRate value\n')
                                Error.close()
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        else:
                            pass
                    else:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('t value or less parameters were expected\n')
                            Error.write('ERROR!! Please check SubstitutionRate value\n')
                            Error.close()
                        print("t value or less parameters were expected")
                        sys.exit('ERROR!! Please check SubstitutionRate value')

        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'beta':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 3:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write( 'Parameter values sampled from a normal distribution of ' + ListSubstitutionRate[1] + ' - ' + ListSubstitutionRate[2] + '\n')
                    Error.close()
                print("Parameter values sampled from a normal distribution of " + ListSubstitutionRate[1] + " - " + ListSubstitutionRate[2])
            elif len(ListSubstitutionRate) < 3:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less substitution rate parameters than was expected\n')
                    Error.write('ERROR!! Please check SubstitutionRate value\n')
                    Error.close()
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                if len(ListSubstitutionRate) < 6:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Incorrect substitution rate parameters than was expected\n')
                        Error.write('ERROR!! Please check SubstitutionRate value\n')
                        Error.close()
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                if len(ListSubstitutionRate) > 6:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Incorrect substitution rate parameters than was expected\n')
                        Error.write('ERROR!! Please check SubstitutionRate value\n')
                        Error.close()
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                else:
                    if ListSubstitutionRate[3] == 't':
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Truncated distribution detected\n')
                            Error.close()
                        print("Truncated distribution detected")
                        if ListSubstitutionRate[4] > ListSubstitutionRate[5]:
                            with open('ProteinModelerABC.out', 'a') as Error:
                                Error.write('Fifth value expected to be higher than the previous one\n')
                                Error.write('ERROR!! Please check SubstitutionRate value\n')
                                Error.close()
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        elif ListSubstitutionRate[4] == ListSubstitutionRate[5]:
                            with open('ProteinModelerABC.out', 'a') as Error:
                                Error.write('Fifth value expected to be higher than the previous one\n')
                                Error.write('ERROR!! Please check SubstitutionRate value\n')
                                Error.close()
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        else:
                            pass
                    else:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('t value or less parameters were expected\n')
                            Error.write('ERROR!! Please check SubstitutionRate value\n')
                            Error.close()
                        print("t value or less parameters were expected")
                        sys.exit('ERROR!! Please check SubstitutionRate value')

        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Substitution rate parameter value error\n')
                Error.write('ERROR!! Please check SubstitutionRate value\n')
                Error.close()
            print("Substitution rate parameter value error")
            sys.exit('ERROR!! Please check SubstitutionRate value')



    ##############################    SUBSTITUTION MODEL    ##############################
    #¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#

    # SubstitutionModel
    for x in bin.Scripts.Variables.SubstitutionModel.split(' '):
        if x == 'Blosum62':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Blosum62 substitution model selected\n')
                Error.close()
            print("Blosum62 substitution model selected")
        elif x == 'CpRev':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('CpRev substitution model selected\n')
                Error.close()
            print("CpRev substitution model selected")
        elif x == 'Dayhoff':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Dayhoff substitution model selected\n')
                Error.close()
            print("Dayhoff substitution model selected")
        elif x == 'DayhoffDCMUT':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('DayhoffDCMUT substitution model selected\n')
                Error.close()
            print("DayhoffDCMUT substitution model selected")
        elif x == 'HIVb':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('HIVb substitution model selected\n')
                Error.close()
            print("HIVb substitution model selected")
        elif x == 'HIVw':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('HIVw substitution model selected\n')
                Error.close()
            print("HIVw substitution model selected")
        elif x == 'JTT':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('JTT substitution model selected\n')
                Error.close()
            print("JTT substitution model selected")
        elif x == 'JonesDCMUT':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('JonesDCMUT substitution model selected\n')
                Error.close()
            print("JonesDCMUT substitution model selected")
        elif x == 'LG':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('LG substitution model selected\n')
                Error.close()
            print("LG substitution model selected")
        elif x == 'Mtart':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Mtart substitution model selected\n')
                Error.close()
            print("Mtart substitution model selected")
        elif x == 'Mtmam':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Mtmam substitution model selected\n')
                Error.close()
            print("Mtmam substitution model selected")
        elif x == 'Mtrev24':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Mtrev24 substitution model selected\n')
                Error.close()
            print("Mtrev24 substitution model selected")
        elif x == 'RtRev':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('RtRev substitution model selected\n')
                Error.close()
            print("RtRev substitution model selected")
        elif x == 'VT':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('VT substitution model selected\n')
                Error.close()
            print("VT substitution model selected")
        elif x == 'WAG':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('WAG substitution model selected\n')
                Error.close()
            print("WAG substitution model selected")
        elif x == 'UserEAAM':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('UserEAAM substitution model selected\n')
                Error.close()
            print("UserEAAM substitution model selected")
        elif x == 'Fitness':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Fitness substitution model selected\n')
                Error.close()
            print("Fitness substitution model selected")
        elif x == 'Neutral':
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Neutral substitution model selected\n')
                Error.close()
            print("Neutral substitution model selected")
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Substitution model selected: ' + str(x) + ' was not recognised\n')
                Error.write('ERROR!! Please check SubstitutionModel value\n')
                Error.close()
            print('Substitution model selected: ' + str(x) + ' was not recognised\n')
            sys.exit('ERROR!! Please check SubstitutionModel value')

            #Empirical substitution model parameters
# """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # AminoacidFrequencies
    if bin.Scripts.Variables.AminoacidFrequencies.split(' ')[0] == 'fix':
        ListAminoacidFrequencies = list(bin.Scripts.Variables.AminoacidFrequencies.split(' '))
        if len(ListAminoacidFrequencies) == 21:
            for AminoacidPosition in range(1, 20):
                if float(ListAminoacidFrequencies[AminoacidPosition]) > 1:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Aminoacid frequency values must be 0-1\n')
                        Error.write('ERROR!! Please check AminoacidFrequencies values\n')
                        Error.close()
                    print('Aminoacid frequency values must be 0-1')
                    sys.exit('ERROR!! Please check AminoacidFrequencies values')
                elif float(ListAminoacidFrequencies[AminoacidPosition]) < 0:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Aminoacid frequency values must be 0-1\n')
                        Error.write('ERROR!! Please check AminoacidFrequencies values\n')
                        Error.close()
                    print('Aminoacid frequency values must be 0-1')
                    sys.exit('ERROR!! Please check AminoacidFrequencies values')
                else:
                    pass
                    #print('Aminoacid frequency values seem correct')
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Aminoacid frequency values seem correct\n')
                Error.close()
            print('Aminoacid frequency values seem correct')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Aminoacid frequencies length seems incorrect\n')
                Error.write('ERROR!! Please check AminoacidFrequencies values\n')
                Error.close()
            print("Aminoacid frequencies length seems incorrect")
            sys.exit('ERROR!! Please check AminoacidFrequencies values')
    elif bin.Scripts.Variables.AminoacidFrequencies.split(' ')[0] == 'dirichlet':
        ListAminoacidFrequencies = list(bin.Scripts.Variables.AminoacidFrequencies.split(' '))
        if len(ListAminoacidFrequencies) == 21:
            for AminoacidPosition in range(1, 20):
                if float(ListAminoacidFrequencies[AminoacidPosition]) > 1:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Aminoacid frequency values must be 0-1\n')
                        Error.write('ERROR!! Please check AminoacidFrequencies values\n')
                        Error.close()
                    print('Aminoacid frequency values must be 0-1')
                    sys.exit('ERROR!! Please check AminoacidFrequencies values')
                elif float(ListAminoacidFrequencies[AminoacidPosition]) < 0:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Aminoacid frequency values must be 0-1\n')
                        Error.write('ERROR!! Please check AminoacidFrequencies values\n')
                        Error.close()
                    print('Aminoacid frequency values must be 0-1')
                    sys.exit('ERROR!! Please check AminoacidFrequencies values')
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Aminoacid frequency values seem correct\n')
                        Error.close()
                    print('Aminoacid frequency values seem correct')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Aminoacid frequencies length seems incorrect\n')
                Error.write('ERROR!! Please check AminoacidFrequencies values\n')
                Error.close()
            print("Aminoacid frequencies length seems incorrect")
            sys.exit('ERROR!! Please check AminoacidFrequencies values')
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Aminoacid frequencies parameter value error\n')
            Error.write('ERROR!! Please check AminoacidFrequencies values\n')
            Error.close()
        print("Aminoacid frequencies parameter value error")
        sys.exit('ERROR!! Please check AminoacidFrequencies values')


    # RateHetSites
    if bin.Scripts.Variables.RateHetSites.split(' ')[0] == '':
        pass
    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'fix':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 2:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Rate of heterogeneity across sites value is fix with a value of ' + ListRateHetSites[1] + '\n')
                Error.close()
            print("Rate of heterogeneity across sites value is fix with a value of " + ListRateHetSites[1])
        elif len(ListRateHetSites) < 2:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less rate of heterogeneity across sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check SubstitutionRate value\n')
                Error.close()
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check SubstitutionRate value')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('More rate of heterogeneity across sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check SubstitutionRate value\n')
                Error.close()
            print("More rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check SubstitutionRate value')

    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'uniform':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Rate of heterogeneity across sites values sampled from a uniform distribution of ' + ListRateHetSites[1] + ' - ' + ListRateHetSites[2] + '\n')
                Error.close()
            print("Rate of heterogeneity across sites values sampled from a uniform distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2])
        elif len(ListRateHetSites) < 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less rate of heterogeneity across sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check RateHetSites value\n')
                Error.close()
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('More rate of heterogeneity across sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check RateHetSites value\n')
                Error.close()
            print("More rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        if ListRateHetSites[1] > ListRateHetSites[2]:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Second value expected to be higher than the previous one\n')
                Error.write('ERROR!! Please check RateHetSites value\n')
                Error.close()
            print("Second value expected to be higher than the previous one\n")
            sys.exit('ERROR!! Please check RateHetSites value')
        elif ListRateHetSites[1] == ListRateHetSites[2]:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Second value expected to be higher than the previous one\n')
                Error.write('ERROR!! Please check RateHetSites value\n')
                Error.close()
            print("Second value expected to be higher than the previous one")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Rate of heterogeneity across sites sampled from a uniform distribution of ' + ListRateHetSites[1] + ' - ' + ListRateHetSites[2] + ' seem correct\n')
                Error.close()
            print("Rate of heterogeneity across sites sampled from a uniform distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2] + " seem correct")

    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'norm':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Rate of heterogeneity across sites values sampled from a normal distribution of ' + ListRateHetSites[1] + ' - ' + ListRateHetSites[2] + '\n')
                Error.close()
            print("Rate of heterogeneity across sites values sampled from a normal distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2])
        elif len(ListRateHetSites) < 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less rate of heterogeneity across sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check RateHetSites value\n')
                Error.close()
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            if len(ListRateHetSites) < 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less rate of heterogeneity across sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check RateHetSites value\n')
                    Error.close()
                print("Less rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            if len(ListRateHetSites) > 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More rate of heterogeneity across sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check RateHetSites value\n')
                    Error.close()
                print("More rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            else:
                if ListRateHetSites[3] == 't':
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Truncated distribution detected\n')
                        Error.close()
                    print("Truncated distribution detected")
                    if ListRateHetSites[4] > ListRateHetSites[5]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check RateHetSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    elif ListRateHetSites[1] == ListRateHetSites[2]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check RateHetSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    else:
                        pass
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('t value or less parameters were expected\n')
                        Error.write('ERROR!! Please check RateHetSites value\n')
                        Error.close()
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check RateHetSites value')

    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'exp':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 2:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Rate of heterogeneity across sites values sampled from a normal distribution of ' + ListRateHetSites[1] + '\n')
                Error.close()
            print("Rate of heterogeneity across sites values sampled from a normal distribution of " + ListRateHetSites[1])
        elif len(ListRateHetSites) < 2:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less rate of heterogeneity across sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check RateHetSites value\n')
                Error.close()
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            if len(ListRateHetSites) < 5:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less rate of heterogeneity across sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check RateHetSites value\n')
                    Error.close()
                print("Less rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            if len(ListRateHetSites) > 5:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More rate of heterogeneity across sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check RateHetSites value\n')
                    Error.close()
                print("More rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            else:
                if ListRateHetSites[2] == 't':
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Truncated distribution detected\n')
                        Error.close()
                    print("Truncated distribution detected")
                    if ListRateHetSites[3] > ListRateHetSites[4]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check RateHetSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    else:
                        pass
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('t value or less parameters were expected\n')
                        Error.write('ERROR!! Please check RateHetSites value\n')
                        Error.close()
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check RateHetSites value')

    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'gamma':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Rate of heterogeneity across sites values sampled from a normal distribution of ' + ListRateHetSites[1] + ' - ' + ListRateHetSites[2] + '\n')
                Error.close()
            print("Rate of heterogeneity across sites values sampled from a normal distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2])
        elif len(ListRateHetSites) < 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less rate of heterogeneity across sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check RateHetSites value\n')
                Error.close()
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            if len(ListRateHetSites) < 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less rate of heterogeneity across sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check RateHetSites value\n')
                    Error.close()
                print("Less rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            if len(ListRateHetSites) > 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More rate of heterogeneity across sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check RateHetSites value\n')
                    Error.close()
                print("More rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            else:
                if ListRateHetSites[3] == 't':
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Truncated distribution detected\n')
                        Error.close()
                    print("Truncated distribution detected")
                    if ListRateHetSites[4] > ListRateHetSites[5]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check RateHetSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    elif ListRateHetSites[1] == ListRateHetSites[2]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check RateHetSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    else:
                        pass
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('t value or less parameters were expected\n')
                        Error.write('ERROR!! Please check RateHetSites value\n')
                        Error.close()
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check RateHetSites value')

    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'beta':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Rate of heterogeneity across sites values sampled from a normal distribution of ' + ListRateHetSites[1] + ' - ' + ListRateHetSites[2] + '\n')
                Error.close()
            print("Rate of heterogeneity across sites values sampled from a normal distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2])
        elif len(ListRateHetSites) < 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less rate of heterogeneity across sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check RateHetSites value\n')
                Error.close()
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            if len(ListRateHetSites) < 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less rate of heterogeneity across sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check RateHetSites value\n')
                    Error.close()
                print("Less rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            if len(ListRateHetSites) > 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More rate of heterogeneity across sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check RateHetSites value\n')
                    Error.close()
                print("More rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            else:
                if ListRateHetSites[3] == 't':
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Truncated distribution detected\n')
                        Error.close()
                    print("Truncated distribution detected")
                    if ListRateHetSites[4] > ListRateHetSites[5]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check RateHetSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    elif ListRateHetSites[1] == ListRateHetSites[2]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check RateHetSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    else:
                        pass
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('t value or less parameters were expected\n')
                        Error.write('ERROR!! Please check RateHetSites value\n')
                        Error.close()
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check RateHetSites value')

    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Rate of heterogeneity across sites parameter value error\n')
            Error.write('ERROR!! Please check RateHetSites value\n')
            Error.close()
        print("Rate of heterogeneity across sites parameter value error")
        sys.exit('ERROR!! Please check RateHetSites value')


    # PropInvSites
    if bin.Scripts.Variables.PropInvSites.split(' ')[0] == '':
        pass
    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'fix':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 2:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Proportion of invariable sites value is fix with a value of ' + ListPropInvSites[1] + '\n')
                Error.close()
            print("Proportion of invariable sites value is fix with a value of " + ListPropInvSites[1])
        elif len(ListPropInvSites) < 2:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less proportion of invariable sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check PropInvSites value\n')
                Error.close()
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('More proportion of invariable sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check PropInvSites value\n')
                Error.close()
            print("More proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')

    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'uniform':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Proportion of invariable sites values sampled from a uniform distribution of ' + ListPropInvSites[1] + ' - ' + ListPropInvSites[2] + '\n')
                Error.close()
            print("Proportion of invariable sites values sampled from a uniform distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2])
        elif len(ListPropInvSites) < 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less proportion of invariable sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check PropInvSites value\n')
                Error.close()
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('More proportion of invariable sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check PropInvSites value\n')
                Error.close()
            print("More proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        if ListPropInvSites[1] > ListPropInvSites[2]:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Second value expected to be higher than the previous one\n')
                Error.write('ERROR!! Please check PropInvSites value\n')
                Error.close()
            print("Second value expected to be higher than the previous one")
            sys.exit('ERROR!! Please check PropInvSites value')
        elif ListPropInvSites[1] == ListPropInvSites[2]:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Second value expected to be higher than the previous one\n')
                Error.write('ERROR!! Please check PropInvSites value\n')
                Error.close()
            print("Second value expected to be higher than the previous one")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Proportion of invariable sites sampled from a uniform distribution of ' + ListPropInvSites[1] + ' - ' + ListPropInvSites[2] + ' seem correct\n')
                Error.close()
            print("Proportion of invariable sites sampled from a uniform distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2] + " seem correct")

    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'norm':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Proportion of invariable sites values sampled from a normal distribution of ' + ListPropInvSites[1] + ' - ' + ListPropInvSites[2] + '\n')
                Error.close()
            print("Proportion of invariable sites values sampled from a normal distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2])
        elif len(ListPropInvSites) < 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less proportion of invariable sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check PropInvSites value\n')
                Error.close()
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            if len(ListPropInvSites) < 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less proportion of invariable sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check PropInvSites value\n')
                    Error.close()
                print("Less proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            if len(ListPropInvSites) > 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More proportion of invariable sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check PropInvSites value\n')
                    Error.close()
                print("More proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            else:
                if ListPropInvSites[3] == 't':
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Truncated distribution detected\n')
                        Error.close()
                    print("Truncated distribution detected")
                    if ListPropInvSites[4] > ListPropInvSites[5]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check PropInvSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    elif ListPropInvSites[1] == ListPropInvSites[2]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check PropInvSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    else:
                        pass
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('t value or less parameters were expected\n')
                        Error.write('ERROR!! Please check PropInvSites value\n')
                        Error.close()
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check PropInvSites value')

    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'exp':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 2:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Proportion of invariable sites values sampled from a normal distribution of ' + ListPropInvSites[1] + '\n')
                Error.close()
            print("Proportion of invariable sites values sampled from a normal distribution of " + ListPropInvSites[1])
        elif len(ListPropInvSites) < 2:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less proportion of invariable sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check PropInvSites value\n')
                Error.close()
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            if len(ListPropInvSites) < 5:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less proportion of invariable sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check PropInvSites value\n')
                    Error.close()
                print("Less proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            if len(ListPropInvSites) > 5:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More proportion of invariable sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check PropInvSites value\n')
                    Error.close()
                print("More proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            else:
                if ListPropInvSites[2] == 't':
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Truncated distribution detected\n')
                        Error.close()
                    print("Truncated distribution detected")
                    if ListPropInvSites[3] > ListPropInvSites[4]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check PropInvSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    else:
                        pass
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('t value or less parameters were expected\n')
                        Error.write('ERROR!! Please check PropInvSites value\n')
                        Error.close()
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check PropInvSites value')

    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'gamma':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Proportion of invariable sites values sampled from a normal distribution of ' + ListPropInvSites[1] + ' - ' + ListPropInvSites[2] + '\n')
                Error.close()
            print("Proportion of invariable sites values sampled from a normal distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2])
        elif len(ListPropInvSites) < 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less proportion of invariable sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check PropInvSites value\n')
                Error.close()
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            if len(ListPropInvSites) < 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less proportion of invariable sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check PropInvSites value\n')
                    Error.close()
                print("Less proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            if len(ListPropInvSites) > 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More proportion of invariable sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check PropInvSites value\n')
                    Error.close()
                print("More proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            else:
                if ListPropInvSites[3] == 't':
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Truncated distribution detected\n')
                        Error.close()
                    print("Truncated distribution detected")
                    if ListPropInvSites[4] > ListPropInvSites[5]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check PropInvSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    elif ListPropInvSites[1] == ListPropInvSites[2]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check PropInvSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    else:
                        pass
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('t value or less parameters were expected\n')
                        Error.write('ERROR!! Please check PropInvSites value\n')
                        Error.close()
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check PropInvSites value')

    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'beta':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Proportion of invariable sites values sampled from a normal distribution of ' + ListPropInvSites[1] + ' - ' + ListPropInvSites[2] + '\n')
                Error.close()
            print("Proportion of invariable sites values sampled from a normal distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2])
        elif len(ListPropInvSites) < 3:
            with open('ProteinModelerABC.out', 'a') as Error:
                Error.write('Less proportion of invariable sites parameters detected than was expected\n')
                Error.write('ERROR!! Please check PropInvSites value\n')
                Error.close()
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            if len(ListPropInvSites) < 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('Less proportion of invariable sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check PropInvSites value\n')
                    Error.close()
                print("Less proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            if len(ListPropInvSites) > 6:
                with open('ProteinModelerABC.out', 'a') as Error:
                    Error.write('More proportion of invariable sites parameters detected than was expected\n')
                    Error.write('ERROR!! Please check PropInvSites value\n')
                    Error.close()
                print("More proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            else:
                if ListPropInvSites[3] == 't':
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('Truncated distribution detected\n')
                        Error.close()
                    print("Truncated distribution detected")
                    if ListPropInvSites[4] > ListPropInvSites[5]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check PropInvSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    elif ListPropInvSites[1] == ListPropInvSites[2]:
                        with open('ProteinModelerABC.out', 'a') as Error:
                            Error.write('Fifth value expected to be higher than the previous one\n')
                            Error.write('ERROR!! Please check PropInvSites value\n')
                            Error.close()
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    else:
                        pass
                else:
                    with open('ProteinModelerABC.out', 'a') as Error:
                        Error.write('t value or less parameters were expected\n')
                        Error.write('ERROR!! Please check PropInvSites value\n')
                        Error.close()
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check PropInvSites value')

    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Proportion of invariable sites parameter value error\n')
            Error.write('ERROR!! Please check PropInvSites value\n')
            Error.close()
        print("Proportion of invariable sites parameter value error")
        sys.exit('ERROR!! Please check PropInvSites value')

        # SCS model parameters
# """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # Thermodynamic temperature
    if isinstance(bin.Scripts.Variables.TEMP, float) == True:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Thermodynamic temperature set at ' + str(bin.Scripts.Variables.TEMP) + '\n')
            Error.close()
        print('Thermodynamic temperature set at ' + str(bin.Scripts.Variables.TEMP))
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Thermodynamic temperature must be an integer\n')
            Error.write('ERROR!! Please check TEMP value\n')
            Error.close()
        print('Thermodynamic temperature must be an integer')
        sys.exit('ERROR!! Please check TEMP value')

    # Configurational entropies per residue
    if isinstance(bin.Scripts.Variables.S0, float) == True:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Unfolded protein configurational entropy per residue set at ' + str(bin.Scripts.Variables.S0) + '\n')
            Error.close()
        print('Unfolded protein configurational entropy per residue set at ' + str(bin.Scripts.Variables.S0))
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Unfolded protein configurational entropy per residue must be a number\n')
            Error.write('ERROR!! Please check S0 value\n')
            Error.close()
        print('Unfolded protein configurational entropy per residue must be a number')
        sys.exit('ERROR!! Please check S0 value')

    if isinstance(bin.Scripts.Variables.SC1,float) == True:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Misfolded protein configurational entropy per residue set at ' + str(bin.Scripts.Variables.SC1) + '\n')
            Error.close()
        print('Misfolded protein configurational entropy per residue set at ' + str(bin.Scripts.Variables.SC1))
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Misfolded protein configurational entropy per residue must be a number\n')
            Error.write('ERROR!! Please check SC1 value\n')
            Error.close()
        print('Misfolded protein configurational entropy per residue must be a number')
        sys.exit('ERROR!! Please check SC1 value')

    if isinstance(bin.Scripts.Variables.SC0, float) == True:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Misfolded protein entropy offset set at ' + str(bin.Scripts.Variables.SC0) + '\n')
            Error.close()
        print('Misfolded protein entropy offset set at ' + str(bin.Scripts.Variables.SC0))
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Misfolded protein entropy offset must be a number\n')
            Error.write('ERROR!! Please check SC0 value\n')
            Error.close()
        print('Misfolded protein entropy offset must be a number')
        sys.exit('ERROR!! Please check SC0 value')

    if isinstance(bin.Scripts.Variables.NPOP, int) == True:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('The population size to simulate molecular evolution under Fitness site-dependent SCS model set at ' + str(bin.Scripts.Variables.NPOP) + '\n')
            Error.close()
        print('The population size to simulate molecular evolution under Fitness site-dependent SCS model set at ' + str(bin.Scripts.Variables.NPOP))
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('The population size to simulate molecular evolution under Fitness site-dependent SCS model must be an integer\n')
            Error.write('ERROR!! Please check NPOP value\n')
            Error.close()
        print('The population size to simulate molecular evolution under Fitness site-dependent SCS model must be an integer')
        sys.exit('ERROR!! Please check NPOP value')


    #######################################
    ## Settings for the estimation phase ##
    #######################################

    # ABCIterations
    if bin.Scripts.Variables.ABCIterations < bin.Scripts.Variables.NumberOfSimulations:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('ABC iterations value of ' + str(bin.Scripts.Variables.ABCIterations) + '\n')
            Error.close()
        print('ABC iterations value of ' + str(bin.Scripts.Variables.ABCIterations))
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('ABC iterations value must be lower than the number of simulations\n')
            Error.write('ERROR!! Please check ABCIterations value\n')
            Error.close()
        print('ABC iterations value must be lower than the number of simulations')
        sys.exit('ERROR!! Please check ABCIterations value')


    # ABCTolerance
    if bin.Scripts.Variables.ABCTolerance < bin.Scripts.Variables.NumberOfSimulations:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('ABC tolerance value of ' + str(bin.Scripts.Variables.ABCTolerance) + '\n')
            Error.close()
        print('ABC tolerance value of ' + str(bin.Scripts.Variables.ABCTolerance))
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('ABC tolerance value must be lower than the number of simulations\n')
            Error.write('ERROR!! Please check ABCTolerance value\n')
            Error.close()
        print('ABC tolerance value must be lower than the number of simulations')
        sys.exit('ERROR!! Please check ABCTolerance value')


    # ABCMethod
    if bin.Scripts.Variables.ABCMethod == 'rejection':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Rejection algorithm was selected\n')
            Error.close()
        print('Rejection algorithm was selected')
    elif bin.Scripts.Variables.ABCMethod == 'mnlogistic':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Mnlogistic algorithm was selected\n')
            Error.close()
        print('Mnlogistic algorithm was selected')
    elif bin.Scripts.Variables.ABCMethod == 'neuralnet':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Neuralnet algorithm was selected\n')
            Error.close()
        print('Neuralnet algorithm was selected')
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('ABC algorithm was not recognised\n')
            Error.write('ERROR!! Please check ABCMethod value\n')
            Error.close()
        print('ABC algorithm was not recognised')
        sys.exit('ERROR!! Please check ABCMethod value')


    # SummaryStatistics
    ListSummaryStatistics = bin.Scripts.Variables.SummaryStatistics.split(' ')
    if len(ListSummaryStatistics) < 1:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Less summary statistics detected than expected\n')
            Error.write('ERROR!! Please check SummaryStatistics value\n')
            Error.close()
        print('Less summary statistics detected than expected')
        sys.exit('ERROR!! Please check SummaryStatistics value')
    elif len(ListSummaryStatistics) > 7:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('More summary statistics detected than expected\n')
            Error.write('ERROR!! Please check SummaryStatistics value\n')
            Error.close()
        print('More summary statistics detected than expected')
        sys.exit('ERROR!! Please check SummaryStatistics value')
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write(str(len(ListSummaryStatistics))  + ' summary statistics detected\n')
            Error.close()
        print(str(len(ListSummaryStatistics))  + ' summary statistics detected')


    # MultiPage
    if bin.Scripts.Variables.MultiPage == 'No':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Plots will be displayed in single-pages files was selected\n')
            Error.close()
        print('Plots will be displayed in single-pages files was selected')
    elif bin.Scripts.Variables.MultiPage == 'Yes':
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Plots will be displayed in multiple-pages files was selected\n')
            Error.close()
        print('Plots will be displayed in multiple-pages files was selected')
    else:
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('Plots displayed option was not recognised\n')
            Error.write('ERROR!! Please check MultiPage value\n')
            Error.close()
        print('Plots displayed option was not recognised')
        sys.exit('ERROR!! Please check MultiPage value')

