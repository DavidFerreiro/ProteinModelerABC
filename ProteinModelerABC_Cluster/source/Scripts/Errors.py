import os
import sys
import bin.Scripts.Variables
from Bio import SeqIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import platform

def err():
    # System detection
    sistema = platform.system()
    if sistema == 'Windows':
        sys.exit('ERROR!!! This program is not executable by Windows\n')
    else:
        print('Working in a ' + str(sistema) + ' machine')


    ################################
    ## Input Data and Information ##
    ################################

    # MPI file exists
    if bin.Scripts.Variables.MPI == '':
        print('Local version selected')
    else:
        if os.path.isfile(os.getcwd() + '/' + bin.Scripts.Variables.MPI):
            print('Cluster version selected')
        else:
            print('MPI file does not exists')
            sys.exit('ERROR!!! Please check the MPI file')

    # NameOfPhylipFile exists
    if os.path.isfile(os.getcwd() + '/' + bin.Scripts.Variables.NameOfPhylipFile):
        print('Alignment file exists')
    else:
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
                    print('Sequence name of line must be 10 characters including spaces before sequences starts')
                    print('Check if alignment file end with and empty line')
                    sys.exit('ERROR!!! Sequence name of line must be 10 characters including spaces before sequences starts')

    except OSError:
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
                        print('Error in alignment file: ' + lines[x][y] + ' character was found in line ' + str(x+1) + ', position ' + str(y))
                        sys.exit('Error in alignment file: ' + lines[x][y] + ' character is not allow')


    # Indels
    if bin.Scripts.Variables.Indels == 'Ignored':
        print("Indels are ignored")
    elif bin.Scripts.Variables.Indels == 'New State':
        print("Indels are considered as a new state")
    else:
        print("Indels value are incorrect")
        sys.exit('ERROR!! Please check the indels value')


    # Template
    if bin.Scripts.Variables.Template != '':
        if bin.Scripts.Variables.Template[-4:] == '.pdb':
            print(bin.Scripts.Variables.Template + ' selected as template')
        else:
            print('You must select a .pdb file as template')
            sys.exit('ERROR!! Please check Template value')
    else:
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
                    print('Chain ' + str(bin.Scripts.Variables.ChainPDB) + ' selected')
                else:
                    print('Chain selected not found in the pdb file')
                    sys.exit('ERROR!! Please check Chain value')
            else:
                print('Chain must be a unique letter')
                sys.exit('ERROR!! Please check Chain value')
        else:
            print('Chain must be describe with letters from A to Z')
            sys.exit('ERROR!! Please check Chain value')
    else:
        print('Chain not selected')
        sys.exit('ERROR!! Please check Chain value')


    #######################################
    ## Settings for the simulation phase ##
    #######################################

    # NumberOfSimulations
    if bin.Scripts.Variables.NumberOfSimulations < 0:
        print("Your number of simulations (" + str(bin.Scripts.Variables.NumberOfSimulations) + ") are not correct")
        sys.exit('ERROR!!! Number of simulations must be to 1-100000000')
    elif bin.Scripts.Variables.NumberOfSimulations > 100000000:
        print("Your number of simulations (" + str(bin.Scripts.Variables.NumberOfSimulations) + ") are not correct")
        sys.exit('ERROR!!! Number of simulations must be to 1-100000000')
    else:
        print("Correct number of simulations")

    # SaveSimulations
    if bin.Scripts.Variables.SaveSimulations == 'No':
        print("Simulated data is not saved")
    elif bin.Scripts.Variables.SaveSimulations == 'Yes':
        print("Simulated data is saved")
    else:
        print("Simulated data value are incorrect")
        sys.exit('ERROR!! Please check simulated data value')


    # ShowInformationScreen
    if bin.Scripts.Variables.ShowInformationScreen == 'No':
        print("Running information will not be displayed on the screen")
    elif bin.Scripts.Variables.ShowInformationScreen == 'Yes':
        print("Running information will be displayed on the screen")
    else:
        print("Running information value are incorrect")
        sys.exit('ERROR!! Please check running information value')

    #############################    EVOLUTIONARY HISTORY    #############################
    #¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#

    # Coalescent or Tree
    if bin.Scripts.Variables.Coalescent != 'Coal' and bin.Scripts.Variables.Coalescent != 'Phylo':
        print("User has to perform coalescent or provide a phylogenetic tree input file")
        sys.exit('ERROR!! Please check CoalescentOrPhylogeny value')


    elif bin.Scripts.Variables.Coalescent == 'Phylo':
        print('User-specified phylogenetic tree selected')

        try:
            with open(bin.Scripts.Variables.Tree) as tree:
                tree.close()

        except OSError:
            print('User-specified phylogenetic tree does not exist')
            sys.exit('ERROR!!! Cannot read ' + str(bin.Scripts.Variables.Tree) + '. Check phylogenetic tree input file')


    elif bin.Scripts.Variables.Coalescent == 'Coal':

        # HaploidDiploid
        if bin.Scripts.Variables.HaploidDiploid == 1:
            print("Haploid data are selected")
        elif bin.Scripts.Variables.HaploidDiploid == 2:
            print("Diploid data are selected")
        else:
            print("Haploid/Diploid information value are incorrect")
            sys.exit('ERROR!! Please check Haploid/Diploid value')


        # PopulationSize
        if bin.Scripts.Variables.PopulationSize > 0:
            print("PopulationSize selected")
        else:
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
                        print("The highest value of the " + str(i) + " interval seems incorrect")
                        sys.exit('ERROR!! Please check DatedTips value')
                for e in range(1, (int(NumOfDatedTips) - 1)):
                    if int(DatedTipsList[(e * 3 - 1)]) < int(DatedTipsList[(e * 3 + 2)]):
                        pass
                    else:
                        print("The lowest value of the " + str(e) + " interval seems incorrect")
                        sys.exit('ERROR!! Please check DatedTips value')
            else:
                print("Number of elements of Dated Tips are incorrect")
                sys.exit('ERROR!! Please check DatedTips value')


        # GenerationTime
        if bin.Scripts.Variables.GenerationTime == '':
            pass
        else:
            if bin.Scripts.Variables.GenerationTime.split(' ')[0] == 'fix':
                ListGenerationTime = list(bin.Scripts.Variables.GenerationTime.split(" "))
                if len(ListGenerationTime) == 2:
                    print("Parameter value is fix with a value of " + ListGenerationTime[1])
                elif len(ListGenerationTime) < 2:
                    print("Less generation time parameters than was expected")
                    sys.exit('ERROR!! Please check GenerationTime value')
                else:
                    print("More generation time parameterss than was expected")
                    sys.exit('ERROR!! Please check GenerationTime value')
            elif bin.Scripts.Variables.GenerationTime.split(' ')[0] == 'uniform':
                ListGenerationTime = list(bin.Scripts.Variables.GenerationTime.split(" "))
                if len(ListGenerationTime) == 3:
                    print("Parameter values sampled from a uniform distribution of " + ListGenerationTime[1] + "-" + ListGenerationTime[2])
                elif len(ListGenerationTime) < 3:
                    print("Less generation time parameters than was expected")
                    sys.exit('ERROR!! Please check GenerationTime value')
                else:
                    print("More generation time parameters than was expected")
                    sys.exit('ERROR!! Please check GenerationTime value')
            else:
                print("Generation time parameter value error")
                sys.exit('ERROR!! Please check GenerationTime value')


        # GrowthRate
        if bin.Scripts.Variables.GrowthRate.split(' ')[0] == 0:
            print("Exponential growth rate selected")
            if len(bin.Scripts.Variables.GrowthRate.split(' ')) == 2:
                print("Exponential growth rate seems correct")
            elif len(bin.Scripts.Variables.GrowthRate.split(' ')) < 2:
                print("Less values than expected in exponential growth rate")
                sys.exit('ERROR!! Please check GrowthRate value')
            elif len(bin.Scripts.Variables.GrowthRate.split(' ')) > 2:
                print("More values than expected in exponential growth rate")
                sys.exit('ERROR!! Please check GrowthRate value')
        elif bin.Scripts.Variables.GrowthRate.split(' ')[0] == 1:
            print("Demographic periods selected")
            if len(bin.Scripts.Variables.GrowthRate.split(' ')) == (3 * int(bin.Scripts.Variables.GrowthRate.split(' ')[1]) + 2):
                print("Demographic periods seems correct")
            elif len(bin.Scripts.Variables.GrowthRate.split(' ')) > (3 * int(bin.Scripts.Variables.GrowthRate.split(' ')[1]) + 2):
                print("Less values than expected in demographic periods")
                sys.exit('ERROR!! Please check GrowthRate value')
            elif len(bin.Scripts.Variables.GrowthRate.split(' ')) < (3 * int(bin.Scripts.Variables.GrowthRate.split(' ')[1]) + 2):
                print("More values than expected in demographic periods")
                sys.exit('ERROR!! Please check GrowthRate value')
        else:
            pass


        # MigrationModel
        if bin.Scripts.Variables.MigrationModel.split(' ')[0] == 1:
            print("Island model selected")
            if len(bin.Scripts.Variables.MigrationModel.split(' ')) == (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                print("Island model seems correct")
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) < (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                print("Less values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) > (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                print("More values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
        elif bin.Scripts.Variables.MigrationModel.split(' ')[0] == 2:
            print("Stepping-stone model selected")
            if len(bin.Scripts.Variables.MigrationModel.split(' ')) == (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                print("Stepping-stone seems correct")
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) < (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                print("Less values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) > (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                print("More values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
        elif bin.Scripts.Variables.MigrationModel.split(' ')[0] == 3:
            print("Continent-island model selected")
            if len(bin.Scripts.Variables.MigrationModel.split(' ')) == (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                print("Continent-island seems correct")
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) < (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                print("Less values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
            elif len(bin.Scripts.Variables.MigrationModel.split(' ')) > (2 + int(bin.Scripts.Variables.MigrationModel.split(' ')[1])):
                print("More values than expected in migration model")
                sys.exit('ERROR!! Please check MigrationModel value')
        else:
            pass


        # MigrationRate
        if bin.Scripts.Variables.MigrationRate == '':
            pass
        else:
            if len(bin.Scripts.Variables.MigrationRate.split(' ')) == (2 * int(bin.Scripts.Variables.MigrationRate.split(' ')[1])):
                print("Migration rate seems correct")
            elif len(bin.Scripts.Variables.MigrationRate.split(' ')) < (2 * int(bin.Scripts.Variables.MigrationRate.split(' ')[1])):
                print("Less values than expected in migration rate")
                sys.exit('ERROR!! Please check MigrationRate value')
            elif len(bin.Scripts.Variables.MigrationRate.split(' ')) > (2 * int(bin.Scripts.Variables.MigrationRate.split(' ')[1])):
                print("More values than expected in migration rate")
                sys.exit('ERROR!! Please check MigrationRate value')


        # ConvergenceDemes
        if bin.Scripts.Variables.ConvergenceDemes == '':
            pass
        else:
            if len(bin.Scripts.Variables.ConvergenceDemes.split(' ')) == (3 * int(bin.Scripts.Variables.ConvergenceDemes(' ')[1]) + 1):
                print("Convergence of demes seems correct")
            elif len(bin.Scripts.Variables.ConvergenceDemes.split(' ')) < (3 * int(bin.Scripts.Variables.ConvergenceDemes.split(' ')[1]) + 1):
                print("Less values than expected in convergence demes")
                sys.exit('ERROR!! Please check ConvergenceDeme value')
            elif len(bin.Scripts.Variables.ConvergenceDemes.split(' ')) > (3 * int(bin.Scripts.Variables.ConvergenceDemes.split(' ')[1]) + 1):
                print("More values than expected in convergence demes")
                sys.exit('ERROR!! Please check ConvergenceDemes value')


        # SubstitutionRate
        if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'fix':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 2:
                print("Parameter value is fix with a value of " + ListSubstitutionRate[1])
            elif len(ListSubstitutionRate) < 2:
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                print("More substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')

        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'uniform':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 3:
                print("Parameter values sampled from a uniform distribution of " + ListSubstitutionRate[1] + " - " + ListSubstitutionRate[2])
            elif len(ListSubstitutionRate) < 3:
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                print("More substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            if ListSubstitutionRate[1] > ListSubstitutionRate[2]:
                print("Second value expected to be higher than the previous one")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            elif ListSubstitutionRate[1] == ListSubstitutionRate[2]:
                print("Second value expected to be higher than the previous one")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                print("Parameter values sampled from a uniform distribution of " + ListSubstitutionRate[1] + " - " + ListSubstitutionRate[2] + " seem correct")

        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'norm':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 3:
                print("Parameter values sampled from a normal distribution of " + ListSubstitutionRate[1] + " - " + ListSubstitutionRate[2])
            elif len(ListSubstitutionRate) < 3:
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                if len(ListSubstitutionRate) < 6:
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                if len(ListSubstitutionRate) > 6:
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                else:
                    if ListSubstitutionRate[3] == 't':
                        print("Truncated distribution detected")
                        if ListSubstitutionRate[4] > ListSubstitutionRate[5]:
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        elif ListSubstitutionRate[4] == ListSubstitutionRate[5]:
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        else:
                            pass
                    else:
                        print("t value or less parameters were expected")
                        sys.exit('ERROR!! Please check SubstitutionRate value')

        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'exp':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 2:
                print("Parameter values sampled from a normal distribution of " + ListSubstitutionRate[1])
            elif len(ListSubstitutionRate) < 2:
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                if len(ListSubstitutionRate) < 5:
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                if len(ListSubstitutionRate) > 5:
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                else:
                    if ListSubstitutionRate[2] == 't':
                        print("Truncated distribution detected")
                        if ListSubstitutionRate[3] > ListSubstitutionRate[4]:
                            print("Fourth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        elif ListSubstitutionRate[3] == ListSubstitutionRate[4]:
                            print("Fourth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        else:
                            pass
                    else:
                        print("t value or less parameters were expected")
                        sys.exit('ERROR!! Please check SubstitutionRate value')

        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 3:
                print("Parameter values sampled from a normal distribution of " + ListSubstitutionRate[1] + "-" + ListSubstitutionRate[2])
            elif len(ListSubstitutionRate) < 3:
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                if len(ListSubstitutionRate) < 6:
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                if len(ListSubstitutionRate) > 6:
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                else:
                    if ListSubstitutionRate[3] == 't':
                        print("Truncated distribution detected")
                        if ListSubstitutionRate[4] > ListSubstitutionRate[5]:
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        elif ListSubstitutionRate[4] == ListSubstitutionRate[5]:
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        else:
                            pass
                    else:
                        print("t value or less parameters were expected")
                        sys.exit('ERROR!! Please check SubstitutionRate value')

        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'beta':
            ListSubstitutionRate = list(bin.Scripts.Variables.SubstitutionRate.split(" "))
            if len(ListSubstitutionRate) == 3:
                print("Parameter values sampled from a normal distribution of " + ListSubstitutionRate[1] + " - " + ListSubstitutionRate[2])
            elif len(ListSubstitutionRate) < 3:
                print("Less substitution rate parameters than was expected")
                sys.exit('ERROR!! Please check SubstitutionRate value')
            else:
                if len(ListSubstitutionRate) < 6:
                    print("Incorrect substitution rate parameters than was expected")
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                if len(ListSubstitutionRate) > 6:
                    print("Incorrect substitution rate parameters than was expected")    
                    sys.exit('ERROR!! Please check SubstitutionRate value')
                else:
                    if ListSubstitutionRate[3] == 't':
                        print("Truncated distribution detected")
                        if ListSubstitutionRate[4] > ListSubstitutionRate[5]:
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        elif ListSubstitutionRate[4] == ListSubstitutionRate[5]:
                            print("Fifth value expected to be higher than the previous one")
                            sys.exit('ERROR!! Please check SubstitutionRate value')
                        else:
                            pass
                    else:
                        print("t value or less parameters were expected")
                        sys.exit('ERROR!! Please check SubstitutionRate value')

        else:
            print("Substitution rate parameter value error")
            sys.exit('ERROR!! Please check SubstitutionRate value')



    ##############################    SUBSTITUTION MODEL    ##############################
    #¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#

    # SubstitutionModel
    for x in bin.Scripts.Variables.SubstitutionModel.split(' '):
        if x == 'Blosum62':
            print("Blosum62 substitution model selected")
        elif x == 'CpRev':
            print("CpRev substitution model selected")
        elif x == 'Dayhoff':
            print("Dayhoff substitution model selected")
        elif x == 'DayhoffDCMUT':
            print("DayhoffDCMUT substitution model selected")
        elif x == 'HIVb':
            print("HIVb substitution model selected")
        elif x == 'HIVw':
            print("HIVw substitution model selected")
        elif x == 'JTT':
            print("JTT substitution model selected")
        elif x == 'JonesDCMUT':
            print("JonesDCMUT substitution model selected")
        elif x == 'LG':
            print("LG substitution model selected")
        elif x == 'Mtart':
            print("Mtart substitution model selected")
        elif x == 'Mtmam':
            print("Mtmam substitution model selected")
        elif x == 'Mtrev24':
            print("Mtrev24 substitution model selected")
        elif x == 'RtRev':
            print("RtRev substitution model selected")
        elif x == 'VT':
            print("VT substitution model selected")
        elif x == 'WAG':
            print("WAG substitution model selected")
        elif x == 'UserEAAM':
            print("UserEAAM substitution model selected")
        elif x == 'Fitness':
            print("Fitness substitution model selected")
        elif x == 'Neutral':
            print("Neutral substitution model selected")
        else:
            print('Substitution model selected: ' + str(x) + ' was not recognised\n')
            sys.exit('ERROR!! Please check SubstitutionModel value')


    # AminoacidFrequencies
    if bin.Scripts.Variables.AminoacidFrequencies.split(' ')[0] == 'fix':
        ListAminoacidFrequencies = list(bin.Scripts.Variables.AminoacidFrequencies.split(' '))
        if len(ListAminoacidFrequencies) == 21:
            for AminoacidPosition in range(1, 20):
                if float(ListAminoacidFrequencies[AminoacidPosition]) > 1:
                    print('Aminoacid frequency values must be 0-1')
                    sys.exit('ERROR!! Please check AminoacidFrequencies values')
                elif float(ListAminoacidFrequencies[AminoacidPosition]) < 0:
                    print('Aminoacid frequency values must be 0-1')
                    sys.exit('ERROR!! Please check AminoacidFrequencies values')
                else:
                    pass
            print('Aminoacid frequency values seem correct')
        else:
            print("Aminoacid frequencies length seems incorrect")
            sys.exit('ERROR!! Please check AminoacidFrequencies values')
    elif bin.Scripts.Variables.AminoacidFrequencies.split(' ')[0] == 'dirichlet':
        ListAminoacidFrequencies = list(bin.Scripts.Variables.AminoacidFrequencies.split(' '))
        if len(ListAminoacidFrequencies) == 21:
            for AminoacidPosition in range(1, 20):
                if float(ListAminoacidFrequencies[AminoacidPosition]) > 1:
                    print('Aminoacid frequency values must be 0-1')
                    sys.exit('ERROR!! Please check AminoacidFrequencies values')
                elif float(ListAminoacidFrequencies[AminoacidPosition]) < 0:
                    print('Aminoacid frequency values must be 0-1')
                    sys.exit('ERROR!! Please check AminoacidFrequencies values')
                else:
                    print('Aminoacid frequency values seem correct')
        else:
            print("Aminoacid frequencies length seems incorrect")
            sys.exit('ERROR!! Please check AminoacidFrequencies values')
    else:
        print("Aminoacid frequencies parameter value error")
        sys.exit('ERROR!! Please check AminoacidFrequencies values')


    # RateHetSites
    if bin.Scripts.Variables.RateHetSites.split(' ')[0] == '':
        pass
    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'fix':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 2:
            print("Rate of heterogeneity across sites value is fix with a value of " + ListRateHetSites[1])
        elif len(ListRateHetSites) < 2:
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check SubstitutionRate value')
        else:
            print("More rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check SubstitutionRate value')

    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'uniform':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 3:
            print("Rate of heterogeneity across sites values sampled from a uniform distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2])
        elif len(ListRateHetSites) < 3:
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            print("More rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        if ListRateHetSites[1] > ListRateHetSites[2]:
            print("Second value expected to be higher than the previous one\n")
            sys.exit('ERROR!! Please check RateHetSites value')
        elif ListRateHetSites[1] == ListRateHetSites[2]:
            print("Second value expected to be higher than the previous one")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            print("Rate of heterogeneity across sites sampled from a uniform distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2] + " seem correct")

    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'norm':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 3:
            print("Rate of heterogeneity across sites values sampled from a normal distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2])
        elif len(ListRateHetSites) < 3:
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            if len(ListRateHetSites) < 6:
                print("Less rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            if len(ListRateHetSites) > 6:
                print("More rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            else:
                if ListRateHetSites[3] == 't':
                    print("Truncated distribution detected")
                    if ListRateHetSites[4] > ListRateHetSites[5]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    elif ListRateHetSites[1] == ListRateHetSites[2]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    else:
                        pass
                else:
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check RateHetSites value')

    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'exp':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 2:
            print("Rate of heterogeneity across sites values sampled from a normal distribution of " + ListRateHetSites[1])
        elif len(ListRateHetSites) < 2:
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            if len(ListRateHetSites) < 5:
                print("Less rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            if len(ListRateHetSites) > 5:
                print("More rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            else:
                if ListRateHetSites[2] == 't':
                    print("Truncated distribution detected")
                    if ListRateHetSites[3] > ListRateHetSites[4]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    else:
                        pass
                else:
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check RateHetSites value')

    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'gamma':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 3:
            print("Rate of heterogeneity across sites values sampled from a normal distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2])
        elif len(ListRateHetSites) < 3:
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            if len(ListRateHetSites) < 6:
                print("Less rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            if len(ListRateHetSites) > 6:
                print("More rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            else:
                if ListRateHetSites[3] == 't':
                    print("Truncated distribution detected")
                    if ListRateHetSites[4] > ListRateHetSites[5]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    elif ListRateHetSites[1] == ListRateHetSites[2]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    else:
                        pass
                else:
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check RateHetSites value')

    elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'beta':
        ListRateHetSites = list(bin.Scripts.Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 3:
            print("Rate of heterogeneity across sites values sampled from a normal distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2])
        elif len(ListRateHetSites) < 3:
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            if len(ListRateHetSites) < 6:
                print("Less rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            if len(ListRateHetSites) > 6:
                print("More rate of heterogeneity across sites parameters detected than was expected")
                sys.exit('ERROR!! Please check RateHetSites value')
            else:
                if ListRateHetSites[3] == 't':
                    print("Truncated distribution detected")
                    if ListRateHetSites[4] > ListRateHetSites[5]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    elif ListRateHetSites[1] == ListRateHetSites[2]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check RateHetSites value')
                    else:
                        pass
                else:
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check RateHetSites value')

    else:
        print("Rate of heterogeneity across sites parameter value error")
        sys.exit('ERROR!! Please check RateHetSites value')


    # PropInvSites
    if bin.Scripts.Variables.PropInvSites.split(' ')[0] == '':
        pass
    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'fix':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 2:
            print("Proportion of invariable sites value is fix with a value of " + ListPropInvSites[1])
        elif len(ListPropInvSites) < 2:
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            print("More proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')

    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'uniform':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 3:
            print("Proportion of invariable sites values sampled from a uniform distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2])
        elif len(ListPropInvSites) < 3:
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            print("More proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        if ListPropInvSites[1] > ListPropInvSites[2]:
            print("Second value expected to be higher than the previous one")
            sys.exit('ERROR!! Please check PropInvSites value')
        elif ListPropInvSites[1] == ListPropInvSites[2]:
            print("Second value expected to be higher than the previous one")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            print("Proportion of invariable sites sampled from a uniform distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2] + " seem correct")

    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'norm':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 3:
            print("Proportion of invariable sites values sampled from a normal distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2])
        elif len(ListPropInvSites) < 3:
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            if len(ListPropInvSites) < 6:
                print("Less proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            if len(ListPropInvSites) > 6:
                print("More proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            else:
                if ListPropInvSites[3] == 't':
                    print("Truncated distribution detected")
                    if ListPropInvSites[4] > ListPropInvSites[5]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    elif ListPropInvSites[1] == ListPropInvSites[2]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    else:
                        pass
                else:
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check PropInvSites value')

    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'exp':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 2:
            print("Proportion of invariable sites values sampled from a normal distribution of " + ListPropInvSites[1])
        elif len(ListPropInvSites) < 2:
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            if len(ListPropInvSites) < 5:
                print("Less proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            if len(ListPropInvSites) > 5:
                print("More proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            else:
                if ListPropInvSites[2] == 't':
                    print("Truncated distribution detected")
                    if ListPropInvSites[3] > ListPropInvSites[4]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    else:
                        pass
                else:
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check PropInvSites value')

    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'gamma':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 3:
            print("Proportion of invariable sites values sampled from a normal distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2])
        elif len(ListPropInvSites) < 3:
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            if len(ListPropInvSites) < 6:
                print("Less proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            if len(ListPropInvSites) > 6:
                print("More proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            else:
                if ListPropInvSites[3] == 't':
                    print("Truncated distribution detected")
                    if ListPropInvSites[4] > ListPropInvSites[5]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    elif ListPropInvSites[1] == ListPropInvSites[2]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    else:
                        pass
                else:
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check PropInvSites value')

    elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'beta':
        ListPropInvSites = list(bin.Scripts.Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 3:
            print("Proportion of invariable sites values sampled from a normal distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2])
        elif len(ListPropInvSites) < 3:
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check PropInvSites value')
        else:
            if len(ListPropInvSites) < 6:
                print("Less proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            if len(ListPropInvSites) > 6:
                print("More proportion of invariable sites parameters detected than was expected")
                sys.exit('ERROR!! Please check PropInvSites value')
            else:
                if ListPropInvSites[3] == 't':
                    print("Truncated distribution detected")
                    if ListPropInvSites[4] > ListPropInvSites[5]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    elif ListPropInvSites[1] == ListPropInvSites[2]:
                        print("Fifth value expected to be higher than the previous one")
                        sys.exit('ERROR!! Please check PropInvSites value')
                    else:
                        pass
                else:
                    print("t value or less parameters were expected")
                    sys.exit('ERROR!! Please check PropInvSites value')

    else:
        print("Proportion of invariable sites parameter value error")
        sys.exit('ERROR!! Please check PropInvSites value')

        # SCS model parameters
        # """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        # Thermodynamic temperature
    if isinstance(bin.Scripts.Variables.TEMP, float) == True:
        print('Thermodynamic temperature set at ' + str(bin.Scripts.Variables.TEMP))
    else:
        print('Thermodynamic temperature must be an integer')
        sys.exit('ERROR!! Please check TEMP value')

    # Configurational entropies per residue
    if isinstance(bin.Scripts.Variables.S0, float) == True:
        print('Unfolded protein configurational entropy per residue set at ' + str(bin.Scripts.Variables.S0))
    else:
        print('Unfolded protein configurational entropy per residue must be a number')
        sys.exit('ERROR!! Please check S0 value')

    if isinstance(bin.Scripts.Variables.SC1, float) == True:
        print('Misfolded protein configurational entropy per residue set at ' + str(bin.Scripts.Variables.SC1))
    else:
        print('Misfolded protein configurational entropy per residue must be a number')
        sys.exit('ERROR!! Please check SC1 value')

    if isinstance(bin.Scripts.Variables.SC0, float) == True:
        print('Misfolded protein entropy offset set at ' + str(bin.Scripts.Variables.SC0))
    else:
        print('Misfolded protein entropy offset must be a number')
        sys.exit('ERROR!! Please check SC0 value')

    if isinstance(bin.Scripts.Variables.NPOP, int) == True:
        print('The population size to simulate molecular evolution under Fitness site-dependent SCS model set at ' + str(bin.Scripts.Variables.NPOP))
    else:
        print('The population size to simulate molecular evolution under Fitness site-dependent SCS model must be an integer')
        sys.exit('ERROR!! Please check NPOP value')


    #######################################
    ## Settings for the estimation phase ##
    #######################################

    # ABCIterations
    if bin.Scripts.Variables.ABCIterations < bin.Scripts.Variables.NumberOfSimulations:
        print('ABC iterations value of ' + str(bin.Scripts.Variables.ABCIterations))
    else:
        print('ABC iterations value must be lower than the number of simulations')
        sys.exit('ERROR!! Please check ABCIterations value')


    # ABCTolerance
    if bin.Scripts.Variables.ABCTolerance < bin.Scripts.Variables.NumberOfSimulations:
        print('ABC tolerance value of ' + str(bin.Scripts.Variables.ABCTolerance))
    else:
        print('ABC tolerance value must be lower than the number of simulations')
        sys.exit('ERROR!! Please check ABCTolerance value')


    # ABCMethod
    if bin.Scripts.Variables.ABCMethod == 'rejection':
        print('Rejection algorithm was selected')
    elif bin.Scripts.Variables.ABCMethod == 'mnlogistic':
        print('Mnlogistic algorithm was selected')
    elif bin.Scripts.Variables.ABCMethod == 'neuralnet':
        print('Neuralnet algorithm was selected')
    else:
        print('ABC algorithm was not recognised')
        sys.exit('ERROR!! Please check ABCMethod value')


    # SummaryStatistics
    ListSummaryStatistics = bin.Scripts.Variables.SummaryStatistics.split(' ')
    if len(ListSummaryStatistics) < 1:
        print('Less summary statistics detected than expected')
        sys.exit('ERROR!! Please check SummaryStatistics value')
    elif len(ListSummaryStatistics) > 7:
        print('More summary statistics detected than expected')
        sys.exit('ERROR!! Please check SummaryStatistics value')
    else:
        print(str(len(ListSummaryStatistics))  + ' summary statistics detected')


    # MultiPage
    if bin.Scripts.Variables.MultiPage == 'No':
        print('Plots will be displayed in single-pages files was selected')
    elif bin.Scripts.Variables.MultiPage == 'Yes':
        print('Plots will be displayed in multiple-pages files was selected')
    else:
        print('Plots displayed option was not recognised')
        sys.exit('ERROR!! Please check MultiPage value')

