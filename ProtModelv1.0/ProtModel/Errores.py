import os
import sys
import Variables
from Bio import SeqIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

def err():
    # Mandatory variables entry data validation
        #NameOfPhylipFile exists
    if os.path.isfile(os.getcwd() + '/' + Variables.NameOfPhylipFile):
        print("Alignment file exists")
    else:
        print("Alignment file does not exist")
        sys.exit('ERROR!!! Please check the alignment file')

        #NameOfPhylipFile format
    try:
        with open(Variables.NameOfPhylipFile) as Filo:
            for line in Filo:
                LongSeq = len(line)
                if LongSeq <= (Variables.NumAA + 11):
                    pass
                else:
                    print('Sequence name of line must be 10 characters including spaces before sequences starts')
                    print('Check if alignment file end with and empty line')
                    sys.exit('ERROR!!! Sequence name of line must be 10 characters including spaces before sequences starts')

    except OSError:
        print('Sequence name of line must be 10 characters including spaces before sequences starts')
        sys.exit('ERROR!!! Sequence name of line must be 10 characters including spaces before sequences starts')


        #NameOfPhylipFile charaacters
        # list to store file lines\n
    lines = []

    with open(Variables.NameOfPhylipFile, "r") as fp:
    # read an store all lines into list
        lines = fp.readlines()
        #if Indels == 'Default':
        for x in range (0,  Variables.NumSeq):
            if len(lines[x]) < 11:
                pass
            else:
                for y in range(10, Variables.NumAA + 8):
                    if lines[x][y] == 'A' or lines[x][y] == 'R' or lines[x][y] == 'N' or lines[x][y] == 'D' or lines[x][y] == 'C' or lines[x][y] == 'Q' or lines[x][y] == 'E' or lines[x][y] == 'G' or lines[x][y] == 'H' or lines[x][y] == 'I' or lines[x][y] == 'L' or lines[x][y] == 'K' or lines[x][y] == 'M' or lines[x][y] == 'F' or lines[x][y] == 'P' or lines[x][y] == 'S' or lines[x][y] == 'T' or lines[x][y] == 'W' or lines[x][y] == 'Y' or lines[x][y] == 'V' or lines[x][y] == '?' or lines[x][y] == '-' or lines[x][y] == 'X':
                        pass
                    else:
                        print('Error in alignment file: ' + lines[x][y] + ' character was found in line ' + str(x) + ', position ' + str(y))
                        sys.exit('Error in alignment file: ' + lines[x][y] + ' character is not allow')

        #NumberOfSimulations
    if Variables.NumberOfSimulations < 0:
        print("Your number of simulations (" + str(Variables.NumberOfSimulations) + ") are not correct")
        sys.exit('ERROR!!! Number of simulations must be to 1-100000')
    elif Variables.NumberOfSimulations > 100000:
        print("Your number of simulations (" + str(Variables.NumberOfSimulations) + ") are not correct")
        sys.exit('ERROR!!! Number of simulations must be to 1-100000')
    else:
        print("Correct number of simulations")

        #Indels
    if Variables.Indels == 'Ignored':
        print("Indels are ignored")
    elif Variables.Indels == 'New State':
        print("Indels are considered as a new state")
    else:
        print("Indels value are incorrect")
        sys.exit('ERROR!! Please check the indels value')

        #SaveSimulations
    if Variables.SaveSimulations == 'No':
        print("Simulated data is not saved")
    elif Variables.SaveSimulations == 'Yes':
        print("Simulated data is saved")
    else:
        print("Simulated data value are incorrect")
        sys.exit('ERROR!! Please check simulated data value')

        #ShowInformationScreen
    if Variables.ShowInformationScreen == 'No':
        print("Running information will not be displayed on the screen")
    elif Variables.ShowInformationScreen == 'Yes':
        print("Running information will be displayed on the screen")
    else:
        print("Running information value are incorrect")
        sys.exit('ERROR!! Please check running information value')

        #HaploidDiploid
    if Variables.HaploidDiploid == 1:
        print("Haploid data are selected")
    elif Variables.HaploidDiploid == 2:
        print("Diploid data are selected")
    else:
        print("Haploid/Diploid information value are incorrect")
        sys.exit('ERROR!! Please check Haploid/Diploid value')

        #PopulationSize
    if Variables.PopulationSize > 0:
        print("PopulationSize selected")

    else:
        print("PopulationSize information value are incorrect")
        sys.exit('ERROR!! Please check PopulationSize value')

        #SubstitutionRate
    if Variables.SubstitutionRate.split(' ')[0] == 'fix':
        ListSubstitutionRate = list(Variables.SubstitutionRate.split(" "))
        if len(ListSubstitutionRate) == 2:
            print("Parameter value is fix with a value of " + ListSubstitutionRate[1])
        elif len(ListSubstitutionRate) < 2:
            print("Less substitution rate parameters than was expected")
            sys.exit('ERROR!! Please check SubstitutionRate value')
        else:
            print("More substitution rate parameters than was expected")
            sys.exit('ERROR!! Please check SubstitutionRate value')

    elif Variables.SubstitutionRate.split(' ')[0] == 'uniform':
        ListSubstitutionRate = list(Variables.SubstitutionRate.split(" "))
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

    elif Variables.SubstitutionRate.split(' ')[0] == 'norm':
        ListSubstitutionRate = list(Variables.SubstitutionRate.split(" "))
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

    elif Variables.SubstitutionRate.split(' ')[0] == 'exp':
        ListSubstitutionRate = list(Variables.SubstitutionRate.split(" "))
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

    elif Variables.SubstitutionRate.split(' ')[0] == 'gamma':
        ListSubstitutionRate = list(Variables.SubstitutionRate.split(" "))
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

    elif Variables.SubstitutionRate.split(' ')[0] == 'beta':
        ListSubstitutionRate = list(Variables.SubstitutionRate.split(" "))
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

    else:
        print("Substitution rate parameter value error")
        sys.exit('ERROR!! Please check SubstitutionRate value')

        #SubstitutionModel
    for x in Variables.SubstitutionModel.split(' '):
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
        elif x == '':
            print("No one substitution model selected")
        else:
            print("Substitution model selected was not recognised")
            sys.exit('ERROR!! Please check SubstitutionModel value')
            
        #StructuralSubstitutionModel
    for x in Variables.StructuralSubstitutionModel.split(' '):
        if x == 'Fitness':
            print("Fitness substitution model selected")
        elif x == 'Neutral':
            print("Neutral substitution model selected")
        else:
            print("Substitution model selected: " + str(x) + " was not recognised")
            sys.exit('ERROR!! Please check StructuralSubstitutionModel value')

        #ABCIterations
    if Variables.ABCIterations == Variables.NumberOfSimulations or Variables.ABCIterations < Variables.NumberOfSimulations:
        print('ABC iterations value of ' + str(Variables.ABCIterations))
    else:
        print('ABC iterations value must be equal or lower than the number of simulations')
        sys.exit('ERROR!! Please check ABCIterations value')

        #ABCTolerance
    if Variables.ABCTolerance < Variables.NumberOfSimulations:
        print('ABC tolerance value of ' + str(Variables.ABCTolerance))
    else:
        print('ABC tolerance value must be lower than the number of simulations')
        sys.exit('ERROR!! Please check ABCTolerance value')

        #ABCMethod
    if Variables.ABCMethod == 'rejection':
        print('Rejection algorithm was selected')
    elif Variables.ABCMethod == 'mnlogistic':
        print('Mnlogistic algorithm was selected')
    elif Variables.ABCMethod == 'neuralnet':
        print('Neuralnet algorithm was selected')
    else:
        print('ABC algorithm was not recognised')
        sys.exit('ERROR!! Please check ABCMethod value')

        #SummaryStatistics
    ListSummaryStatistics = Variables.SummaryStatistics.split(' ')
    if len(ListSummaryStatistics) < 1:
        print('Less summary statistics detected than expected')
        sys.exit('ERROR!! Please check SummaryStatisticsvalue')
    elif len(ListSummaryStatistics) > 7:
        print('More summary statistics detected than expected')
        sys.exit('ERROR!! Please check SummaryStatisticsvalue')
    else:
        print(str(len(ListSummaryStatistics))  + ' summary statistics detected')

        #MultiPage
    if Variables.MultiPage == 'No':
        print('Plots will be displayed in single-pages files was selected')
    elif Variables.MultiPage == 'Yes':
        print('Plots will be displayed in multiple-pages files was selected')
    else:
        print('Plots displayed option was not recognised')
        sys.exit('ERROR!! Please check MultiPage value')



    # Optional variables

        #NumberOfProcessors
    #MacOSInformationCommand = 'system_profiler SPHardwareDataType'    #Command to check the computer number of cores
    #MacOSInformation2 = os.popen(MacOSInformationCommand).read()       #Read system information
    #NumProcesList = list(MacOSInformation2.split("Cores: "))   #Find Processors
    #NumProcesList1 = NumProcesList[1]   #Take rhe right part from cores word
    #NumProces = NumProcesList1[0]      #Take the first word = Number of Processors
    ComputerProcessors = os.cpu_count()

    if Variables.NumberOfProcessors == 'Default':
        NumberOfProcessors = ComputerProcessors
    else:
        if int(Variables.NumberOfProcessors) < 1:
            print("Number of processor must have a value of 1 at least")
            sys.exit('ERROR!! Please check NumberOfProcessors value')
        #elif int(NumberOfProcessors) == int(NumProces) or int(NumberOfProcessors) < int(NumProces):
        elif int(Variables.NumberOfProcessors) <= int(ComputerProcessors):
            print(str(Variables.NumberOfProcessors) + " processors are selected")
        else:
            print("Number of processors introduce is higher than yours computer number of processors")
            print("Please check your computer number of processors")
            sys.exit('ERROR!! Please check NumberOfProcessors value')

        #DatedTips
    if Variables.DatedTips == '':
        pass
    else:
        NumOfDatedTips = Variables.DatedTips[0]
        DatedTipsList = list(Variables.DatedTips.split(" "))
        NumOfTotalDatedTips = int(NumOfDatedTips) * 3 + 1
        if len(DatedTipsList) == int(NumOfTotalDatedTips):
            for i in range(1, int(NumOfDatedTips)):
                if int(DatedTipsList[(i * 3)]) < Variables.NumSeq:
                    pass
                else:
                    print("The highest value of the " + str(i) + " intervRal seems incorrect")
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

        #GenerationTime
    if Variables.GenerationTime == '':
        pass
    else:
        if Variables.GenerationTime.split(' ')[0] == 'fix':
            ListGenerationTime = list(Variables.GenerationTime.split(" "))
            if len(ListGenerationTime) == 2:
                print("Parameter value is fix with a value of " + ListGenerationTime[1])
            elif len(ListGenerationTime) < 2:
                print("Less generation time paraments than was expected")
                sys.exit('ERROR!! Please check GenerationTime value')
            else:
                print("More generation time paraments than was expected")
                sys.exit('ERROR!! Please check GenerationTime value')
        elif Variables.GenerationTime.split(' ')[0] == 'uniform':
            ListGenerationTime = list(Variables.GenerationTime.split(" "))
            if len(ListGenerationTime) == 3:
                print("Parameter values sampled from a uniform distribution of " + ListGenerationTime[1] + "-" + ListGenerationTime[2])
            elif len(ListGenerationTime) < 3:
                print("Less generation time paraments than was expected")
                sys.exit('ERROR!! Please check GenerationTime value')
            else:
                print("More generation time paraments than was expected")
                sys.exit('ERROR!! Please check GenerationTime value')
        else:
            print("Generation time parameter value error")
            sys.exit('ERROR!! Please check GenerationTime value')

        #AminoacidFrequencies
    if Variables.AminoacidFrequencies.split(' ')[0] == 'Default':
        pass
    elif Variables.AminoacidFrequencies.split(' ')[0] == 'fix':
        ListAminoacidFrequencies = list(Variables.AminoacidFrequencies.split(' '))
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
                    #print('Aminoacid frequency values seem correct')
            print('Aminoacid frequency values seem correct')
        else:
            print("Aminoacid frequencies length seems incorrect")
            sys.exit('ERROR!! Please check AminoacidFrequencies values')
    elif Variables.AminoacidFrequencies.split(' ')[0] == 'dirichlet':
        ListAminoacidFrequencies = list(Variables.AminoacidFrequencies.split(' '))
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

        #RateHetSites
    if Variables.RateHetSites.split(' ')[0] == '':
        pass
    elif Variables.RateHetSites.split(' ')[0] == 'fix':
        ListRateHetSites = list(Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 2:
            print("Rate of heterogeneity across sites value is fix with a value of " + ListRateHetSites[1])
        elif len(ListRateHetSites) < 2:
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check SubstitutionRate value')
        else:
            print("More rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check SubstitutionRate value')

    elif Variables.RateHetSites.split(' ')[0] == 'uniform':
        ListRateHetSites = list(Variables.RateHetSites.split(" "))
        if len(ListRateHetSites) == 3:
            print("Rate of heterogeneity across sites values sampled from a uniform distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2])
        elif len(ListRateHetSites) < 3:
            print("Less rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            print("More rate of heterogeneity across sites parameters detected than was expected")
            sys.exit('ERROR!! Please check RateHetSites value')
        if ListRateHetSites[1] > ListRateHetSites[2]:
            print("Second value expected to be higher than the previous one")
            sys.exit('ERROR!! Please check RateHetSites value')
        elif ListRateHetSites[1] == ListRateHetSites[2]:
            print("Second value expected to be higher than the previous one")
            sys.exit('ERROR!! Please check RateHetSites value')
        else:
            print("Rate of heterogeneity across sites sampled from a uniform distribution of " + ListRateHetSites[1] + "-" + ListRateHetSites[2] + " seem correct")

    elif Variables.RateHetSites.split(' ')[0] == 'norm':
        ListRateHetSites = list(Variables.RateHetSites.split(" "))
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

    elif Variables.RateHetSites.split(' ')[0] == 'exp':
        ListRateHetSites = list(Variables.RateHetSites.split(" "))
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

    elif Variables.RateHetSites.split(' ')[0] == 'gamma':
        ListRateHetSites = list(Variables.RateHetSites.split(" "))
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

    elif Variables.RateHetSites.split(' ')[0] == 'beta':
        ListRateHetSites = list(Variables.RateHetSites.split(" "))
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

        #PropInvSites
    if Variables.PropInvSites.split(' ')[0] == '':
        pass
    elif Variables.PropInvSites.split(' ')[0] == 'fix':
        ListPropInvSites = list(Variables.PropInvSites.split(" "))
        if len(ListPropInvSites) == 2:
            print("Proportion of invariable sites value is fix with a value of " + ListPropInvSites[1])
        elif len(ListPropInvSites) < 2:
            print("Less proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check SubstitutionRate value')
        else:
            print("More proportion of invariable sites parameters detected than was expected")
            sys.exit('ERROR!! Please check SubstitutionRate value')

    elif Variables.PropInvSites.split(' ')[0] == 'uniform':
        ListPropInvSites = list(Variables.PropInvSites.split(" "))
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
            print("proportion of invariable sites sampled from a uniform distribution of " + ListPropInvSites[1] + "-" + ListPropInvSites[2] + " seem correct")

    elif Variables.PropInvSites.split(' ')[0] == 'norm':
        ListPropInvSites = list(Variables.PropInvSites.split(" "))
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

    elif Variables.PropInvSites.split(' ')[0] == 'exp':
        ListPropInvSites = list(Variables.PropInvSites.split(" "))
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

    elif Variables.PropInvSites.split(' ')[0] == 'gamma':
        ListPropInvSites = list(Variables.PropInvSites.split(" "))
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

    elif Variables.PropInvSites.split(' ')[0] == 'beta':
        ListPropInvSites = list(Variables.PropInvSites.split(" "))
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

    if Variables.GrowthRate.split(' ')[0] == 0:
        print("Exponential growth rate selected")
        if len(Variables.GrowthRate.split(' ')) == 2:
            print("Exponential growth rate seems correct")
        elif len(Variables.GrowthRate.split(' ')) < 2:
            print("Less values than expected in exponential growth rate")
            sys.exit('ERROR!! Please check GrowthRate value')
        elif len(Variables.GrowthRate.split(' ')) > 2:
            print("More values than expected in exponential growth rate")
            sys.exit('ERROR!! Please check GrowthRate value')
    elif Variables.GrowthRate.split(' ')[0] == 1:
        print("Demographic periods selected")
        if len(Variables.GrowthRate.split(' ')) == (3 * int(Variables.GrowthRate.split(' ')[1]) + 2):
            print("Demographic periods seems correct")
        elif len(Variables.GrowthRate.split(' ')) > (3 * int(Variables.GrowthRate.split(' ')[1]) + 2):
            print("Less values than expected in demographic periods")
            sys.exit('ERROR!! Please check GrowthRate value')
        elif len(Variables.GrowthRate.split(' ')) < (3 * int(Variables.GrowthRate.split(' ')[1]) + 2):
            print("More values than expected in demographic periods")
            sys.exit('ERROR!! Please check GrowthRate value')
    else:
        pass

    if Variables.MigrationModel.split(' ')[0] == 1:
        print("Island model selected")
        if len(Variables.MigrationModel.split(' ')) == (2 + int(Variables.MigrationModel.split(' ')[1])):
            print("Island model seems correct")
        elif len(Variables.MigrationModel.split(' ')) > (2 + int(Variables.MigrationModel.split(' ')[1])):
            print("Less values than expected in migration model")
            sys.exit('ERROR!! Please check MigrationModel value')
        elif len(Variables.MigrationModel.split(' ')) < (2 + int(Variables.MigrationModel.split(' ')[1])):
            print("More values than expected in migration model")
            sys.exit('ERROR!! Please check MigrationModel value')
    elif Variables.MigrationModel.split(' ')[0] == 2:
        print("Stepping-stone model selected")
        if len(Variables.MigrationModel.split(' ')) == (2 + int(Variables.MigrationModel.split(' ')[1])):
            print("Stepping-stone seems correct")
        elif len(Variables.MigrationModel.split(' ')) > (2 + int(Variables.MigrationModel.split(' ')[1])):
            print("Less values than expected in migration model")
            sys.exit('ERROR!! Please check MigrationModel value')
        elif len(Variables.MigrationModel.split(' ')) < (2 + int(Variables.MigrationModel.split(' ')[1])):
            print("More values than expected in migration model")
            sys.exit('ERROR!! Please check MigrationModel value')
    elif Variables.MigrationModel.split(' ')[0] == 3:
        print("Continent-island model selected")
        if len(Variables.MigrationModel.split(' ')) == (2 + int(Variables.MigrationModel.split(' ')[1])):
            print("Continent-island seems correct")
        elif len(Variables.MigrationModel.split(' ')) > (2 + int(Variables.MigrationModel.split(' ')[1])):
            print("Less values than expected in migration model")
            sys.exit('ERROR!! Please check MigrationModel value')
        elif len(Variables.MigrationModel.split(' ')) < (2 + int(Variables.MigrationModel.split(' ')[1])):
            print("More values than expected in migration model")
            sys.exit('ERROR!! Please check MigrationModel value')
    else:
        pass

    if Variables.MigrationRate == '':
        pass
    else:
        if len(Variables.MigrationRate.split(' ')) == (2 * int(Variables.MigrationRate.split(' ')[1])):
            print("Migration rate seems correct")
        elif len(Variables.MigrationRate.split(' ')) == (2 * int(Variables.MigrationRate.split(' ')[1])):
            print("Less values than expected in migration rate")
            sys.exit('ERROR!! Please check MigrationRate value')
        elif len(Variables.MigrationRate.split(' ')) == (2 * int(Variables.MigrationRate.split(' ')[1])):
            print("Less values than expected in migration rate")
            sys.exit('ERROR!! Please check MigrationRate value')

    if Variables.ConvergenceDemes == '':
        pass
    else:
        if len(Variables.ConvergenceDeme.split(' ')) == (3 * int(Variables.ConvergenceDeme(' ')[1]) + 1):
            print("Convergence of demes seems correct")
        elif len(Variables.ConvergenceDeme.split(' ')) > (3 * int(Variables.ConvergenceDeme.split(' ')[1]) + 1):
            print("Less values than expected in convergence demes")
            sys.exit('ERROR!! Please check ConvergenceDeme value')
        elif len(Variables.ConvergenceDeme.split(' ')) < (3 * int(Variables.ConvergenceDeme.split(' ')[1]) + 1):
            print("More values than expected in convergence demes")
            sys.exit('ERROR!! Please check ConvergenceDeme value')

        # Template
    if Variables.Template != '':
        if Variables.Template[-4:] == '.pdb':
            print(Variables.Template + ' selected as template')
        else:
            print('You must select a .pdb file as template')
            sys.exit('ERROR!! Please check Template value')
    else:
        print('Template not selected')

        #Chain
    if Variables.ChainPDB != '':
        if type(Variables.ChainPDB)==str:
            if len(Variables.ChainPDB) == 1:
                warnings.simplefilter('ignore', PDBConstructionWarning)
                PDBFile = './' + Variables.Template
                check = 0
                with open(PDBFile, 'r') as pdb_file:
                    for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                        if record.annotations["chain"] == Variables.ChainPDB:
                            check = 1
                        else:
                            pass
                if check == 1:
                    print('Chain ' + str(Variables.ChainPDB) + ' selected')
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



        # GMRCA
    lines = []
    if Variables.GMRCA != '':
        with open(Variables.GMRCA, "r") as fp:
            # read an store all lines into list
            lines = fp.readlines()
            for y in range(0, Variables.NumAA):
                if lines[x][y] == 'A' or lines[x][y] == 'R' or lines[x][y] == 'N' or lines[x][y] == 'D' or lines[x][y] == 'C' or lines[x][y] == 'Q' or lines[x][y] == 'E' or lines[x][y] == 'G' or lines[x][y] == 'H' or lines[x][y] == 'I' or lines[x][y] == 'L' or lines[x][y] == 'K' or lines[x][y] == 'M' or lines[x][y] == 'F' or lines[x][y] == 'P' or lines[x][y] == 'S' or lines[x][y] == 'T' or lines[x][y] == 'W' or lines[x][y] == 'Y' or lines[x][y] == 'V':
                    pass
                else:
                    print('Error in GMRCA: ' + lines[y] + ' character was found in position ' + str(y))
                    sys.exit('Error in GMRCA: ' + lines[y] + ' character is not allow')
    else:
        print('GMRCA not selected')
