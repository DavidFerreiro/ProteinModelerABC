import sys
import re
import bin.Scripts.Variables


def leer():
    file = 'Settings.txt'
    # Reading Settings and creating each variable
    try:
        with open(file) as sett:
            for line in sett:
                line = line.split('=')
                line[0] = line[0].strip()
                if line[0] == '*MPI':
                    MPI = re.sub(r'','', line[1].strip()).split()
                if line[0] == '*NameOfPhylipFile':
                    NameOfPhylipFile = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '*NumberOfSimulations':
                    NumberOfSimulations = float(line[1])
                elif line[0] == '*Indels':
                    Indels = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '*NumberOfProcessors':
                    try:
                        NumberOfProcessors = float(line[1])
                    except:
                        NumberOfProcessors = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '*SaveSimulations':
                    SaveSimulations = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '*ShowInformationScreen':
                    ShowInformationScreen = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '|*Haploid/Diploid':
                    try:
                        HaploidDiploid = float(line[1])
                    except:
                        HaploidDiploid = ''
                elif line[0] == '|*PopulationSize':
                    try:
                        PopulationSize = float(line[1])
                    except:
                        PopulationSize = ''
                elif line[0] == '|DatedTips':
                    DatedTips = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '|GenerationTime':
                    GenerationTime = re.sub(r'','', line[1].strip()).split()
                #elif line[0] == '*RecombinationRate':
                    #RecombinationRate = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '|*SubstitutionRate':
                    SubstitutionRate = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '*SubstitutionModel':
                    SubstitutionModel = re.sub(r'','', line[1].strip()).split()
                #elif line[0] == '*StructuralSubstitutionModel':
                    #StructuralSubstitutionModel = re.sub(r'', '', line[1].strip()).split()
                elif line[0] == '*Template':
                    Template = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '*Chain':
                    ChainPDB = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '|*AminoacidFrequencies':
                    AminoacidFrequencies = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '|RateHetSites':
                    RateHetSites = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '|PropInvSites':
                    PropInvSites = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '|*TEMP':
                    TEMP = float(line[1])
                elif line[0] == '|*S0':
                    S0 = float(line[1])
                elif line[0] == '|*SC1':
                    SC1 = float(line[1])
                elif line[0] == '|*SC0':
                    SC0 = float(line[1])
                elif line[0] == '|*REM3':
                    REM3 = int(line[1])
                elif line[0] == '|*NPOP':
                    NPOP = int(line[1])
                elif line[0] == '*ABCIterations':
                    ABCIterations = float(line[1])
                elif line[0] == '*ABCTolerance':
                    ABCTolerance = float(line[1])
                elif line[0] == '*ABCMethod':
                    ABCMethod = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '*SummaryStatistics':
                    SummaryStatistics = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '*MultiPage':
                    MultiPage = re.sub(r'','', line[1].strip()).split()
                elif line[0] == '|GrowthRate':
                    GrowthRate = re.sub(r'', '', line[1].strip()).split()
                elif line[0] == '|MigrationModel':
                    MigrationModel = re.sub(r'', '', line[1].strip()).split()
                elif line[0] == '|MigrationRate':
                    MigrationRate = re.sub(r'', '', line[1].strip()).split()
                elif line[0] == '|ConvergenceDemes':
                    ConvergenceDemes = re.sub(r'', '', line[1].strip()).split()
                elif line[0] == '*CoalescentOrPhylogeny':
                    Coalescent = re.sub(r'', '', line[1].strip()).split()
                elif line[0] == '|*Tree':
                    Tree = re.sub(r'', '', line[1].strip()).split()


    except OSError:
        sys.exit('ERROR!!! Cannot open ' + file + '\n')


    #Correct Variables
    if 'mi_variable' in locals():
        bin.Scripts.Variables.MPI = ' '.join(MPI)
    bin.Scripts.Variables.NameOfPhylipFile = ' '.join(NameOfPhylipFile)
    bin.Scripts.Variables.NumberOfSimulations = int(NumberOfSimulations)
    bin.Scripts.Variables.Indels = ' '.join(Indels)
    try:
        bin.Scripts.Variables.NumberOfProcessors = int(NumberOfProcessors)
    except:
        bin.Scripts.Variables.NumberOfProcessors = ' '.join(NumberOfProcessors)
    bin.Scripts.Variables.SaveSimulations = ' '.join(SaveSimulations)
    bin.Scripts.Variables.ShowInformationScreen = ' '.join(ShowInformationScreen)
    try:
        bin.Scripts.Variables.HaploidDiploid = int(HaploidDiploid)
    except:
        bin.Scripts.Variables.HaploidDiploid = ''
    try:
        bin.Scripts.Variables.PopulationSize = int(PopulationSize)
    except:
        bin.Scripts.Variables.PopulationSize = ''
    bin.Scripts.Variables.DatedTips = ' '.join(DatedTips)
    bin.Scripts.Variables.GenerationTime = ' '.join(GenerationTime)
    #RecombinationRate = ' '.join(RecombinationRate)
    bin.Scripts.Variables.SubstitutionRate = ' '.join(SubstitutionRate)
    bin.Scripts.Variables.SubstitutionModel = ' '.join(SubstitutionModel)
    #bin.Scripts.Variables.StructuralSubstitutionModel = ' '.join(StructuralSubstitutionModel)
    bin.Scripts.Variables.Template = ' '.join(Template)
    bin.Scripts.Variables.ChainPDB= ' '.join(ChainPDB)
    bin.Scripts.Variables.AminoacidFrequencies = ' '.join(AminoacidFrequencies)
    bin.Scripts.Variables.RateHetSites = ' '.join(RateHetSites)
    bin.Scripts.Variables.PropInvSites = ' '.join(PropInvSites)
    bin.Scripts.Variables.TEMP = float(TEMP)
    bin.Scripts.Variables.S0 = float(S0)
    bin.Scripts.Variables.SC1 = float(SC1)
    bin.Scripts.Variables.SC0 = float(SC0)
    bin.Scripts.Variables.REM3 = int(REM3)
    bin.Scripts.Variables.NPOP = int(NPOP)
    bin.Scripts.Variables.ABCIterations = int(ABCIterations)
    bin.Scripts.Variables.ABCTolerance = float(ABCTolerance)
    bin.Scripts.Variables.ABCMethod = ' '.join(ABCMethod)
    bin.Scripts.Variables.SummaryStatistics = ' '.join(SummaryStatistics)
    bin.Scripts.Variables.MultiPage = ' '.join(MultiPage)
    bin.Scripts.Variables.GrowthRate = ' '.join(GrowthRate)
    bin.Scripts.Variables.MigrationModel = ' '.join(MigrationModel)
    bin.Scripts.Variables.MigrationRate = ' '.join(MigrationRate)
    bin.Scripts.Variables.ConvergenceDemes = ' '.join(ConvergenceDemes)
    bin.Scripts.Variables.Coalescent = ' '.join(Coalescent)
    bin.Scripts.Variables.Tree = ' '.join(Tree)
    
    
    #Number of sequences and aminoacids
    bin.Scripts.Variables.NumSeq=int()
    bin.Scripts.Variables.NumAA=int()


    try:
        with open(bin.Scripts.Variables.NameOfPhylipFile) as Filo:
            bin.Scripts.Variables.NumSeq = sum(1 for line in Filo)
            bin.Scripts.Variables.NumSeq = bin.Scripts.Variables.NumSeq - 1
    except OSError:
        sys.exit('ERROR!!! Cannot read ' + bin.Scripts.Variables.NameOfPhylipFile + '. Check the alignment file\n')



    try:
        with open(bin.Scripts.Variables.NameOfPhylipFile) as Filo:
            for line in Filo:
                bin.Scripts.Variables.NumAA = line.count('',12)


    except OSError:
        sys.exit('ERROR!!! Cannot read ' + bin.Scripts.Variables.NameOfPhylipFile + '. Check the alignment file\n')


    bin.Scripts.Variables.LongSM = 0
    for x in bin.Scripts.Variables.SubstitutionModel.split(' '):
        if x == 'Fitness' or x == 'Neutral':
            bin.Scripts.Variables.LongSM = bin.Scripts.Variables.LongSM + 1

    if bin.Scripts.Variables.LongSM == 0:
        sys.exit('ERROR!!! User must select at least one structural constrained substitution (SCS) model . Check Settings.txt file\n')


    if bin.Scripts.Variables.SubstitutionModel == '':
        bin.Scripts.Variables.LongEM = 0
    else:
        bin.Scripts.Variables.LongEM = len(bin.Scripts.Variables.SubstitutionModel.split(' ')) - bin.Scripts.Variables.LongSM


    bin.Scripts.Variables.Total_simu = bin.Scripts.Variables.NumberOfSimulations * (bin.Scripts.Variables.LongSM + bin.Scripts.Variables.LongEM)
    for x in range(1, bin.Scripts.Variables.Total_simu + 1):  # name of simulation files
        if x < 10:
            bin.Scripts.Variables.Simulations_Name.append('sequences00000000' + str(x))
        elif x < 100:
            bin.Scripts.Variables.Simulations_Name.append('sequences0000000' + str(x))
        elif x < 1000:
            bin.Scripts.Variables.Simulations_Name.append('sequences000000' + str(x))
        elif x < 10000:
            bin.Scripts.Variables.Simulations_Name.append('sequences00000' + str(x))
        elif x < 100000:
            bin.Scripts.Variables.Simulations_Name.append('sequences0000' + str(x))
        elif x < 1000000:
            bin.Scripts.Variables.Simulations_Name.append('sequences000' + str(x))
        elif x < 10000000:
            bin.Scripts.Variables.Simulations_Name.append('sequences00' + str(x))
        elif x < 100000000:
            bin.Scripts.Variables.Simulations_Name.append('sequences0' + str(x))
        elif x < 1000000000:
            bin.Scripts.Variables.Simulations_Name.append('sequences' + str(x))

    bin.Scripts.Variables.Total_E_simu = bin.Scripts.Variables.NumberOfSimulations * (bin.Scripts.Variables.LongEM)
     
    bin.Scripts.Variables.Total_S_simu = bin.Scripts.Variables.NumberOfSimulations * (bin.Scripts.Variables.LongSM)
