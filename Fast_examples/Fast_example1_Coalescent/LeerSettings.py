import sys
import re
import Variables


def leer():
    #print(sistema + ' system detected\n')
    file = 'Settings.txt'
    # Reading Settings and creating each variable
    try:
        with open(file) as sett:
            for line in sett:
                line = line.split('=')
                line[0] = line[0].strip()
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
        with open ('ProteinModelerABC.out', 'a') as Error:
            Error.write('ERROR!!! Cannot open ' + file + '\n')
            Error.close()
        sys.exit('ERROR!!! Cannot open ' + file + '\n')

    #Print some information
    #print("> Input file uploaded: " + file + "\n\n")



    #Correct Variables
    Variables.NameOfPhylipFile = ' '.join(NameOfPhylipFile)
    Variables.NumberOfSimulations = int(NumberOfSimulations)
    Variables.Indels = ' '.join(Indels)
    try:
        Variables.NumberOfProcessors = int(NumberOfProcessors)
    except:
        Variables.NumberOfProcessors = ' '.join(NumberOfProcessors)
    Variables.SaveSimulations = ' '.join(SaveSimulations)
    Variables.ShowInformationScreen = ' '.join(ShowInformationScreen)
    try:
        Variables.HaploidDiploid = int(HaploidDiploid)
    except:
        Variables.HaploidDiploid = ''
    try:
        Variables.PopulationSize = int(PopulationSize)
    except:
        Variables.PopulationSize = ''
    Variables.DatedTips = ' '.join(DatedTips)
    Variables.GenerationTime = ' '.join(GenerationTime)
    #RecombinationRate = ' '.join(RecombinationRate)
    Variables.SubstitutionRate = ' '.join(SubstitutionRate)
    Variables.SubstitutionModel = ' '.join(SubstitutionModel)
    #Variables.StructuralSubstitutionModel = ' '.join(StructuralSubstitutionModel)
    Variables.Template = ' '.join(Template)
    Variables.ChainPDB= ' '.join(ChainPDB)
    Variables.AminoacidFrequencies = ' '.join(AminoacidFrequencies)
    Variables.RateHetSites = ' '.join(RateHetSites)
    Variables.PropInvSites = ' '.join(PropInvSites)
    Variables.TEMP = float(TEMP)
    Variables.S0 = float(S0)
    Variables.SC1 = float(SC1)
    Variables.SC0 = float(SC0)
    Variables.REM3 = int(REM3)
    Variables.NPOP = int(NPOP)
    Variables.ABCIterations = int(ABCIterations)
    Variables.ABCTolerance = float(ABCTolerance)
    Variables.ABCMethod = ' '.join(ABCMethod)
    Variables.SummaryStatistics = ' '.join(SummaryStatistics)
    Variables.MultiPage = ' '.join(MultiPage)
    Variables.GrowthRate = ' '.join(GrowthRate)
    Variables.MigrationModel = ' '.join(MigrationModel)
    Variables.MigrationRate = ' '.join(MigrationRate)
    Variables.ConvergenceDemes = ' '.join(ConvergenceDemes)
    Variables.Coalescent = ' '.join(Coalescent)
    Variables.Tree = ' '.join(Tree)
    
    
    #Number of sequences and aminoacids
    Variables.NumSeq=int()
    Variables.NumAA=int()


    try:
        with open(Variables.NameOfPhylipFile) as Filo:
            Variables.NumSeq = sum(1 for line in Filo)
            Variables.NumSeq = Variables.NumSeq - 1
    except OSError:
        with open ('ProteinModelerABC.out', 'a') as Error:
            Error.write('ERROR!!! Cannot read ' + Variables.NameOfPhylipFile + '. Check the alignment file\n')
            Error.close()
        sys.exit('ERROR!!! Cannot read ' + Variables.NameOfPhylipFile + '. Check the alignment file\n')



    try:
        with open(Variables.NameOfPhylipFile) as Filo:
            for line in Filo:
                Variables.NumAA = line.count('',12)


    except OSError:
        with open ('ProteinModelerABC.out', 'a') as Error:
            Error.write('ERROR!!! Cannot read ' + Variables.NameOfPhylipFile + '. Check the alignment file\n')
            Error.close()
        sys.exit('ERROR!!! Cannot read ' + Variables.NameOfPhylipFile + '. Check the alignment file\n')


    Variables.LongSM = 0
    for x in Variables.SubstitutionModel.split(' '):
        if x == 'Fitness' or x == 'Neutral':
            Variables.LongSM = Variables.LongSM + 1

    if Variables.LongSM == 0:
        sys.exit('ERROR!!! User must select at least one structural constrained substitution (SCS) model . Check Settings.txt file\n')


    if Variables.SubstitutionModel == '':
        Variables.LongEM = 0
    else:
        Variables.LongEM = len(Variables.SubstitutionModel.split(' ')) - Variables.LongSM


    Variables.Total_simu = Variables.NumberOfSimulations * (Variables.LongSM + Variables.LongEM)
    for x in range(1, Variables.Total_simu + 1):  # name of simulation files
        if x < 10:
            Variables.Simulations_Name.append('sequences00000000' + str(x))
        elif x < 100:
            Variables.Simulations_Name.append('sequences0000000' + str(x))
        elif x < 1000:
            Variables.Simulations_Name.append('sequences000000' + str(x))
        elif x < 10000:
            Variables.Simulations_Name.append('sequences00000' + str(x))
        elif x < 100000:
            Variables.Simulations_Name.append('sequences0000' + str(x))
        elif x < 1000000:
            Variables.Simulations_Name.append('sequences000' + str(x))
        elif x < 10000000:
            Variables.Simulations_Name.append('sequences00' + str(x))
        elif x < 100000000:
            Variables.Simulations_Name.append('sequences0' + str(x))
        elif x < 1000000000:
            Variables.Simulations_Name.append('sequences' + str(x))

    Variables.Total_E_simu = Variables.NumberOfSimulations * (Variables.LongEM)
     
    Variables.Total_S_simu = Variables.NumberOfSimulations * (Variables.LongSM)