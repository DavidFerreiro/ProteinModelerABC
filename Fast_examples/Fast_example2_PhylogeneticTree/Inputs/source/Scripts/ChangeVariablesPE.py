import bin.Scripts.Variables
import os
import sys
from Bio import SeqIO
from random import randint, uniform
import numpy as np
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning


# Change variables to PE
def chavar():

    # Coalescent or Tree
    if bin.Scripts.Variables.Coalescent == 'Phylo':
        Treefile_PE = '-p' + str(bin.Scripts.Variables.Tree) + ' '

    else:
        # Change variables to PE
        if bin.Scripts.Variables.NumSeq != '':
            if bin.Scripts.Variables.NumAA != '':
                SampleSize_SequenceLength_PE = '-s' + str(bin.Scripts.Variables.NumSeq) + ' ' + str(bin.Scripts.Variables.NumAA) + ' '
            else:
                SampleSize_SequenceLength_PE = '-s' + str(bin.Scripts.Variables.NumSeq) + ' '
        else:
            SampleSize_SequenceLength_PE = ''


        # Population size
        if bin.Scripts.Variables.PopulationSize != '':
            if bin.Scripts.Variables.HaploidDiploid != '':
                PopSize_Hap_Dip_PE = '-e' + str(bin.Scripts.Variables.PopulationSize) + ' ' + str(bin.Scripts.Variables.HaploidDiploid) + ' '
            else:
                PopSize_Hap_Dip_PE = '-e' + str(bin.Scripts.Variables.PopulationSize) + ' '
        else:
            PopSize_Hap_Dip_PE = ''


        # Dated tips
        if bin.Scripts.Variables.DatedTips != '':
            DatedTips_PE = '-=' + str(bin.Scripts.Variables.DatedTips) + ' '
        else:
            DatedTips_PE = ''


        # Generation time
        if bin.Scripts.Variables.GenerationTime != '':
            if bin.Scripts.Variables.GenerationTime[0] == 'fix':
                GenTime_PE = '-/' + str(bin.Scripts.Variables.GenerationTime[1]) + ' '
            elif bin.Scripts.Variables.GenerationTime[0] == 'uniform':
                GenTime_PE = '-/' + str(randint(bin.Scripts.Variables.GenerationTime[1], bin.Scripts.Variables.GenerationTime[2])) + ' '
        else:
            GenTime_PE = ''


        # Geowth rate
        if bin.Scripts.Variables.GrowthRate != '':
            GrowthRate_PE = '-g' + bin.Scripts.Variables.GrowthRate + ' '
        else:
            GrowthRate_PE = ''


        # Migration model
        if bin.Scripts.Variables.MigrationModel != '':
            MigrationModel_PE = '-q' + bin.Scripts.Variables.MigrationModel + ' '
        else:
            MigrationModel_PE = ''


        # Migration rate
        if bin.Scripts.Variables.MigrationRate != '':
            MigrationRate_PE = '-t' + bin.Scripts.Variables.MigrationRate + ' '
        else:
            MigrationRate_PE = ''


        # Convergence demes
        if bin.Scripts.Variables.ConvergenceDemes != '':
            ConvergenceDemes_PE = '-%' + bin.Scripts.Variables.ConvergenceDemes + ' '
        else:
            ConvergenceDemes_PE = ''


        # Substitution rate
        if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] != '':
            if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'fix':
                SubstitutionRate_PE = '-u' + str(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]) + ' '
            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                SubstitutionRate_PE = '-u' + str(np.random.uniform(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'normal':
                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'exp':
                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(
                        np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'beta':
                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
        else:
            SubstitutionRate_PE = ''


    #Amino acid frequencies
    if bin.Scripts.Variables.AminoacidFrequencies != '':
        if bin.Scripts.Variables.AminoacidFrequencies == "Default":
            AminoacidFrequencies_PE = '-f20 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 '
        if bin.Scripts.Variables.AminoacidFrequencies.split(' ')[0] == "fix":
            AminoacidFrequencies_PE = '-f20 ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[1]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[2]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[3]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[4]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[5]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[6]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[7]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[8]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[9]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[10]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[11]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[12]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[13]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[14]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[15]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[16]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[17]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[18]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[19]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[20]) + ' '
        if bin.Scripts.Variables.AminoacidFrequencies.split(' ')[0] == "dirichlet":
            AminoacidFrequencies_PE = '-f20 ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[1]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[2]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[3]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[4]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[5]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[6]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[7]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[8]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[9]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[10]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[11]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[12]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[13]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[14]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[15]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[16]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[17]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[18]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[19]) + ' ' + str(bin.Scripts.Variables.AminoacidFrequencies.split(' ')[20]) + ' '
    else:
        AminoacidFrequencies_PE = ''


    # Rate heterogenous sites
    if bin.Scripts.Variables.RateHetSites.split(' ')[0] != '':
        if bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'fix':
            RateHetSites_PE = '-a' + str(bin.Scripts.Variables.RateHetSites.split(' ')[1]) + ' '
        elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'uniform':
            RateHetSites_PE = '-a' + str(np.random.uniform(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1))[1:-1] + ' '
        elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'normal':
            if bin.Scripts.Variables.RateHetSites.split(' ')[3] == 't':
                val = np.random.normal(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(bin.Scripts.Variables.RateHetSites.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.RateHetSites.split(' ')[4]):
                    val = np.random.normal(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1)[1:-1] + ' '
                RateHetSites_PE = '-a' + str(val) + ' '
            else:
                RateHetSites_PE = '-a' + str(np.random.normal(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1))[1:-1] + ' '
        elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'exp':
            if bin.Scripts.Variables.RateHetSites.split(' ')[3] == 't':
                val = np.random.exponential(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(bin.Scripts.Variables.RateHetSites.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.RateHetSites.split(' ')[4]):
                    val = np.random.exponential(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1)[1:-1] + ' '
                RateHetSites_PE = '-a' + str(val) + ' '
            else:
                RateHetSites_PE = '-a' + str(np.random.exponential(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1))[1:-1] + ' '
        elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'gamma':
            if bin.Scripts.Variables.RateHetSites.split(' ')[3] == 't':
                val = np.random.gamma(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(bin.Scripts.Variables.RateHetSites.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.RateHetSites.split(' ')[4]):
                    val = np.random.gamma(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]),float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1)[1:-1] + ' '
                RateHetSites_PE = '-a' + str(val) + ' '
            else:
                RateHetSites_PE = '-a' + str(np.random.gamma(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]),1))[1:-1] + ' '
        elif bin.Scripts.Variables.RateHetSites.split(' ')[0] == 'beta':
            if bin.Scripts.Variables.RateHetSites.split(' ')[3] == 't':
                val = np.random.exponential(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(bin.Scripts.Variables.RateHetSites.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.RateHetSites.split(' ')[4]):
                    val = np.random.exponential(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]), 1)[1:-1] + ' '
                RateHetSites_PE = '-a' + str(val) + ' '
            else:
                RateHetSites_PE = '-a' + str(np.random.exponential(float(bin.Scripts.Variables.RateHetSites.split(' ')[1]), float(bin.Scripts.Variables.RateHetSites.split(' ')[2]),1))[1:-1] + ' '
    else:
        RateHetSites_PE = ''


    # Proportion of invariables sites
    if bin.Scripts.Variables.PropInvSites.split(' ')[0] != '':
        if bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'fix':
            PropInvSites_PE = '-i' + str(bin.Scripts.Variables.PropInvSites.split(' ')[1]) + ' '
        elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'uniform':
            bin.Scripts.Variables.PropInvSites_PE = '-i' + str(np.random.uniform(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1)) + ' '
        elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'normal':
            if bin.Scripts.Variables.PropInvSites.split(' ')[3] == 't':
                val = np.random.normal(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(bin.Scripts.Variables.PropInvSites.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.PropInvSites.split(' ')[4]):
                    val = np.random.normal(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1)[1:-1] + ' '
                PropInvSites_PE = '-i' + str(val) + ' '
            else:
                PropInvSites_PE = '-i' + str(np.random.normal(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1))[1:-1] + ' '
        elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'exp':
            if bin.Scripts.Variables.PropInvSites.split(' ')[3] == 't':
                val = np.random.exponential(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(bin.Scripts.Variables.PropInvSites.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.PropInvSites.split(' ')[4]):
                    val = np.random.exponential(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1)[1:-1] + ' '
                PropInvSites_PE = '-i' + str(val) + ' '
            else:
                PropInvSites_PE = '-i' + str(np.random.exponential(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1))[1:-1] + ' '
        elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'gamma':
            if bin.Scripts.Variables.PropInvSites.split(' ')[3] == 't':
                val = np.random.gamma(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(bin.Scripts.Variables.PropInvSites.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.PropInvSites.split(' ')[4]):
                    val = np.random.gamma(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]),float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1)[1:-1] + ' '
                PropInvSites_PE = '-i' + str(val) + ' '
            else:
                PropInvSites_PE = '-i' + str(np.random.gamma(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]),1))[1:-1] + ' '
        elif bin.Scripts.Variables.PropInvSites.split(' ')[0] == 'beta':
            if bin.Scripts.Variables.PropInvSites.split(' ')[3] == 't':
                val = np.random.exponential(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(bin.Scripts.Variables.PropInvSites.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.PropInvSites.split(' ')[4]):
                    val = np.random.exponential(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]), 1)[1:-1] + ' '
                PropInvSites_PE = '-i' + str(val) + ' '
            else:
                PropInvSites_PE = '-i' + str(np.random.exponential(float(bin.Scripts.Variables.PropInvSites.split(' ')[1]), float(bin.Scripts.Variables.PropInvSites.split(' ')[2]),1))[1:-1] + ' '
    else:
        PropInvSites_PE = ''


    # Template
    if bin.Scripts.Variables.Template != '':
        warnings.simplefilter('ignore', PDBConstructionWarning)
        Template_PE = '-zPop_evol.in '
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

        GMRCA_Name = 'seqGMRCA'
        fp = open(GMRCA_Name, 'w')
        fp.write(str(PDB_Seq))
        fp.close()
        GMRCA_PE = '-x' + GMRCA_Name + ' '
        with open('ProteinModelerABC.out', 'a') as Error:
            Error.write('PDB sequence will be used as GMRCA\n')
            Error.close()
        print('PDB sequence will be used as GMRCA\n')
    else:
        Template_PE = ''
        if bin.Scripts.Variables.GMRCA != '':
            GMRCA_PE = '-x' + bin.Scripts.Variables.GMRCA + ' '
        else:
            GMRCA_PE = ''


    # Information on the screen
    if bin.Scripts.Variables.ShowInformationScreen == 'No':
        ShowInformationScreen_PE = 0
    else:
        ShowInformationScreen_PE = 1


    # Create variables for simulations with PE
    if len(bin.Scripts.Variables.SubstitutionModel.split(' ')) < 2:
        sys.exit('Error!! You must select at least 2 substitution models')
    else:
        SubModel_PE=[]
        SubstitutionRate_List=[]
        SubModel_PE_SCS = []

        for x in range(0, len(bin.Scripts.Variables.SubstitutionModel.split(' '))):
            if bin.Scripts.Variables.SubstitutionModel.split(' ')[x] != 'Fitness' and bin.Scripts.Variables.SubstitutionModel.split(' ')[x] != 'Neutral':
                for y in range(1, (bin.Scripts.Variables.NumberOfSimulations + 1)):
                    if x == 0:
                        if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] != '':
                            if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'fix':
                                SubstitutionRate_PE = '-u' + str(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]) + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                                SubstitutionRate_PE = '-u' + str(np.random.uniform(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'normal':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                    val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                        val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                    SubstitutionRate_PE = '-u' + str(val) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                else:
                                    SubstitutionRate_PE = '-u' + str(np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'exp':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                    val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                    SubstitutionRate_PE = '-u' + str(val) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                else:
                                    SubstitutionRate_PE = '-u' + str(np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                    val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                        val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                    SubstitutionRate_PE = '-u' + str(val) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                else:
                                    SubstitutionRate_PE = '-u' + str(np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'beta':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                    val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                    SubstitutionRate_PE = '-u' + str(val) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                else:
                                    SubstitutionRate_PE = '-u' + str(np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                        else:
                            SubstitutionRate_PE = ''
                            SubstitutionRate_List.append(SubstitutionRate_PE)
                        if bin.Scripts.Variables.Tree != '':
                            SubModel_PE_used = '-n1 ' + str(Treefile_PE) + str(AminoacidFrequencies_PE) + '-@' + str(bin.Scripts.Variables.SubstitutionModel.split(' ')[x]) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(y)
                            SubModel_PE.append(str(SubModel_PE_used))
                        else:
                            SubModel_PE_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@' + str(bin.Scripts.Variables.SubstitutionModel.split(' ')[x]) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(y)
                            SubModel_PE.append(str(SubModel_PE_used))
                    elif x != 0:
                        w = y + (x * bin.Scripts.Variables.NumberOfSimulations)
                        if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] != '':
                            if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'fix':
                                SubstitutionRate_PE = '-u' + str(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]) + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                                SubstitutionRate_PE = '-u' + str(np.random.uniform(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'normal':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                    val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                    SubstitutionRate_PE = '-u' + str(val) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                else:
                                    SubstitutionRate_PE = '-u' + str(np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'exp':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                    val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                    SubstitutionRate_PE = '-u' + str(val) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                else:
                                    SubstitutionRate_PE = '-u' + str(np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                    val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                        val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                    SubstitutionRate_PE = '-u' + str(val) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                else:
                                    SubstitutionRate_PE = '-u' + str(np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'beta':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                    val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                    SubstitutionRate_PE = '-u' + str(val) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                else:
                                    SubstitutionRate_PE = '-u' + str(np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                        else:
                            SubstitutionRate_PE = ''
                            SubstitutionRate_List.append(SubstitutionRate_PE)

                        if bin.Scripts.Variables.Tree != '':
                            SubModel_PE_used = '-n1 ' + str(Treefile_PE) + str(AminoacidFrequencies_PE) + '-@' + str(bin.Scripts.Variables.SubstitutionModel.split(' ')[x]) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                            SubModel_PE.append(str(SubModel_PE_used))
                        else:
                            SubModel_PE_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@' + str(bin.Scripts.Variables.SubstitutionModel.split(' ')[x]) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                            SubModel_PE.append(str(SubModel_PE_used))

            else:
                for y in range(1, (bin.Scripts.Variables.NumberOfSimulations + 1)):
                    w = y + (x * bin.Scripts.Variables.NumberOfSimulations)
                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] != '':
                        if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'fix':
                            SubstitutionRate_PE = '-u' + str(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]) + ' '
                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                            SubstitutionRate_PE = '-u' + str(np.random.uniform(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'normal':
                            if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1)[1:-1]
                                while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                    val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                SubstitutionRate_PE = '-u' + str(val) + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = '-u' + str(np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'exp':
                            if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]), float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                    val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                SubstitutionRate_PE = '-u' + str(val) + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = '-u' + str(np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                            if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1)[1:-1]
                                while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                    val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                SubstitutionRate_PE = '-u' + str(val) + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = '-u' + str(np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                        elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'beta':
                            if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                    val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                SubstitutionRate_PE = '-u' + str(val) + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = '-u' + str(np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                    else:
                        SubstitutionRate_PE = ''
                        SubstitutionRate_List.append(SubstitutionRate_PE)

                    if bin.Scripts.Variables.Tree != '':
                        SubModel_PE_SCS_used = '-n1 ' + str(Treefile_PE) + str(AminoacidFrequencies_PE) + '-@JTT ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(Template_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                        SubModel_PE_SCS.append(str(SubModel_PE_SCS_used))

                    else:
                        SubModel_PE_SCS_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@JTT ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(Template_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                        SubModel_PE_SCS.append(str(SubModel_PE_SCS_used))

    # if file exists, delate
    if os.path.exists('ProteinEvolver_Arguments.txt') == True:
        os.system('rm ProteinEvolver_Arguments.txt')
    # Next, write it
    with open('ProteinEvolver_Arguments.txt', 'a') as file:
        for emp in SubModel_PE:
            file.write(emp + '\n')
        for scs in SubModel_PE_SCS:
            file.write(scs + '\n')
        file.close()
    #bin.Scripts.Variables.PE_used = SubModel_PE
    #bin.Scripts.Variables.PE_SCS_used = SubModel_PE_SCS

    if bin.Scripts.Variables.Coalescent == 'Phylo':
        pass
    else:
        bin.Scripts.Variables.SubstitutionRate_L = SubstitutionRate_List

        for x in range(0, len(SubstitutionRate_List)):
            bin.Scripts.Variables.ThetaRate_L.append(2 * bin.Scripts.Variables.HaploidDiploid * bin.Scripts.Variables.PopulationSize * float(SubstitutionRate_List[x]) * bin.Scripts.Variables.NumAA)
        #bin.Scripts.Variables.AminoacidFrequencies = AminoacidFrequencies_PE
        bin.Scripts.Functions.plot()

    #Write Pop_evol.in file needed for simulations under SCS models
    with open('Pop_evol.in', 'w') as fp:
        fp.write('PDB= ' + bin.Scripts.Variables.Template + ' 		# pdb file\n')
        fp.write('CHAIN=  ' + str(bin.Scripts.Variables.ChainPDB) + '			# Chain of pdb file\n')
        fp.write('FILE_STR=	structures.in 	# List of contact matrices\n')
        fp.write('TEMP=	' + str(bin.Scripts.Variables.TEMP) + '			# Temperature\n')
        fp.write('S0=	' + str(bin.Scripts.Variables.S0) + '			# s0, configurational entropy per residue (unfolded)\n')
        fp.write('SC1=    ' + str(bin.Scripts.Variables.SC1) + '			# configurational entropy per residue (misfolded)\n')
        fp.write('SC0=    ' + str(bin.Scripts.Variables.SC0) + '			# configurational entropy offset (misfolded)\n')
        fp.write('REM3=	' + str(bin.Scripts.Variables.REM3) + '			# Third cumulant in REM calculation\n')
        fp.write('NEUTRAL=	 0		# If 1, Neutral landscape, otherwise population size dependent selection \n')
        fp.write('NPOP=	' + str(bin.Scripts.Variables.NPOP) + '			# N_pop, population size\n')
        fp.write('TYPE_BL=	2		# Branch type: "1" branches by mutations, "2" branches by substitutions\n')
        fp.write('OUTPUT_LEVEL=	0		# Print all output files "2"; Print only "final" output files "1"; Do not print them (default) "0"\n')
        fp.close()
