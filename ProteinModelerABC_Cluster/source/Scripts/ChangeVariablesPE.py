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
                SubstitutionRate_PE = '-u' + str(np.random.uniform(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                   float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[
                                             1:-1] + ' '
            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'normal':
                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                           float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(
                            bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                               float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                      float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),
                                                                      1))[1:-1] + ' '
            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'exp':
                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(
                            bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                    float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(
                        np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                              float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                          float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(
                            bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                              float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                     float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),
                                                                     1))[1:-1] + ' '
            elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'beta':
                if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(
                            bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                    float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(
                        np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                              float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
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
        #with open('ProteinModelerABC.out', 'a') as Error:
            #Error.write('PDB sequence will be used as GMRCA\n')
            #Error.close()
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


    # Create documents for simulations with PE
    if len(bin.Scripts.Variables.SubstitutionModel.split(' ')) < 2:
        sys.exit('Error!! You must select at least 2 substitution models')
    else:
        SubModel_PE=[]
        SubstitutionRate_List=[]
        SubModel_PE_SCS = []

        for x in bin.Scripts.Variables.SubstitutionModel.split(' '):
            if x != 'Fitness' and x != 'Neutral':
                with open('ProteinModelerABC_arguments.txt', 'a') as Prot_a:
                    for y in range(1, (bin.Scripts.Variables.NumberOfSimulations + 1)):
                        if bin.Scripts.Variables.SubstitutionModel.split(' ').index(x) == 0:
                            if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] != '':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'fix':
                                    SubstitutionRate_PE = '-u' + str(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                                    SubstitutionRate_PE = '-u' + str(
                                        np.random.uniform(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                          float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'normal':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                               float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                   float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                             float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                     1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'exp':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                    float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                              1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                        float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                  float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                     1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                              float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                  float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                            float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'beta':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                    float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                              1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                        float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                  float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                     1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = ''
                                SubstitutionRate_List.append(SubstitutionRate_PE)
                            if bin.Scripts.Variables.Tree != '':
                                SubModel_PE_used = '-n1 ' + str(Treefile_PE) + str(AminoacidFrequencies_PE) + '-@' + str(x) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(y)
                                SubModel_PE.append(str(SubModel_PE_used))
                                Prot_a.write(SubModel_PE_used + '\n')
                            else:
                                SubModel_PE_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@' + str(x) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(y)
                                SubModel_PE.append(str(SubModel_PE_used))
                                Prot_a.write(SubModel_PE_used + '\n')
                        elif bin.Scripts.Variables.SubstitutionModel.split(' ').index(x) != 0:
                            w = y + ((int(bin.Scripts.Variables.SubstitutionModel.split(' ').index(x))) * bin.Scripts.Variables.NumberOfSimulations)
                            if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] != '':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'fix':
                                    SubstitutionRate_PE = '-u' + str(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                                    SubstitutionRate_PE = '-u' + str(
                                        np.random.uniform(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                          float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'normal':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                               float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                   float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                             float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                     1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'exp':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                    float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                              1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                        float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                  float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                     1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                              float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                  float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                            float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'beta':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                    float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                              1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                                        float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = ''
                                SubstitutionRate_List.append(SubstitutionRate_PE)

                            if bin.Scripts.Variables.Tree != '':
                                SubModel_PE_used = '-n1 ' + str(Treefile_PE) + str(AminoacidFrequencies_PE) + '-@' + str(x) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                                SubModel_PE.append(str(SubModel_PE_used))
                                Prot_a.write(SubModel_PE_used + '\n')
                            else:
                                SubModel_PE_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@' + str(x) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                                SubModel_PE.append(str(SubModel_PE_used))
                                Prot_a.write(SubModel_PE_used + '\n')
                    Prot_a.close()
            else:
                comprobar1 = 'Fitness' in bin.Scripts.Variables.SubstitutionModel.split(' ')
                comprobar2 = 'Neutral' in bin.Scripts.Variables.SubstitutionModel.split(' ')
                if comprobar1 == True and comprobar2 == True:
                    for k in range(0,2):
                        with open('ProteinModelerABC_S_arguments' + str(k) + '.txt', 'w') as Prot_Sa:
                            for y in range(1, (bin.Scripts.Variables.NumberOfSimulations + 1)):
                                w = y + ((len(bin.Scripts.Variables.SubstitutionModel.split(' ')) - 2 + k) * bin.Scripts.Variables.NumberOfSimulations)
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
                                            val = np.random.exponential(
                                                float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),
                                                float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                            while float(val) > float(
                                                    bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                                val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                            SubstitutionRate_PE = '-u' + str(val) + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                        else:
                                            SubstitutionRate_PE = '-u' + str(
                                                np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                        if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                            val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1)[1:-1]
                                            while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                                val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                            SubstitutionRate_PE = '-u' + str(val) + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                        else:
                                            SubstitutionRate_PE = '-u' + str(
                                                np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
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
                                    Prot_Sa.write(SubModel_PE_SCS_used + '\n')
                                else:
                                    SubModel_PE_SCS_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@JTT ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(Template_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                                    SubModel_PE_SCS.append(str(SubModel_PE_SCS_used))
                                    Prot_Sa.write(SubModel_PE_SCS_used + '\n')
                            Prot_Sa.close()
                elif comprobar1 == True and comprobar2 == False:
                    with open('ProteinModelerABC_S_arguments.txt', 'w') as Prot_Sa:
                        for y in range(1, (bin.Scripts.Variables.NumberOfSimulations + 1)):
                            w = y + ((len(bin.Scripts.Variables.SubstitutionModel.split(' ')) -1) * bin.Scripts.Variables.NumberOfSimulations)
                            if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] != '':
                                if bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'fix':
                                    SubstitutionRate_PE = '-u' + str(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                                    SubstitutionRate_PE = '-u' + str(np.random.uniform(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'normal':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.normal(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1)[1:-1] + ' '
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
                                        SubstitutionRate_PE = '-u' + str(np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1)[1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.gamma(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif bin.Scripts.Variables.SubstitutionRate.split(' ')[0] == 'beta':
                                    if bin.Scripts.Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1)[1:-1]
                                        while float(val) > float(bin.Scripts.Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(bin.Scripts.Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.exponential(float(bin.Scripts.Variables.SubstitutionRate.split(' ')[1]),float(bin.Scripts.Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = ''
                                SubstitutionRate_List.append(SubstitutionRate_PE)

                            if bin.Scripts.Variables.Tree != '':
                                SubModel_PE_SCS_used = '-n1 ' + str(Treefile_PE) + str(AminoacidFrequencies_PE) + '-@JTT ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(Template_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                                SubModel_PE_SCS.append(str(SubModel_PE_SCS_used))
                                Prot_Sa.write(SubModel_PE_SCS_used + '\n')
                            else:
                                SubModel_PE_SCS_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@JTT ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(Template_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                                SubModel_PE_SCS.append(str(SubModel_PE_SCS_used))
                                Prot_Sa.write(SubModel_PE_SCS_used + '\n')
                        Prot_Sa.close()
                elif comprobar1 == False and comprobar2 == True:
                    with open('ProteinModelerABC_S_arguments.txt', 'w') as Prot_Sa:
                        for y in range(1, (bin.Scripts.Variables.NumberOfSimulations + 1)):
                            w = y + ((len(bin.Scripts.Variables.SubstitutionModel.split(' ')) -1) * bin.Scripts.Variables.NumberOfSimulations)
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
                                Prot_Sa.write(SubModel_PE_SCS_used + '\n')
                            else:
                                SubModel_PE_SCS_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@JTT ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(Template_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                                SubModel_PE_SCS.append(str(SubModel_PE_SCS_used))
                                Prot_Sa.write(SubModel_PE_SCS_used + '\n')
                        Prot_Sa.close()
                else:
                    pass


    if bin.Scripts.Variables.Coalescent == 'Phylo':
        pass
    else:
        with open('PSimulations.txt', 'w') as thetafile:
            thetafile.write('SubstitutionRate,Theta\n')
            for x in range(0, len(SubstitutionRate_List)):
                thetafile.write(str(SubstitutionRate_List[x]) + ',' + str(2 * bin.Scripts.Variables.HaploidDiploid * bin.Scripts.Variables.PopulationSize * float(SubstitutionRate_List[x]) * bin.Scripts.Variables.NumAA) + '\n')
            thetafile.close()

