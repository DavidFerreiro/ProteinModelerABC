import Variables
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
    if Variables.Coalescent == 'Phylo':
        Treefile_PE = '-p' + str(Variables.Tree) + ' '

    else:
        # Change variables to PE
        if Variables.NumSeq != '':
            if Variables.NumAA != '':
                SampleSize_SequenceLength_PE = '-s' + str(Variables.NumSeq) + ' ' + str(Variables.NumAA) + ' '
            else:
                SampleSize_SequenceLength_PE = '-s' + str(Variables.NumSeq) + ' '
        else:
            SampleSize_SequenceLength_PE = ''


        # Population size
        if Variables.PopulationSize != '':
            if Variables.HaploidDiploid != '':
                PopSize_Hap_Dip_PE = '-e' + str(Variables.PopulationSize) + ' ' + str(Variables.HaploidDiploid) + ' '
            else:
                PopSize_Hap_Dip_PE = '-e' + str(Variables.PopulationSize) + ' '
        else:
            PopSize_Hap_Dip_PE = ''


        # Dated tips
        if Variables.DatedTips != '':
            DatedTips_PE = '-=' + str(Variables.DatedTips) + ' '
        else:
            DatedTips_PE = ''


        # Generation time
        if Variables.GenerationTime != '':
            if Variables.GenerationTime[0] == 'fix':
                GenTime_PE = '-/' + str(Variables.GenerationTime[1]) + ' '
            elif Variables.GenerationTime[0] == 'uniform':
                GenTime_PE = '-/' + str(randint(Variables.GenerationTime[1], Variables.GenerationTime[2])) + ' '
        else:
            GenTime_PE = ''


        # Geowth rate
        if Variables.GrowthRate != '':
            GrowthRate_PE = '-g' + Variables.GrowthRate + ' '
        else:
            GrowthRate_PE = ''


        # Migration model
        if Variables.MigrationModel != '':
            MigrationModel_PE = '-q' + Variables.MigrationModel + ' '
        else:
            MigrationModel_PE = ''


        # Migration rate
        if Variables.MigrationRate != '':
            MigrationRate_PE = '-t' + Variables.MigrationRate + ' '
        else:
            MigrationRate_PE = ''


        # Convergence demes
        if Variables.ConvergenceDemes != '':
            ConvergenceDemes_PE = '-%' + Variables.ConvergenceDemes + ' '
        else:
            ConvergenceDemes_PE = ''


        # Substitution rate
        if Variables.SubstitutionRate.split(' ')[0] != '':
            if Variables.SubstitutionRate.split(' ')[0] == 'fix':
                SubstitutionRate_PE = '-u' + str(Variables.SubstitutionRate.split(' ')[1]) + ' '
            elif Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                SubstitutionRate_PE = '-u' + str(np.random.uniform(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                   float(Variables.SubstitutionRate.split(' ')[2]), 1))[
                                             1:-1] + ' '
            elif Variables.SubstitutionRate.split(' ')[0] == 'normal':
                if Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),
                                           float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(
                            Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),
                                               float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                      float(Variables.SubstitutionRate.split(' ')[2]),
                                                                      1))[1:-1] + ' '
            elif Variables.SubstitutionRate.split(' ')[0] == 'exp':
                if Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(
                            Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                    float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(
                        np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                              float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
            elif Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                if Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),
                                          float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(
                            Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),
                                              float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                     float(Variables.SubstitutionRate.split(' ')[2]),
                                                                     1))[1:-1] + ' '
            elif Variables.SubstitutionRate.split(' ')[0] == 'beta':
                if Variables.SubstitutionRate.split(' ')[3] == 't':
                    val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                    while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(
                            Variables.SubstitutionRate.split(' ')[4]):
                        val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                    float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                    SubstitutionRate_PE = '-u' + str(val) + ' '
                else:
                    SubstitutionRate_PE = '-u' + str(
                        np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                              float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
        else:
            SubstitutionRate_PE = ''


    #Amino acid frequencies
    if Variables.AminoacidFrequencies != '':
        if Variables.AminoacidFrequencies == "Default":
            AminoacidFrequencies_PE = '-f20 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 '
        if Variables.AminoacidFrequencies.split(' ')[0] == "fix":
            AminoacidFrequencies_PE = '-f20 ' + str(Variables.AminoacidFrequencies.split(' ')[1]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[2]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[3]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[4]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[5]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[6]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[7]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[8]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[9]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[10]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[11]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[12]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[13]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[14]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[15]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[16]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[17]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[18]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[19]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[20]) + ' '
        if Variables.AminoacidFrequencies.split(' ')[0] == "dirichlet":
            AminoacidFrequencies_PE = '-f20 ' + str(Variables.AminoacidFrequencies.split(' ')[1]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[2]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[3]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[4]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[5]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[6]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[7]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[8]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[9]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[10]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[11]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[12]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[13]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[14]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[15]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[16]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[17]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[18]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[19]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[20]) + ' '
    else:
        AminoacidFrequencies_PE = ''


    # Rate heterogenous sites
    if Variables.RateHetSites.split(' ')[0] != '':
        if Variables.RateHetSites.split(' ')[0] == 'fix':
            RateHetSites_PE = '-a' + str(Variables.RateHetSites.split(' ')[1]) + ' '
        elif Variables.RateHetSites.split(' ')[0] == 'uniform':
            RateHetSites_PE = '-a' + str(np.random.uniform(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]), 1))[1:-1] + ' '
        elif Variables.RateHetSites.split(' ')[0] == 'normal':
            if Variables.RateHetSites.split(' ')[3] == 't':
                val = np.random.normal(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.RateHetSites.split(' ')[5]) or float(val) < float(Variables.RateHetSites.split(' ')[4]):
                    val = np.random.normal(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]), 1)[1:-1] + ' '
                RateHetSites_PE = '-a' + str(val) + ' '
            else:
                RateHetSites_PE = '-a' + str(np.random.normal(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]), 1))[1:-1] + ' '
        elif Variables.RateHetSites.split(' ')[0] == 'exp':
            if Variables.RateHetSites.split(' ')[3] == 't':
                val = np.random.exponential(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.RateHetSites.split(' ')[5]) or float(val) < float(Variables.RateHetSites.split(' ')[4]):
                    val = np.random.exponential(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]), 1)[1:-1] + ' '
                RateHetSites_PE = '-a' + str(val) + ' '
            else:
                RateHetSites_PE = '-a' + str(np.random.exponential(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]), 1))[1:-1] + ' '
        elif Variables.RateHetSites.split(' ')[0] == 'gamma':
            if Variables.RateHetSites.split(' ')[3] == 't':
                val = np.random.gamma(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.RateHetSites.split(' ')[5]) or float(val) < float(Variables.RateHetSites.split(' ')[4]):
                    val = np.random.gamma(float(Variables.RateHetSites.split(' ')[1]),float(Variables.RateHetSites.split(' ')[2]), 1)[1:-1] + ' '
                RateHetSites_PE = '-a' + str(val) + ' '
            else:
                RateHetSites_PE = '-a' + str(np.random.gamma(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]),1))[1:-1] + ' '
        elif Variables.RateHetSites.split(' ')[0] == 'beta':
            if Variables.RateHetSites.split(' ')[3] == 't':
                val = np.random.exponential(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.RateHetSites.split(' ')[5]) or float(val) < float(Variables.RateHetSites.split(' ')[4]):
                    val = np.random.exponential(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]), 1)[1:-1] + ' '
                RateHetSites_PE = '-a' + str(val) + ' '
            else:
                RateHetSites_PE = '-a' + str(np.random.exponential(float(Variables.RateHetSites.split(' ')[1]), float(Variables.RateHetSites.split(' ')[2]),1))[1:-1] + ' '
    else:
        RateHetSites_PE = ''


    # Proportion of invariables sites
    if Variables.PropInvSites.split(' ')[0] != '':
        if Variables.PropInvSites.split(' ')[0] == 'fix':
            PropInvSites_PE = '-i' + str(Variables.PropInvSites.split(' ')[1]) + ' '
        elif Variables.PropInvSites.split(' ')[0] == 'uniform':
            Variables.PropInvSites_PE = '-i' + str(np.random.uniform(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]), 1)) + ' '
        elif Variables.PropInvSites.split(' ')[0] == 'normal':
            if Variables.PropInvSites.split(' ')[3] == 't':
                val = np.random.normal(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.PropInvSites.split(' ')[5]) or float(val) < float(Variables.PropInvSites.split(' ')[4]):
                    val = np.random.normal(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]), 1)[1:-1] + ' '
                PropInvSites_PE = '-i' + str(val) + ' '
            else:
                PropInvSites_PE = '-i' + str(np.random.normal(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]), 1))[1:-1] + ' '
        elif Variables.PropInvSites.split(' ')[0] == 'exp':
            if Variables.PropInvSites.split(' ')[3] == 't':
                val = np.random.exponential(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.PropInvSites.split(' ')[5]) or float(val) < float(Variables.PropInvSites.split(' ')[4]):
                    val = np.random.exponential(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]), 1)[1:-1] + ' '
                PropInvSites_PE = '-i' + str(val) + ' '
            else:
                PropInvSites_PE = '-i' + str(np.random.exponential(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]), 1))[1:-1] + ' '
        elif Variables.PropInvSites.split(' ')[0] == 'gamma':
            if Variables.PropInvSites.split(' ')[3] == 't':
                val = np.random.gamma(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.PropInvSites.split(' ')[5]) or float(val) < float(Variables.PropInvSites.split(' ')[4]):
                    val = np.random.gamma(float(Variables.PropInvSites.split(' ')[1]),float(Variables.PropInvSites.split(' ')[2]), 1)[1:-1] + ' '
                PropInvSites_PE = '-i' + str(val) + ' '
            else:
                PropInvSites_PE = '-i' + str(np.random.gamma(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]),1))[1:-1] + ' '
        elif Variables.PropInvSites.split(' ')[0] == 'beta':
            if Variables.PropInvSites.split(' ')[3] == 't':
                val = np.random.exponential(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.PropInvSites.split(' ')[5]) or float(val) < float(Variables.PropInvSites.split(' ')[4]):
                    val = np.random.exponential(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]), 1)[1:-1] + ' '
                PropInvSites_PE = '-i' + str(val) + ' '
            else:
                PropInvSites_PE = '-i' + str(np.random.exponential(float(Variables.PropInvSites.split(' ')[1]), float(Variables.PropInvSites.split(' ')[2]),1))[1:-1] + ' '
    else:
        PropInvSites_PE = ''


    # Template
    if Variables.Template != '':
        warnings.simplefilter('ignore', PDBConstructionWarning)
        Template_PE = '-zPop_evol.in '
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
        if Variables.GMRCA != '':
            GMRCA_PE = '-x' + Variables.GMRCA + ' '
        else:
            GMRCA_PE = ''


    # Information on the screen
    if Variables.ShowInformationScreen == 'No':
        ShowInformationScreen_PE = 0
    else:
        ShowInformationScreen_PE = 1


    # Create documents for simulations with PE
    if len(Variables.SubstitutionModel.split(' ')) < 2:
        sys.exit('Error!! You must select at least 2 substitution models')
    else:
        SubModel_PE=[]
        SubstitutionRate_List=[]
        SubModel_PE_SCS = []

        for x in Variables.SubstitutionModel.split(' '):
            if x != 'Fitness' and x != 'Neutral':
                with open('ProteinModelerABC_arguments.txt', 'a') as Prot_a:
                    for y in range(1, (Variables.NumberOfSimulations + 1)):
                        if Variables.SubstitutionModel.split(' ').index(x) == 0:
                            if Variables.SubstitutionRate.split(' ')[0] != '':
                                if Variables.SubstitutionRate.split(' ')[0] == 'fix':
                                    SubstitutionRate_PE = '-u' + str(Variables.SubstitutionRate.split(' ')[1]) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                                    SubstitutionRate_PE = '-u' + str(
                                        np.random.uniform(float(Variables.SubstitutionRate.split(' ')[1]),
                                                          float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'normal':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),
                                                               float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                   float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),
                                                             float(Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                     1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'exp':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                    float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                              1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                        float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                  float(Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                     1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),
                                                              float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                  float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),
                                                            float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'beta':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                    float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                              1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                        float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                  float(Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                     1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = ''
                                SubstitutionRate_List.append(SubstitutionRate_PE)
                            if Variables.Tree != '':
                                SubModel_PE_used = '-n1 ' + str(Treefile_PE) + str(AminoacidFrequencies_PE) + '-@' + str(x) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(y)
                                SubModel_PE.append(str(SubModel_PE_used))
                                Prot_a.write(SubModel_PE_used + '\n')
                            else:
                                SubModel_PE_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@' + str(x) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(y)
                                SubModel_PE.append(str(SubModel_PE_used))
                                Prot_a.write(SubModel_PE_used + '\n')
                        elif Variables.SubstitutionModel.split(' ').index(x) != 0:
                            w = y + ((int(Variables.SubstitutionModel.split(' ').index(x))) * Variables.NumberOfSimulations)
                            if Variables.SubstitutionRate.split(' ')[0] != '':
                                if Variables.SubstitutionRate.split(' ')[0] == 'fix':
                                    SubstitutionRate_PE = '-u' + str(Variables.SubstitutionRate.split(' ')[1]) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                                    SubstitutionRate_PE = '-u' + str(
                                        np.random.uniform(float(Variables.SubstitutionRate.split(' ')[1]),
                                                          float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'normal':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),
                                                               float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                   float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),
                                                             float(Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                     1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'exp':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                    float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                              1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                        float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                  float(Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                     1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),
                                                              float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                  float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(
                                            np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),
                                                            float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'beta':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                    float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                              1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(
                                                val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                                        float(Variables.SubstitutionRate.split(' ')[2]), 1)[
                                                  1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = ''
                                SubstitutionRate_List.append(SubstitutionRate_PE)

                            if Variables.Tree != '':
                                SubModel_PE_used = '-n1 ' + str(Treefile_PE) + str(AminoacidFrequencies_PE) + '-@' + str(x) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                                SubModel_PE.append(str(SubModel_PE_used))
                                Prot_a.write(SubModel_PE_used + '\n')
                            else:
                                SubModel_PE_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@' + str(x) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                                SubModel_PE.append(str(SubModel_PE_used))
                                Prot_a.write(SubModel_PE_used + '\n')
                    Prot_a.close()
            else:
                comprobar1 = 'Fitness' in Variables.SubstitutionModel.split(' ')
                comprobar2 = 'Neutral' in Variables.SubstitutionModel.split(' ')
                if comprobar1 == True and comprobar2 == True:
                    for k in range(0,2):
                        with open('ProteinModelerABC_S_arguments' + str(k) + '.txt', 'w') as Prot_Sa:
                            for y in range(1, (Variables.NumberOfSimulations + 1)):
                                w = y + ((len(Variables.SubstitutionModel.split(' ')) - 2 + k) * Variables.NumberOfSimulations)
                                if Variables.SubstitutionRate.split(' ')[0] != '':
                                    if Variables.SubstitutionRate.split(' ')[0] == 'fix':
                                        SubstitutionRate_PE = '-u' + str(Variables.SubstitutionRate.split(' ')[1]) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    elif Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                                        SubstitutionRate_PE = '-u' + str(np.random.uniform(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    elif Variables.SubstitutionRate.split(' ')[0] == 'normal':
                                        if Variables.SubstitutionRate.split(' ')[3] == 't':
                                            val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1)[1:-1]
                                            while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                                val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                            SubstitutionRate_PE = '-u' + str(val) + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                        else:
                                            SubstitutionRate_PE = '-u' + str(np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    elif Variables.SubstitutionRate.split(' ')[0] == 'exp':
                                        if Variables.SubstitutionRate.split(' ')[3] == 't':
                                            val = np.random.exponential(
                                                float(Variables.SubstitutionRate.split(' ')[1]),
                                                float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                            while float(val) > float(
                                                    Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                                val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                            SubstitutionRate_PE = '-u' + str(val) + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                        else:
                                            SubstitutionRate_PE = '-u' + str(
                                                np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    elif Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                        if Variables.SubstitutionRate.split(' ')[3] == 't':
                                            val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1)[1:-1]
                                            while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                                val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                            SubstitutionRate_PE = '-u' + str(val) + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                        else:
                                            SubstitutionRate_PE = '-u' + str(
                                                np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    elif Variables.SubstitutionRate.split(' ')[0] == 'beta':
                                        if Variables.SubstitutionRate.split(' ')[3] == 't':
                                            val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                            while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                                val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                            SubstitutionRate_PE = '-u' + str(val) + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                        else:
                                            SubstitutionRate_PE = '-u' + str(np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                            SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                else:
                                    SubstitutionRate_PE = ''
                                    SubstitutionRate_List.append(SubstitutionRate_PE)

                                if Variables.Tree != '':
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
                        for y in range(1, (Variables.NumberOfSimulations + 1)):
                            w = y + ((len(Variables.SubstitutionModel.split(' ')) -1) * Variables.NumberOfSimulations)
                            if Variables.SubstitutionRate.split(' ')[0] != '':
                                if Variables.SubstitutionRate.split(' ')[0] == 'fix':
                                    SubstitutionRate_PE = '-u' + str(Variables.SubstitutionRate.split(' ')[1]) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                                    SubstitutionRate_PE = '-u' + str(np.random.uniform(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'normal':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1)[1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'exp':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1)[1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'beta':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = ''
                                SubstitutionRate_List.append(SubstitutionRate_PE)

                            if Variables.Tree != '':
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
                        for y in range(1, (Variables.NumberOfSimulations + 1)):
                            w = y + ((len(Variables.SubstitutionModel.split(' ')) -1) * Variables.NumberOfSimulations)
                            if Variables.SubstitutionRate.split(' ')[0] != '':
                                if Variables.SubstitutionRate.split(' ')[0] == 'fix':
                                    SubstitutionRate_PE = '-u' + str(Variables.SubstitutionRate.split(' ')[1]) + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'uniform':
                                    SubstitutionRate_PE = '-u' + str(np.random.uniform(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'normal':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'exp':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'gamma':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                elif Variables.SubstitutionRate.split(' ')[0] == 'beta':
                                    if Variables.SubstitutionRate.split(' ')[3] == 't':
                                        val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                                        while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                                            val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                                        SubstitutionRate_PE = '-u' + str(val) + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                                    else:
                                        SubstitutionRate_PE = '-u' + str(np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
                                        SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                            else:
                                SubstitutionRate_PE = ''
                                SubstitutionRate_List.append(SubstitutionRate_PE)

                            if Variables.Tree != '':
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


    if Variables.Coalescent == 'Phylo':
        pass
    else:
        with open('PSimulations.txt', 'w') as thetafile:
            thetafile.write('SubstitutionRate,Theta\n')
            for x in range(0, len(SubstitutionRate_List)):
                thetafile.write(str(SubstitutionRate_List[x]) + ',' + str(2 * Variables.HaploidDiploid * Variables.PopulationSize * float(SubstitutionRate_List[x]) * Variables.NumAA) + '\n')
            thetafile.close()

