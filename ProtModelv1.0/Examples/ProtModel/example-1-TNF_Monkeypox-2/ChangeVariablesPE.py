import Variables
import os
import sys
from Bio import SeqIO
from random import randint, uniform
import numpy as np
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

def chavar():
    #Change variables to PE
    if Variables.NumSeq != '':
        if Variables.NumAA != '':
            SampleSize_SequenceLength_PE = '-s' + str(Variables.NumSeq) + ' ' + str(Variables.NumAA) + ' '
        else:
            SampleSize_SequenceLength_PE = '-s' + str(Variables.NumSeq) + ' '
    else:
        SampleSize_SequenceLength_PE = ''

    if Variables.PopulationSize != '':
        if Variables.HaploidDiploid!= '':
            PopSize_Hap_Dip_PE = '-e' + str(Variables.PopulationSize) + ' ' + str(Variables.HaploidDiploid) + ' '
        else:
            PopSize_Hap_Dip_PE = '-e' + str(Variables.PopulationSize) + ' '
    else:
        PopSize_Hap_Dip_PE = ''

    if Variables.DatedTips != '':
        DatedTips_PE = '-=' + str(Variables.DatedTips) + ' '
    else:
        DatedTips_PE = ''

    if Variables.GenerationTime != '':
        if Variables.GenerationTime[0] == 'fix':
            GenTime_PE = '-/' + str(Variables.GenerationTime[1]) + ' '
        elif Variables.GenerationTime[0] == 'uniform':
            GenTime_PE = '-/' + str(randint(Variables.GenerationTime[1],Variables.GenerationTime[2])) + ' '
    else:
        GenTime_PE = ''

    if Variables.GrowthRate != '':
        GrowthRate_PE = '-g' + Variables.GrowthRate + ' '
    else:
        GrowthRate_PE = ''

    if Variables.MigrationModel != '':
        MigrationModel_PE = '-q' + Variables.MigrationModel + ' '
    else:
        MigrationModel_PE = ''

    if Variables.MigrationRate != '':
        MigrationRate_PE = '-t' + Variables.MigrationRate + ' '
    else:
        MigrationRate_PE = ''

    if Variables.ConvergenceDemes != '':
        ConvergenceDemes_PE = '-%' + Variables.ConvergenceDemes + ' '
    else:
        ConvergenceDemes_PE = ''

    if Variables.SubstitutionRate.split(' ')[0] != '':
        if Variables.SubstitutionRate.split(' ')[0] == 'fix':
            SubstitutionRate_PE = '-u' + str(Variables.SubstitutionRate.split(' ')[1]) + ' '
        elif Variables.SubstitutionRate.split(' ')[0] == 'uniform':
            SubstitutionRate_PE = '-u' + str(np.random.uniform(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
        elif Variables.SubstitutionRate.split(' ')[0] == 'normal':
            if Variables.SubstitutionRate.split(' ')[3] == 't':
                val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                    val = np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                SubstitutionRate_PE = '-u' + str(val) + ' '
            else:
                SubstitutionRate_PE = '-u' + str(np.random.normal(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
        elif Variables.SubstitutionRate.split(' ')[0] == 'exp':
            if Variables.SubstitutionRate.split(' ')[3] == 't':
                val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                    val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                SubstitutionRate_PE = '-u' + str(val) + ' '
            else:
                SubstitutionRate_PE = '-u' + str(np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]), 1))[1:-1] + ' '
        elif Variables.SubstitutionRate.split(' ')[0] == 'gamma':
            if Variables.SubstitutionRate.split(' ')[3] == 't':
                val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                    val = np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]),float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                SubstitutionRate_PE = '-u' + str(val) + ' '
            else:
                SubstitutionRate_PE = '-u' + str(np.random.gamma(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
        elif Variables.SubstitutionRate.split(' ')[0] == 'beta':
            if Variables.SubstitutionRate.split(' ')[3] == 't':
                val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1]
                while float(val) > float(Variables.SubstitutionRate.split(' ')[5]) or float(val) < float(Variables.SubstitutionRate.split(' ')[4]):
                    val = np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]), 1)[1:-1] + ' '
                SubstitutionRate_PE = '-u' + str(val) + ' '
            else:
                SubstitutionRate_PE = '-u' + str(np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]), float(Variables.SubstitutionRate.split(' ')[2]),1))[1:-1] + ' '
    else:
        SubstitutionRate_PE = ''

    if Variables.AminoacidFrequencies != '':
        if Variables.AminoacidFrequencies == "Default":
            AminoacidFrequencies_PE = '-f20 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 '
        if Variables.AminoacidFrequencies.split(' ')[0] == "fix":
            AminoacidFrequencies_PE = '-f20 ' + str(Variables.AminoacidFrequencies.split(' ')[1]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[2]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[3]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[4]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[5]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[6]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[7]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[8]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[9]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[10]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[11]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[12]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[13]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[14]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[15]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[16]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[17]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[18]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[19]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[20]) + ' '
        if Variables.AminoacidFrequencies.split(' ')[0] == "dirichlet":
            AminoacidFrequencies_PE = '-f20 ' + str(Variables.AminoacidFrequencies.split(' ')[1]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[2]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[3]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[4]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[5]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[6]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[7]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[8]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[9]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[10]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[11]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[12]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[13]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[14]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[15]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[16]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[17]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[18]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[19]) + ' ' + str(Variables.AminoacidFrequencies.split(' ')[20]) + ' '
    else:
        AminoacidFrequencies_PE = ''

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


        if Variables.GMRCA != '':
            with open('ProtModel.out', 'a') as Error:
                Error.write('ERROR!! You can not select a template and GMRCA together\n')
                Error.close()
            sys.exit('ERROR!! You can not select a template and GMRCA together')
        else:
            GMRCA_Name = 'seqGMRCA'
            fp = open(GMRCA_Name, 'w')
            fp.write(str(PDB_Seq))
            fp.close()
            GMRCA_PE = '-x' + GMRCA_Name + ' '
            with open('ProtModel.out', 'a') as Error:
                Error.write('PDB sequence will be used as GMRCA instead your input sequences in SCS models\n\n')
                Error.close()
            print('PDB sequence will be used as GMRCA instead your input sequences in SCS models\n')
    else:
        Template_PE = ''
        if Variables.GMRCA != '':
            GMRCA_PE = '-x' + Variables.GMRCA + ' '
        else:
            GMRCA_PE = ''

    if Variables.ShowInformationScreen == 'No':
        ShowInformationScreen_PE = 0
    else:
        ShowInformationScreen_PE = 1

    with open('ProtModel_arguments.txt', 'w') as Prot_a:
        if Variables.SubstitutionModel == '':
            NSubsMo = 0
        else:
            w=0
            SubModel_PE=[]
            SubstitutionRate_List=[]
            for x in Variables.SubstitutionModel.split(' '):
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
                    #globals()['Model_%s' % x + '_Simu_%s' % y] = 'n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@' + str(x) + ' ' + '-bsequences -c1 1 0 -y0 -:' + str(y)
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
                                    SubstitutionRate_PE = '-u' + str(
                                        np.random.exponential(float(Variables.SubstitutionRate.split(' ')[1]),
                                                              float(Variables.SubstitutionRate.split(' ')[2]), 1))[
                                                                 1:-1] + ' '
                                    SubstitutionRate_List.append(SubstitutionRate_PE[2:])
                        else:
                            SubstitutionRate_PE = ''
                            SubstitutionRate_List.append(SubstitutionRate_PE)

                        SubModel_PE_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@' + str(x) + ' ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                        SubModel_PE.append(str(SubModel_PE_used))
                        Prot_a.write(SubModel_PE_used + '\n')
        Prot_a.close()

    with open('ProtModel_S_arguments.txt', 'w') as Prot_Sa:
        w=0
        SubModel_PE_SCS=[]
        if Variables.StructuralSubstitutionModel == '':
            Prot_Sa.close()
        if len(Variables.StructuralSubstitutionModel.split(' ')) == 1:
            Prot_Sa.close()
            os.system('rm -r ProtModel_S_arguments.txt')
            with open('ProtModel_S_arguments0' + '.txt', 'w') as Prot_Sa:
                for x in range(0, 1):
                    for y in range(1, (Variables.NumberOfSimulations + 1)):
                        w = y + ((x + NSubsMo) * Variables.NumberOfSimulations)
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

                        SubModel_PE_SCS_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@JTT ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(Template_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                        SubModel_PE_SCS.append(str(SubModel_PE_SCS_used))
                        Prot_Sa.write(SubModel_PE_SCS_used + '\n')
                Prot_Sa.close()
        if len(Variables.StructuralSubstitutionModel.split(' ')) == 2:
            Prot_Sa.close()
            os.system('rm -r ProtModel_S_arguments.txt')
            for x in range(0, 2):
                with open('ProtModel_S_arguments' + str(x) + '.txt', 'w') as Prot_Sa:
                    for y in range(1, (Variables.NumberOfSimulations + 1)):
                        w = y + ((x + len(Variables.SubstitutionModel.split(' '))) * Variables.NumberOfSimulations)
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

                        SubModel_PE_SCS_used = '-n1 ' + str(SampleSize_SequenceLength_PE) + str(PopSize_Hap_Dip_PE) + str(DatedTips_PE) + str(GenTime_PE) + str(GrowthRate_PE) + str(MigrationModel_PE) + str(MigrationRate_PE) + str(ConvergenceDemes_PE) + str(SubstitutionRate_PE) + str(AminoacidFrequencies_PE) + '-@JTT ' + str(RateHetSites_PE) + str(PropInvSites_PE) + str(Template_PE) + str(GMRCA_PE) + '-bsequences -c1 1 0 -y' + str(ShowInformationScreen_PE) + ' -:' + str(w)
                        SubModel_PE_SCS.append(str(SubModel_PE_SCS_used))
                        Prot_Sa.write(SubModel_PE_SCS_used + '\n')
        Prot_Sa.close()

    with open('PSimulations.txt', 'w') as thetafile:
        thetafile.write('SubstitutionRate,Theta\n')
        for x in range(0, len(SubstitutionRate_List)):
            thetafile.write(str(SubstitutionRate_List[x]) + ',' + str(2 * Variables.HaploidDiploid * Variables.PopulationSize * float(SubstitutionRate_List[x]) * Variables.NumAA) + '\n')
        thetafile.close()