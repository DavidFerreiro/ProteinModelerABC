#!/usr/bin/env python
#import libraries
import os
import sys
import multiprocessing
import Variables
import LeerSettings
import Errores
import ChangeVariablesPE
import Functions
import time


if __name__ == "__main__":
    start_time = time.time()
    # Call variables
    Variables.init()

    LeerSettings.leer()

    Errores.err()

    ChangeVariablesPE.chavar()
    
    print('Starting simulations based on ' + str(Variables.SubstitutionModel) + ' amino acid substitution model')
    
    e_pool = multiprocessing.Pool(int(Variables.NumberOfProcessors))

    result = e_pool.map(Functions.E_Simu, range(int(Variables.Total_E_simu)), chunksize=1)

    e_pool.close()

    e_pool.join()
    
    print('Simulations based on ' + str(Variables.SubstitutionModel) + ' amino acid substitution models ended\n')
    
    if Variables.LongSM == 0:
        pass
    elif Variables.LongSM == 1:
        print('Starting simulations based on ' + str(Variables.StructuralSubstitutionModel) + ' amino acid substitution model')
    
        o_pool = multiprocessing.Pool(int(Variables.NumberOfProcessors))

        result = o_pool.map(Functions.S_Simu_1, range(int(Variables.NumberOfSimulations)), chunksize=1)

        o_pool.close()

        o_pool.join()
        
    elif Variables.LongSM == 2:
        print('Starting simulations based on ' + str(Variables.StructuralSubstitutionModel.split(' ')[0]) + ' amino acid substitution model')
        o_pool = multiprocessing.Pool(int(Variables.NumberOfProcessors))

        result = o_pool.map(Functions.S_Simu_F, range(int(Variables.NumberOfSimulations)), chunksize=1)

        o_pool.close()

        o_pool.join()

        print('Simulations based on ' + str(Variables.StructuralSubstitutionModel.split(' ')[0]) + ' amino acid substitution models ended\n')

        print('Starting simulations based on ' + str(Variables.StructuralSubstitutionModel.split(' ')[1]) + ' amino acid substitution model')

        o_pool = multiprocessing.Pool(int(Variables.NumberOfProcessors))

        result = o_pool.map(Functions.S_Simu_N, range(int(Variables.NumberOfSimulations)), chunksize=1)

        o_pool.close()

        o_pool.join()
    
    print('Simulations based on ' + str(Variables.StructuralSubstitutionModel.split(' ')[1]) + ' amino acid substitution models ended\n')

    Functions.CheckSimu()

    Functions.RealSummarySta(Variables.NameOfPhylipFile)

    Functions.Head()

    #Sumarios.SummarySta(0, 'Simu', Variables.SummaryStatistics)
    a_pool = multiprocessing.Pool(int(Variables.NumberOfProcessors))

    result = a_pool.map(Functions.SimuSummarySta, range(int(Variables.Total_simu)), chunksize=1)
    
    #elapsed_time = time.time() - start_time

    a_pool.close()

    a_pool.join()

    #print('\n' + str(elapsed_time))

    if Variables.SaveSimulations == 'Yes':
        os.system('tar -czf Simulations.tar.gz ./Simulations/ ')
        os.system('rm -r ./Simulations')
        os.system('rm -r ./Results')
        os.system('rm -r ./__pycache__')
        os.system('rm -r newalignment.fasta')
        os.system('rm -r Local_interactions.dat')
        os.system('rm -r REM.txt')
        os.system('rm -r E_loc_0.txt')
        os.system('rm -r seqGMRCA')
        os.system('rm -r ProtModel_arguments.txt')
        os.system('rm -r ProtModel_S_arguments0.txt')
        os.system('rm -r ProtModel_S_arguments1.txt')
        if os.path.exists('ProtModel_E-Repeat_arguments.txt') == True:
            os.system('rm -r ProtModel_E-Repeat_arguments.txt')
        if os.path.exists('ProtModel_F-Repeat_arguments.txt') == True:
            os.system('rm -r ProtModel_F-Repeat_arguments.txt')
        if os.path.exists('ProtModel_N-Repeat_arguments.txt') == True:
            os.system('rm -r ProtModel_N-Repeat_arguments.txt')
        os.system('rm -r Pop_evol.in')

    else:
        os.system('rm -r newalignment.fasta')
        os.system('rm -r ./Results')
        os.system('rm -r ./__pycache__')
        os.system('rm -r Local_interactions.dat')
        os.system('rm -r REM.txt')
        os.system('rm -r E_loc_0.txt')
        os.system('rm -r seqGMRCA')
        os.system('rm -r ProtModel_arguments.txt')
        os.system('rm -r ProtModel_S_arguments0.txt')
        os.system('rm -r ProtModel_S_arguments1.txt')
        if os.path.exists('ProtModel_E-Repeat_arguments.txt') == True:
            os.system('rm -r ProtModel_E-Repeat_arguments.txt')
        if os.path.exists('ProtModel_F-Repeat_arguments.txt') == True:
            os.system('rm -r ProtModel_F-Repeat_arguments.txt')
        if os.path.exists('ProtModel_N-Repeat_arguments.txt') == True:
            os.system('rm -r ProtModel_N-Repeat_arguments.txt')
        os.system('rm -r Pop_evol.in')

    with open("ABCAnalysis.r", "w") as P:
        P.write('############################\n')
        P.write('suppressPackageStartupMessages(library(abc))\n')
        P.write('############################\n\n')
        P.write('#####################################################\n')
        P.write('################### ABC VARIABLES ###################\n')
        P.write('#####################################################\n')
        P.write('ABC_Method <- "' + str(Variables.ABCMethod) + '"\n')
        P.write('ABC_Tolerance <- ' + str(Variables.ABCTolerance) + '\n')
        P.write('ABC_N_Iterations <- ' + str(Variables.ABCIterations) + '\n')
        P.write('#####################################################\n')
        P.write('#####################################################\n')
        P.write('#####################################################\n\n\n')
        P.write('#Path\n')
        P.write('address<-paste("' + os.getcwd() + '", sep="")\n')
        P.write('setwd(address)\n')
        P.write('############################\n\n')
        P.write('#Load Priors\n')
        P.write('matrix <- paste("PSimulations.txt", sep="")\n')
        P.write('ThisName <- paste("PSimulations",sep="")\n')
        P.write('FullPriormatrix <- read.table(matrix, header=TRUE, sep=",")\n')
        P.write('############################\n\n')
        P.write('#Load Simulations SS\n')
        P.write('matrix <- paste("SSSimulations.csv", sep="")\n')
        P.write('ThisName <- paste("SSSimulations",sep="")\n')
        P.write('FullSSmatrix <- read.csv(matrix, header=TRUE, sep=",")\n')
        P.write('col_names_matrixSS <- dimnames(FullSSmatrix)\n')
        P.write('col_names_matrixSS <- col_names_matrixSS[[2]]\n')
        P.write('ncol_matrixSS <- length(FullSSmatrix[1,])\n')
        P.write('nrow_matrixSS <- length(FullSSmatrix[,1])\n')
        P.write('############################\n\n')
        P.write('#Load Real Data SS\n')
        P.write('matrixReal <- paste("SSRealData.csv", sep="")\n')
        P.write('ThisNameReal <- paste("SSReal",sep="")\n')
        P.write('FullRealmatrix <- read.csv(matrixReal, header=TRUE, sep=",")\n')
        P.write('col_names_matrixRD <- dimnames(FullRealmatrix)\n')
        P.write('col_names_matrixRD <- col_names_matrixRD[[2]]\n')
        P.write('ncol_matrixRD <- length(FullSSmatrix[1,])\n')
        P.write('nrow_matrixRD <- length(FullSSmatrix[,1])\n')
        P.write('############################\n\n')
        P.write('# make vector with assignemt models\n')
        P.write('ModelsVector <- vector(,nrow_matrixSS)\n')
        P.write('#Assing models in the vector\n')
        P.write('j <- 0\n')
        Models = []
        if Variables.LongEM != 0: 
            for x in Variables.SubstitutionModel.split(' '):
                Models.append(x)
        if Variables.LongSM != 0: 
            for x in Variables.StructuralSubstitutionModel.split(' '):
                Models.append(x)
        for x in range(1, len(Models) + 1):
            P.write('for (j in ' + str((x - 1) * Variables.NumberOfSimulations + 1) + ':' + str(x * Variables.NumberOfSimulations) + ')\n')
            P.write('    {\n')
            P.write('    ModelsVector[j] <- "' + Models[x - 1] + '"\n')
            P.write('    }\n')

        SS_list = []
        for x in Variables.SummaryStatistics.split(' '):
            if x == '1':
                SS_list.append('DGREM_Mean')
            elif x == '2':
                SS_list.append('DGREM_sd')
            elif x == '3':
                SS_list.append('SegSites')
            elif x == '4':
                SS_list.append('Grantham_mean_Position')
            elif x == '5':
                SS_list.append('Grantham_sd_Position')
            elif x == '6':
                SS_list.append('Grantham_sk_Position')
            elif x == '7':
                SS_list.append('Grantham_ku_Position')

        P.write('############################\n\n')
        P.write('figureName1<-paste("./Histogram_Priors.pdf",sep="")\n')
        P.write('pdf(figureName1)\n')
        P.write('if (nrow_matrixSS < 1000)	{\n')
        P.write('	hist (FullPriormatrix[,1], breaks=nrow_matrixSS, main="Substitution Rate Prior", xlab="Substitution Rate")\n')
        P.write('	} else {\n')
        P.write('	hist (FullPriormatrix[,1], breaks=1000, main="Substitution Rate Prior", xlab="Substitution Rate")\n')
        P.write('	}\n')
        #P.write('invisible(dev.off())\n')
        #P.write('############################\n\n')
        #P.write('figureName1<-paste("./Histogram_PriorTheta.pdf",sep="")\n')
        #P.write('pdf(figureName1)\n')
        P.write('if (nrow_matrixSS < 1000)	{\n')
        P.write('	hist (FullPriormatrix[,2], breaks=nrow_matrixSS, main="Theta Prior", xlab="Theta")\n')
        P.write('	} else {\n')
        P.write('	hist (FullPriormatrix[,2], breaks=1000, main="Theta Prior", xlab="Theta")\n')
        P.write('	}\n')
        P.write('invisible(dev.off())\n')
        P.write('############################\n\n')

        if len(SS_list) == 7:
            #P.write('figureName<-paste("Results_SS_1",ThisName,".jpeg",sep="")\n')
            #P.write('jpeg(figureName, w=1200, h=900)\n')
            P.write('figureName1<-paste("./Results_SS_Energy.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(1,2))\n')
            else:
                pass
            for x in range(0,2):
                P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(SS_list[x]) + '")\n')
                P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(SS_list[x]) + '")\n')
            P.write('invisible(dev.off())\n\n')
            #P.write('figureName<-paste("Results_SS_2",ThisName,".jpeg",sep="")\n')
            #P.write('jpeg(figureName, w=1200, h=900)\n')
            P.write('figureName1<-paste("./Results_SS_AAReplacements.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(2,3))\n')
            else:
                pass
            for x in range(2,7):
                P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(SS_list[x]) + '")\n')
                P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(SS_list[x]) + '")\n')
            P.write('invisible(dev.off())\n')
        elif len(SS_list) < 7 and len(SS_list) > 4:
            P.write('figureName1<-paste("./Results_SS.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(2,3))\n')
            else:
                pass
            for x in range(0,len(SS_list)):
                P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(SS_list[x]) + '")\n')
                P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(SS_list[x]) + '")\n')
            P.write('invisible(dev.off())\n')
        elif len(SS_list) < 5 and len(SS_list) > 2:
            P.write('figureName1<-paste("./Results_SS.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(1,4))\n')
            else:
                pass
            for x in range(0,len(SS_list)):
                P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(SS_list[x]) + '")\n')
                P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(SS_list[x]) + '")\n')
            P.write('invisible(dev.off())\n')
        else:
            P.write('figureName1<-paste("./Results_SS.pdf",sep="")\n')
            P.write('pdf(figureName1)\n')
            if Variables.MultiPage == 'Yes':
                P.write('par(mfrow=c(1,2))\n')
            else:
                pass
            for x in range(0,len(SS_list)):
                P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(SS_list[x]) + '")\n')
                P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(SS_list[x]) + '")\n')
            P.write('invisible(dev.off())\n')
        P.write('############################\n\n')

        P.write('#Create output file\n')
        P.write('unlink("Results_text.txt", recursive = FALSE)\n')
        P.write('############################\n\n')
        #P.write('cv.modsel <- NULL\n')
        #P.write('attempt <- 1\n')
        #P.write('while(is.null(cv.modsel) && attempt <= 10){\n')
        #P.write('   attempt <- attempt + 1\n')
        #P.write('   try(cv.modsel <- suppressWarnings(cv4postpr(ModelsVector, FullSSmatrix, nval=ABC_N_Iterations, tol=ABC_Tolerance, method=ABC_Method)))\n')
        #P.write('}\n')
        P.write('cv.modsel <- suppressWarnings(cv4postpr(ModelsVector, FullSSmatrix, nval=ABC_N_Iterations, tol=ABC_Tolerance, method=ABC_Method))\n')
        P.write('write(paste("Model selection with abc - cross validation - based on ' + str(Variables.ABCIterations) + ' samples",sep=""),"Results_text.txt",append=T)\n')
        P.write('try(capture.output(summary(cv.modsel), file = "Results_text.txt", append = TRUE), silent=TRUE)\n')
        P.write('write(paste("\\n\\n",sep=""),"Results_text.txt",append=T)\n')
        P.write('figureName1<-paste("./Results_ConfusionMatrix_' + str(Variables.ABCIterations) + 'samp_",ThisName,".pdf",sep="")\n')
        P.write('pdf(figureName1)\n')
        #P.write('figureName<-paste("Results_ConfusionMatrix_' + str(Variables.ABCIterations) + 'samp_",ThisName,".jpeg",sep="")\n')
        #P.write('jpeg(figureName, w=1200, h=900)\n')
        Models.sort()
        Str_Models = ''
        for x in Models:
            Str_Models = Str_Models + '"' + str(x) + '"' + ', '
        Str_Models = Str_Models[:-2]
        P.write('try(plot(cv.modsel, names.arg=c(' + Str_Models + ')), silent=TRUE)\n')
        P.write('invisible(dev.off())\n')
        P.write('############################\n\n')

        P.write('#Posterior probabilities of Real data\n')
        P.write('modsel.results <- postpr(FullRealmatrix, ModelsVector, FullSSmatrix, tol=ABC_Tolerance, method=ABC_Method)\n')
        P.write('write(paste("Model selection with abc - Real data",sep=""),"Results_text.txt",append=T)\n')
        P.write('capture.output(summary(modsel.results), file = "Results_text.txt", append = TRUE)\n')
        P.write('write(paste("\\n\\n",sep=""),"Results_text.txt",append=T)\n')
        P.write('############################\n\n')

        P.write('#Goodness-of-fit of Real data\n')
        P.write('write(paste("Goodness of fit of Real data",sep=""),"Results_text.txt",append=T)\n')
        P.write('figureName1<-paste("./Histogram_GoodnessOfFit.pdf",sep="")\n')
        P.write('pdf(figureName1)\n')
        for x in Models:
            P.write('res.gfit.' + str(x) + '=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="' + str(x) + '",], statistic=mean, nb.replicate=' + str(Variables.ABCIterations) + '))\n')
            P.write('capture.output(summary(res.gfit.' + str(x) + '), file = "Results_text.txt", append = TRUE)\n')
            P.write('write(paste("\\n\\n",sep=""),"Results_text.txt",append=T)\n')
            #P.write('figureName<-paste("Hist_' + str(x) + '",ThisName,".jpeg",sep="")\n')
            #P.write('jpeg(figureName, w=1200, h=900)\n')
            P.write('plot(res.gfit.' + str(x) + ', main="Histogram under ' + str(x) + '")\n')
        P.write('invisible(dev.off())\n\n')
        P.write('############################\n\n')
        P.write('gfitpca = function(target, sumstat, index, cprob=0.1, xlim=NULL, ylim=NULL, ...){\n')
        P.write('    loc2plot = function(x, y, cprob, ...)\n')
        P.write('    {\n')
        P.write('        fit = locfit(~lp(x, y, nn=.2), maxk=200, mint=100, maxit=100)\n')
        P.write('        lev = sort(fitted(fit))[floor(cprob * length(x))]\n')
        P.write('        plot(fit, lev=lev, m=100, drawlabels=FALSE, ...)\n')
        P.write('        return (list(fit=fit, lev=lev))\n')
        P.write('    }\n\n')
        P.write('    # when target is a vector\n')
        P.write('    if (is.vector(target)){\n')
        P.write('        target = t( as.data.frame(target))\n')
        P.write('        if (is.data.frame(sumstat)){colnames(target)=names(sumstat)}\n')
        P.write('        if (is.matrix(sumstat)){colnames(target)=colnames(sumstat)}\n')
        P.write('    }\n\n')
        P.write('    # acp\n')
        P.write('    res.prcomp = prcomp(sumstat, scale=T, center=T)\n\n')
        P.write('    nmod = length(table(index))\n')
        P.write('    theindex = names(table(index))\n\n')
        P.write('    # plot\n')
        P.write('    figureName <- paste("PCA.pdf", sep="")\n')
        P.write('    pdf(figureName)\n')
        P.write('    par(mfrow=c(1, 1))\n')
        P.write('    if (is.null(xlim)){xlim = ylim}\n')
        P.write('    if (is.null(ylim)){ylim=xlim}\n\n')
        P.write('    if (! is.null(xlim)){\n')
        P.write('        plot(0, type="n", xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2")\n')
        P.write('    }\n')
        P.write('    for (i in 1:nmod){\n')
        P.write('        ind = index == theindex[i]\n')
        P.write('        if ((i == 1) & (is.null(xlim))){add=FALSE} else {add=TRUE}\n')
        P.write('        loc2plot(res.prcomp$x[ind, 1], res.prcomp$x[ind, 2], cprob, col = i, lty = 1, lwd = 2, add = add, ...)\n')
        P.write('    }\n\n')
        P.write('    # observed data\n')
        P.write('    points(predict(res.prcomp, target)[1], predict(res.prcomp, target)[2], col=1, cex=2, pch=3, lwd=2)\n\n')
        P.write('    # legend\n')
        P.write('    legend("topright", legend=theindex, cex=1.5, col=c(1: nmod), lty = 0, pch = 15)\n\n')
        P.write('    # saveplot\n')
        P.write('    invisible(dev.off())\n')
        P.write('}\n\n')

        P.write('############################\n')
        P.write('pca <- gfitpca(target=FullRealmatrix, sumstat=FullSSmatrix, index=ModelsVector, cprob=' + str(Variables.ABCTolerance) + ')\n')
        #P.write('figureName1<-paste("./PCA.pdf",sep="")\n')
        #P.write('pdf(figureName1)\n')
        #P.write('#figureName<-paste("PCA.jpeg",sep="")\n')
        #P.write('#jpeg(figureName, w=1200, h=900)\n')
        P.write('#plot(pca, main="PCA")\n')
        P.write('#invisible(dev.off())\n')
        P.write('############################\n\n')

        P.write('############################\n')
        P.write('SS_All <- matrix(nrow=' + str(Variables.NumberOfSimulations) + ')\n')
        P.write('for (i in 1:' + str(len(Variables.SummaryStatistics.split(' '))) + '){\n')
        P.write('  for(p in 1:3){\n')
        P.write('    SS <- FullSSmatrix[,i][(' + str(Variables.NumberOfSimulations) + ' * (p-1) + 1):(' + str(Variables.NumberOfSimulations) + ' * p)]\n')
        P.write('    SS_All <- cbind(SS_All, SS)\n')
        P.write('  }\n')
        P.write('}\n\n')
        P.write('SS_All <- SS_All[,-1]\n\n')
        P.write('pdf(file=paste("Histograms_SStats.pdf",sep=""))\n')
        P.write('for (i in 1:' + str(len(Variables.SummaryStatistics.split(' '))) + '){\n')
        P.write('  for(p in 1:3){\n')
        P.write('    hist(SS_All[,(i * 3 + (-3 + p))], main=ModelsVector[' + str(Variables.NumberOfSimulations) + ' * p], xlab=col_names_matrixSS[i])\n')
        P.write('    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
        P.write('  }\n')
        P.write('}\n')
        P.write('invisible(dev.off())\n\n')
        P.write('MinMax <- c()\n')
        P.write('for (i in 1:' + str(len(Variables.SummaryStatistics.split(' '))) + '){\n')
        P.write('  h <- min(FullSSmatrix[,i])\n')
        P.write('  MinMax <- c(MinMax, h)\n')
        P.write('  g <- max(FullSSmatrix[,i])\n')
        P.write('  MinMax <- c(MinMax, g)\n')
        P.write('}\n\n')
        P.write('SS_Abs <- matrix(nrow=' + str(Variables.NumberOfSimulations * len(Models)) + ')\n')
        P.write('SS_Tol <- matrix(nrow=ABC_N_Iterations)\n\n')
        P.write('for (i in 1:' + str(len(Variables.SummaryStatistics.split(' '))) + '){\n')
        P.write('  vec <- c()\n')
        P.write('  Number <- 0\n')
        P.write('  for(p in 1:' + str(Variables.NumberOfSimulations * len(Models)) + '){\n')
        P.write('    Number <- abs(FullSSmatrix[p,i] - FullRealmatrix[1,i])\n')
        P.write('    vec = c(vec, Number)\n')
        P.write('  }\n')
        P.write('  SS_Abs <- cbind(SS_Abs, FullSSmatrix[,i])\n')
        P.write('  SS_Abs <- cbind(SS_Abs, Abs=vec)\n')
        P.write('  SS_Abs <- cbind(SS_Abs, FullPriormatrix[,2])\n')
        P.write('}\n')
        P.write('SS_Abs <- SS_Abs[,-1]\n\n')
        P.write('for (i in 1:' + str(len(Variables.SummaryStatistics.split(' '))) + '){\n')
        P.write('  for(p in 1:3){\n')
        P.write('    SS <- SS_Abs[,(i * 3 - 2):(i * 3)][(' + str(Variables.NumberOfSimulations) + ' * (p-1) + 1):(' + str(Variables.NumberOfSimulations) + ' * p),]\n')
        P.write('    SS <- SS[order(SS[,2]),]\n')
        P.write('    SS_Tol <- cbind(SS_Tol, SS[1:ABC_N_Iterations, 1])\n')
        P.write('    SS_Tol <- cbind(SS_Tol, SS[1:ABC_N_Iterations, 3])\n')
        P.write('  }\n')
        P.write('}\n')
        P.write('SS_Tol <- SS_Tol[,-1]\n\n')
        #P.write('pdf(file=paste("SummaryStatistics2.pdf",sep=""))\n')
        #if Variables.MultiPage == 'Yes':
            #P.write('par(mfrow=c(1,3))\n')
        #else:
            #pass
        #P.write('for (i in 1:' + str(len(Variables.SummaryStatistics.split(' '))) + '){\n')
        #P.write('  for(p in 1:3){\n')
        #P.write('    plot(y=SS_Tol[,(i * 3 * 2 + (-7 + p*2))], x=SS_Tol[,(i * 3 * 2 + (-6 + p*2))], main=ModelsVector[' + str(Variables.NumberOfSimulations) + ' * p], ylab=col_names_matrixSS[i], xlab="Theta")\n')
        #P.write('    abline(h = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
        #P.write('  }\n')
        #P.write('}\n')
        #P.write('invisible(dev.off())\n\n')
        P.write('pdf(file=paste("Scaterplots_SStatsVSParams.pdf",sep=""))\n')
        if Variables.MultiPage == 'Yes':
            P.write('par(mfrow=c(1,3))\n')
        else:
            pass
        P.write('for (i in 1:' + str(len(Variables.SummaryStatistics.split(' '))) + '){\n')
        P.write('  for(p in 1:3){\n')
        P.write('    plot(y=SS_Tol[,(i * 3 * 2 + (-7 + p*2))], x=SS_Tol[,(i * 3 * 2 + (-6 + p*2))], main=ModelsVector[' + str(Variables.NumberOfSimulations) + ' * p], ylab=col_names_matrixSS[i], xlab="Theta", ylim=c(MinMax[i*2-1], MinMax[i*2]))\n')
        P.write('    abline(h = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
        P.write('  }\n')
        P.write('}\n')
        P.write('invisible(dev.off())\n')
        P.close()

    print("\nExecuting ABC analysis\n")

    #Execute R script
    #attempt = 0
    #while attempt <= 10:
    os.system('Rscript --vanilla ABCAnalysis.r')
    #attempt = attempt + 1
    
    if os.path.exists('ABCOutputs') == True:
        os.system('rm -r ABCOutputs')
    if os.path.exists('SimulationsOutputs') == True:
        os.system('rm -r SimulationsOutputs')
        
    os.system('mkdir ABCOutputs')
    os.system('mkdir SimulationsOutputs')
    
    os.system('mv Histograms_SStats.pdf ABCOutputs/')
    os.system('mv Scaterplots_SStatsVSParams.pdf ABCOutputs/')
    os.system('mv PCA.pdf ABCOutputs/')
    os.system('mv Histogram_GoodnessOfFit.pdf ABCOutputs/')
    os.system('mv Results_text.txt ABCOutputs/')
    os.system('mv Results_ConfusionMatrix_' + str(Variables.ABCIterations) + 'samp_SSSimulations.pdf ABCOutputs/')
    
    if os.path.exists('Results_SS_Energy.pdf') == True:
        os.system('mv Results_SS_Energy.pdf ABCOutputs/')
    if os.path.exists('Results_SS_AAReplacements.pdf') == True:
        os.system('mv Results_SS_AAReplacements.pdf ABCOutputs/')
    if os.path.exists('Results_SS.pdf') == True:
        os.system('mv Results_SS.pdf ABCOutputs/')
        
    os.system('mv Histogram_Priors.pdf ABCOutputs/')
    os.system('mv ABCAnalysis.r ABCOutputs/')
    os.system('cp PSimulations.txt ABCOutputs/')
    os.system('cp SSSimulations.csv ABCOutputs/')
    os.system('cp SSRealData.csv ABCOutputs/')
    
    if os.path.exists('Simulations.tar.gz') == True:
        os.system('mv Simulations.tar.gz SimulationsOutputs/')
    os.system('mv PSimulations.txt SimulationsOutputs/')
    os.system('mv SSSimulations.csv SimulationsOutputs/')
    os.system('mv SSRealData.csv SimulationsOutputs/')

    rut = str(os.getcwd())
    rut2 = rut + '/ABCOutputs'

    with open('./ABCOutputs/Results_text.txt') as f:
            with open('./ABCOutputs/ABCAnalysis.r', 'r') as a:
                with open('./ABCOutputs/ABCAnalysis2.r', 'w') as n:
                    for line in a:
                        line = str(line)
                        line = line.replace(rut, rut2)
                        n.write(line)
            os.system('rm -r ./ABCOutputs/ABCAnalysis.r')
            os.system('mv ./ABCOutputs/ABCAnalysis2.r ./ABCOutputs/ABCAnalysis.r')
            f.close()
            a.close()
            n.close()

    attempt = 1

    for x in range(0, 11):
        if attempt < 10:
            with open('./ABCOutputs/Results_text.txt') as f:
                if not 'Mean model posterior probabilities' in f.read():
                    attempt = attempt + 1
                    print('ABC estimation failed\n')
                    print('Starting another try: ' + str(x +1))
                    os.system('Rscript --vanilla ./ABCOutputs/ABCAnalysis.r')

                else:
                    if x == 0:
                        print("ABC analysis  has finished\n")
                        break
                    else:
                        print('Try ' + str(x) + '  succeeded')
                        print("ABC analysis  has finished\n")
                        break
        else:
            sys.exit("\n"
                     " ______                     _ \n"
                     "|  ____|                   | |\n"
                     "| |__   _ __ _ __ ___  _ __| |\n"
                     "|  __| | '__| '__/ _ \| '__| |\n"
                     "| |____| |  | | | (_) | |  |_|\n"
                     "|______|_|  |_|  \___/|_|  (_)\n\n"


                     "But dont panic!! You dont have to repeat everything, just go to the ./ABCOutputs and re execute ABCAnalysis.r file increasing the tolerance or changing the ABC method to rejection. Check manual for further information\n")
        # else:
        # print("ABC analysis  has finished/n")

    print('PROTMODEL has finished!!\n')

    elapsed_time = time.time() - start_time
    print('\n' + str(elapsed_time) + ' seconds')




