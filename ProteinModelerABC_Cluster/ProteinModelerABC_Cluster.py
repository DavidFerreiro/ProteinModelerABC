import os
import sys

#Bin PATH
sys.path.append('/Users/DavidFerreiro/Desktop/ProteinModelerABC.v1.1/Cluster_F')
ROOT = "'/Users/DavidFerreiro/Desktop/ProteinModelerABC.v1.1/Cluster_F'"

import bin.Scripts.Variables
import bin.Scripts.ReadSettings
import bin.Scripts.Errors
import bin.Scripts.ChangeVariablesPE
import bin.Scripts.Functions


# Call variables
bin.Scripts.ReadSettings.leer()
bin.Scripts.Errors.err()
bin.Scripts.ChangeVariablesPE.chavar()

#Read cluster information file to create the executable file
with open('Cluster_info.txt', 'r') as f:
    data = f.readlines()

sche = ''
nod = ''
mem = ''
tim = ''
mod = ''
for line in data:
    if line.startswith('SCHE'):
        sche = line.split('=')[1]
        sche = sche.split()[0]
    elif line.startswith('NODES'):
        nod = line.split('=')[1]
        nod = nod.split()[0]
        if nod != '#':
            nod = ' -N ' + str(nod)
        else:
            nod = ''
    elif line.startswith('MEM'):
        mem = line.split('=')[1]
        mem = mem.split()[0]
        if mem != '#':
            mem = ' --mem=' + str(mem)
        else:
            mem = ''
    elif line.startswith('TIME'):
        tim = line.split('=')[1]
        tim = tim.split()[0]
        if tim != '#':
            tim = ' -t ' + str(tim)
        else:
            tim = ''
    elif line.startswith('MOD'):
        mod = line.split('=')[1]
        mod = mod.split('#')[0].strip()

if sche == 'Slurm':
    sis = True
    mpi_w = 'srun'
    with open('launch_PMABC.sh', 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH' + str(mem) + str(nod) + ' -n ' + str(bin.Scripts.Variables.NumberOfProcessors) + str(tim) + '\n\n')
        f.write(str(mod) + '\n\n')
        f.write('STARTTIME_F=$(date +%s)\n')
        f.write('echo -e "\\n"\n\n')
        for g in range(0, len(bin.Scripts.Variables.SubstitutionModel.split(' '))):
            model = "'" + str(bin.Scripts.Variables.SubstitutionModel.split(' ')[g]) + "'"
            f.write('echo "_______________' + model[1:-1] + ' model_______________"\n')
            f.write('echo "Starting simulations based on ' + model + ' amino acid substitution model"\n')
            f.write(mpi_w + ' python -c "import sys; sys.path.append(' + ROOT + '); import bin.Scripts.Functions; bin.Scripts.Functions.Cluster_Simu(' + str(g * bin.Scripts.Variables.NumberOfSimulations) + ', ' + str((g + 1) * bin.Scripts.Variables.NumberOfSimulations) + ', ' + model + ')"\n')
            f.write(mpi_w + ' python -c "import sys; sys.path.append(' + ROOT + '); import bin.Scripts.Functions; bin.Scripts.Functions.CheckSimu(' + str(g * bin.Scripts.Variables.NumberOfSimulations) + ', ' + str((g + 1) * bin.Scripts.Variables.NumberOfSimulations) + ', ' + model + ')"\n')
            f.write('echo "Simulations based on ' + model[1:-1] + ' amino acid substitution model ended"\n')
            f.write('echo -e "\\n"\n\n')
            f.write('echo "Calculating summary statistics of ' + model[1:-1] + ' amino acid substitution model simulations"\n')
            f.write(mpi_w + ' python -c "import sys; sys.path.append(' + ROOT + '); import bin.Scripts.Functions; bin.Scripts.Functions.Cluster_SS(' + str(g * bin.Scripts.Variables.NumberOfSimulations) + ', ' + str((g + 1) * bin.Scripts.Variables.NumberOfSimulations) + ')"\n')
            f.write('echo "Summary statistics calculations of ' + model[1:-1] + ' amino acid substitution model simulations ended"\n')
            f.write('echo -e "\\n\\n"\n\n')
        f.write('echo "_____________________Executing ABC analysis_____________________"\n')
        f.write('python -c "import sys; sys.path.append(' + ROOT + '); import bin.Scripts.Functions; bin.Scripts.Functions.Order()"\n\n')
        f.write('ENDTIME_F=$(date +%s)\n')
        f.write('echo "ProteinModelerABC has finished!!"\n\n')
        f.write('echo "ProteinModelerABC finished after $(($ENDTIME_F - $STARTTIME_F)) seconds"\n\n')
        f.close()
elif sche == 'PBS Pro':
    sis = True
    mpi_w = 'mpirun -np ' + str(bin.Scripts.Variables.NumberOfProcessors)
    with open('launch_PMABC.pbs', 'w') as f:
        f.write('#PBS -l ' + str(mem)[2:] + '\n')
        f.write('#PBS -l nodes=' + str(bin.Scripts.Variables.NumberOfProcessors) + '\n')
        f.write('#PBS -l walltime= ' + str(tim)[3:] + '\n\n')
        f.write(str(mod) + '\n\n')
        for g in range(0, len(bin.Scripts.Variables.SubstitutionModel.split(' '))):
            model = "'" + str(bin.Scripts.Variables.SubstitutionModel.split(' ')[g]) + "'"
            f.write('echo "_______________' + model + ' model_______________"\n')
            f.write('echo "Starting simulations based on ' + model + ' amino acid substitution model"\n')
            f.write(mpi_w + ' python -c "import sys; sys.path.append(' + ROOT + '); import bin.Scripts.Functions; bin.Scripts.Functions.Cluster_Simu(' + str(g * bin.Scripts.Variables.NumberOfSimulations) + ', ' + str((g + 1) * bin.Scripts.Variables.NumberOfSimulations) + ', ' + model + ')"\n')
            f.write(mpi_w + ' python -c "import sys; sys.path.append(' + ROOT + '); import bin.Scripts.Functions; bin.Scripts.Functions.CheckSimu(' + str(g * bin.Scripts.Variables.NumberOfSimulations) + ', ' + str((g + 1) * bin.Scripts.Variables.NumberOfSimulations) + ', ' + model + ')"\n')
            f.write('echo "Simulations based on ' + model + ' amino acid substitution model ended"\n')
            f.write('echo -e "\\n"\n\n')
            f.write('echo "Calculating summary statistics of ' + model + ' amino acid substitution model simulations"\n')
            f.write(mpi_w + ' python -c "import sys; sys.path.append(' + ROOT + '); import bin.Scripts.Functions; bin.Scripts.Functions.Cluster_SS(' + str(g * bin.Scripts.Variables.NumberOfSimulations) + ', ' + str((g + 1) * bin.Scripts.Variables.NumberOfSimulations) + ')"\n')
            f.write('echo "Summary statistics calculations of ' + model + ' amino acid substitution model simulations ended"\n')
            f.write('echo -e "\\n\\n"\n\n')
        f.write('echo "_____________________Executing ABC analysis_____________________"\n')
        f.write('python -c "import sys; sys.path.append(' + ROOT + '); import bin.Scripts.Functions; bin.Scripts.Functions.Order()"\n\n')
        f.write('ENDTIME=$(date +%s)\n')
        f.write('echo "ProteinModelerABC has finished!!"\n\n')
        f.write('echo "ProteinModelerABC finished after $(($ENDTIME - $STARTTIME)) seconds"\n\n')
        f.close()

#Create ABC.r script
with open("ABCAnalysis.r", "w") as P:
    P.write('############################\n')
    P.write('suppressPackageStartupMessages(library(abc))\n')
    P.write('############################\n\n')
    P.write('#####################################################\n')
    P.write('################### ABC VARIABLES ###################\n')
    P.write('#####################################################\n')
    P.write('ABC_Method <- "' + str(bin.Scripts.Variables.ABCMethod) + '"\n')
    P.write('ABC_Tolerance <- ' + str(bin.Scripts.Variables.ABCTolerance) + '\n')
    P.write('ABC_N_Iterations <- ' + str(bin.Scripts.Variables.ABCIterations) + '\n')
    P.write('#####################################################\n')
    P.write('#####################################################\n')
    P.write('#####################################################\n\n\n')
    P.write('#Path\n')
    P.write('address<-paste("' + os.getcwd() + '/ABCOutputs", sep="")\n')
    P.write('setwd(address)\n')
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
    for x in bin.Scripts.Variables.SubstitutionModel.split(' '):
        Models.append(x)
    for x in range(1, len(Models) + 1):
        P.write('for (j in ' + str((x - 1) * bin.Scripts.Variables.NumberOfSimulations + 1) + ':' + str(
            x * bin.Scripts.Variables.NumberOfSimulations) + ')\n')
        P.write('    {\n')
        P.write('    ModelsVector[j] <- "' + Models[x - 1] + '"\n')
        P.write('    }\n')

    SS_list = []
    for x in bin.Scripts.Variables.SummaryStatistics.split(' '):
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

    if len(SS_list) == 7:
        # P.write('figureName<-paste("Results_SS_1",ThisName,".jpeg",sep="")\n')
        # P.write('jpeg(figureName, w=1200, h=900)\n')
        P.write('figureName1<-paste("./Results_SS_Energy.pdf",sep="")\n')
        P.write('pdf(figureName1)\n')
        if bin.Scripts.Variables.MultiPage == 'Yes':
            P.write('par(mfrow=c(1,2))\n')
        else:
            pass
        for x in range(0, 2):
            P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(
                SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
            P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(
                SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
        P.write('invisible(dev.off())\n\n')
        P.write('figureName1<-paste("./Results_SS_AAReplacements.pdf",sep="")\n')
        P.write('pdf(figureName1)\n')
        if bin.Scripts.Variables.MultiPage == 'Yes':
            P.write('par(mfrow=c(2,3))\n')
        else:
            pass
        for x in range(2, 7):
            P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(
                SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
            P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(
                SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
        P.write('invisible(dev.off())\n')
    elif len(SS_list) < 7 and len(SS_list) > 4:
        P.write('figureName1<-paste("./Results_SS.pdf",sep="")\n')
        P.write('pdf(figureName1)\n')
        if bin.Scripts.Variables.MultiPage == 'Yes':
            P.write('par(mfrow=c(2,3))\n')
        else:
            pass
        for x in range(0, len(SS_list)):
            P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(
                SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
            P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(
                SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
        P.write('invisible(dev.off())\n')
    elif len(SS_list) < 5 and len(SS_list) > 2:
        P.write('figureName1<-paste("./Results_SS.pdf",sep="")\n')
        P.write('pdf(figureName1)\n')
        if bin.Scripts.Variables.MultiPage == 'Yes':
            P.write('par(mfrow=c(2,2))\n')
        else:
            pass
        for x in range(0, len(SS_list)):
            P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(
                SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
            P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(
                SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
        P.write('invisible(dev.off())\n')
    else:
        P.write('figureName1<-paste("./Results_SS.pdf",sep="")\n')
        P.write('pdf(figureName1)\n')
        if bin.Scripts.Variables.MultiPage == 'Yes':
            P.write('par(mfrow=c(1,2))\n')
        else:
            pass
        for x in range(0, len(SS_list)):
            P.write('boxplot(FullSSmatrix[,"' + str(SS_list[x]) + '"]~ModelsVector, main="' + str(
                SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
            P.write('abline( h = FullRealmatrix[,"' + str(SS_list[x]) + '"], col = "blue", main="' + str(
                SS_list[x]) + '", xlab="Substitution Models", ylab="' + str(SS_list[x]) + '")\n')
        P.write('invisible(dev.off())\n')
    P.write('############################\n\n')

    P.write('#Create output file\n')
    P.write('unlink("Results_text.txt", recursive = FALSE)\n')
    P.write('############################\n\n')
    P.write(
        'cv.modsel <- suppressWarnings(cv4postpr(ModelsVector, FullSSmatrix, nval=ABC_N_Iterations, tol=ABC_Tolerance, method=ABC_Method))\n')
    P.write('write(paste("Model selection with abc - cross validation - based on ' + str(
        bin.Scripts.Variables.ABCIterations) + ' samples",sep=""),"Results_text.txt",append=T)\n')
    P.write('try(capture.output(summary(cv.modsel), file = "Results_text.txt", append = TRUE), silent=TRUE)\n')
    P.write('write(paste("\\n\\n",sep=""),"Results_text.txt",append=T)\n')
    P.write('figureName1<-paste("./Results_ConfusionMatrix_' + str(
        bin.Scripts.Variables.ABCIterations) + 'samp_",ThisName,".pdf",sep="")\n')
    P.write('pdf(figureName1)\n')
    Models.sort()
    Str_Models = ''
    for x in Models:
        Str_Models = Str_Models + '"' + str(x) + '"' + ', '
    Str_Models = Str_Models[:-2]
    P.write('try(plot(cv.modsel, names.arg=c(' + Str_Models + ')), silent=TRUE)\n')
    col = []
    for y in range(1, len(Models) + 1):
        col.append(int(100 / (len(Models) + 1) * y))
    colors = []
    for z in col:
        colors.append("grey" + str(z))
    P.write(
        'try(legend(x="topright", inset = c(-0.05, -0.16), legend = c(' + Str_Models + '), xpd = TRUE, bty = "n", fill = c(' + str(
            colors)[1:-1] + ')), silent=TRUE)\n')
    P.write('invisible(dev.off())\n')
    P.write('############################\n\n')
    P.write('#Posterior probabilities of Real data\n')
    P.write(
        'modsel.results <- postpr(FullRealmatrix, ModelsVector, FullSSmatrix, tol=ABC_Tolerance, method=ABC_Method)\n')
    P.write('if (ABC_Method == "rejection") {')
    P.write(
        '  write(paste("Model selection with abc - Real data\\nEstimation – Posterior probability of every substitution model (rejection):",sep=""),"Results_text.txt",append=T)\n')
    P.write('  summary_output <- capture.output(summary(modsel.results))\n')
    P.write('  write(paste(summary_output[12:13], sep=""), "Results_text.txt",append=T)\n')
    P.write('  write(paste("\\n\\n",sep=""),"Results_text.txt",append=T)\n')
    P.write('}else{\n')
    P.write('  write(paste("Model selection with abc - Real data",sep=""),"Results_text.txt",append=T)\n')
    P.write('  summary_output <- capture.output(summary(modsel.results))\n')
    P.write('  write(paste(summary_output[11:14], sep=""), "Results_text.txt",append=T)\n')
    P.write(
        '  write(paste("Estimation – Posterior probability of every substitution model (", ABC_Method,"):",sep=""), "Results_text.txt",append=T)\n')
    P.write('  write(paste(summary_output[23:24], sep=""), "Results_text.txt",append=T)\n')
    P.write('  write(paste("\\n\\n",sep=""),"Results_text.txt",append=T)\n')
    P.write('}\n')
    P.write('############################\n\n')
    P.write('#Goodness-of-fit of Real data\n')
    P.write('write(paste("Goodness of fit of Real data\\n",sep=""),"Results_text.txt",append=T)\n')
    P.write('figureName1<-paste("./Histogram_GoodnessOfFit.pdf",sep="")\n')
    P.write('pdf(figureName1)\n')
    if bin.Scripts.Variables.MultiPage == 'Yes':
        P.write('par(mfrow=c(1,' + str(len(Models)) + '))\n')
        for x in Models:
            P.write('res.gfit.' + str(
                x) + '=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="' + str(
                x) + '",], statistic=mean, nb.replicate=' + str(bin.Scripts.Variables.ABCIterations) + '))\n')
            P.write('write(paste("   -' + str(x) + ':",sep=""),"Results_text.txt",append=T)\n')
            P.write('capture.output(summary(res.gfit.' + str(x) + '), file = "Results_text.txt", append = TRUE)\n')
            P.write('plot(res.gfit.' + str(x) + ', main="Histogram under ' + str(x) + '")\n')
    else:
        for x in Models:
            P.write('res.gfit.' + str(
                x) + '=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="' + str(
                x) + '",], statistic=mean, nb.replicate=' + str(bin.Scripts.Variables.ABCIterations) + '))\n')
            P.write('write(paste("   -' + str(x) + ':",sep=""),"Results_text.txt",append=T)\n')
            P.write('capture.output(summary(res.gfit.' + str(x) + '), file = "Results_text.txt", append = TRUE)\n')
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
    P.write(
        '        loc2plot(res.prcomp$x[ind, 1], res.prcomp$x[ind, 2], cprob, col = i, lty = 1, lwd = 2, add = add, ...)\n')
    P.write('    }\n\n')
    P.write('    # observed data\n')
    P.write(
        '    points(predict(res.prcomp, target)[1], predict(res.prcomp, target)[2], col=1, cex=2, pch=3, lwd=2)\n\n')
    P.write('    # legend\n')
    P.write('    legend("topright", legend=theindex, cex=1.5, col=c(1: nmod), lty = 0, pch = 15)\n\n')
    P.write('    # saveplot\n')
    P.write('    invisible(dev.off())\n')
    P.write('}\n\n')

    P.write('############################\n')
    P.write('pca <- gfitpca(target=FullRealmatrix, sumstat=FullSSmatrix, index=ModelsVector, cprob=' + str(
        bin.Scripts.Variables.ABCTolerance) + ')\n')
    P.write('#plot(pca, main="PCA")\n')
    P.write('#invisible(dev.off())\n')
    P.write('############################\n\n')

    P.write('############################\n')
    P.write('SS_All <- matrix(nrow=' + str(bin.Scripts.Variables.NumberOfSimulations) + ')\n')
    P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
    P.write('  for(p in 1:' + str(len(Models)) + '){\n')
    P.write('    SS <- FullSSmatrix[,i][(' + str(bin.Scripts.Variables.NumberOfSimulations) + ' * (p-1) + 1):(' + str(
        bin.Scripts.Variables.NumberOfSimulations) + ' * p)]\n')
    P.write('    SS_All <- cbind(SS_All, SS)\n')
    P.write('  }\n')
    P.write('}\n\n')
    P.write('SS_All <- SS_All[,-1]\n\n')
    P.write('pdf(file=paste("Histograms_SStats.pdf",sep=""))\n')
    if bin.Scripts.Variables.MultiPage == 'Yes':
        if len(SS_list) == 7:
            P.write('par(mfrow=c(4,3))\n')
            P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
            P.write('  for(p in 1:' + str(len(Models)) + '){\n')
            P.write('    hist(SS_All[,(i * ' + str(len(Models)) + ' + (-' + str(
                len(Models)) + ' + p))], main=ModelsVector[' + str(
                bin.Scripts.Variables.NumberOfSimulations) + ' * p], xlab=col_names_matrixSS[i])\n')
            P.write('    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
            P.write('  }\n')
            P.write('}\n')
        elif len(SS_list) < 7 and len(SS_list) > 4:
            P.write('par(mfrow=c(3,3))\n')
            P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
            P.write('  for(p in 1:' + str(len(Models)) + '){\n')
            P.write('    hist(SS_All[,(i * ' + str(len(Models)) + ' + (-' + str(
                len(Models)) + ' + p))], main=ModelsVector[' + str(
                bin.Scripts.Variables.NumberOfSimulations) + ' * p], xlab=col_names_matrixSS[i])\n')
            P.write('    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
            P.write('  }\n')
            P.write('}\n')
        else:
            P.write('par(mfrow=c(2,3))\n')
            P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
            P.write('  for(p in 1:' + str(len(Models)) + '){\n')
            P.write('    hist(SS_All[,(i * ' + str(len(Models)) + ' + (-' + str(
                len(Models)) + ' + p))], main=ModelsVector[' + str(
                bin.Scripts.Variables.NumberOfSimulations) + ' * p], xlab=col_names_matrixSS[i])\n')
            P.write('    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
            P.write('  }\n')
            P.write('}\n')
    else:
        P.write('for (i in 1:' + str(len(bin.Scripts.Variables.SummaryStatistics.split(' '))) + '){\n')
        P.write('  for(p in 1:' + str(len(Models)) + '){\n')
        P.write('    hist(SS_All[,(i * ' + str(len(Models)) + ' + (-' + str(
            len(Models)) + ' + p))], main=ModelsVector[' + str(
            bin.Scripts.Variables.NumberOfSimulations) + ' * p], xlab=col_names_matrixSS[i])\n')
        P.write('    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])\n')
        P.write('  }\n')
        P.write('}\n')
    P.write('invisible(dev.off())\n\n')
    P.close()

# Calculate summary statistics of real data
print('Real data summary statistics calculation')
print(' Calculating ' + str(bin.Scripts.Variables.NameOfPhylipFile)[:-4] + ' summary statistics')
bin.Scripts.Functions.SScalc(bin.Scripts.Variables.NameOfPhylipFile, 'real')
print('Calculation of real data summary statistics ended\n\n')
