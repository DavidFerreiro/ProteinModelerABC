############################
suppressPackageStartupMessages(library(abc))
############################

#####################################################
################### ABC VARIABLES ###################
#####################################################
ABC_Method <- "rejection"
ABC_Tolerance <- 0.005
ABC_N_Iterations <- 100
#####################################################
#####################################################
#####################################################


#Path
address<-paste("Route to example9-SE/ABCOutputs", sep="")
setwd(address)
############################

#Load Simulations SS
matrix <- paste("SSSimulations.csv", sep="")
ThisName <- paste("SSSimulations",sep="")
FullSSmatrix <- read.csv(matrix, header=TRUE, sep=",")
col_names_matrixSS <- dimnames(FullSSmatrix)
col_names_matrixSS <- col_names_matrixSS[[2]]
ncol_matrixSS <- length(FullSSmatrix[1,])
nrow_matrixSS <- length(FullSSmatrix[,1])
############################

#Load Real Data SS
matrixReal <- paste("SSRealData.csv", sep="")
ThisNameReal <- paste("SSReal",sep="")
FullRealmatrix <- read.csv(matrixReal, header=TRUE, sep=",")
col_names_matrixRD <- dimnames(FullRealmatrix)
col_names_matrixRD <- col_names_matrixRD[[2]]
ncol_matrixRD <- length(FullSSmatrix[1,])
nrow_matrixRD <- length(FullSSmatrix[,1])
############################

# make vector with assignemt models
ModelsVector <- vector(,nrow_matrixSS)
#Assing models in the vector
j <- 0
for (j in 1:10000)
    {
    ModelsVector[j] <- "WAG"
    }
for (j in 10001:20000)
    {
    ModelsVector[j] <- "Fitness"
    }
for (j in 20001:30000)
    {
    ModelsVector[j] <- "Neutral"
    }
############################

figureName1<-paste("./Results_SS_Energy.pdf",sep="")
pdf(figureName1)
par(mfrow=c(1,2))
boxplot(FullSSmatrix[,"DGREM_Mean"]~ModelsVector, main="DGREM_Mean", xlab="Substitution Models", ylab="DGREM_Mean")
abline( h = FullRealmatrix[,"DGREM_Mean"], col = "blue", main="DGREM_Mean", xlab="Substitution Models", ylab="DGREM_Mean")
boxplot(FullSSmatrix[,"DGREM_sd"]~ModelsVector, main="DGREM_sd", xlab="Substitution Models", ylab="DGREM_sd")
abline( h = FullRealmatrix[,"DGREM_sd"], col = "blue", main="DGREM_sd", xlab="Substitution Models", ylab="DGREM_sd")
invisible(dev.off())

figureName1<-paste("./Results_SS_AAReplacements.pdf",sep="")
pdf(figureName1)
par(mfrow=c(2,3))
boxplot(FullSSmatrix[,"SegSites"]~ModelsVector, main="SegSites", xlab="Substitution Models", ylab="SegSites")
abline( h = FullRealmatrix[,"SegSites"], col = "blue", main="SegSites", xlab="Substitution Models", ylab="SegSites")
boxplot(FullSSmatrix[,"Grantham_mean_Position"]~ModelsVector, main="Grantham_mean_Position", xlab="Substitution Models", ylab="Grantham_mean_Position")
abline( h = FullRealmatrix[,"Grantham_mean_Position"], col = "blue", main="Grantham_mean_Position", xlab="Substitution Models", ylab="Grantham_mean_Position")
boxplot(FullSSmatrix[,"Grantham_sd_Position"]~ModelsVector, main="Grantham_sd_Position", xlab="Substitution Models", ylab="Grantham_sd_Position")
abline( h = FullRealmatrix[,"Grantham_sd_Position"], col = "blue", main="Grantham_sd_Position", xlab="Substitution Models", ylab="Grantham_sd_Position")
boxplot(FullSSmatrix[,"Grantham_sk_Position"]~ModelsVector, main="Grantham_sk_Position", xlab="Substitution Models", ylab="Grantham_sk_Position")
abline( h = FullRealmatrix[,"Grantham_sk_Position"], col = "blue", main="Grantham_sk_Position", xlab="Substitution Models", ylab="Grantham_sk_Position")
boxplot(FullSSmatrix[,"Grantham_ku_Position"]~ModelsVector, main="Grantham_ku_Position", xlab="Substitution Models", ylab="Grantham_ku_Position")
abline( h = FullRealmatrix[,"Grantham_ku_Position"], col = "blue", main="Grantham_ku_Position", xlab="Substitution Models", ylab="Grantham_ku_Position")
invisible(dev.off())
############################

#Create output file
unlink("Results_text.txt", recursive = FALSE)
############################

cv.modsel <- suppressWarnings(cv4postpr(ModelsVector, FullSSmatrix, nval=ABC_N_Iterations, tol=ABC_Tolerance, method=ABC_Method))
write(paste("Model selection with abc - cross validation - based on 100 samples",sep=""),"Results_text.txt",append=T)
try(capture.output(summary(cv.modsel), file = "Results_text.txt", append = TRUE), silent=TRUE)
write(paste("\n\n",sep=""),"Results_text.txt",append=T)
figureName1<-paste("./Results_ConfusionMatrix_100samp_",ThisName,".pdf",sep="")
pdf(figureName1)
try(plot(cv.modsel, names.arg=c("Fitness", "Neutral", "WAG")), silent=TRUE)
try(legend(x="topright", inset = c(-0.05, -0.16), legend = c("Fitness", "Neutral", "WAG"), xpd = TRUE, bty = "n", fill = c('grey25', 'grey50', 'grey75')), silent=TRUE)
invisible(dev.off())
############################

#Posterior probabilities of Real data
modsel.results <- postpr(FullRealmatrix, ModelsVector, FullSSmatrix, tol=ABC_Tolerance, method=ABC_Method)
if (ABC_Method == "rejection") {  write(paste("Model selection with abc - Real data\nEstimation – Posterior probability of every substitution model (rejection):",sep=""),"Results_text.txt",append=T)
  summary_output <- capture.output(summary(modsel.results))
  write(paste(summary_output[12:13], sep=""), "Results_text.txt",append=T)
  write(paste("\n\n",sep=""),"Results_text.txt",append=T)
}else{
  write(paste("Model selection with abc - Real data",sep=""),"Results_text.txt",append=T)
  summary_output <- capture.output(summary(modsel.results))
  write(paste(summary_output[11:14], sep=""), "Results_text.txt",append=T)
  write(paste("Estimation – Posterior probability of every substitution model (", ABC_Method,"):",sep=""), "Results_text.txt",append=T)
  write(paste(summary_output[23:24], sep=""), "Results_text.txt",append=T)
  write(paste("\n\n",sep=""),"Results_text.txt",append=T)
}
############################

#Goodness-of-fit of Real data
write(paste("Goodness of fit of Real data\n",sep=""),"Results_text.txt",append=T)
figureName1<-paste("./Histogram_GoodnessOfFit.pdf",sep="")
pdf(figureName1)
par(mfrow=c(1,3))
res.gfit.Fitness=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="Fitness",], statistic=mean, nb.replicate=100))
write(paste("   -Fitness:",sep=""),"Results_text.txt",append=T)
capture.output(summary(res.gfit.Fitness), file = "Results_text.txt", append = TRUE)
plot(res.gfit.Fitness, main="Histogram under Fitness")
res.gfit.Neutral=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="Neutral",], statistic=mean, nb.replicate=100))
write(paste("   -Neutral:",sep=""),"Results_text.txt",append=T)
capture.output(summary(res.gfit.Neutral), file = "Results_text.txt", append = TRUE)
plot(res.gfit.Neutral, main="Histogram under Neutral")
res.gfit.WAG=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="WAG",], statistic=mean, nb.replicate=100))
write(paste("   -WAG:",sep=""),"Results_text.txt",append=T)
capture.output(summary(res.gfit.WAG), file = "Results_text.txt", append = TRUE)
plot(res.gfit.WAG, main="Histogram under WAG")
invisible(dev.off())

############################

gfitpca = function(target, sumstat, index, cprob=0.1, xlim=NULL, ylim=NULL, ...){
    loc2plot = function(x, y, cprob, ...)
    {
        fit = locfit(~lp(x, y, nn=.2), maxk=200, mint=100, maxit=100)
        lev = sort(fitted(fit))[floor(cprob * length(x))]
        plot(fit, lev=lev, m=100, drawlabels=FALSE, ...)
        return (list(fit=fit, lev=lev))
    }

    # when target is a vector
    if (is.vector(target)){
        target = t( as.data.frame(target))
        if (is.data.frame(sumstat)){colnames(target)=names(sumstat)}
        if (is.matrix(sumstat)){colnames(target)=colnames(sumstat)}
    }

    # acp
    res.prcomp = prcomp(sumstat, scale=T, center=T)

    nmod = length(table(index))
    theindex = names(table(index))

    # plot
    figureName <- paste("PCA.pdf", sep="")
    pdf(figureName)
    par(mfrow=c(1, 1))
    if (is.null(xlim)){xlim = ylim}
    if (is.null(ylim)){ylim=xlim}

    if (! is.null(xlim)){
        plot(0, type="n", xlim=xlim, ylim=ylim, xlab="PC1", ylab="PC2")
    }
    for (i in 1:nmod){
        ind = index == theindex[i]
        if ((i == 1) & (is.null(xlim))){add=FALSE} else {add=TRUE}
        loc2plot(res.prcomp$x[ind, 1], res.prcomp$x[ind, 2], cprob, col = i, lty = 1, lwd = 2, add = add, ...)
    }

    # observed data
    points(predict(res.prcomp, target)[1], predict(res.prcomp, target)[2], col=1, cex=2, pch=3, lwd=2)

    # legend
    legend("topright", legend=theindex, cex=1.5, col=c(1: nmod), lty = 0, pch = 15)

    # saveplot
    invisible(dev.off())
}

############################
pca <- gfitpca(target=FullRealmatrix, sumstat=FullSSmatrix, index=ModelsVector, cprob=0.005)
#plot(pca, main="PCA")
#invisible(dev.off())
############################

############################
SS_All <- matrix(nrow=10000)
for (i in 1:7){
  for(p in 1:3){
    SS <- FullSSmatrix[,i][(10000 * (p-1) + 1):(10000 * p)]
    SS_All <- cbind(SS_All, SS)
  }
}

SS_All <- SS_All[,-1]

pdf(file=paste("Histograms_SStats.pdf",sep=""))
par(mfrow=c(4,3))
for (i in 1:7){
  for(p in 1:3){
    hist(SS_All[,(i * 3 + (-3 + p))], main=ModelsVector[10000 * p], xlab=col_names_matrixSS[i])
    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])
  }
}
invisible(dev.off())

