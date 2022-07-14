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
address<-paste("/Users route to /ABCOutputs", sep="")
setwd(address)
############################

#Load Priors
matrix <- paste("PSimulations.txt", sep="")
ThisName <- paste("PSimulations",sep="")
FullPriormatrix <- read.table(matrix, header=TRUE, sep=",")
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
    ModelsVector[j] <- "JTT"
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

figureName1<-paste("./Histogram_Priors.pdf",sep="")
pdf(figureName1)
if (nrow_matrixSS < 1000)	{
	hist (FullPriormatrix[,1], breaks=nrow_matrixSS, main="Substitution Rate Prior", xlab="Substitution Rate")
	} else {
	hist (FullPriormatrix[,1], breaks=1000, main="Substitution Rate Prior", xlab="Substitution Rate")
	}
if (nrow_matrixSS < 1000)	{
	hist (FullPriormatrix[,2], breaks=nrow_matrixSS, main="Theta Prior", xlab="Theta")
	} else {
	hist (FullPriormatrix[,2], breaks=1000, main="Theta Prior", xlab="Theta")
	}
invisible(dev.off())
############################

figureName1<-paste("./Results_SS_Energy.pdf",sep="")
pdf(figureName1)
par(mfrow=c(1,2))
boxplot(FullSSmatrix[,"DGREM_Mean"]~ModelsVector, main="DGREM_Mean")
abline( h = FullRealmatrix[,"DGREM_Mean"], col = "blue", main="DGREM_Mean")
boxplot(FullSSmatrix[,"DGREM_sd"]~ModelsVector, main="DGREM_sd")
abline( h = FullRealmatrix[,"DGREM_sd"], col = "blue", main="DGREM_sd")
invisible(dev.off())

figureName1<-paste("./Results_SS_AAReplacements.pdf",sep="")
pdf(figureName1)
par(mfrow=c(2,3))
boxplot(FullSSmatrix[,"SegSites"]~ModelsVector, main="SegSites")
abline( h = FullRealmatrix[,"SegSites"], col = "blue", main="SegSites")
boxplot(FullSSmatrix[,"Grantham_mean_Position"]~ModelsVector, main="Grantham_mean_Position")
abline( h = FullRealmatrix[,"Grantham_mean_Position"], col = "blue", main="Grantham_mean_Position")
boxplot(FullSSmatrix[,"Grantham_sd_Position"]~ModelsVector, main="Grantham_sd_Position")
abline( h = FullRealmatrix[,"Grantham_sd_Position"], col = "blue", main="Grantham_sd_Position")
boxplot(FullSSmatrix[,"Grantham_sk_Position"]~ModelsVector, main="Grantham_sk_Position")
abline( h = FullRealmatrix[,"Grantham_sk_Position"], col = "blue", main="Grantham_sk_Position")
boxplot(FullSSmatrix[,"Grantham_ku_Position"]~ModelsVector, main="Grantham_ku_Position")
abline( h = FullRealmatrix[,"Grantham_ku_Position"], col = "blue", main="Grantham_ku_Position")
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
try(plot(cv.modsel, names.arg=c("Fitness", "JTT", "Neutral")), silent=TRUE)
invisible(dev.off())
############################

#Posterior probabilities of Real data
modsel.results <- postpr(FullRealmatrix, ModelsVector, FullSSmatrix, tol=ABC_Tolerance, method=ABC_Method)
write(paste("Model selection with abc - Real data",sep=""),"Results_text.txt",append=T)
capture.output(summary(modsel.results), file = "Results_text.txt", append = TRUE)
write(paste("\n\n",sep=""),"Results_text.txt",append=T)
############################

#Goodness-of-fit of Real data
write(paste("Goodness of fit of Real data\n",sep=""),"Results_text.txt",append=T)
figureName1<-paste("./Histogram_GoodnessOfFit.pdf",sep="")
pdf(figureName1)
res.gfit.Fitness=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="Fitness",], statistic=mean, nb.replicate=100))
write(paste("   -Fitness:",sep=""),"Results_text.txt",append=T)
capture.output(summary(res.gfit.Fitness), file = "Results_text.txt", append = TRUE)
plot(res.gfit.Fitness, main="Histogram under Fitness")
res.gfit.JTT=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="JTT",], statistic=mean, nb.replicate=100))
write(paste("   -JTT:",sep=""),"Results_text.txt",append=T)
capture.output(summary(res.gfit.JTT), file = "Results_text.txt", append = TRUE)
plot(res.gfit.JTT, main="Histogram under JTT")
res.gfit.Neutral=suppressWarnings(gfit(target=FullRealmatrix, sumstat=FullSSmatrix[ModelsVector=="Neutral",], statistic=mean, nb.replicate=100))
write(paste("   -Neutral:",sep=""),"Results_text.txt",append=T)
capture.output(summary(res.gfit.Neutral), file = "Results_text.txt", append = TRUE)
plot(res.gfit.Neutral, main="Histogram under Neutral")
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
for (i in 1:7){
  for(p in 1:3){
    hist(SS_All[,(i * 3 + (-3 + p))], main=ModelsVector[10000 * p], xlab=col_names_matrixSS[i])
    abline(v = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])
  }
}
invisible(dev.off())

MinMax <- c()
for (i in 1:7){
  h <- min(FullSSmatrix[,i])
  MinMax <- c(MinMax, h)
  g <- max(FullSSmatrix[,i])
  MinMax <- c(MinMax, g)
}

SS_Abs <- matrix(nrow=30000)
SS_Tol <- matrix(nrow=ABC_N_Iterations)

for (i in 1:7){
  vec <- c()
  Number <- 0
  for(p in 1:30000){
    Number <- abs(FullSSmatrix[p,i] - FullRealmatrix[1,i])
    vec = c(vec, Number)
  }
  SS_Abs <- cbind(SS_Abs, FullSSmatrix[,i])
  SS_Abs <- cbind(SS_Abs, Abs=vec)
  SS_Abs <- cbind(SS_Abs, FullPriormatrix[,2])
}
SS_Abs <- SS_Abs[,-1]

for (i in 1:7){
  for(p in 1:3){
    SS <- SS_Abs[,(i * 3 - 2):(i * 3)][(10000 * (p-1) + 1):(10000 * p),]
    SS <- SS[order(SS[,2]),]
    SS_Tol <- cbind(SS_Tol, SS[1:ABC_N_Iterations, 1])
    SS_Tol <- cbind(SS_Tol, SS[1:ABC_N_Iterations, 3])
  }
}
SS_Tol <- SS_Tol[,-1]

pdf(file=paste("Scaterplots_SStatsVSParams.pdf",sep=""))
par(mfrow=c(1,3))
for (i in 1:7){
  for(p in 1:3){
    plot(y=SS_Tol[,(i * 3 * 2 + (-7 + p*2))], x=SS_Tol[,(i * 3 * 2 + (-6 + p*2))], main=ModelsVector[10000 * p], ylab=col_names_matrixSS[i], xlab="Theta", ylim=c(MinMax[i*2-1], MinMax[i*2]))
    abline(h = FullRealmatrix[,i], col = "blue", main=col_names_matrixSS[i])
  }
}
invisible(dev.off())
