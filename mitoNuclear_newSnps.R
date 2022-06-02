rm(list=ls())

#setwd("/Volumes/thirdeyesqueegee/Dobslab Dropbox/Adam Dobson/mitoNuclear")
setwd("Dobslab Dropbox/Adam Dobson/mitoNuclear/")
# library(remotes)

 #install_github("nyiuab/NBZIMM", force=T, build_vignettes=F)
library(ggsci)
library(wesanderson)
library(gplots)
library(superheat)
library(NBZIMM)
library(irr)
library(pscl)
library(performance)
library(rsq)
library(boot)
library(survival)
library(insight)
library(ggplot2)
library(dplyr)
library(MASS)
#library(compute.es)
library(MuMIn)
library(lme4)
library(car)
library(emmeans)
library(piecewiseSEM)
#library(glmmADMB)
library(fitdistrplus)
#library(logspline)
#library(R2admb)
library(blmeco)
library(pheatmap)
library(effectsize)
library(glmmTMB)
library(see)
#library(tidyverse)
library(viridis)
library(scales)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(coxme)
library(snow)

	#set na action & contrasts
options(na.action="na.fail", contrasts=c("contr.sum", "contr.poly"))

	#make cluster
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 2), type = clusterType))

	#conversion factor for inches to mm
mmToInch <- 0.03937008

plotEtas_fixed <- function(x){
	ahnova <- anova(x)
	leng <- nrow(ahnova)
	effs <- joint_tests(x)[,4]
	dEffs <- ahnova$Df[2:leng]
	dEffError <- ahnova[2:leng,3]
	eta2 <- F_to_eta2(f=effs, df=dEffs, df_error=dEffError)
	eta2$Parameter <- rownames(ahnova)[2:leng]
	p <- plot(eta2)
	p + 
		scale_color_manual(values=c("hotpink", "hotpink")) + 
		theme_bw() +
		xlim(0,1) +
		theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90))

}

plotEtas_mixed <- function(x){
	ahnova <- joint_tests(x)
	effs <- ahnova[,4]
	dEffs <- ahnova[,2]
	eta2 <- F_to_eta2(f=effs, df=dEffs, df_error=df.residual(x))
	eta2$Parameter <- ahnova[,1]
	p <- plot(eta2)
	p + 
	scale_color_manual(values=c("hotpink", "hotpink")) + 
	theme_bw() +
	xlim(0,1)	+
	theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90))
}

#plotEtas_mixed <- function(x){
#	ahnova <- joint_tests(x)
#	effs <- ahnova[,4]
#	dEffs <- ahnova[,2]
#	eta2 <- F_to_eta2(f=effs, df=dEffs, df_error=df.residual(x))
#	eta2$Parameter <- ahnova[,1]
#	plot(eta2)
#}

bubbleFC_process <- function(x, sigsOnly, withinDiet, withinGenoOnly, vsSYAonly, pThresh, meanCentre, eaaVSmar, doubleDip){
	
	#create a response index: 
	#a signed, logged version of the absolute FCs
signs <- sign(x$fc)
absolutes <- abs(x$fc)
logs <- log(absolutes+1)
x$RI <- logs * signs
	
x$dietTraitWhen <- with(x, factor(paste(dietY, trait, when, sep="")))
	#mean-centre / trait / diet
if(meanCentre){mns <- aggregate(RI ~ dietTraitWhen, x, mean)
colnames(mns)[2] <- "RImean"
x <- merge(x, mns, by="dietTraitWhen")
x$RInorm <- x$RI - x$RImean}else{x$RInorm <- x$RI}

	#normalise to max / trait / diet
maxes <- aggregate(RInorm ~ dietTraitWhen, x, function(xx){max(abs(xx))})
colnames(maxes)[2] <- "RImax"
x <- merge(x, maxes, by="dietTraitWhen")
x$RInorm <- x$RInorm / x$RImax
	
		#if within-diet comparisons are not of interest
dim(x)
if(withinDiet==F){x <- subset(x, dietX!=dietY)}
dim(x)

	#EAA:MAR comparison is not informative
if(eaaVSmar==F){
x <- subset(x, !(dietX=="eaa" & dietY=="mar"))
dim(x)
x <- subset(x, !(dietX=="mar" & dietY=="eaa"))
dim(x)}

	#if we want to compare only within genotype (ie AA to AA, BB to BB)
if(withinGenoOnly==T){x <- subset(x, substr(genoY, 1, 2)==substr(genoX, 1, 2))}
dim(x)

	#if we want to not double-dip out comparisons
if(doubleDip==F){
x <- subset(x, !(dietY=="sya" & dietX %in% c("eaa", "mar")))
x <- subset(x, !(dietY=="mar" & dietX=="eaa"))
}

	#if we only comparison of SYA to experimental diets
dim(x)
if(vsSYAonly==T){x <- subset(x, dietX=="sya")}
dim(x)

x <- droplevels(x)

	#correct for multiple testing
	#and take -log10 of adjusted values
x$padj <- p.adjust(x$p, "fdr")
x$neglog <- -log10(x$padj)
dim(x)

	#add a variable for p < 0.05
x$sig <- x$padj <= pThresh

	#remove non-significant comparisons,if desired
if(sigsOnly){x <- droplevels(subset(x, padj<=pThresh))}

x
}

bubbleFC <- function(CELLS, trait, when, exponentiate, invLogit, mirror){

		#matrix of adjusted P-values, adjusted only for number of non-redundant comparisons
	pvs <- pwpm(CELLS, diffs=F, adjust="none", flip=T)
	diag(pvs) <- NA
	pvs <- matrix(as.numeric(gsub("<", "0", pvs)), nrow=nrow(pvs), dimnames=dimnames(pvs))	#so entries of 0.0001 are now p≤0.0001
	pvs <- matrix(pvs, nrow=nrow(pvs), ncol=ncol(pvs), dimnames=dimnames(pvs))
	if(mirror){pvs[upper.tri(pvs)] <- t(pvs)[upper.tri(t(pvs))]}	#mirror

		#calculate fold-changes
		#ugly loop alert - strange behaviour emerged with apply so playing it safe.
	cells <- as.data.frame(CELLS)	
	cells$ind <- paste(cells[,1], cells[,2])
	if(exponentiate){cells$emmean <- exp(cells$emmean)}
	if(invLogit){cells$emmean <- inv.logit(cells$emmean)}
	fcs <- matrix(nrow=nrow(pvs), ncol=ncol(pvs), dimnames=dimnames(pvs))
		for(i in colnames(fcs)){
			for(ii in rownames(fcs)){
				if(i != ii){
				effsee <- (cells$emmean[cells$ind==ii] - cells$emmean[cells$ind==i]) / cells$emmean[cells$ind==i]
				fcs[i,ii] <- effsee}
			}
		}
		
	if(!mirror){fcs[upper.tri(fcs)] <- NA}	

	xx <- data.frame(fc=as.numeric(fcs), p=as.numeric(pvs), x=rownames(fcs)[row(fcs)], y=colnames(fcs)[col(fcs)])

	xx <- subset(xx, complete.cases(xx))
	xx$trait <- rep(trait, nrow(xx))
	xx$when <- rep(when, nrow(xx))
	xx$dietX <- substr(xx$x, 1, 3)
	xx$dietY <- substr(xx$y, 1, 3)
	xx$genoX <- substr(xx$x, 5, 7)
	xx$genoY <- substr(xx$y, 5, 7)
	xx
	}

colChoose <- function(x){
	ifelse(x=="eaa", "#E52521", ifelse(x=="mar", "#3CB7EB", "#010101"))
	}


plotDat <- function(y, x, z, datCex, meanLwd, meanCex, xoffset, jitterFactor, datsAlpha, seed, hgt, wdt, fn, meanLineCol, meanLineLty, wideVar, exponent, invlogit, ...){

	set.seed(seed)

	if(exponent==T){z <- transform(z, yvar=exp(yvar), LCL=exp(LCL), UCL=exp(UCL))}
	if(invlogit==T){z <- transform(z, yvar=inv.logit(yvar), LCL=inv.logit(LCL), UCL=inv.logit(UCL))}
	zLines <- reshape(z[,1:3], idvar=wideVar, timevar="diet", direction="wide")
	rownames(zLines) <- zLines[,1]
	zLines <- t(zLines[,2:4])

	LCL <- reshape(z[,c(1,2,6)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(LCL) <- LCL[,1]
	LCL <- t(LCL[,2:4])

	UCL <- reshape(z[,c(1,2,7)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(UCL) <- UCL[,1]
	UCL <- t(UCL[,2:4])

	durts <- reshape(z[,c(1,2,1)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(durts) <- durts[,1]
	durts <- t(durts[,2:4])

	xvals <- t(matrix(c(rep(1:3, 3), 
						rep(4:6, 3), 
						rep(7:9, 3), 
						rep(10:12, 3)), 
						ncol=3, byrow=T))+xoffset

	xvals <- jitter(xvals, jitterFactor)

	par(bty="n", las=2, lwd=1)
	plot(y ~ x, xlab="", ...)

	matlines(y=zLines, x=xvals, col=meanLineCol, lty=meanLineLty, lwd=meanLwd)

	arrows(x0=xvals, x1=xvals, y0=LCL, y1=UCL, length=0, col=colChoose(durts), lwd=meanLwd)
	points(y=zLines, x=xvals, bg=colChoose(durts), col="white", pch=21, cex=meanCex, lwd=meanCex)
	points(jitter(as.numeric(x), factor=jitterFactor)-xoffset, y, bg=alpha(colChoose(substr(x, 5, 7)), datsAlpha), col="white", pch=21, cex=datCex)
	mtext(side=1, text=levels(factor(substr(levels(x), 1, 3))), at=seq(2, 11, 3), cex=0.8, line=2)
}

plotDatGenoSurv <- function(z, meanlwd, meancex, seed, datCex, meanLineCol, meanLineLty, meanLwd, meanCex, wideVar, CIs, ...){

	set.seed(seed)

#	z <- z[order(z$geno),]

	zLines <- reshape(z[,1:3], idvar=wideVar, timevar="diet", direction="wide")
	rownames(zLines) <- zLines[,1]
	zLines <- t(zLines[,2:4])
	zLines <- zLines[c(2,1,3),]

	LCL <- reshape(z[,c(1,2,6)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(LCL) <- LCL[,1]
	LCL <- t(LCL[,2:4])
	LCL <- LCL[c(2,1,3),]
	
	UCL <- reshape(z[,c(1,2,7)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(UCL) <- UCL[,1]
	UCL <- t(UCL[,2:4])
	UCL <- UCL[c(2,1,3),]

	durts <- reshape(z[,c(1,2,1)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(durts) <- durts[,1]
	durts <- t(durts[,2:4])
	durts <- durts[c(2,1,3),]
	
	xvals <- t(matrix(1:length(zLines), 
						ncol=3, byrow=T))#+xoffset

	if(CIs==T){
	theMax <- max(UCL, na.rm=T)
	theMin <- min(LCL, na.rm=T)
		}else{
	theMax <- max(zLines, na.rm=T)
	theMin <- min(zLines, na.rm=T)
		}
	
	par(bty="n", las=2, lwd=1)
	plot(y=c(theMax, theMin), x=c(1,length(zLines)), xlab="", col=0, ...)

#	matlines(y=zLines, x=xvals, col="white", lty=1, lwd=meanlwd*1.5)
	matlines(y=zLines, x=xvals, col=meanLineCol, lty=meanLineLty, lwd=meanLwd)

#	arrows(x0=xvals, x1=xvals, y0=LCL, y1=UCL, length=0, col="white", lwd=meanLwd)
	if(CIs==T){arrows(x0=xvals, x1=xvals, y0=LCL, y1=UCL, length=0, col=colChoose(durts), lwd=meanLwd)}
	points(y=zLines, x=xvals, bg=colChoose(durts), col="white", pch=21, cex=meanCex, lwd=meanCex)
	mtext(side=1, text=colnames(zLines), at=seq(2, 24, 3), cex=0.5, line=2)
}

plotDatGeno <- 
function(y, x, z, datCex, meanLwd, meanCex, xoffset, jitterFactor, datsAlpha, seed, hgt, wdt, fn, meanLineCol, meanLineLty, wideVar, exponent, invlogit, ...){
	set.seed(seed)

	if(exponent==T){z <- transform(z, yvar=exp(yvar), LCL=exp(LCL), UCL=exp(UCL))}
	if(invlogit==T){z <- transform(z, yvar=inv.logit(yvar), LCL=inv.logit(LCL), UCL=inv.logit(UCL))}
	
	zLines <- reshape(z[,1:3], idvar=wideVar, timevar="diet", direction="wide")
	rownames(zLines) <- zLines[,1]
	zLines <- t(zLines[,2:4])

	LCL <- reshape(z[,c(1,2,6)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(LCL) <- LCL[,1]
	LCL <- t(LCL[,2:4])

	UCL <- reshape(z[,c(1,2,7)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(UCL) <- UCL[,1]
	UCL <- t(UCL[,2:4])

	durts <- reshape(z[,c(1,2,1)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(durts) <- durts[,1]
	durts <- t(durts[,2:4])

	xvals <- t(matrix(1:length(zLines), 
						ncol=3, byrow=T))+xoffset

	xvals <- jitter(xvals, jitterFactor)

	par(bty="n", las=2, lwd=1)
	plot(y ~ x, xlab="", ...)
	points(jitter(as.numeric(x), factor=jitterFactor)-xoffset, y, bg=alpha(colChoose(substr(as.character(x), 3, 5)), datsAlpha), col="white", pch=21, cex=datCex)
	matlines(y=zLines, x=xvals, col=meanLineCol, lty=meanLineLty, lwd=meanLwd)

	arrows(x0=xvals, x1=xvals, y0=LCL, y1=UCL, length=0, col=colChoose(durts), lwd=meanLwd)
	points(y=zLines, x=xvals, bg=colChoose(durts), col="white", pch=21, cex=meanCex, lwd=meanCex)
	mtext(side=1, text=colnames(zLines), at=seq(2, 24, 3), cex=0.5, line=2)
}

plotDatSurv <- function(z, meanlwd, meancex, seed, datCex, meanLineCol, meanLineLty, meanLwd, meanCex, wideVar, ...){

	set.seed(seed)

#	z <- z[order(z$geno),]

	zLines <- reshape(z[,1:3], idvar=wideVar, timevar="diet", direction="wide")
	rownames(zLines) <- zLines[,1]
	zLines <- t(zLines[,2:4])
#	zLines <- zLines[c(1,3,2),]

	LCL <- reshape(z[,c(1,2,6)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(LCL) <- LCL[,1]
	LCL <- t(LCL[,2:4])
#	LCL <- LCL[c(1,3,2),]
	
	UCL <- reshape(z[,c(1,2,7)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(UCL) <- UCL[,1]
	UCL <- t(UCL[,2:4])
#	UCL <- UCL[c(1,3,2),]

	durts <- reshape(z[,c(1,2,1)], idvar=wideVar, timevar="diet", direction="wide")
	rownames(durts) <- durts[,1]
	durts <- t(durts[,2:4])
#	durts <- durts[c(1,3,2),]
	
	xvals <- t(matrix(	c(rep(1:3, 3), 
						rep(4:6, 3), 
						rep(7:9, 3), 
						rep(10:12, 3)), 
						ncol=3, byrow=T))

		#UGLY CODE ALERT:
		#this will truncate xvals if aa3 is absent.
		#it will need further tuning to work correctly if other genos are excluded
	if(length(xvals) > length(zLines)){
		warning("fewer genotypes than expected, aa3 is assumed absent. Revise this function if otherwise")
		xvals <- xvals[,2:ncol(xvals)]
	}

	theMax <- max(UCL, na.rm=T)
	theMin <- min(LCL, na.rm=T)
	
	par(bty="n", las=2, lwd=1)
	plot(y=c(theMax, theMin), x=c(1,12), xlab="", col=0, ...)

#	matlines(y=zLines, x=xvals, col="white", lty=1, lwd=meanlwd*1.5)
	matlines(y=zLines, x=xvals, col=meanLineCol, lty=meanLineLty, lwd=meanLwd)

#	arrows(x0=xvals, x1=xvals, y0=LCL, y1=UCL, length=0, col="white", lwd=meanLwd)
	arrows(x0=xvals, x1=xvals, y0=LCL, y1=UCL, length=0, col=colChoose(durts), lwd=meanLwd)
	points(y=zLines, x=xvals, bg=colChoose(durts), col="white", pch=21, cex=meanCex, lwd=meanCex)
	mtext(side=1, text=levels(factor(substr(colnames(durts), 1, 2))), at=seq(2, 11, 3), cex=0.8, line=2)
}

##################################################
##	impact of margarine on dahomey egg laying	##
##################################################

dah <- read.csv("RI1_dahomeyEggLaying/dah.csv", stringsAsFactors=T)
dah$eggsPerFly <- with(dah, Neggs/Nflies)
pdf("figsFromR/MarEggsRita.pdf")
par(bty="n")
plot(eggsPerFly ~ diet, dah)
dev.off()
with(dah, points(y=eggsPerFly, x=jitter(as.numeric(diet))))
with(dah, t.test(eggsPerFly[diet=="sya"], eggsPerFly[diet=="mar"], "greater"))
with(dah, tapply(eggsPerFly, diet, mean))

percChange <- 100 * ((7.412037 - 5.921384) / 7.412037)
percChange

##########################
##	load fecundity data	##
##########################
	
dee <- read.csv("20190730_fecundity/fecundity.csv", stringsAsFactors=T)
str(dee)
levels(dee$diet)

	#set levels so comparisons are made relative to SYA
dee$diet <- factor(dee$diet, levels=levels(dee$diet)[c(3,1,2)])

	#add caloric estimates
dee$calories <- NA
dee$calories[dee$diet=="sya"] <- 58.3
dee$calories[dee$diet=="eaa"] <- 71.1
dee$calories[dee$diet=="mar"] <- 166.3

	#make levels of plotting variables sane
levels(dee$cond2)
dee$cond2 <- factor(dee$cond2, levels=levels(dee$cond2)[c(1,3,2,4,6,5,7,9,8,10,12,11)])
levels(dee$cond2)

##########################
## plot eggs ~ calories	##	
##########################

plot(eggs ~ calories, dee, log="y")
dee$mn <- with(dee, factor(paste(mito, nuclear, sep="")))

pdf("FigsFromR/experiment1_eggsVsCalories.pdf", h=3, w=4)
ggplot(dee, aes(log(calories), log(eggs + 1))) + 
	geom_smooth() + 
	geom_point(aes(col=diet), alpha=0.5) + 
	labs(title="eggs ~ calories") + 
	facet_grid(cols=vars(mito), rows=vars(nuclear)) +
	scale_color_manual(values=colChoose(c("sya", "eaa", "mar"))) +
	theme_bw() +
	theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90))
dev.off()

##########################
##	exploratory plot	##
##########################

plot(eggs ~ cond2, data=dee, las=2, bty="n", log="y", xlab="")

##################################
##	inter-experiment correlation	##
##################################

experimeans <- data.frame(tapply(dee$eggs, list(dee$cond, dee$experi), mean))
experrors <- data.frame(tapply(dee$eggs, list(dee$cond, dee$experi), function(x){sd(x, na.rm=T)/sqrt(length(x))}))

colnames(experimeans) <- paste("Experi.", 1:3)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, method="spearman"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex=3)
}

par(las=2)
pairs(experimeans, lower.panel=panel.cor, col=colChoose(substr(rownames(experimeans), 5, 7)), pch=16)

##############################
## ditribution of the data?	##
##############################

#dev.new()
plot(density(dee$eggs))
	
	#where do data sit on a skewness-kurtosis plot?
descdist(dee$eggs, discrete=T)

	#fit distribs and compare
fitNbinom <- fitdist(dee$eggs, "nbinom")
fitPoisson <- fitdist(dee$eggs, "pois")
plot(fitNbinom)
plot(fitPoisson)

fitNbinom$aic
fitPoisson$aic
	#neg. binom is superior

##############################
##	model - experiment 1	##
##############################

m1nb <- glmer.nb(eggs ~ mito * nuclear * diet + (1|geno) + experi, data=dee, na.action = "na.fail")

Anova(m1nb, "III")
joint_tests(m1nb)
m1nb_diets <- as.data.frame(joint_tests(m1nb, by="diet"))
m1nb_diets <- m1nb_diets[order(m1nb_diets[,2], m1nb_diets[,1]),] 
m1nb_diets
rsqm1nb <- r.squaredGLMM(m1nb)
rsqm1nb
pdf("figsFromR/experiment1_etaSquared.pdf", h=3, w=6)
plotEtas_mixed(m1nb)
dev.off()

##########################################
##	reaction norms of experiments 1+2	##
##########################################

m3nb <- glm.nb(eggs ~ geno * diet + experi, dee)
Anova(m3nb, type="III")

#pdf("figsFromR/experiment1_reactionNorms.pdf", useDingbats=F, h=110*mmToInch, w=165*mmToInch)
with(dee, plotDat(
y = eggs+1,
x = cond2,
z = emmip(m3nb, ~ diet | geno, CIs=T, plotit=F)[c(13:24, 1:12, 25:36),], #simple way to order EMM levels 
exponent = T, invlogit = F, xoffset = 0.25, jitterFactor = 0.25, datsAlpha = 1, datCex = 0.75, meanCex = 1, meanLineLty = 1, meanLwd = 1, meanLineCol = 1, seed = 1984, log="y", col=0, border="lightgrey", xaxt="n", ylab="Eggs+1", outline=F, wideVar="geno"))
text(8, 10, paste("R^2 =", round(as.numeric(rsqm1nb),2)))
#dev.off()

######################################
##	bubble plot of experiments 1+2	##
######################################

experi1_FCdata <- bubbleFC_process(bubbleFC(emmeans(m3nb, ~ diet * geno), "", "", exponentiate=T, invLogit=F, mirror=T), sigsOnly=F, withinDiet=T, withinGenoOnly=F, vsSYAonly=F, pThresh=0.05, meanCentre=F, eaaVSmar=T, doubleDip=T)

#pdf("figsFromR/experiment1_FCs.pdf", useDingbats=F, h=100*mmToInch, w=125*mmToInch)
ggplot(NULL, aes(x,y)) +
	geom_point(data=experi1_FCdata, aes(x=x, y=y, size=neglog, fill= RInorm, col=sig), alpha=0.75, shape=21) +
	scale_fill_gradientn(colors=c("darkblue","blue", "cyan", "white", "orange", "red", "darkred"), limits=c(-1,1)) +
	theme_minimal() +
	xlab("") +
	ylab("") +
	theme(axis.text.x = element_text(angle = 90)) +
	scale_size(range = c(0.1, 5), name="-log10 FDR") +
	labs(size="-log10 FDR", fill="Response index", color="Sig. diff.") +
	scale_color_manual(values=c("grey","black")) +
	labs(title="SYA vs enriched diets") + 
	theme (plot.title = element_text (hjust = 0.5))
#dev.off()

##################################################
##	full reproductive metrics (experiment 3)	##
##################################################

d <- read.csv("20190815_fertility+transgenerational/reproductive_metrics.csv", stringsAsFactors=T)
str(d)

	#rm columns that we don't want
d <- d[,!colnames(d) %in% c("inputIndex", "adultsPre", "adultsPost", "pupaePre", "pupaePost", "pupaePreMedian", "pupaePostMedian")]

	#intuitive colnames
colnames(d)[which(colnames(d) == "eggs")] <- "eggsPre"
colnames(d)[which(colnames(d) == "eggsD2")] <- "eggsPost"
colnames(d)[which(colnames(d) == "adultsPreMedian")] <- "progenyPreMedian"
colnames(d)[which(colnames(d) == "adultsPostMedian")] <- "progenyPostMedian"
colnames(d)[which(colnames(d) == "progenyPreProp")] <- "fertilityPre"
colnames(d)[which(colnames(d) == "progenyPostProp")] <- "fertilityPost"
colnames(d)
head(d)

	#resequence levels
levels(d$diet)
d$diet <- factor(d$diet, levels=levels(d$diet)[c(1,3,2)])
levels(d$diet)

levels(d$cond2)
d$cond2 <- factor(d$cond2, levels=levels(d$cond2)[c(1,3,2,4,6,5,7,9,8,10,12,11)])
levels(d$cond2)

rownames(d) <- d$vial

pheno <- d[,colnames(d) %in% c("eggsPre", "progenyPre", "fertilityPre", "progenyPreMedian", "eggsPost", "progenyPost", "fertilityPost", "progenyPostMedian")]
phenoInd <- d[,colnames(d) %in% c("vial", "mito", "nuclear", "diet", "geno", "geno2", "cond2")]

	#add caloric estimates
d$calories <- NA
d$calories[d$diet=="sya"] <- 58.3
d$calories[d$diet=="eaa"] <- 71.1
d$calories[d$diet=="mar"] <- 166.3

########################
##	exploratory plots ##
########################

par(mfcol=c(2,4))
for(i in 1:ncol(pheno)){
	plot(pheno[,i] ~ phenoInd$cond2, ylab=colnames(pheno)[i], xlab="", las=2, horizontal=T)
}

pdf("FigsFromR/experiment2_eggsPreVsCalories.pdf", h=3, w=4)
ggplot(d, aes(log(calories), log(eggsPre + 1))) + 
	geom_smooth() + 
	geom_point(aes(col=diet), alpha=0.5) + 
	labs(title="eggs ~ calories (chronic)") + 
	facet_grid(cols=vars(mito), rows=vars(nuclear)) +
	scale_color_manual(values=colChoose(c("eaa", "sya", "mar"))) +
	theme_bw() +
	theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90)) +
	ylab("log(eggs +1)")
dev.off()

pdf("FigsFromR/experiment2_progenyPreVsCalories.pdf", h=3, w=4)
ggplot(d, aes(log(calories), log(progenyPre + 1))) + 
	geom_smooth() + 
	geom_point(aes(col=diet), alpha=0.5) + 
	labs(title="progeny ~ calories (chronic)") + 
	facet_grid(cols=vars(mito), rows=vars(nuclear)) +
	scale_color_manual(values=colChoose(c("eaa", "sya", "mar"))) +
	theme_bw() +
	theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90)) +
	ylab("log(progeny +1)")
dev.off()

pdf("FigsFromR/experiment2_fertilityPreVsCalories.pdf", h=3, w=4)
ggplot(d, aes(log(calories), fertilityPre)) + 
	geom_smooth() + 
	geom_point(aes(col=diet), alpha=0.5) + 
	labs(title="fertility ~ calories (chronic)") + 
	facet_grid(cols=vars(mito), rows=vars(nuclear)) +
	scale_color_manual(values=colChoose(c("eaa", "sya", "mar"))) +
	theme_bw() +
	theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90)) +
	ylab("fertility")
dev.off()


pdf("FigsFromR/experiment2_eggsPostVsCalories.pdf", h=3, w=4)
ggplot(d, aes(log(calories), log(eggsPost + 1))) + 
	geom_smooth() + 
	geom_point(aes(col=diet), alpha=0.5) + 
	labs(title="eggs ~ calories (parental)") + 
	facet_grid(cols=vars(mito), rows=vars(nuclear)) +
	scale_color_manual(values=colChoose(c("eaa", "sya", "mar"))) +
	theme_bw() +
	theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90)) +
	ylab("log(eggs +1)")
dev.off()

pdf("FigsFromR/experiment2_progenyPostVsCalories.pdf", h=3, w=4)
ggplot(d, aes(log(calories), log(progenyPost + 1))) + 
	geom_smooth() + 
	geom_point(aes(col=diet), alpha=0.5) + 
	labs(title="progeny ~ calories (parental)") + 
	facet_grid(cols=vars(mito), rows=vars(nuclear)) +
	scale_color_manual(values=colChoose(c("eaa", "sya", "mar"))) +
	theme_bw() +
	theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90)) +
	ylab("log(progeny +1)")
dev.off()

pdf("FigsFromR/experiment2_fertilityPostVsCalories.pdf", h=3, w=4)
ggplot(d, aes(log(calories), fertilityPost)) + 
	geom_smooth() + 
	geom_point(aes(col=diet), alpha=0.5) + 
	labs(title="fertility ~ calories (parental)") + 
	facet_grid(cols=vars(mito), rows=vars(nuclear)) +
	scale_color_manual(values=colChoose(c("eaa", "sya", "mar"))) +
	theme_bw() +
	theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90)) +
	ylab("fertility")
dev.off()

##############
##	models	##
##############

resDens <- function(x){
	plot(density(resid(x)))
}

eggsPreM1 <- glmer.nb(eggsPre ~ mito * nuclear * diet + (1|geno), d, na.action = "na.fail")
resDens(eggsPreM1)

progenyPreM1 <- update(eggsPreM1, progenyPre ~ .)
resDens(progenyPreM1)

fertilityPreM1 <- glmer(cbind(progenyPre, eggsPre) ~ mito * nuclear * diet + (1|geno), family="binomial", na.action = "na.fail", data=d)
dispersion_glmer(fertilityPreM1)	#likely overdispersed
resDens(fertilityPreM1)	
	#fit an observation-level random effect (ie. vial): https://r-sig-mixed-models.r-project.narkive.com/fi8Nu155/r-sig-me-quasi-binomial-family-in-lme4
fertilityPreM2 <- update(fertilityPreM1, . ~ . + (1|vial))
dispersion_glmer(fertilityPreM2)

timePreM1 <- lmer(progenyPreMedian ~ mito * nuclear * diet + (1|geno), subset(d, !is.na(progenyPreMedian)))	#this is a numerical ordinal with equal distances between levels.
resDens(timePreM1)

eggsPostM1 <- update(eggsPreM1, eggsPost ~ .)
resDens(eggsPostM1)

progenyPostM1 <- update(progenyPreM1, progenyPost ~ .)
resDens(progenyPostM1)

fertilityPostM1 <- with(d, update(fertilityPreM2, cbind(progenyPost, eggsPost) ~ .))
resDens(fertilityPostM1)
dispersion_glmer(fertilityPostM1)

timePostM1 <- update(timePreM1, progenyPostMedian ~ ., data=subset(d, !is.na(progenyPostMedian)))
resDens(timePostM1)

##############################
##	open development data	##
##############################

adultsPre <- read.csv("20190815_fertility+transgenerational/adultsPre.csv", stringsAsFactors=T)
adultsPost <- read.csv("20190815_fertility+transgenerational/adultsPost.csv", stringsAsFactors=T)

adultsPre$diet <- factor(adultsPre$diet, levels=c("eaa", "sya", "mar"))
adultsPost$diet <- factor(adultsPost$diet, levels=c("eaa", "sya", "mar"))

	#create survival objects
adultsPreS <- with(adultsPre, Surv(time, censor))
adultsPostS <- with(adultsPost, Surv(time, censor))

##########################################################
##	plot development data								##
##	y axis will need to be flipped in graphics editor	##
##########################################################

par(mfrow=c(2,4))

for(i in levels(adultsPre$mitoNuclear)){
	incl <- adultsPre$mitoNuclear==i
	aa <- droplevels(adultsPre[incl,])
	SS <- adultsPreS[incl]
	with(aa, plot(survfit(SS ~ cond), col=colChoose(substr(levels(cond), 5, 7)), main=i, xlab="days", ylab="prop. eclosed", conf.int=F, axes=F, ylim=c(1,0), xlim=c(9,13)))
	axis(1, at=c(9,11,13), labels=c(9,11,13))
	axis(2, at=c(0,0.5,1), labels=c(1,0.5,0))
	}

for(i in levels(adultsPost$mitoNuclear)){
	incl <- adultsPost$mitoNuclear==i
	aa <- droplevels(adultsPost[incl,])
	SS <- adultsPostS[incl]
	with(aa, plot(survfit(SS ~ cond), col=colChoose(substr(levels(cond), 5, 7)), main=i, xlab="days", ylab="prop. eclosed", conf.int=F, axes=F, ylim=c(1,0), xlim=c(9,13)))
	axis(1, at=c(9,11,13), labels=c(9,11,13))
	axis(2, at=c(0,0.5,1), labels=c(1,0.5,0))
}

	#plot mean/genotype
for(i in levels(adultsPre$mitoNuclear)){
	incl <- adultsPre$mitoNuclear==i
	aa <- droplevels(adultsPre[incl,])
	SS <- adultsPreS[incl]
	with(aa, plot(survfit(SS ~ factor(paste(mitoNuclear, diet))), col=alpha(colChoose(levels(diet)), 0.75), main=i, xlab="days", ylab="prop. eclosed", conf.int=F, axes=F, ylim=c(1,0), xlim=c(9,13)))
	axis(1, at=c(9,11,13), labels=c(9,11,13))
	axis(2, at=c(0,0.5,1), labels=c(1,0.5,0))
	}

for(i in levels(adultsPost$mitoNuclear)){
	incl <- adultsPost$mitoNuclear==i
	aa <- droplevels(adultsPost[incl,])
	SS <- adultsPostS[incl]
	with(aa, plot(survfit(SS ~ factor(paste(mitoNuclear, diet))), col=alpha(colChoose(levels(diet)), 0.75), main=i, xlab="days", ylab="prop. eclosed", conf.int=F, axes=F, ylim=c(1,0), xlim=c(9,13)))
	axis(1, at=c(9,11,13), labels=c(9,11,13))
	axis(2, at=c(0,0.5,1), labels=c(1,0.5,0))
}

##########################
##	Cox/frailty models	##
##########################

adultsPreM <- coxme(Surv(time, censor) ~ mito * nuclear * diet * sex2 + eggs + (1|geno), data=adultsPre, na.action="na.fail")
adultsPostM <- coxme(Surv(time, censor) ~ mito * nuclear * diet * sex2 + eggsD2 + (1|geno), data=adultsPost, na.action="na.fail")
	#note - including vial as nested random effect i.e. + (1|geno/vial) yields sig. 4-way interaction terms. This was included in initial pre-print. However this result is pretty unintelligible, without much visible in the data to support the statistics. vial term is also largely redundant with egg counts. A result contingent on including a highly granular random effect is not easily generalizable. Therefore we accept the simpler model.

survMods <- list(pre=adultsPreM, post=adultsPostM)

##############
##	testing	##
##############

theModsRepro <- list(
		eggsPre=eggsPreM1, 
		progenyPre=progenyPreM1, 
		fertilityPre=fertilityPreM2, 
#		timePre=timePreM1, 
		eggsPost=eggsPostM1, 
		progenyPost=progenyPostM1, 
		fertilityPost=fertilityPostM1#, 
#		timePost=timePostM1
		)

lapply(theModsRepro, function(x){Anova(x, type="III")})
lapply(survMods, function(x){Anova(x, type="III")})

rsqs <- lapply(theModsRepro, r.squaredGLMM)
rsqs

mods_tests <- lapply(theModsRepro, joint_tests)
mods_tests
survMods_tests <- lapply(survMods, joint_tests)
survMods_tests

mods_tests_diets <- lapply(theModsRepro, joint_tests, by="diet")
mods_tests_diets <- lapply(mods_tests_diets, function(x){
	x <- as.data.frame(x)
	x <- x[order(x[,2], x[,1]),]
	})
mods_tests_diets
	
survMods_tests_diets <- lapply(survMods, joint_tests, by="diet")
survMods_tests_diets <- lapply(survMods_tests_diets, function(x){
	x <- as.data.frame(x)
	x <- x[order(x[,2], x[,1]),]
	})
survMods_tests_diets


survMods_tests_sex <- lapply(survMods, joint_tests, by="sex2")
survMods_tests_sex

survMods_tests_diets <- lapply(survMods_tests_diets, function(x){
	x <- as.data.frame(x)
	x <- x[order(x[,2], x[,1]),]
	})
survMods_tests_diets

######################
##	effect sizes		##
######################

lapply(theModsRepro, plotEtas_mixed)

for(i in names(theModsRepro)){
	pdf(paste("figsFromR/experiment2_etaSquared", i, "all.pdf", sep="_"), h=3, w=6)
	print(plotEtas_mixed(theModsRepro[i][[1]]))	#plenty of warnings so definitely doing it right.
	dev.off()
	pdf(paste("figsFromR/experiment2_etaSquared", i, "interactions.pdf", sep="_"), h=3, w=6)
	print(plotEtas_mixed(theModsRepro[i][[1]]))
	dev.off()
}	
	
	#calculate survival eta2 using integrated view of DF
adultsPreM_eta2 <- F_to_eta2(f=survMods_tests$pre$F.ratio, df=survMods_tests$pre$df1, df_error=adultsPreM$df[1])
adultsPreM_eta2$Parameter <- survMods_tests$pre[,1]
adultsPostM_eta2 <- F_to_eta2(f=survMods_tests$post$F.ratio, df=survMods_tests$post$df1, df_error=adultsPostM$df[1])
adultsPostM_eta2$Parameter <- survMods_tests$post[,1]

	#plot
pdf("figsFromR/experiment2_etaSquared_developmentPre_all.pdf", h=3 * 2.285714, w=6)
plot(adultsPreM_eta2) + 
		scale_color_manual(values=c("hotpink", "hotpink")) + 
		theme_bw() +
		xlim(0,1) +
		theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90))
dev.off()
pdf("figsFromR/experiment2_etaSquared_developmentPre_interactions.pdf", h=3 * 2.285714, w=6)
plot(adultsPreM_eta2[grepl(":", adultsPreM_eta2$Parameter),]) + 
		scale_color_manual(values=c("hotpink", "hotpink")) + 
		theme_bw() +
		xlim(0,1) +
		theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90))
dev.off()
pdf("figsFromR/experiment2_etaSquared_developmentPost_all.pdf", h=3 * 2.285714, w=6)
plot(adultsPostM_eta2) + 
		scale_color_manual(values=c("hotpink", "hotpink")) + 
		theme_bw() +
		xlim(0,1) +
		theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90))
dev.off()
pdf("figsFromR/experiment2_etaSquared_developmentPost_interactions.pdf", h=3 * 2.285714, w=6)
plot(adultsPostM_eta2[grepl(":", adultsPostM_eta2$Parameter),]) + 
		scale_color_manual(values=c("hotpink", "hotpink")) + 
		theme_bw() +
		xlim(0,1) +
		theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90))
dev.off()

	#development ~ sex
pdf("figsFromR/eta2_males_pre.pdf")
plot(F_to_eta2(
	subset(survMods_tests_sex$pre, sex2=="m")$F.ratio, 
	subset(survMods_tests_sex$pre, sex2=="m")$df1, 
	malesPreM$df[1]
	))
dev.off()
pdf("figsFromR/eta2_females_pre.pdf")
plot(F_to_eta2(
	subset(survMods_tests_sex$pre, sex2=="f")$F.ratio, 
	subset(survMods_tests_sex$pre, sex2=="f")$df1, 
	femalesPreM$df[1]
	))
dev.off()
pdf("figsFromR/eta2_males_post.pdf")
plot(F_to_eta2(
	subset(survMods_tests_sex$post, sex2=="m")$F.ratio, 
	subset(survMods_tests_sex$post, sex2=="m")$df1, 
	malesPostM$df[1]
	))
dev.off()
pdf("figsFromR/eta2_females_post.pdf")
plot(F_to_eta2(
	subset(survMods_tests_sex$post, sex2=="f")$F.ratio, 
	subset(survMods_tests_sex$post, sex2=="f")$df1, 
	femalesPostM$df[1]
	))
dev.off()

######################################
##	simplified genotype model		##
##	(ignoring time)					##
######################################

repsMods <- list(
					eggsPre=glm.nb(eggsPre+1 ~ geno * diet, d),
					eggsPost=glm.nb(eggsPost+1 ~ geno * diet, d),
					progenyPre=glm.nb(progenyPre+1 ~ geno * diet, d),
					progenyPost=glm.nb(progenyPost+1 ~ geno * diet, d),
					fertilityPre=glm(cbind(progenyPre, eggsPre) ~ geno * diet, family="binomial", d),
					fertilityPost=glm(cbind(progenyPost, eggsPost) ~ geno * diet, family="binomial", d))

	#model replicates with cph (because can extract resids and df for eff size calculation)
	#EAA as baseline level for plotting
adultsPreM2.1 <- coxph(Surv(time, censor) ~ geno * diet + sex2 + eggs, adultsPre)	#due to zero counts in geno aa3 ?
adultsPre_noAA3 <- droplevels(subset(adultsPre, geno!="aa3"))
adultsPreS_noaa3 <- with(adultsPre_noAA3, Surv(time, censor))	#dredge requires a variable rather than function value
adultsPreM2 <- coxph(adultsPreS_noaa3 ~ geno * diet + sex2 + eggs, adultsPre_noAA3)
adultsPostM2 <- coxph(adultsPostS ~  geno * diet + sex2 + eggsD2, adultsPost)

##############################
##	test genotype models		##
##############################

#lapply(repsMods, summary)
#lapply(repsMods, function(x){Anova(x, type="III")})
#lapply(repsMods, joint_tests)
#summary(adultsPreM3.1)
#summary(adultsPostM3)
#Anova(adultsPreM3.1, type="III")
#Anova(adultsPostM3, type="III")
	#yup

	#update eggs post
#repsMods$eggsPost <- update(repsMods$eggsPost, . ~ geno + diet)
#joint_tests(repsMods$eggsPost)
#Anova(repsMods$eggsPost, type="III")

##########
##	AIC	##
##########

theModsRepro_AICR2 <- lapply(theModsRepro2, function(x){
	drd <- dredge(x, m.min=1, trace=2, evaluate=T)
	daModz <- get.models(drd, subset=T)
	rs <- unlist(lapply(daModz, function(x){r.squaredGLMM(x)[1,]}))
	drd <- cbind(drd, data.frame(
		r2_marginal=rs[seq(1, length(rs)-1, 2)],
		r2_conditional=rs[seq(2, length(rs), 2)]))
#	drd <- as.data.frame(drd)
	drd$formula <- as.character(lapply(daModz, formula))
	tests <- anova(daModz[18][[1]], daModz[17][[1]], daModz[16][[1]], daModz[15][[1]], daModz[14][[1]], daModz[13][[1]], daModz[12][[1]], daModz[11][[1]], daModz[10][[1]], daModz[9][[1]], daModz[8][[1]], daModz[7][[1]], daModz[6][[1]], daModz[5][[1]], daModz[4][[1]], daModz[3][[1]], daModz[2][[1]], daModz[1][[1]], test="Chisq")
	drd <- drd[order(drd[,2], drd[,3], drd[,4], drd[,5], drd[,6], drd[,7], drd[,8]),]
	drd$mostVariance <- ifelse(drd$r2_marginal==max(drd$r2_marginal), "*", "")
	drd$weightiest <- ifelse(drd$weight==max(drd$weight), "*", "")
	list(anovas=tests, summary=drd)
})

for(i in names(theModsRepro_AICR2)){
	write.csv(theModsRepro_AICR2[i][[1]]$anovas, row.names=F, na="", file=paste("figsFromR/dredge_anovas", i, ".csv", sep=""))		
	write.csv(theModsRepro_AICR2[i][[1]]$summary, row.names=F, na="", file=paste("figsFromR/dredge_summary", i, ".csv", sep=""))
}


	#now for surv mods
survMods_AIC <- list(
	pre=dredge(survMods$pre, fixed=c("eggs")),
	post=dredge(survMods$post, fixed=c("eggsD2")))
lapply(survMods_AIC, function(x){head(as.data.frame(x))})
write.csv(survMods_AIC$pre, row.names=F, "figsFromR/survPre_AIC.csv", na="")
write.csv(survMods_AIC$post, row.names=F, "figsFromR/survPost_AIC.csv", na="")
								
survMods_AICmodels <- lapply(survMods_AIC, get.models, subset=T)
survMods_AICmodels

#survMods_AICtests <- with(survMods_AICmodels, list(
#	pre=anova(pre[3][[1]], pre[2][[1]], pre[5][[1]], pre[4][[1]], pre[1][[1]]),
#	post=anova(post[4][[1]], post[5][[1]], post[3][[1]], post[2][[1]], post[1][[1]])))
#survMods_AIC <- lapply(survMods_AIC, function(x){
#	x <- as.data.frame(x)
#	x[order(x[,5], x[,3], x[,1]),]
#})
#survMods_AIC$pre$P <- rev(survMods_AICtests$pre[,4])
#survMods_AIC$post$P <- rev(survMods_AICtests$post[,4])
#survMods_AIC

##############################################
##	plots with replicated reaction norms	##
##############################################

pdf("figsFromR/experiment2_reactionNorms.pdf", useDingbats=F, h=110*mmToInch, w=165*mmToInch)
par(mfrow=c(4,2), mar=c(2,4,2,2))

with(d, plotDat(
y = eggsPre+1,
x = cond2,
z = emmip(repsMods$eggsPre, ~ diet | geno, CIs=T, plotit=F),
exponent = T,
invlogit = F, 
xoffset = 0.2,
jitterFactor = 0.25,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984, col=0, border="lightgrey", xaxt="n", ylab="Eggs+1", outline=F, wideVar="geno"))

with(d, plotDat(
y = eggsPost+1,
x = cond2,
z = emmip(repsMods$eggsPost, ~ diet | geno, CIs=T, plotit=F),
exponent = T,
invlogit = F, 
xoffset = 0.2,
jitterFactor = 0.25,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984, col=0, border="lightgrey", xaxt="n", ylab="Eggs+1", outline=F, wideVar="geno"))

with(d, plotDat(
y = progenyPre+1,
x = cond2,
z = emmip(repsMods$progenyPre, ~ diet | geno, CIs=T, plotit=F),
exponent = T,
invlogit = F, 
xoffset = 0.2,
jitterFactor = 0.25,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984, col=0, border="lightgrey", xaxt="n", ylab="Progeny+1", outline=F, wideVar="geno"))

with(d, plotDat(
y = progenyPost+1,
x = cond2,
z = emmip(repsMods$progenyPost, ~ diet | geno, CIs=T, plotit=F),
exponent = T,
invlogit = F, 
xoffset = 0.2,
jitterFactor = 0.25,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984, col=0, border="lightgrey", xaxt="n", ylab="Progeny+1", outline=F, wideVar="geno"))

with(d, plotDat(
y = fertilityPre,
x = cond2,
z = emmip(repsMods$fertilityPre, ~ diet | geno, CIs=T, plotit=F),
exponent = T,
invlogit = F, 
xoffset = 0.2,
jitterFactor = 0.25,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984, col=0, border="lightgrey", xaxt="n", ylab="Fertility", outline=F, wideVar="geno"))

with(d, plotDat(
y = fertilityPost,
x = cond2,
z = emmip(repsMods$fertilityPost, ~ diet | geno, CIs=T, plotit=F),
exponent = T,
invlogit = F, 
xoffset = 0.2,
jitterFactor = 0.25,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984, col=0, border="lightgrey", xaxt="n", ylab="Fertility", outline=F, wideVar="geno"))


plotDatSurv(z=emmip(adultsPreM2, ~ diet | geno, CIs = TRUE, plotit=F), meanlwd=1, meancex=1, seed=1984, datCex=1, meanLineCol=1, meanLineLty=1, meanLwd=1, meanCex=1, xaxt="n", ylab="development index", ylim=c(-6,2), wideVar="geno")
#plotDatSurv(z=emmip(adultsPreM2.1, ~ diet | geno, CIs = TRUE, plotit=F), meanlwd=1, meancex=1, seed=1984, datCex=1, meanLineCol=1, meanLineLty=1, meanLwd=1, meanCex=1, xaxt="n", ylab="development index", ylim=c(-20,11), wideVar="geno")
plotDatSurv(z=emmip(adultsPostM2, ~ diet | geno, CIs = TRUE, plotit=F), meanlwd=1, meancex=1, seed=1984, datCex=1, meanLineCol=1, meanLineLty=1, meanLwd=1, meanCex=1, xaxt="n", ylab="development index", wideVar="geno")
dev.off()

######################################
## bubble plots of EMM fold-change	##
######################################

	#FDR threshold
PT <- 0.01

FCdata <- rbind(
bubbleFC_process(bubbleFC(emmeans(repsMods$eggsPre, ~ diet * geno), trait="eggs", when="pre", exponentiate=T, invLogit=F, mirror=T), sigsOnly=F, withinDiet=T, withinGenoOnly=F, vsSYAonly=F, pThresh=PT, meanCentre=F, eaaVSmar=T, doubleDip=T),
bubbleFC_process(bubbleFC(emmeans(repsMods$eggsPost, ~ diet * geno), trait="eggs", when="post", exponentiate=T, invLogit=F, mirror=T),sigsOnly=F, withinDiet=T, withinGenoOnly=F, vsSYAonly=F, pThresh=PT, meanCentre=F, eaaVSmar=T, doubleDip=T),
bubbleFC_process(bubbleFC(emmeans(repsMods$progenyPre, ~ diet * geno), trait="progeny", when="pre", exponentiate=T, invLogit=F, mirror=T),sigsOnly=F, withinDiet=T, withinGenoOnly=F, vsSYAonly=F, pThresh=PT, meanCentre=F, eaaVSmar=T, doubleDip=T),
bubbleFC_process(bubbleFC(emmeans(repsMods$progenyPost, ~ diet * geno), trait="progeny", when="post", exponentiate=T, invLogit=F, mirror=T),sigsOnly=F, withinDiet=T, withinGenoOnly=F, vsSYAonly=F, pThresh=PT, meanCentre=F, eaaVSmar=T, doubleDip=T),
bubbleFC_process(bubbleFC(emmeans(repsMods$fertilityPre, ~ diet * geno), trait="fertility", when="pre", exponentiate=F, invLogit=T, mirror=T),sigsOnly=F, withinDiet=T, withinGenoOnly=F, vsSYAonly=F, pThresh=PT, meanCentre=F, eaaVSmar=T, doubleDip=T),
bubbleFC_process(bubbleFC(emmeans(repsMods$fertilityPost, ~ diet * geno), trait="fertility", when="post", exponentiate=F, invLogit=T, mirror=T),sigsOnly=F, withinDiet=T, withinGenoOnly=F, vsSYAonly=F, pThresh=PT, meanCentre=F, eaaVSmar=T, doubleDip=T),
bubbleFC_process(bubbleFC(emmeans(adultsPreM2, ~ diet * geno), trait="development", when="pre", exponentiate=T, invLogit=F, mirror=T),sigsOnly=F, withinDiet=T, withinGenoOnly=F, vsSYAonly=F, pThresh=PT, meanCentre=F, eaaVSmar=T, doubleDip=T),
bubbleFC_process(bubbleFC(emmeans(adultsPostM2, ~ diet * geno), trait="development", when="post", exponentiate=T, invLogit=F, mirror=T), sigsOnly=F, withinDiet=T, withinGenoOnly=F, vsSYAonly=F, pThresh=PT, meanCentre=F, eaaVSmar=T, doubleDip=T))

	#arrange levels for faceting
FCdata$when <- factor(FCdata$when, levels=c("pre", "post"))
FCdata$trait <- factor(FCdata$trait, levels=c("eggs", "progeny", "fertility", "development"))
FCdata$y <- factor(FCdata$y, levels=rev(levels(factor(FCdata$y))))
FCdata$x <- factor(FCdata$x, levels=rev(levels(factor(FCdata$x))))

pdf("figsFromR/experiment2_FCs.pdf", useDingbats=F, h=125*mmToInch*1.5, w=250*mmToInch*1.5)
ggplot(NULL, aes(x,y)) +
	geom_point(data=FCdata, aes(x=x, y=y, size=neglog, fill= RInorm, col=sig), alpha=0.75, shape=21) +
	scale_fill_gradientn(colors=c("darkblue","blue", "cyan", "white", "orange", "red", "darkred"), limits=c(-1,1)) +
	theme_minimal() +
	xlab("") +
	ylab("") +
	theme(axis.text.x = element_text(angle = 90)) +
	facet_grid(cols=vars(trait), rows=vars(when)) +
	scale_size(range = c(0.01, 2.5), name="-log10 FDR") +
	labs(size="-log10 FDR", fill="Response index", color="FDR ≤ 0.05") +
	scale_color_manual(values=c("grey","black")) +
	labs(title="SYA vs enriched diets") + 
	theme (plot.title = element_text (hjust = 0.5))
dev.off()

##################################
##	architecture of phenotype	##
##################################

emmsList <- list(
eggsPre = emmip(repsMods$eggsPre, ~ diet | geno, CIs=T, plotit=F),
eggsPost=emmip(repsMods$eggsPost, ~ diet | geno, CIs=T, plotit=F),
progenyPre=emmip(repsMods$progenyPre, ~ diet | geno, CIs=T, plotit=F),
progenyPost=emmip(repsMods$progenyPost, ~ diet | geno, CIs=T, plotit=F),
fertilityPre=emmip(repsMods$fertilityPre, ~ diet | geno, CIs=T, plotit=F),
fertilityPost=emmip(repsMods$fertilityPost, ~ diet | geno, CIs=T, plotit=F),
survPre=emmip(adultsPreM2.1, ~ diet | geno, CIs = TRUE, plotit=F),
survPost=emmip(adultsPostM2, ~ diet | geno, CIs = TRUE, plotit=F)
)

emmsList <- lapply(emmsList, function(x){
	x <- x[order(x$diet,x$geno),]
	x
})

emmsDF <- data.frame(
	eggsPre = emmsList$eggsPre$yvar,
	eggsPost = emmsList$eggsPost$yvar,
	progenyPre = emmsList$progenyPre$yvar,
	progenyPost = emmsList$progenyPost$yvar,
	fertilityPre = emmsList$fertilityPre$yvar,
	fertilityPost = emmsList$fertilityPost$yvar,
	survPre = emmsList$survPre$yvar,
	survPost = emmsList$survPost$yvar, 
	row.names=with(emmsList$eggsPre, paste(geno, diet)))

emm_rowVars <- data.frame(
	nuclear=toupper(substr(emmsList$eggsPre$geno,2,2)),
	mito=toupper(substr(emmsList$eggsPre$geno,1,1)),
	diet=emmsList$eggsPre$diet
)

emm_annColours <- list(
	mito=c(A="#F49622", B="#774997"),
	nuclear=c(A="#F49622", B="#774997"),
	diet=c(sya=colChoose("sya"), eaa=colChoose("eaa"), mar=colChoose("mar"))
)

	#scale
emms <- apply(emmsDF,2,scale)
dimnames(emms) <- dimnames(emmsDF)
emms <- t(emms)

topval <- max(abs(as.matrix(emms)),na.rm=T)
brks <- seq(-topval,topval,0.1)

#brks <- seq(-6,6,0.1)

myColsScaled <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(brks)-1)

pheatmap(emms, annotation_col=emm_rowVars, annotation_colors=emm_annColours, breaks=brks, col=myColsScaled,cellwidth=10,cellheight=10, border_color="white",na_col="black")

cor(emms)

########################
##	mapping mito SNPs ##
########################

mitoSnps <- read.csv("mitoSnps/susi/ABonly_mtAmtBsigdiff_major.csv")
head(mitoSnps)
dim(mitoSnps)

rownames(mitoSnps) <- mitoSnps$X
mitoSnps <- mitoSnps[,2:ncol(mitoSnps)]
table(as.matrix(mitoSnps))
mitoSnps <- data.frame(t(mitoSnps))
rownames(mitoSnps) <- gsub("X", "", rownames(mitoSnps))

	#are all SNPs biallelic?
apply(mitoSnps, 1, table)
		#yes

	#decode SNP to gene
snpAnnot <- read.table("mitoSNPs/snp-gene.txt", sep="\t", header=T)
snpAnnot$substitution[which(snpAnnot$substitution=="-")] <- NA

loci <- sort(unique(snpAnnot$gene))
nSnpColours <- length(loci)
snpAnnotCols <- brewer.pal(nSnpColours, "Paired")
names(snpAnnotCols) <- loci

lineAnnot <- data.frame(	nuclear=substr(colnames(mitoSnps), 2, 2),
mito=substr(colnames(mitoSnps), 1, 1),
				row.names=colnames(mitoSnps))

snpAnnotDf <- snpAnnot[,c("substitution", "gene")]
rownames(snpAnnotDf) <- snpAnnot$position

	#colour coding
ann_colors <- list(
	mito=c(A="#F49622", B="#774997"),
	nuclear=c(A="#F49622", B="#774997"),
	locus=snpAnnotCols,
	substitution=c(synonymous="white", nonsynonymous="black"))

mitoSnps <- as.matrix(mitoSnps)

	#bases as numbers
mitoSnpsNum <- matrix(as.numeric(as.factor(mitoSnps)), nrow=nrow(mitoSnps))
dimnames(mitoSnpsNum) <- dimnames(mitoSnps)

snpClust <- hclust(dist(t(mitoSnpsNum)))

pheatmap(t(mitoSnpsNum), col=pal_aaas()(4), display_numbers=t(mitoSnps), number_color="white", border_color="white", annotation_colors=ann_colors, annotation_row=lineAnnot, annotation_col=snpAnnotDf, cluster_cols=F, cluster_rows=snpClust)

	#5 clusters of mito genotypes	
	#8 clusters of mito-nuclear genotypes

lineAnnot2 <- lineAnnot
lineAnnot2$mitoNuclear <-as.numeric(factor(paste(cutree(snpClust, 5), lineAnnot$nuclear, sep="")))
lineAnnot2$mitoNuclear <- paste("mn", lineAnnot2$mitoNuclear, sep="")

colnames(lineAnnot2)[2] <- "mito_geographic"
bp <- brewer.pal(8, "Paired")
ann_colors2 <- list(
	mito_geographic=c(A="#F49622", B="#774997"),
	nuclear=c(A="#F49622", B="#774997"),
	locus=snpAnnotCols,
	substitution=c(synonymous="white", nonsynonymous="black"), 
	mitoNuclear = c(mn1=bp[1], mn2=bp[2], mn3=bp[3], mn4=bp[4], mn5=bp[5], mn6=bp[6], mn7=bp[7], mn8=bp[8]))

set.seed(2)
pheatmap(t(mitoSnpsNum), col=pal_aaas()(10), display_numbers=t(mitoSnps), number_color="white", border_color="white", annotation_colors=ann_colors2, annotation_row=lineAnnot2, annotation_col=snpAnnotDf, cluster_cols=F, cluster_rows=snpClust)

##################################
##	model mito as one variable	##
##################################

rename_mitoGenos <- function(x){
	daDf <- data.frame(
			geog=factor(c("ab3", "ab2", "ab1", "aa1", "aa2", "aa3", "ba1", "bb2", "bb3", "bb1", "ba2", "ba3")),
	geno=factor(c(1,2,2,3,3,4,5,6,7,7,8,8)))
	daDf$geno[match(x, daDf$geog)]
}

d$mitoNuclear <- rename_mitoGenos(d$geno)
adultsPre$mitoNuclear <- rename_mitoGenos(adultsPre$geno)
adultsPost$mitoNuclear <- rename_mitoGenos(adultsPost$geno)
adultsPre_noAA3$mitoNuclear <- droplevels(rename_mitoGenos(adultsPre_noAA3$geno))
adultsPre <- droplevels(adultsPre)
adultsPost <- droplevels(adultsPost)
adultsPre_noAA3 <- droplevels(adultsPre_noAA3)

adultsPre$diet <- factor(adultsPre$diet, levels=c("sya", "eaa", "mar"))
adultsPost$diet <- factor(adultsPost$diet, levels=c("sya", "eaa", "mar"))
d$diet <- factor(d$diet, levels=c("sya", "eaa", "mar"))

	#update models
theModsRepro_mitoGeno <- list(
		eggsPre=glm.nb(eggsPre ~ diet * mitoNuclear, d), 
		progenyPre=glm.nb(progenyPre ~ diet * mitoNuclear, d),  
		fertilityPre=glm(cbind(progenyPre, eggsPre) ~ diet * mitoNuclear, d, family="binomial"), 
		eggsPost=glm.nb(eggsPost ~ diet * mitoNuclear, d), 
		progenyPost=glm.nb(progenyPost ~ diet * mitoNuclear, d), 
		fertilityPost=glm(cbind(progenyPost, eggsPost) ~ diet * mitoNuclear, d, family="binomial")
		)

lapply(theModsRepro_mitoGeno, Anova, type="III")
lapply(theModsRepro_mitoGeno, summary)

survMods_mitoGeno <- list(
	pre=coxph(Surv(time, censor) ~ mitoNuclear * diet * sex2 + eggs, adultsPre),
	post=coxph(Surv(time, censor) ~ mitoNuclear * diet * sex2 + eggsD2, adultsPost)
#	pre=coxph(Surv(time, censor) ~ mitoNuclear * diet + sex2, adultsPre), 
#	preForPlot=coxme(Surv(time, censor) ~ mitoNuclear * diet  + (1|geno), adultsPre),
#	preNoAA3=coxph(Surv(time, censor) ~ mitoNuclear * diet + sex2, adultsPre_noAA3),
#	pre_male=coxph(Surv(time, censor) ~ mitoNuclear * diet, subset(adultsPre, sex2=="m")), 
#	preForPlot_male=coxme(Surv(time, censor) ~ mitoNuclear * diet  + (1|geno), subset(adultsPre, sex2=="m")),
#	preNoAA3_male=coxph(Surv(time, censor) ~ mitoNuclear * diet, subset(adultsPre_noAA3, sex2=="m")),
#	pre_female=coxph(Surv(time, censor) ~ mitoNuclear * diet, subset(adultsPre, sex2=="f")), 
#	preForPlot_female=coxme(Surv(time, censor) ~ mitoNuclear * diet + (1|geno), subset(adultsPre, sex2=="f")),
#	preNoAA3_female=coxph(Surv(time, censor) ~ mitoNuclear * diet, subset(adultsPre_noAA3, sex2=="f")),
#	post=coxph(Surv(time, censor) ~ mitoNuclear * diet + sex2, adultsPost),
#	post_male=coxph(Surv(time, censor) ~ mitoNuclear * diet, subset(adultsPost, sex2=="m")),
#	post_female=coxph(Surv(time, censor) ~ mitoNuclear * diet, subset(adultsPost, sex2=="f"))
	)

lapply(survMods_mitoGeno, Anova, type="III")

lapply(survMods_mitoGeno, joint_tests)
lapply(survMods_mitoGeno, summary)

lapply(theModsRepro_mitoGeno, function(x){
	plot(density(resid(x)))
})


mods_tests <- lapply(theModsRepro_mitoGeno, joint_tests, by="mitoNuclear")
mods_tests
survMods_mitoGeno_tests <- lapply(survMods_mitoGeno, joint_tests, by="mitoNuclear")
survMods_mitoGeno_tests

##########################
##	plot reaction norms	##
##########################

d$mitoNuclearDiet <- with(d, factor(paste(mitoNuclear, diet)))
d$mitoNuclearDiet <- factor(d$mitoNuclearDiet, levels=levels(d$mitoNuclearDiet)[c(1,3,2,4,6,5,7,9,8,10,12,11,13,15,14,16,18,17,19,21,20,22,24,23)])

pdf("figsFromR/mitoNuclear_rxns.pdf", h=7,w=7)
par(mfcol=c(5,2), mar=c(5,5,3,1))

with(d, plotDatGeno(
y = eggsPre,
x = mitoNuclearDiet,
z = emmip(theModsRepro_mitoGeno$eggsPre, ~ diet | mitoNuclear, CIs=T, plotit=F)[c(9:16, 1:8, 17:24),],
exponent = T,
invlogit = F,
xoffset = 0.25,
jitterFactor = 0.5,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984,
col=0,
border=0,
xaxt="n",
ylab="Eggs",
outline=F,
wideVar="mitoNuclear",
main="eggs mito nuclear pre",
las=2
))

with(d, plotDatGeno(
y = progenyPre,
x = mitoNuclearDiet,
z = emmip(theModsRepro_mitoGeno$progenyPre, ~ diet | mitoNuclear, CIs=T, plotit=F)[c(9:16, 1:8, 17:24),],
exponent = T,
invlogit = F,
xoffset = 0.25,
jitterFactor = 0.5,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984,
col=0,
border=0,
xaxt="n",
ylab="Progeny",
outline=F,
wideVar="mitoNuclear",
main="Progeny mito nuclear pre"
))

with(d, plotDatGeno(
y = fertilityPre,
x = mitoNuclearDiet,
z = emmip(theModsRepro_mitoGeno$fertilityPre, ~ diet | mitoNuclear, CIs=T, plotit=F)[c(9:16, 1:8, 17:24),],
exponent = T,
invlogit = F,
xoffset = 0.25,
jitterFactor = 0.5,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984,
col=0,
border=0,
xaxt="n",
ylab="Fertility",
outline=F,
wideVar="mitoNuclear",
main="Fertility mito nuclear pre"
))

plotDatGenoSurv(z=emmip(survMods_mitoGeno$pre, ~ diet | mitoNuclear, CIs = TRUE, plotit=F), meanlwd=1, meancex=1, seed=1984, datCex=1, meanLineCol=1, meanLineLty=1, meanLwd=1, meanCex=1, xaxt="n", ylab="development index", wideVar="mitoNuclear", main="Development index pre", CIs=F)

with(d, plotDatGeno(
y = eggsPost,
x = mitoNuclearDiet,
z = emmip(theModsRepro_mitoGeno$eggsPost, ~ diet | mitoNuclear, CIs=T, plotit=F)[c(9:16, 1:8, 17:24),],
exponent = T,
invlogit = F,
xoffset = 0.25,
jitterFactor = 0.5,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984,
col=0,
border=0,
xaxt="n",
ylab="Eggs",
outline=F,
wideVar="mitoNuclear",
main="eggs mito nuclear Post"
))

with(d, plotDatGeno(
y = progenyPost,
x = mitoNuclearDiet,
z = emmip(theModsRepro_mitoGeno$progenyPost, ~ diet | mitoNuclear, CIs=T, plotit=F)[c(9:16, 1:8, 17:24),],
exponent = T,
invlogit = F,
xoffset = 0.25,
jitterFactor = 0.5,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984,
col=0,
border=0,
xaxt="n",
ylab="Progeny",
outline=F,
wideVar="mitoNuclear",
main="Progeny mito nuclear Post"
))


with(d, plotDatGeno(
y = fertilityPost,
x = mitoNuclearDiet,
z = emmip(theModsRepro_mitoGeno$fertilityPost, ~ diet | mitoNuclear, CIs=T, plotit=F)[c(9:16, 1:8, 17:24),],
exponent = T,
invlogit = F,
xoffset = 0.25,
jitterFactor = 0.5,
datsAlpha = 1,
datCex = 0.75,
meanCex = 1,
meanLineLty = 1,
meanLwd = 1,
meanLineCol = 1,
seed = 1984,
col=0,
border=0,
xaxt="n",
ylab="Fertility",
outline=F,
wideVar="mitoNuclear",
main="Fertility mito nuclear Post"
))

plotDatGenoSurv(z=emmip(survMods_mitoGeno$post, ~ diet | mitoNuclear, CIs = TRUE, plotit=F), meanlwd=1, meancex=1, seed=1984, datCex=1, meanLineCol=1, meanLineLty=1, meanLwd=1, meanCex=1, xaxt="n", ylab="development index", wideVar="mitoNuclear", main="Development index post", CIs=T)

dev.off()

	surv without AA3 EAA
plotDatGenoSurv(z=emmip(survMods_mitoGeno$preNoAA3, ~ diet | mitoNuclear, CIs = TRUE, plotit=F), meanlwd=1, meancex=1, seed=1984, datCex=1, meanLineCol=1, meanLineLty=1, meanLwd=1, meanCex=1, xaxt="n", ylab="development index", wideVar="mitoNuclear", main="Development index pre",CIs=T)
plotDatGenoSurv(z=emmip(survMods_mitoGeno$preNoAA3_female, ~ diet | mitoNuclear, CIs = TRUE, plotit=F), meanlwd=1, meancex=1, seed=1984, datCex=1, meanLineCol=1, meanLineLty=1, meanLwd=1, meanCex=1, xaxt="n", ylab="development index", wideVar="mitoNuclear", main="Development index pre female no aa3",CIs=T, ylim=c(-6,3))
plotDatGenoSurv(z=emmip(survMods_mitoGeno$preNoAA3_male, ~ diet | mitoNuclear, CIs = TRUE, plotit=F), meanlwd=1, meancex=1, seed=1984, datCex=1, meanLineCol=1, meanLineLty=1, meanLwd=1, meanCex=1, xaxt="n", ylab="development index", wideVar="mitoNuclear", main="Development index pre male no aa3",CIs=T, ylim=c(-6,3))


##############################################
##	testing effects of particular mito SNPs	##
##############################################

	#subset of lines differentiated only by nuclear genotype and lnRNA SNP
lnRNA <- droplevels(subset(d, mitoNuclear %in% c(5,6,7,8)))
dim(lnRNA)
lnRNA$lnRNA <- ifelse(lnRNA$mitoNuclear %in% c(5,6), "T", "C")

lnRNA_eggsPre <- glm.nb(eggsPre ~ lnRNA * nuclear * diet, lnRNA)
Anova(lnRNA_eggsPre, type="III")
joint_tests(lnRNA_eggsPre, by="diet")
lnRNA_eggsPre_eaa <- update(lnRNA_eggsPre, data=droplevels(subset(lnRNA, diet!="mar")))
Anova(lnRNA_eggsPre_eaa, "III")

lnRNA_progenyPre <- update(lnRNA_eggsPre, progenyPre ~ .)
Anova(lnRNA_progenyPre, "III")
joint_tests(lnRNA_progenyPre, by="diet")

lnRNA_fertilityPre <- glm(cbind(progenyPre, eggsPre) ~ lnRNA * nuclear * diet, data=lnRNA, family="binomial")
Anova(lnRNA_fertilityPre, type="III")
joint_tests(lnRNA_fertilityPre, by="diet")

lnRNA_eggsPost <- glm.nb(eggsPost ~ lnRNA * nuclear * diet, lnRNA)
Anova(lnRNA_eggsPost, type="III")
joint_tests(lnRNA_eggsPost, by="diet")

lnRNA_progenyPost <- update(lnRNA_eggsPost, progenyPost ~ .)
Anova(lnRNA_progenyPost, "III")
joint_tests(lnRNA_progenyPost, by="diet")

lnRNA_fertilityPost <- glm(cbind(progenyPost, eggsPost) ~ lnRNA * nuclear * diet, data=lnRNA, family="binomial")
Anova(lnRNA_fertilityPost, type="III")

adultsPre_lnRNA <- droplevels(subset(adultsPre, mitoNuclear %in% c(5,6,7,8)))
dim(adultsPre_lnRNA)
adultsPre_lnRNA$lnRNA <- ifelse(adultsPre_lnRNA$mitoNuclear %in% c(5,6), "T", "C")
adultsPre_lnRNA_m1 <- coxme(Surv(time, censor) ~ lnRNA  * nuclear * diet * sex2 + eggs + (1|geno), adultsPre_lnRNA)
Anova(adultsPre_lnRNA_m1, type="III")

	#calculate survival eta2 using integrated view of DF
adultsPre_lnRNA_m1_jt <- joint_tests(adultsPre_lnRNA_m1)
adultsPre_lnRNA_m1_eta2 <- F_to_eta2(f=adultsPre_lnRNA_m1_jt$F.ratio, df=adultsPre_lnRNA_m1_jt$df1, df_error=adultsPre_lnRNA_m1$df[1])
adultsPre_lnRNA_m1_eta2$Parameter <- adultsPre_lnRNA_m1_jt[,1]

adultsPost_lnRNA <- droplevels(subset(adultsPost, mitoNuclear %in% c(5,6,7,8)))
dim(adultsPost_lnRNA)
adultsPost_lnRNA$lnRNA <- ifelse(adultsPost_lnRNA$mitoNuclear %in% c(5,6), "T", "C")
adultsPost_lnRNA_m1 <- coxme(Surv(time, censor) ~ lnRNA  * nuclear * diet * sex2 + eggsD2 + (1|geno), adultsPost_lnRNA)
Anova(adultsPost_lnRNA_m1, type="III")

adultsPost_lnRNA_m1_jt <- joint_tests(adultsPost_lnRNA_m1)
adultsPost_lnRNA_m1_eta2 <- F_to_eta2(f=adultsPost_lnRNA_m1_jt$F.ratio, df=adultsPost_lnRNA_m1_jt$df1, df_error=adultsPost_lnRNA_m1$df[1])
adultsPost_lnRNA_m1_eta2$Parameter <- adultsPost_lnRNA_m1_jt[,1]



Anova(lnRNA_eggsPre, type="III")
Anova(lnRNA_progenyPre, "III")
Anova(lnRNA_fertilityPre, type="III")
Anova(adultsPre_lnRNA_m1, type="III")
Anova(lnRNA_eggsPost, type="III")
Anova(lnRNA_progenyPost, "III")
Anova(lnRNA_fertilityPost, type="III")
Anova(adultsPost_lnRNA_m1, type="III")



pdf("figsFromR/ETA2_lnRNA_eggsPre.pdf", h=3, w=6)
plotEtas_fixed(lnRNA_eggsPre)
dev.off()
pdf("figsFromR/ETA2_lnRNA_progenyPre.pdf", h=3, w=6)
plotEtas_fixed(lnRNA_progenyPre)
dev.off()
pdf("figsFromR/ETA2_lnRNA_fertilityPre.pdf", h=3, w=6)
plotEtas_fixed(lnRNA_fertilityPre)
dev.off()
pdf("figsFromR/ETA2_lnRNA_adultsPre_m1.pdf", h=3 * 2.285714, w=6)
plot(adultsPre_lnRNA_m1_eta2)+ 
		scale_color_manual(values=c("hotpink", "hotpink")) + 
		theme_bw() +
		xlim(0,1) +
		theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90))
dev.off()

pdf("figsFromR/ETA2_lnRNA_eggsPost.pdf", h=3, w=6)
plotEtas_fixed(lnRNA_eggsPost)
dev.off()
pdf("figsFromR/ETA2_lnRNA_progenyPost.pdf", h=3, w=6)
plotEtas_fixed(lnRNA_progenyPost)
dev.off()
pdf("figsFromR/ETA2_lnRNA_fertilityPost.pdf", h=3, w=6)
plotEtas_fixed(lnRNA_fertilityPost)
dev.off()
pdf("figsFromR/ETA2_lnRNA_adultsPost_m1.pdf", h=3 * 2.285714, w=6)
plot(adultsPost_lnRNA_m1_eta2)+ 
		scale_color_manual(values=c("hotpink", "hotpink")) + 
		theme_bw() +
		xlim(0,1) +
		theme(axis.text=element_text(size=10), axis.title=element_text(size=12, face="bold"), axis.text.x = element_text(angle = 90))
dev.off()

##############################################
##	repeatability of mito-nuclear effects	##
##############################################

b <- prcomp(t(emms), scale = F, center = T)
summary(b)
df <- as.data.frame(predict(b)[,1:2])
df$diet <- factor(substr(colnames(emms), 5, 7), levels=c("sya", "eaa", "mar"))
df$mitoGeno <- rename_mitoGenos(substr(colnames(emms), 1, 3))
df$mito <- substr(colnames(emms), 1, 1)
df$nuclear <- substr(colnames(emms), 2, 2)

ggplot(df, aes(PC1, PC2, aes(mito, nuclear, diet)))   + 
  geom_point(aes(colour=diet)) +
  theme_bw() +
  facet_grid(cols=vars(mito), rows=vars(nuclear))+
  scale_color_manual(values=alpha(colChoose(c("sya", "eaa", "mar")), 0.75))

  

ggplot(df, aes(PC1, PC2, aes(mitoGeno, diet)))  +
  geom_point(aes(colour=diet)) +
  theme_bw() +
	facet_grid(cols=vars(mitoGeno))+
	scale_color_manual(values=colChoose(c("sya", "eaa", "mar")))

ggplot(subset(df, mitoGeno!="3"), aes(PC1, PC2, aes(mitoGeno, diet)))  +
  geom_point(aes(colour=diet)) +
  theme_bw() +
	facet_grid(cols=vars(mitoGeno))+
	scale_color_manual(values=colChoose(c("sya", "eaa", "mar")))


	#just lines that are replicate genotypes
emmVars <- emm_rowVars
emmVars$geno <- substr(colnames(emms), 1, 3)
emmVars$mitoNuclear <- paste("x", rename_mitoGenos(emmVars$geno), sep="")
emmVars <- emmVars[with(emmVars, order(diet, geno)),]
rownames(emmVars) <- with(emmVars, paste(geno, diet))
all(rownames(emmVars) == colnames(emms))
emmReps <- emms[,emmVars$mitoNuclear %in% paste("x",c(2,3,7,8),sep="")]
emmRepsVars <- droplevels(subset(emmVars, mitoNuclear %in% paste("x",c(2,3,7,8),sep="")))
all(rownames(emmRepsVars)==colnames(emmReps))
emmRepsVars <- emmRepsVars[,c(1,2,3,5)]
emmRepsVars$mito <- factor(emmRepsVars$mito)
emmRepsVars$nuclear <- factor(emmRepsVars$nuclear)
emmRepsVars$mitoNuclear <- factor(emmRepsVars$mitoNuclear)

emm_annColours3 <- list(
	mito=c(A="#F49622", B="#774997"),
	nuclear=c(A="#F49622", B="#774997"),
	diet=c(sya=colChoose("sya"), eaa=colChoose("eaa"), mar=colChoose("mar")),
	mitoNuclear=c(x2="#15AC8A", x3="#E2740B", x7="#B58825", x8="#7A797A"))

pheatmap(cor(emmReps, method="kendall"), annotation_col=emmRepsVars, annotation_row=emmRepsVars, annotation_colors=emm_annColours3, border_color="lightgrey")