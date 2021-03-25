#!/usr/bin/Rscript

options(stringsAsFactors=FALSE)

if (!("yogitools" %in% installed.packages()[,1])) {
	stop("This script requires the package 'yogitools'\n=> https://github.com/jweile/yogitools")
}
if (!("hgvsParseR" %in% installed.packages()[,1])) {
	stop("This script requires the package 'hgvsParseR'\n=> https://github.com/VariantEffect/hgvsParseR")
}
if (!("yogiroc" %in% installed.packages()[,1])) {
	stop("This script requires the package 'yogiroc'\n=> https://github.com/jweile/yogiroc")
}
if (!("maveLLR" %in% installed.packages()[,1])) {
	stop("This script requires the package 'maveLLR'\n=> https://github.com/jweile/maveLLR")
}

library(pbmcapply)
library(beeswarm)
library(hgvsParseR)
library(yogitools)
library(yogiroc)
library(maveLLR)


###########################
# Read age-of-onset table #
###########################
aooTable <- read.csv("reference_data/age_of_onset_cleaned_corrected.csv")
colnames(aooTable) <- c(
	"cellLine","aooMonths","aooLabel","caseActivity","ctrlActivity","relActivity",
	"nucl1","prot1","type1","nucl2","prot2","type2","homhet",
	"comments","a222v","e429a","thermolabile","source1","source2","source3"
)


###############################################
# Extract disease missense variant candidates #
# (a.k.a positive reference set)
###############################################


prsPreCandidates <- unique(with(aooTable,do.call(c,lapply(which(aooLabel%in%c("early","late") & (type1 == "missense" | type2 == "missense")), function(i) {
	if (type1[[i]] == "missense" && type2[[i]] == "missense") c(prot1[[i]],prot2[[i]])
	else if (type1[[i]] == "missense") prot1[[i]]
	else prot2[[i]]
}))))

subtable <- with(aooTable,aooTable[which(aooLabel%in%c("early","late") & (type1 == "missense" | type2 == "missense")),])
preAoo <- sapply(prsPreCandidates,function(pc) {
	evl <- table(with(subtable,aooLabel[which(prot1 == pc | prot2 == pc)]))
	if (length(evl)==2 && evl[[1]]==evl[[2]]) return("none")
	else names(which.max(evl))
})
preAoo <- preAoo[!grepl("=|Ter$",names(preAoo))]
prsCandidates <- names(which(preAoo=="early"))


################################
# and for late onset variants
################################
prsLateCandidates <- names(which(preAoo=="late"))

#parse the PRS HGVS strings
prsCanDetail <- parseHGVS(prsCandidates)
prsLateCanDetail <- parseHGVS(prsLateCandidates)


#######################################################
# Build a negative reference set from Gnomad variants #
#######################################################
gnomad <- read.csv("reference_data/gnomad_mthfr.csv")
nrsCandidates <- gnomad[which(with(gnomad,Homozygote.Count > 0 | Allele.Frequency > 1e-4)),"Protein.Consequence"]

#correct the MVN that is entered incorrectly in gnomad
nrsCandidates <- c(setdiff(nrsCandidates,c("p.Glu470Val","p.Glu470Gln")),"p.Glu470Leu")

#if they were seen in symptomatic patients, they are probably not benign
nrsCandidates <- setdiff(nrsCandidates,c(prsCandidates,prsLateCandidates))
nrsCanDetail <- parseHGVS(nrsCandidates)

prsnrs <- rbind(
	cbind(prsCanDetail,reference="positive"),
	cbind(nrsCanDetail,reference="negative")
)
write.csv(prsnrs,"results/mthfr_prsNrs2.csv",row.names=FALSE)

#########################################
# Load fitted models and filter by logL #
##########################################
fitted <- read.csv("results/folate_response_model5.csv")
rownames(fitted) <- fitted$hgvs
fitted <- fitted[which(fitted$w.logl > -10),]

#contains edmCorr(sm.fitness,conc) function
load("results/edm_correction.Rdata")
concs <- c(12,25,100,200)
avsm <- function(conc) fitted["p.Ala222Val","w.fitness"]+fitted["p.Ala222Val","w.remediation"]*conc
modelFuns <- lapply(1:nrow(fitted), function(i) {
	#select single mutant model
	if (!is.na(fitted[i,"w.post"]) && fitted[i,"w.post"] > 0.5) {
		sm <- function(conc) fitted[i,"w.fitness"]+fitted[i,"w.remediation"]*conc
	} else {
		meanf <- mean(unlist(fitted[i,sprintf("w%d.score",concs)]),na.rm=TRUE)
		sm <- function(conc) rep(meanf,length(conc))
	}
	#expected double-mutant fitness
	edm <- function(conc) sm(conc) * avsm(conc) + edmCorr(fitted[i,"w.fitness"],conc)
	#select double mutant model
	if (!is.na(fitted[i,"e.post.r"]) && fitted[i,"e.post.r"] > 0.5) {
		#dynamic interaction
		dm <- function(conc) edm(conc) + fitted[i,"e.b"] + fitted[i,"e.r"]*conc
	} else if (!is.na(fitted[i,"e.post.b"]) && fitted[i,"e.post.b"] > 0.5) {
		#static interaction
		dm <- function(conc) edm(conc) + fitted[i,"e.b.static"]
	} else if (any(is.na(fitted[i,c("e.post.r","e.post.b")]))) {
		#no data
		dm <- function(conc) rep(NA,length(conc))
	} else {
		#no interaction
		dm <- edm
	}
	list(sm=sm,edm=edm,dm=dm)
})
names(modelFuns) <- fitted$hgvs






##########################################################
# Define function for domain-specific LLR transformation #
##########################################################

#inter-domain border
idb <- 337
cataPRS <- prsCanDetail$hgvs[prsCanDetail$start < idb]
cataNRS <- nrsCanDetail$hgvs[nrsCanDetail$start < idb]
reguPRS <- prsCanDetail$hgvs[prsCanDetail$start >= idb]
reguNRS <- nrsCanDetail$hgvs[nrsCanDetail$start >= idb]

#Load Polyphen-2 and PROVEAN predictions
b <- new.hgvs.builder.p(aacode=3)
insilico <- read.csv("reference_data/insilico.csv")
rownames(insilico) <- sapply(1:nrow(insilico), function(i) with(insilico[i,],{
	b$substitution(pos,ref,alt)
}))

#Load CADD predictions
cadd <- read.csv("reference_data/cadd_mthfr.csv")
rownames(cadd) <- cadd$hgvp
insilico$cadd <- cadd[rownames(insilico),"cadd.mean"]

#Load REVEL predictions
revel <- read.csv("reference_data/revel_mthfr.csv")
rownames(revel) <- revel$hgvsp
insilico$revel <- revel[rownames(insilico),"revel.mean"]

#Load Deogen2 predictions
deogen2 <- read.csv("reference_data/deogen2_mthfr.csv")
rownames(deogen2) <- deogen2$hgvsp
insilico$deogen2 <- deogen2[rownames(insilico),"deogen"]

#Load Varity predictions
varity <- read.csv("reference_data/varity_mthfr.csv")
rownames(varity) <- varity$hgvsp
insilico$varity <- varity[rownames(insilico),"varity"]

#Load Varity predictions
snap2 <- read.csv("reference_data/snap2_mthfr.csv")
rownames(snap2) <- snap2$hgvsp
insilico$snap2 <- snap2[rownames(insilico),"snap2"]

#################################################
# Draw PRC curves for all different maps and models  
################################################

evaluatePRCs <- function(fitted,prs,nrs,title) {
	prsNrs <- c(prs,nrs)
	truth <- c(rep(TRUE,length(prs)),rep(FALSE,length(nrs)))

	wtAverage <- rowMeans(fitted[prsNrs,sprintf("w%d.score",c(12,25,100,200))],na.rm=TRUE)
	wtMedian <- apply(fitted[prsNrs,sprintf("w%d.score",c(12,25,100,200))],1,median,na.rm=TRUE)
	mutAverage <- rowMeans(fitted[prsNrs,sprintf("m%d.score",c(12,25,100,200))],na.rm=TRUE)
	mutMedian <- apply(fitted[prsNrs,sprintf("m%d.score",c(12,25,100,200))],1,median,na.rm=TRUE)
	allAverage <- rowMeans(fitted[prsNrs,
		do.call(c,lapply(c("w%d.score","m%d.score"),function(x)sprintf(x,c(12,25,100,200)))),
	],na.rm=TRUE)

	yrAverages <- yr2(truth,cbind(
		`SM average`=wtAverage,
		`DM average`=mutAverage,
		`All average`=allAverage
	),high=FALSE)
	avAUC <- auprc(yrAverages,balanced=TRUE)

	folateConcs <- 1:220

	smVirtual <- do.call(cbind,setNames(lapply(folateConcs, function(folateConc) {
		sapply(prsNrs, function(mut) {
			if (mut %in% names(modelFuns)) {
				modelFuns[[mut]]$sm(folateConc)
			} else NA
		})
	}),paste0("smVirtual",folateConcs)))
	
	dmVirtual <- do.call(cbind,setNames(lapply(folateConcs, function(folateConc) {
		sapply(prsNrs, function(mut) {
			if (mut %in% names(modelFuns)) {
				modelFuns[[mut]]$dm(folateConc)
			} else NA
		})
	}),paste0("dmVirtual",folateConcs)))

	yrVirt <- yr2(truth,cbind(smVirtual,dmVirtual),high=FALSE)
	virtAUC <- auprc(yrVirt,balanced=TRUE)

	smExperim <- do.call(cbind,setNames(lapply(c(12,25,100,200), function(folateConc) {
		fitted[prsNrs,sprintf("w%d.score",folateConc)]
	}),paste0("smExp",c(12,25,100,200))))
	dmExperim <- do.call(cbind,setNames(lapply(c(12,25,100,200), function(folateConc) {
		fitted[prsNrs,sprintf("m%d.score",folateConc)]
	}),paste0("dmExp",c(12,25,100,200))))

	yrExperim <- yr2(truth,cbind(smExperim,dmExperim),high=FALSE)
	experimAUC <- auprc(yrExperim,balanced=TRUE)

	isScores <- insilico[prsNrs,
		c("provean","pp2div","cadd","revel","deogen2","varity","snap2")
	]
	yrInsilico <- yr2(truth,as.matrix(isScores),
		high=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
	)
	insilicoAUC <- auprc(yrInsilico,balanced=TRUE)

	yrJoint <- structure(
		c(yrExperim[c("smExp100","dmExp25")],yrAverages,yrInsilico),
		class="yr2"
	)

	layout(cbind(1,2))
	op <- par(mar=c(5,4,1,1)+.1,oma=c(0,0,3,0))
	#PRC curves
	yrColors <- c("steelblue3","chartreuse3","steelblue4",
				  "chartreuse4","gold","gray30","gray40","gray50","gray60","gray70","gray80","gray90")
	draw.prc(yrJoint,col=yrColors,lwd=2,balanced=TRUE)
	# draw.prc.CI(yrJoint,col=yrColors,lwd=2)
	abline(h=90,lty="dotted",col="gray80")
	#AUC tracks
	plot(folateConcs,virtAUC[paste0("smVirtual",folateConcs)],
		type="l",lwd=2,col="steelblue3",
		xlab="[folate] (ug/ml)",ylab="area under balanced PRC curve",
		# ylim=range(c(virtAUC,experimAUC,avAUC,insilicoAUC))
		ylim=c(0.5,1)
	)
	lines(folateConcs,virtAUC[paste0("dmVirtual",folateConcs)],lwd=2,col="chartreuse3")
	points(c(12.5,25,100,200),experimAUC[paste0("smExp",c(12,25,100,200))],col="steelblue3",pch=20)
	points(c(12.5,25,100,200),experimAUC[paste0("dmExp",c(12,25,100,200))],col="chartreuse3",pch=20)
	abline(h=avAUC[["SM average"]],col="steelblue3",lty="dashed")
	abline(h=avAUC[["DM average"]],col="chartreuse3",lty="dashed")
	abline(h=avAUC[["All average"]],col="gold",lty="dashed")
	abline(h=insilicoAUC[["provean"]],col="gray50",lty="dashed")
	legend("bottomleft",
		c(
			"DM folate response model","DM map mean",#"A222V map median",
			"raw A222V maps",
			"SM folate response model","SM map mean",#"WT map median",
			"raw SM maps","All maps average","PROVEAN"
		), cex=0.7,
		lwd=c(rep(c(2,1,NA),2),1,1),pch=c(rep(c(NA,NA,20),2),NA,NA),
		lty=c(rep(c("solid","dashed",NA),2),"dashed","dashed"),
		col=c(rep("chartreuse3",3),rep("steelblue3",3),"gold","gray50"),
		bg="white"
	)
	mtext(title,outer=TRUE)
	par(op)

}

pdf("vis/PRC_catalytic06.pdf",10,5)
evaluatePRCs(fitted,cataPRS,cataNRS,"Catalytic domain") 
invisible(dev.off())

pdf("vis/PRC_regulatory06.pdf",10,5)
evaluatePRCs(fitted,reguPRS,reguNRS,"Regulatory domain") 
invisible(dev.off())




#####################################################
# Plot the spatial distribution of disease variants #
# across the protein sequence and across the 
# score spectrum #
######################################################

pdf("vis/prsNrsDistr.pdf",15,5)
domains <- rbind(
	all=c(1,656),
	ser=c(1,35),
	catalytic=c(48,337),
	regulatory=c(363,644)
)
op <- par(mfrow=c(2,1))
plot(NA, type="n",xlim=c(0,720),ylim=c(0,4),xlab="AA position",ylab="",axes=FALSE)
axis(1)
abline(v=prsCanDetail$start,col="firebrick3",lwd=2)
abline(v=prsLateCanDetail$start,col="firebrick3",lty=2)
abline(v=nrsCanDetail$start,col="chartreuse3",lwd=2)
rect(domains["all",1],0,domains["all",2],1,col="white")
rect(domains[-1,1],0,domains[-1,2],1,col="gray80")
text(rowMeans(domains[-1,]),0.5,rownames(domains[-1,]))
legend("right",
	c("Early onset","Late onset","Gnomad homozygous"),
	col=c("firebrick3","firebrick3","chartreuse3"),
	lty=c(1,2,1),cex=0.7,bg="white"
)

hist(
	with(fitted,m25.score[!grepl("=|Ter",hgvs)]),
	breaks=100,col="gray",border=NA,
	xlim=c(0,2),main="",xlab="fitness for A222V at 25ug/ml folate"
)
prsColors <- sapply(prsCanDetail$start <= 337,ifelse,"firebrick3","orange")
prsLateColors <- sapply(prsLateCanDetail$start <= 337,ifelse,"firebrick3","orange")
nrsColors <- sapply(nrsCanDetail$start <= 337,ifelse,"chartreuse4","chartreuse3")
abline(v=fitted[prsCandidates,"m25.score"],col=prsColors,lwd=2)
abline(v=fitted[prsLateCandidates,"m25.score"],col=prsLateColors,lty=2)
abline(v=fitted[nrsCandidates,"m25.score"],col=nrsColors,lwd=2)

legend("right",
	c(
		"Early onset catalytic","Late onset catalytic",
		"Early onset regulatory","Late onset regulatory",
		"Gnomad homozygous"
	),
	col=c("firebrick3","firebrick3","orange","orange","chartreuse3"),
	lty=c(1,2,1,2,1),cex=0.7,bg="white"
)
par(op)
dev.off()


pdf("vis/prsNrsSpatialVFitness.pdf",10,4)
domains <- rbind(
	all=c(1,656),
	ser=c(1,35),
	catalytic=c(48,337),
	regulatory=c(363,644)
)
plot(NA, type="n",xlim=c(0,720),ylim=c(-.15,1.3),bty="n",
	xlab="AA position",ylab="fitness for A222Va at 25ug/ml folinic acid"
)
segments(prsCanDetail$start,-0.1,prsCanDetail$start,fitted[prsCandidates,"m25.score"],col="firebrick3",lwd=2)
segments(prsLateCanDetail$start,-0.1,prsLateCanDetail$start,fitted[prsLateCandidates,"m25.score"],col="firebrick3",lty=2)
segments(nrsCanDetail$start,-0.1,nrsCanDetail$start,fitted[nrsCandidates,"m25.score"],col="chartreuse3",lwd=2)
abline(h=0)
rect(domains["all",1],-0.1,domains["all",2],-0.2,col="white")
rect(domains[-1,1],-0.1,domains[-1,2],-0.2,col="gray80")
text(rowMeans(domains[-1,]),-0.15,rownames(domains[-1,]))
legend("topleft",
	c("Early onset","Late onset","Gnomad homozygous"),
	col=c("firebrick3","firebrick3","chartreuse3"),
	lty=c(1,2,1),cex=0.7,bg="white"
)
dev.off()



##########################
# beeswarm plot of scores
##########################

beedata <- list(
	earlyCatal=fitted[prsCandidates[prsCanDetail$start <= 337],"m25.score"],
	lateCatal=fitted[prsLateCandidates[prsLateCanDetail$start <= 337],"m25.score"],
	gnomadCatal=fitted[nrsCandidates[nrsCanDetail$start <= 337],"m25.score"],
	earlyRegul=fitted[prsCandidates[prsCanDetail$start > 337],"m25.score"],
	lateRegul=fitted[prsLateCandidates[prsLateCanDetail$start > 337],"m25.score"],
	gnomadRegul=fitted[nrsCandidates[nrsCanDetail$start > 337],"m25.score"]
)

drawPval <- function(leftScores,rightScores,xLeft,xRight,alternative="less",height=1.6) {
	p <- wilcox.test(leftScores,rightScores,alternative=alternative)$p.value
	lines(c(xLeft,xLeft,xRight,xRight),c(height-0.05,height,height,height-0.05))
	if (p > 0.01) {
		pLabel <- sprintf("p = %.02f",p)
	} else {
		expo <- floor(1-log10(p))
		lead <- format(p*10^expo,digits=3)
		pLabel <- bquote(p == .(lead) %*% 10 ^ -.(expo))
	}
	text(mean(c(xLeft,xRight)),height+0.1,pLabel,cex=0.8)
}

pdf("vis/prsNrsBeeswarm.pdf",5,5)
a222vSyn <- fitted["p.Ala222Val","w25.score"]
op <- par(las=2,mar=c(8,4,1,1)+1)
beeswarm(
	beedata,pch=20,ylim=c(-0.2,2.5),
	labels=rep(c("Early onset","Late onset","Gnomad"),2),
	ylab="fitness score\nin cis with A222V at 25ug/ml folate",
	col=c(
		"firebrick3","firebrick2","chartreuse3",
		"firebrick3","firebrick2","chartreuse4"
	)
)
abline(h=0:1,lty="dashed",col=c("firebrick3","chartreuse3"))
abline(h=a222vSyn,lty="dashed",col="chartreuse4")
bxplot(beedata,add=TRUE,col="gray30")
with(beedata,drawPval(earlyCatal,lateCatal,1,2))
with(beedata,drawPval(lateCatal,gnomadCatal,2,3))
with(beedata,drawPval(earlyCatal,gnomadCatal,1,3,height=2))
with(beedata,drawPval(earlyRegul,lateRegul,4,5))
with(beedata,drawPval(earlyRegul,gnomadRegul,4,6,height=2))
par(op)
mtext(c("Catalytic","Regulatory"),side=1,at=c(2,5),line=3)
dev.off()


##########################
#calculate LLR funtions
########################
m25cataLLR <- buildLLR.kernel(
	na.omit(fitted[cataPRS,"m25.score"]),
	na.omit(fitted[cataNRS,"m25.score"]),
	bw=0.1
)
m25reguLLR <- buildLLR.kernel(
	na.omit(fitted[reguPRS,"m25.score"]),
	na.omit(fitted[reguNRS,"m25.score"]),
	bw=0.1
)
proveanLLR <- buildLLR.kernel(
	na.omit(insilico[prsCandidates,"provean"]),
	na.omit(insilico[nrsCandidates,"provean"]),
	bw="nrd"
)
m25cataFun <- buildLLR.kernel(na.omit(fitted[cataPRS,"m25.score"]),na.omit(fitted[cataNRS,"m25.score"]))$llr
m25reguFun <- buildLLR.kernel(na.omit(fitted[reguPRS,"m25.score"]),na.omit(fitted[reguNRS,"m25.score"]))$llr
proveanFun <-  buildLLR.kernel(
	na.omit(insilico[prsCandidates,"provean"]),
	na.omit(insilico[nrsCandidates,"provean"]),
	bw=0.2
)$llr


#Visualize LLR functions
pdf("vis/LLRfun_density_cata.pdf",5,5)
drawDensityLLR(
	fitted$m25.score,
	m25cataLLR$llr,m25cataLLR$posDens,
	m25cataLLR$negDens,
	na.omit(fitted[cataPRS,"m25.score"]),
	na.omit(fitted[cataNRS,"m25.score"])
)
dev.off()

pdf("vis/LLRfun_density_regu.pdf",5,5)
drawDensityLLR(
	fitted$m25.score,
	m25reguLLR$llr,m25reguLLR$posDens,
	m25reguLLR$negDens,
	na.omit(fitted[reguPRS,"m25.score"]),
	na.omit(fitted[reguNRS,"m25.score"])
)
dev.off()

pdf("vis/LLRfun_density_provean.pdf",5,5)
drawDensityLLR(
	insilico$provean,
	proveanLLR$llr,proveanLLR$posDens,
	proveanLLR$negDens,
	na.omit(insilico[prsCandidates,"provean"]),
	na.omit(insilico[nrsCandidates,"provean"])
)
dev.off()


pdf("vis/LLRfun_density.pdf",5,5)
curve(m25cataFun,
	xlim=c(-.1,1.1),ylim=c(-20,20),lwd=2,
	xlab="fitness score\nin cis with A222V at 25ug/ml folate",
	ylab="log likelihood ratio of harmfulness"
)
curve(m25reguFun,add=TRUE,col="royalblue3",lwd=2)
# abline(v=0:1,col=c(2,3))
abline(h=0,lty="dotted")
drawScoreBars <- function(scores,col,bottom=TRUE) {
	top <- (-1)^((!bottom)+1)
	segments(scores,top*100,scores,top*10,col=col)
}
drawScoreBars(fitted[cataPRS,"m25.score"],"firebrick3",bottom=TRUE)
drawScoreBars(fitted[cataNRS,"m25.score"],"darkolivegreen3",bottom=TRUE)
drawScoreBars(fitted[reguPRS,"m25.score"],"firebrick2",bottom=FALSE)
drawScoreBars(fitted[reguNRS,"m25.score"],"darkolivegreen2",bottom=FALSE)
text(0.5,8,"regulatory domain",col="royalblue3")
text(0.5,-8,"catalytic domain")
dev.off()


##############################
# Beeswarm plot of LLR values
##############################

beedata <- list(
	earlyCatal=m25cataFun(fitted[prsCandidates[prsCanDetail$start <= 337],"m25.score"]),
	lateCatal=m25cataFun(fitted[prsLateCandidates[prsLateCanDetail$start <= 337],"m25.score"]),
	gnomadCatal=m25cataFun(fitted[nrsCandidates[nrsCanDetail$start <= 337],"m25.score"]),
	earlyRegul=m25reguFun(fitted[prsCandidates[prsCanDetail$start > 337],"m25.score"]),
	lateRegul=m25reguFun(fitted[prsLateCandidates[prsLateCanDetail$start > 337],"m25.score"]),
	gnomadRegul=m25reguFun(fitted[nrsCandidates[nrsCanDetail$start > 337],"m25.score"])
)

pdf("vis/prsNrsBeeswarmLLR.pdf",5,5)
op <- par(las=2,mar=c(8,4,1,1)+1)
beeswarm(
	beedata,pch=20,ylim=c(-3,5),
	labels=rep(c("Early onset","Late onset","Gnomad"),2),
	ylab="Pathogenicity LLR",
	col=c(
		"firebrick3","firebrick2","chartreuse3",
		"firebrick3","firebrick2","chartreuse4"
	)
)
# abline(h=0:1,lty="dashed",col=c("firebrick3","chartreuse3"))
abline(h=0,lty="dashed",col="gray")
bxplot(beedata,add=TRUE,col="gray30")
with(beedata,drawPval(lateCatal,earlyCatal,1,2,height=2))
with(beedata,drawPval(gnomadCatal,lateCatal,2,3,height=2))
with(beedata,drawPval(gnomadCatal,earlyCatal,1,3,height=3))
with(beedata,drawPval(lateRegul,earlyRegul,4,5,height=3))
with(beedata,drawPval(gnomadRegul,earlyRegul,4,6,height=4))
par(op)
mtext(c("Catalytic","Regulatory"),side=1,at=c(2,5),line=3)
dev.off()




################################
#Plot LLRs across folate space #
################################

llrFun25 <- function(fit,start) if (start < idb) m25cataFun(fit) else m25reguFun(fit)

outLLR <- with(fitted,mapply(llrFun25,m25.score,start))
outLLRciL <- with(fitted,mapply(llrFun25,qnorm(0.025,m25.score,m25.se),start))
outLLRciH <- with(fitted,mapply(llrFun25,qnorm(0.975,m25.score,m25.se),start))
outLLRci <- sprintf("[%.02f;%.02f]",outLLRciH,outLLRciL)
outLLRTable <- with(fitted,data.frame(
	hgvs,type,start,ancestral,variant,
	m25.score,m25.se,llr=outLLR,llrCI=outLLRci
))
outLLRTable <- outLLRTable[order(outLLRTable$llr,decreasing=TRUE),]
write.csv(outLLRTable,"results/m25LLRs.csv",row.names=FALSE)

folateConcs <- seq(1,200,2)
bins <- c(-Inf,seq(log10(1/1000),log10(1000),.2),Inf)
cmap <- colmap(
	c(log10(1/1000),log10(1/100),log10(1/2),log10(2),log10(100),log10(1000)),
	colStops=c("chartreuse4","chartreuse3","gold","gold","firebrick3","firebrick4")
)
# mids <- NULL

bincounts <- pbmclapply(folateConcs, function(folateConc) {
	mutfit <- sapply(modelFuns, function(funs) funs$dm(folateConc))
	wtfit <- sapply(modelFuns, function(funs) funs$sm(folateConc))
	mutllr <- mapply(llrFun25,mutfit,fitted$start)
	wtllr <- mapply(llrFun25,wtfit,fitted$start)
	muthist <- hist(mutllr,breaks=bins,plot=FALSE)
	wthist <- hist(wtllr,breaks=bins,plot=FALSE)
	return(list(mids=muthist$mids,mutcounts=muthist$counts,wtcounts=wthist$counts))
},mc.cores=8)

mids <- bincounts[[1]]$mids


mutcounts <- do.call(rbind,lapply(bincounts,`[[`,"mutcounts"))
mutfreqs <- mutcounts/max(rowSums(mutcounts))
mutcumus <- cbind(0,t(apply(mutfreqs,1,function(xs) sapply(1:length(xs),function(i)sum(xs[1:i],na.rm=TRUE)))))
wtcounts <- do.call(rbind,lapply(bincounts,`[[`,"wtcounts"))
wtfreqs <- wtcounts/max(rowSums(wtcounts))
wtcumus <- cbind(0,t(apply(wtfreqs,1,function(xs) sapply(1:length(xs),function(i)sum(xs[1:i],na.rm=TRUE)))))


pdf("vis/llr_folate_space_V3.pdf",10,5)
layout(cbind(1,2,3),widths=c(1,1,.2))
plot(
	NA,type="n",xlim=c(0,201),ylim=c(0,1),
	xlab="[folinic acid] (ug/ml)",ylab="proportion of variants",
	main="WT background"
)
invisible(lapply(1:length(mids), function(i) {
	polygon(
		c(folateConcs,rev(folateConcs)),c(wtcumus[,i],rev(wtcumus[,i+1])),
		col=cmap(mids[[i]]),border=NA
	)
}))
plot(
	NA,type="n",xlim=c(0,201),ylim=c(0,1),
	xlab="[folinic acid] (ug/ml)",ylab="proportion of variants",
	main="p.Ala222Val background"
)
invisible(lapply(1:length(mids), function(i) {
	polygon(
		c(folateConcs,rev(folateConcs)),c(mutcumus[,i],rev(mutcumus[,i+1])),
		col=cmap(mids[[i]]),border=NA
	)
}))
op <- par(mar=c(5,0,4,4))
plot(NA,type="n",xlim=c(0,1),ylim=c(log10(1/1000),log10(1000)),axes=FALSE,xlab="",ylab="")
axis(4)
mtext(expression(log[10]("Likelihood Ratio")),side=4,line=3,cex=.7)
legendRange <- seq(log10(1/1000),log10(1000),0.5)
rect(0,legendRange-0.25,1,legendRange+0.25,col=cmap(legendRange),border=NA)
par(op)
invisible(dev.off())


#########################################
# Evalute LLR functions for WT maps
#########################################

getModelScores <- function(vars,conc) {
	vars <- intersect(vars,names(modelFuns))
	setNames(sapply(modelFuns[vars],function(funs)funs$sm(conc)),vars)
}

w120cataFuns <- buildLLR.kernel(
	getModelScores(cataPRS,120),
	getModelScores(cataNRS,120),
	bw=0.2
)
pdf("vis/llrW120_cata.pdf",5,5)
drawDensityLLR(
	getModelScores(names(modelFuns),120),
	w120cataFuns$llr,w120cataFuns$posDens,
	w120cataFuns$negDens,
	getModelScores(cataPRS,120),
	getModelScores(cataNRS,120)
)
invisible(dev.off())

w120reguFuns <- buildLLR.kernel(
	getModelScores(reguPRS,120),
	getModelScores(reguNRS,120),
	bw=0.2
)
pdf("vis/llrW120_regu.pdf",5,5)
drawDensityLLR(
	getModelScores(names(modelFuns),120),
	w120reguFuns$llr,w120reguFuns$posDens,
	w120reguFuns$negDens,
	getModelScores(reguPRS,120),
	getModelScores(reguNRS,120)
)
invisible(dev.off())

#############################################
# Derive both WT and A222V LLRs functions in preparation for diploids
#############################################

m25cataFun <- buildLLR.kernel(
	na.omit(fitted[cataPRS,"m25.score"]),
	na.omit(fitted[cataNRS,"m25.score"]))$llr
m25reguFun <- buildLLR.kernel(
	na.omit(fitted[reguPRS,"m25.score"]),
	na.omit(fitted[reguNRS,"m25.score"]))$llr

w120cataFun <- buildLLR.kernel(
	getModelScores(cataPRS,120),
	getModelScores(cataNRS,120),
	bw=0.2
)$llr
w120reguFun <- buildLLR.kernel(
	getModelScores(reguPRS,120),
	getModelScores(reguNRS,120),
	bw=0.2
)$llr

llrFun25 <- function(fit,start) if (start < idb) m25cataFun(fit) else m25reguFun(fit)
llrFun120 <- function(fit,start) if (start < idb) w120cataFun(fit) else w120reguFun(fit)


##################################
# Prepare diploid reference data
##################################

aooTable <- read.csv("reference_data/age_of_onset+1000g.csv")
colnames(aooTable) <- c(
	"cellLine","aooMonths","aooLabel","caseActivity","ctrlActivity","relActivity",
	"nucl1","prot1","type1","nucl2","prot2","type2","homhet",
	"comments","a222v","e429a","thermolabile","source1","source2","source3"
)

aooTable$meanActivity <- sapply(strsplit(aooTable$relActivity,","),function(xs)mean(as.numeric(trimws(xs))))

#filter down to usable cases (with all required data)
filteredAoo <- aooTable[with(aooTable,
	((type1 == "missense" & type2 %in% c("","missense")) |
	(type1 %in% c("","missense") & type2 == "missense")) & 
	a222v %in% c("AA","AV","VA","VV","AV/VA") &
	prot1 %in% c(fitted$hgvs,"") & prot2 %in% c(fitted$hgvs,"")
),]

#explicitly label WT genotypes from 1000g
filteredAoo$prot1[which(filteredAoo$prot1=="")] <- "WT"
filteredAoo$prot2[which(filteredAoo$prot2=="")] <- "WT"


#calculate the applicable scores for all in-cis scenarios (and models)
# scoreVersions <- with(fitted,data.frame(
# 	row.names=c(hgvs,"WT"),
# 	ae=c(w25.score,1),
# 	aa=c(w25.score,1)*fitted["p.Glu429Ala","w25.score"],
# 	ve=c(m25.score,fitted["p.Ala222Val","w25.score"]),
# 	ve.edm=c(w25.score,1)*fitted["p.Ala222Val","w25.score"],
# 	va=c(
# 		m25.score*fitted["p.Glu429Ala","m25.score"]/fitted["p.Ala222Val","w25.score"],
# 		fitted["p.Glu429Ala","m25.score"]
# 	)
# ))

getEDMScores <- function(vars,conc) {
	vars <- intersect(vars,names(modelFuns))
	setNames(sapply(modelFuns[vars],function(funs)funs$edm(conc)),vars)
}
getDMScores <- function(vars,conc) {
	vars <- intersect(vars,names(modelFuns))
	setNames(sapply(modelFuns[vars],function(funs)funs$dm(conc)),vars)
}

scoreVersions <- with(fitted,data.frame(
	row.names=c(hgvs,"WT"),
	ae=c(getModelScores(hgvs,120),1),
	aa=c(getModelScores(hgvs,120),1)*getModelScores("p.Glu429Ala",120),
	ve=c(getDMScores(hgvs,120),getModelScores("p.Ala222Val",120)),
	ve.edm=c(getEDMScores(hgvs,120),getModelScores("p.Ala222Val",120)),
	va=c(
		getDMScores(hgvs,120)*getDMScores("p.Glu429Ala",120)/getModelScores("p.Ala222Val",120),
		getDMScores("p.Glu429Ala",120)
	)
))

#do the same for provean
reScaledProvean <- setNames((insilico$provean+15)/15,rownames(insilico))
proveanVersions <- with(insilico,data.frame(
	row.names=c(rownames(insilico),"WT"),
	ae=c(reScaledProvean,1),
	aa=c(reScaledProvean,1)*reScaledProvean[["p.Glu429Ala"]],
	ve=c(reScaledProvean,1)*reScaledProvean[["p.Ala222Val"]],
	va=c(reScaledProvean,1)*reScaledProvean[["p.Ala222Val"]]*reScaledProvean[["p.Glu429Ala"]]
))

#convert the scores to LLRs
llrVersions <- do.call(data.frame,lapply(colnames(scoreVersions), function(cn) {
	mapply(llrFun120,
		scoreVersions[,cn],
		c(as.integer(gsub("\\D+","",rownames(scoreVersions)[-nrow(scoreVersions)])),1)
	)
}))
dimnames(llrVersions) <- dimnames(scoreVersions)

llrProveanVersions <- do.call(data.frame,lapply(colnames(proveanVersions), function(cn) {
	sapply(proveanVersions[,cn],proveanFun)
}))
dimnames(llrProveanVersions) <- dimnames(proveanVersions)

#################################
#define diploid modeling functions
###################################

#epistasis model
epiLLR <- function(mut,av,lr=1,edm=FALSE) {
	avs <- substr(strsplit(av,"/")[[1]],lr,lr)
	results <- sapply(avs, switch, 
		A={
			llrVersions[mut,"ae"]
		},
		V={
			if (edm) {
				llrVersions[mut,"ve.edm"]
			} else {
				llrVersions[mut,"ve"]
			}
		}
	)
	paste(signif(results,digits=3),collapse=";")
}

#full diploid model
fullLLR <- function(mut,av,ea,lr=1) {
	if (ea==""){
		ea <- "EE"
	}
	avs <- substr(strsplit(av,"/")[[1]],lr,lr)
	eas <- substr(strsplit(ea,"/")[[1]],lr,lr)
	results <- apply(expand.grid(av=avs,ea=eas),1,function(avea) {
		if (avea[["av"]]=="A") {
			if (avea[["ea"]]=="E") {
				llrVersions[mut,"ae"]
			} else {
				llrVersions[mut,"aa"]
			}
		} else {
			if (avea[["ea"]]=="E") {
				llrVersions[mut,"ve"]
			} else {
				llrVersions[mut,"va"]
			}
		}
	})
	paste(signif(results,digits=3),collapse=";")
}

proveanFullLLR <- function(mut,av,ea,lr=1) {
	if (ea==""){
		ea <- "EE"
	}
	avs <- substr(strsplit(av,"/")[[1]],lr,lr)
	eas <- substr(strsplit(ea,"/")[[1]],lr,lr)
	results <- apply(expand.grid(av=avs,ea=eas),1,function(avea) {
		if (avea[["av"]]=="A") {
			if (avea[["ea"]]=="E") {
				llrProveanVersions[mut,"ae"]
			} else {
				llrProveanVersions[mut,"aa"]
			}
		} else {
			if (avea[["ea"]]=="E") {
				llrProveanVersions[mut,"ve"]
			} else {
				llrProveanVersions[mut,"va"]
			}
		}
	})
	paste(signif(results,digits=3),collapse=";")
}

#calculate the results for all models
diploidModels <- as.df(lapply(1:nrow(filteredAoo), function(i) with(filteredAoo[i,],{
	list(
		allele1=prot1,
		allele2=prot2,
		a222v=a222v,
		e429a=e429a,
		aoo=aooLabel,
		months=aooMonths,
		activity=meanActivity,
		simpleLLR1=llrVersions[prot1,"ae"],
		simpleLLR2=llrVersions[prot2,"ae"],
		edmLLR1=epiLLR(prot1,a222v,1,edm=TRUE),
		edmLLR2=epiLLR(prot2,a222v,2,edm=TRUE),
		epiLLR1=epiLLR(prot1,a222v,1),
		epiLLR2=epiLLR(prot2,a222v,2),
		fullLLR1=fullLLR(prot1,a222v,e429a,1),
		fullLLR2=fullLLR(prot2,a222v,e429a,2),
		simpleProvean1=llrProveanVersions[prot1,"ae"],
		simpleProvean2=llrProveanVersions[prot2,"ae"],
		fullProvean1=proveanFullLLR(prot1,a222v,e429a,1),
		fullProvean2=proveanFullLLR(prot2,a222v,e429a,2)
	)
})))

#helper function to apply sum,mean,median,min etc to comma-separate lists
splitFuns <- function(left,right,fun) {
	lsplit <- if (is.character(left)) {
		lapply(strsplit(left,";"),as.numeric)
	} else left
	rsplit <- if (is.character(right)) {
		lapply(strsplit(right,";"),as.numeric)
	} else right
	result <- mapply(function(lelem,relem) mapply(fun,lelem,relem),lsplit,rsplit)
	sapply(result,paste,collapse=";")
}
diploidModels$simpleSum <- with(diploidModels,splitFuns(simpleLLR1,simpleLLR2,sum))
diploidModels$simpleMin <- with(diploidModels,splitFuns(simpleLLR1,simpleLLR2,min))
diploidModels$edmSum <- with(diploidModels,splitFuns(edmLLR1,edmLLR2,sum))
diploidModels$edmMin <- with(diploidModels,splitFuns(edmLLR1,edmLLR2,min))
diploidModels$epiSum <- with(diploidModels,splitFuns(epiLLR1,epiLLR2,sum))
diploidModels$epiMin <- with(diploidModels,splitFuns(epiLLR1,epiLLR2,min))
diploidModels$fullSum <- with(diploidModels,splitFuns(fullLLR1,fullLLR2,sum))
diploidModels$fullMin <- with(diploidModels,splitFuns(fullLLR1,fullLLR2,min))
diploidModels$provSum <- with(diploidModels,splitFuns(simpleProvean1,simpleProvean2,sum))
diploidModels$provMin <- with(diploidModels,splitFuns(simpleProvean1,simpleProvean2,min))
diploidModels$provFullSum <- with(diploidModels,splitFuns(fullProvean1,fullProvean2,sum))
diploidModels$provFullMin <- with(diploidModels,splitFuns(fullProvean1,fullProvean2,min))


write.csv(diploidModels,"results/diploidModels3.csv")

m <- function(xs) {
	if (is.character(xs)) {
		mean(as.numeric(strsplit(xs,";")[[1]]))
	} else {
		xs
	}
}

meanDMs <- apply(diploidModels[,c("simpleMin","edmMin","epiMin","fullMin","provMin","provFullMin")],1:2,m)
subpops <- apply(meanDMs,2,function(submodel) tapply(submodel,diploidModels$aoo,c))

######################################
#precision-recall analysis of diploid models
######################################
rocs <- structure(do.call(c,lapply(names(subpops),function(spname) {
	subpop <- subpops[[spname]]
	yogiroc::yr2(
		truth=c(rep(TRUE,length(subpop$early)),rep(FALSE,length(subpop$asymp))),
		scores=t(t(c(subpop$early,subpop$asymp))),names=spname
	)
})),class="yr2")
aucs <- auprc(rocs,balanced=TRUE)
r90ps <- recall.at.prec(rocs,.9,balanced=TRUE)
rocCols <- c(paste0("chartreuse",1:4),paste0("gray",c(60,80)))

rocs2 <- structure(do.call(c,lapply(names(subpops),function(spname) {
	subpop <- subpops[[spname]]
	yogiroc::yr2(
		truth=c(rep(TRUE,length(subpop$early)),rep(FALSE,length(subpop$late))),
		scores=t(t(c(subpop$early,subpop$late))),names=spname
	)
})),class="yr2")
aucs2 <- auprc(rocs2,balanced=TRUE)
r90ps2 <- recall.at.prec(rocs2,.9,balanced=TRUE)

pdf("vis/w120llr_diploid_PRC.pdf",10,5)
op <- par(mfrow=c(1,2))
draw.prc(rocs,balanced=TRUE,col=rocCols,lwd=2)
mtext("Early onset vs 1000genomes",line=1)
abline(h=90,col="gray",lty="dashed")
draw.prc(rocs2,balanced=TRUE,col=rocCols,lwd=2)
mtext("Early onset vs Late onset",line=1)
abline(h=90,col="gray",lty="dashed")
par(op)
dev.off()

#test AUPRC for significance
auprc.signif(rocs)
auprc.signif(rocs2)



medianDiffs <- do.call(rbind,lapply(subpops,function(spop) {
	c(
		md.early.1000=median(spop$early)-median(spop$asymp),
		md.late.1000=median(spop$late)-median(spop$asymp),
		md.early.late=median(spop$early)-median(spop$late)
	)
}))
pvals <- do.call(rbind,lapply(subpops,function(spop) {
	c(
		p.early.1000=wilcox.test(spop$early,spop$asymp)$p.value,
		p.late.1000=wilcox.test(spop$late,spop$asymp)$p.value,
		p.early.late=wilcox.test(spop$early,spop$late)$p.value
	)
}))


pdf("vis/diploidModels.pdf",8,8)
op <- par(mar=c(5,4.5,1,1))
modelNames <- c("WT","multiplicative","epistatic","full")
plot(NA,xlim=c(1,5.5),ylim=c(-3,3),axes=FALSE,xlab="model",
	ylab=expression(min(log[10]("likelihood ratio")))
)
axis(2)
axis(1,at=1:4,modelNames)
rect(-1,-1,4,1,col="gray90",border=NA)
rect(-1,log10(1/5),4,log10(5),col="gray80",border=NA)
rect(-1,log10(1/2),4,log10(2),col="gray70",border=NA)
linecol <- c(early="firebrick3",late="gold",asymp="steelblue3")
invisible(lapply(1:nrow(diploidModels), function(i) with(diploidModels[i,],{
	lines(1:4,c(m(simpleMin),m(edmMin),m(epiMin),m(fullMin)),col=linecol[[aoo]])
	points(1:4,c(m(simpleMin),m(edmMin),m(epiMin),m(fullMin)),col=linecol[[aoo]],pch=20)
})))
legend("right",
	c("early onset","late onset","1000 genomes",
		"LR factor <2","LR factor <5","LR factor <10"),
	col=c(linecol,NA,NA,NA),fill=c(NA,NA,NA,"gray70","gray80","gray90"),
	lty=c(1,1,1,NA,NA,NA),pch=c(20,20,20,NA,NA,NA),border=NA
)
dev.off()



