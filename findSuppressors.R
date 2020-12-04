#!/usr/bin/Rscript
options(stringsAsFactors=FALSE)

if (!("yogitools" %in% installed.packages()[,1])) {
	stop("This script requires the package 'yogitools'\n=> https://github.com/jweile/yogitools")
}

fitted <- read.csv("results/folate_response_model5.csv")
rownames(fitted) <- fitted$hgvs
fitted$variant[fitted$type=="synonymous"] <- fitted$ancestral[fitted$type=="synonymous"]

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



synDists <-  with(fitted,list(
	syn12 = m12.score[type=="synonymous" & e.logl > -8],
	syn25 = m25.score[type=="synonymous" & e.logl > -8],
	syn100 = m100.score[type=="synonymous" & e.logl > -8],
	syn200 = m200.score[type=="synonymous" & e.logl > -8]
))
# layout(cbind(1:4))
# lapply(synDists, function(syns) {
# 	hist(syns,breaks=seq(0,1.5,0.05))
# 	abline(v=syns,col="gray90")
# })
# fitted[which(fitted$type=="synonymous" & fitted$m.logl > -8),c("hgvs","m12.score")]
avQuants <- do.call(rbind,lapply(synDists, function(syns) {
	sapply(c(0.025,0.25,0.5,0.75,0.975), function(q) quantile(syns,q,na.rm=TRUE))
}))



diffs <- do.call(rbind,lapply(1:nrow(fitted), function(i) {
	#discard rows with missing data
	if (is.na(fitted[i,"e.post.b"])) return(rep(NA,3))
	#is it a good fit?
	goodFit <- with(fitted[i,],{
		w.logl > -8 && (e.logl > -8 || e.logl.static > -8)
	})
	if (!goodFit) return(rep(NA,3))
	#is positive interactor?
	posIA <- with(fitted[i,],{
		e.post.b > 0.95 && e.b.static > 0 || e.post.r > 0.95 && e.r > 0
	})
	if (!posIA) return(rep(NA,3))

	#check at low, mid and high concentrations whether the double mutant is 
	#fitter than 95% quantile of A222V+synonymous
	sapply(c(12,100,200),function(conc) {
		# modelFuns[[i]]$dm(conc) - avsm(conc)
		modelFuns[[i]]$dm(conc) - avQuants[paste0("syn",conc),"97.5%"]
	})
}))
dimnames(diffs) <- list(fitted$hgvs,c(12,100,200))
suppr <- apply(diffs,c(1,2),`>`,0)
supprVars <- fitted$hgvs[which(rowSums(suppr)>0)]

# colSums(suppr,na.rm=TRUE)

low <- names(which(suppr[,1] & !suppr[,3]))
high <- names(which(!suppr[,1] & suppr[,3]))
both <- names(which(suppr[,1] & suppr[,3]))

#order by magnitude of difference
low <- low[order(diffs[low,1],decreasing=TRUE)]
high <- high[order(diffs[high,3],decreasing=TRUE)]
both <- both[order(rowSums(diffs[both,c(1,3)]),decreasing=TRUE)]

#positional enrichments?
bothpos <- as.integer(yogitools::extract.groups(both,"(\\d+)"))
pdf("vis/suppressors_both_histo.pdf",7,3)
op <- par(mar=c(5,4,1,1)+.1)
hist(bothpos,breaks=45,col="steelblue3",border=NA,main="",xlab="AA position",ylab="#Suppressors")
par(op)
dev.off()

drawDoubleTrajectory <- function(i) {
	condNames <- c(w="WT",m="A222V")
	condCols <- c(w="black",m="steelblue3")
	#prepare plot frame
	plot(NA,type="n",
		xlim=c(0,200),ylim=c(-.1,2),axes=FALSE,
		main=fitted[i,"hgvs"],
		xlab="folate (ug/ml)",ylab="score"
	)
	axis(1,at=concs)
	axis(2)
	#draw A222V-SM quantiles
	polygon(c(concs,rev(concs)),c(avQuants[,"2.5%"],rev(avQuants[,"97.5%"])),border=NA,col="gray90")
	polygon(c(concs,rev(concs)),c(avQuants[,"25%"],rev(avQuants[,"75%"])),border=NA,col="gray80")
	#draw datapoints with error bars
	for (cond in c("w","m")) {
		score <- unlist(fitted[i,sprintf("%s%d.score",cond,concs)])
		se <- unlist(fitted[i,sprintf("%s%d.se",cond,concs)])
		with(fitted[i,],segments(concs-2,score,concs+2,score,col=condCols[[cond]]))
		with(fitted,arrows(concs,score-se,concs,score+se,code=3,length=0.01,angle=90,col=condCols[[cond]]))
		abline(h=0:1,col=c("firebrick3","chartreuse3"),lty="dotted")
	}
	sm <- modelFuns[[i]]$sm
	edm <- modelFuns[[i]]$edm
	dm <- modelFuns[[i]]$dm
	#draw curves
	curve(sm, from=0,to=200,col=condCols[[1]],add=TRUE)
	curve(avsm, from=0,to=200,col="gray",lty="dashed",add=TRUE)
	curve(edm,from=0,to=200,col="orange",add=TRUE)
	curve(dm,from=0,to=200,col=condCols[[2]],add=TRUE)
	return(invisible(NULL))
}

pdf("vis/suppressors_2020-04.pdf",6,5)
layout(t(matrix(1:12,ncol=3)))
invisible(lapply(both[1:4],drawDoubleTrajectory))
invisible(lapply(low[1:4],drawDoubleTrajectory))
invisible(lapply(high[1:4],drawDoubleTrajectory))
dev.off()

bothTable <- fitted[both,]
lowTable <- fitted[low,]
highTable <- fitted[high,]

write.csv(bothTable,"results/suppressors_both.csv",row.names=FALSE)
write.csv(lowTable,"results/suppressors_low.csv",row.names=FALSE)
write.csv(highTable,"results/suppressors_high.csv",row.names=FALSE)
