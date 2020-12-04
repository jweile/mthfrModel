#!/usr/bin/Rscript

options(stringsAsFactors=FALSE)

if (!("yogitools" %in% installed.packages()[,1])) {
	stop("This script requires the package 'yogitools'\n=> https://github.com/jweile/yogitools")
}
library(yogitools)


fitted <- read.csv("results/folate_response_model5.csv")
rownames(fitted) <- fitted$hgvs

#fit error model
avstderr.w <- apply(fitted[,c("w12.se","w25.se","w100.se","w200.se")],1,median,na.rm=TRUE)
loglBins <- seq(-15,-1)
bin.medians <- sapply(loglBins, function(bin) median(avstderr.w[which(fitted$w.logl >= bin & fitted$w.logl < bin+1)],na.rm=TRUE))
z <- lm(log(bin.medians) ~ loglBins)
co <- coefficients(z)
logl2stderr <- function(x) {
	exp(co[[1]] + co[[2]] * x)
}

fitted$w.fitness.sd <- logl2stderr(fitted$w.logl)*sqrt(4)
# fitted$m.fitness.sd <- logl2stderr(fitted$m.logl)*sqrt(4)

qualFilter <- 0.3

syn.distr <- with(fitted,w.fitness[which(type=="synonymous" & w.fitness.sd/sqrt(4) < qualFilter)])

syn.sd <- sd(syn.distr) 

welch.test <- function(m1,m2,sd1,sd2,n1,n2,tail=c("lower","upper","both")) {
	tail <- match.arg(tail,c("lower","upper","both"))
	t <- (m1-m2)/sqrt((sd1^2)/n1 + (sd2^2)/n2)
	nu <- (((sd1^2)/n1 + (sd2^2)/n2)^2) / ((sd1^4)/((n1-1)*n1^2) + (sd2^4)/((n2-1)*n2^2))
	p <- switch(tail,
		lower=pt(t,nu),
		upper=pt(t,nu,lower.tail=FALSE),
		both=2*pt(-abs(t),nu)
	)
	return(c(t=t,p=p))
}

candidates <- with(fitted,fitted[which(w.fitness > 1 & w.fitness.sd/sqrt(4) < qualFilter),])

wtest <- do.call(rbind,lapply(1:nrow(candidates), function(i) with(candidates[i,], {
	welch.test(1,w.fitness,syn.sd,w.fitness.sd,length(syn.distr),4,tail="lower")
})))
wtest <- cbind(candidates[,c("w.fitness","w.fitness.sd")],wtest)
wtest$q <- p.adjust(wtest$p,method="fdr")

hits <- wtest[which(wtest$q < 0.05),]
hits <- hits[order(hits$w.fitness,decreasing=TRUE),]

write.csv(hits,"results/hypercomplementers.csv")


###########################################################
# Are hypercomplementers enriched in the regulatory domain?
###########################################################

hit.resi <- fitted[rownames(hits),"start"]
catalytic <- 48:337
regulatory <- 363:644

catalytic.vars <- length(catalytic)*19
regulatory.vars <- length(regulatory)*19

catalytic.hits <- sum(hit.resi %in% catalytic)
regulatory.hits <- sum(hit.resi %in% regulatory)

regu.test <- fisher.test(rbind(
	regulatory=c(
		hit=regulatory.hits,
		nothit=regulatory.vars - regulatory.hits
	),
	catalytic=c(
		hit=catalytic.hits,
		nothit=catalytic.vars - catalytic.hits
	)
))

cat("Test if hypercomplementers are enriched in regulatory domain\n")
print(regu.test)

###########################################
# Are they enriched in the SAM interface? #
###########################################

strucfeats <- read.csv("reference_data/MTHFR_structural_features.csv")

sam.resi <- with(strucfeats,Position[SAM.SAH=="resi"])
non.resi <- setdiff(regulatory,sam.resi)

all.sam.hits <- table(hit.resi[which(hit.resi %in% sam.resi)])

sam.hits <- sum(hit.resi %in% sam.resi)
non.hits <- sum(hit.resi %in% non.resi)

sam.test <- fisher.test(rbind(
	sam=c(
		hit=sam.hits,
		nothit=length(sam.resi)*19 - sam.hits
	),
	rest=c(
		hit=non.hits,
		nothit=length(non.resi)*19 - non.hits
	)
))
cat("Test if hypercomplementars are enriched in SAM interface\n")
print(sam.test)


########################################
#Draw a map of variant effects near interface
#######################################

submap <- fitted[which(fitted$start %in% names(all.sam.hits)),]
aas <- toChars("AVLIMFYWHKRDENQSTGCP")
aanames <- c("Ala","Val","Leu","Ile","Met","Phe","Tyr","Trp","His","Lys",
			"Arg","Asp","Glu","Asn","Gln","Ser","Thr","Gly","Cys","Pro")
fitvals <- lapply(names(all.sam.hits), function(pos) {
	sapply(aas, function(aa) with(submap,{
		if (sum(start==pos & variant==aa,na.rm=TRUE)==1){
			w.fitness[which(start==pos & variant==aa)]
		} else NA
	}))
})
logls <- lapply(names(all.sam.hits), function(pos) {
	sapply(aas, function(aa) with(submap,{
		if (sum(start==pos & variant==aa,na.rm=TRUE)==1){
			w.logl[which(start==pos & variant==aa)]
		} else NA
	}))
})

pdf("vis/hypercomplementers_near_SAM.pdf",5,5)
op <- par(bg="gray",xaxs="i",las=2,mar=c(5,4,1,1))
plot(NA,type="n",
	xlim=c(0.5,length(fitvals)+.5),ylim=c(0.5,20.5),
	axes=FALSE,xlab="Position",ylab="Variant"
)
axis(1,1:length(fitvals),names(all.sam.hits))
axis(2,20:1,aanames)
cm <- yogitools::colmap()
x <- do.call(c,lapply(1:length(fitvals),rep,20))
y <- rep(20:1,length(fitvals))
rect(x-.4,y-.5,x+.4,y+.5,col=cm(do.call(c,fitvals)),border=NA)
slash <- sapply(do.call(c,logls),function(logl) {
	if (is.na(logl) || logl < -10) .3 else -.3*logl/10
})
backslash <- sapply(do.call(c,logls),function(logl) {
	if (is.na(logl) || logl < -12) .3 else 0
})
segments(x-slash,y-slash,x+slash,y+slash)
segments(x-backslash,y+backslash,x+backslash,y-backslash)
par(op)
dev.off()

######################
# draw enrichment Figure
######################

pdf("vis/hypercomp_fig.pdf",9,3)
layout(cbind(1,2,3))

op <- par(mar=c(5,4,1,1)+.1)
with(wtest,plot(w.fitness,-log10(q),
	pch=".",xlab="model fitness",ylab=expression(-log[10](q))
))
abline(h=-log10(0.05),col="red",lty="dashed")
text(2.25,-log10(0.05),"q=0.05",col="red",pos=1)

bars <- c(
	catalytic=100*catalytic.hits/catalytic.vars,
	regulatory=100*regulatory.hits/regulatory.vars
)
barplot(bars,ylab="% significant hypercomplementers",xlab="domain",
	col=c("darkolivegreen3","steelblue3"),border=NA
)
label <- with(regu.test,sprintf("OR=%.02f; p=%f",estimate,p.value))
text(mean(par("usr")[1:2]),mean(par("usr")[3:4]),label)

bars <- c(
	`SAM binding`=100*sam.hits/(length(sam.resi)*19),
	`other regulatory`=100*non.hits/(length(non.resi)*19)
)
barplot(bars,ylab="% significant hypercomplementers",xlab="residues",
	col=c("firebrick3","gold"),border=NA
)
label <- with(sam.test,sprintf("OR=%.02f; p=%f",estimate,p.value))
text(mean(par("usr")[1:2]),mean(par("usr")[3:4]),label)

par(op)
dev.off()

