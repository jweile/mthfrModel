#!/usr/bin/Rscript

if (!("yogitools" %in% installed.packages()[,1])) {
	stop("This script requires the package 'yogitools'\n=> https://github.com/jweile/yogitools")
}
if (!("hgvsParseR" %in% installed.packages()[,1])) {
	stop("This script requires the package 'hgvsParseR'\n=> https://github.com/VariantEffect/hgvsParseR")
}

options(stringsAsFactors=FALSE)
library(optimization)
library(pbmcapply)
library(yogitools)
library(hgvsParseR)

cat("Reading input data...\n")

#import experimentally measured map data
maps <- list(
	w12 = read.csv("map_data/WT12.csv"),
	w25 = read.csv("map_data/WT25.csv"),
	w100 = read.csv("map_data/WT100.csv"),
	w200 = read.csv("map_data/WT200.csv"),
	m12 = read.csv("map_data/V12.csv"),
	m25 = read.csv("map_data/V25.csv"),
	m100 = read.csv("map_data/V100.csv"),
	m200 = read.csv("map_data/V200.csv")
)
for (i in 1:length(maps)) {
	rownames(maps[[i]]) <- maps[[i]]$hgvs_pro
}

allMuts <- Reduce(union,lapply(maps, function(x) x[,1]))
mutants <- parseHGVS(allMuts,aacode=1)
mutants$type[mutants$variant=="Ter"] <- "nonsense"
mutants <- cbind(mutants,do.call(cbind,lapply(maps, `[`, allMuts, c("score","se") )))
rownames(mutants) <- mutants$hgvs

#normalize A222V double mutants to overall WT fitness (instead of A222V single mutant fitness)
concs <- c(12,25,100,200)
for (conc in concs) {
	wscore <- sprintf("w%d.score",conc)
	mscore <- sprintf("m%d.score",conc)
	mse <- sprintf("m%d.se",conc)
	avfactor <- mutants["p.Ala222Val",wscore]
	mutants[,mscore] <- mutants[,mscore] * avfactor
	mutants[,mse] <- mutants[,mse] * avfactor
}

###############
# Run linear fit to folate levels
###############

cat("Fitting WT models...\n")

runFit <- function(allele="w") {
	do.call(rbind,pbmclapply(1:nrow(mutants), function(i) {
		data <- data.frame(
			concentration=concs,
			score=unlist(mutants[i,sprintf("%s%d.score",allele,concs)]),
			se=unlist(mutants[i,sprintf("%s%d.se",allele,concs)])
		)
		if (sum(is.na(data$score)) > 2) {
			#if more than two datapoints are missing we can't do anything
			return(c(fitness=NA,remediation=NA,logl=NA,null.logl=NA))
		} else if (all(data$score <= 0,na.rm=TRUE)) {
			#if everything is dead, the model is simple
			logl <- sum(dnorm(0,0,data$se,log=TRUE),na.rm=TRUE)
			return(c(fitness=0,remediation=0,logl=logl,null.logl=logl))
		}
		calcLogl <- function(x) {
			base <- x[[1]]
			remed <- x[[2]]
			sum(apply(data,1,function(row) {
				pred <- base+row["concentration"]*remed
				dnorm(pred,row["score"],row["se"],log=TRUE) 
				# - dnorm(row["score"],row["score"],row["se"],log=TRUE)
			}),na.rm=TRUE)
		}
		#the null hypothesis is that there is no remediation
		null.logl <- calcLogl(c(mean(data$score,na.rm=TRUE),0))
		#default case
		fitbase <- fitremed <- logl <- NA
		try({
			start <- c(0.5,0)
			optim <- optim_nm(calcLogl,start=start,maximum=TRUE)
			stopifnot(!identical(optim$par,start)) 
			fitbase <- optim$par[[1]]
			fitremed <- optim$par[[2]]
			logl <- optim$function_value
		},silent=TRUE)
		return(c(fitness=fitbase,remediation=fitremed,logl=logl,null.logl=null.logl))
	},mc.cores=8))
}

#define assumed prior probability of folate responsiveness
priorProb <- 0.01
priorLogOdds <- log(priorProb/(1-priorProb))

#define logistic function
logistic <- function(x) exp(x)/(1+exp(x))

#run model fitting
fitting.w <- runFit("w")
colnames(fitting.w) <- paste0("w.",colnames(fitting.w))
#normalize log likelihoods
loglmax <- max(fitting.w[,"w.logl"],na.rm=TRUE)
fitting.w[,"w.logl"] <- fitting.w[,"w.logl"] - loglmax
fitting.w[,"w.null.logl"] <- fitting.w[,"w.null.logl"] - loglmax
#calculate log odds by adding prior to LLRs and convert to posterior probability
fitting.w <- cbind(fitting.w, w.lod=fitting.w[,"w.logl"] - fitting.w[,"w.null.logl"] + priorLogOdds)
fitting.w <- cbind(fitting.w, w.post=sapply(fitting.w[,"w.lod"], logistic))

# fitted <- cbind(mutants,fitting.w,fitting.m)
fitted <- cbind(mutants,fitting.w)
rownames(fitted) <- fitted$hgvs

###################################################
# Run preliminary fitting for genetic interactions #
###################################################
cat("First round of epistasis fitting\n")
fitting.e <- do.call(rbind,pbmclapply(1:nrow(fitted), function(i) {
	expected <- function(conc) {
		if (fitted[i,"w.post"] > 0.5) {
			(fitted[i,"w.fitness"] + conc*fitted[i,"w.remediation"]) * 
			(fitted["p.Ala222Val","w.fitness"] + conc*fitted["p.Ala222Val","w.remediation"])
		} else {
			mean(unlist(mutants[i,sprintf("w%d.score",concs)]),na.rm=TRUE) * 
			(fitted["p.Ala222Val","w.fitness"] + conc*fitted["p.Ala222Val","w.remediation"])
		}
	}
	data <- data.frame(
		concentration=concs,
		score=unlist(fitted[i,sprintf("m%d.score",concs)]),
		se=unlist(fitted[i,sprintf("m%d.se",concs)])
	)
	if (sum(is.na(data$score)) > 2 || is.na(fitted[i,"w.fitness"])) {
		#if more than two datapoints are missing we can't do anything
		return(c(e.b=NA,e.r=NA,e.b.static=NA,e.logl=NA,e.logl.static=NA,e.logl.null=NA))
	} 
	calcLogl <- function(x) {
		b <- x[[1]]
		r <- if (length(x) > 1) x[[2]] else 0
		sum(apply(data,1,function(row) {
			pred <- b+row["concentration"]*r + expected(row["concentration"])
			dnorm(pred,row["score"],row["se"],log=TRUE) 
		}),na.rm=TRUE)
	}
	null.logl <- calcLogl(0)
	#default case in case of failure
	b <- r <- logl <- b.static <- logl.static <- NA
	try({
		start <- c(0,0)
		optim <- optim_nm(calcLogl,k=2,start=start,maximum=TRUE)
		stopifnot(!identical(optim$par,start)) 
		b <- optim$par[[1]]
		r <- optim$par[[2]]
		logl <- optim$function_value
	},silent=TRUE)
	try({
		optim <- optimize(calcLogl,c(-1,1),maximum=TRUE)
		b.static <- optim$maximum
		logl.static <- optim$objective
	},silent=TRUE)
	return(c(e.b=b,e.r=r,e.b.static=b.static,e.logl=logl,e.logl.static=logl.static,e.logl.null=null.logl))
},mc.cores=8))

#normalize log likelihoods
loglmax <- max(fitting.e[,c("e.logl","e.logl.static","e.logl.null")],na.rm=TRUE)
fitting.e[,"e.logl"] <- fitting.e[,"e.logl"] - loglmax
fitting.e[,"e.logl.null"] <- fitting.e[,"e.logl.null"] - loglmax
fitting.e[,"e.logl.static"] <- fitting.e[,"e.logl.static"] - loglmax
#calculate log odds by adding prior to LLRs and convert to posterior probability
fitting.e <- cbind(fitting.e, e.lod.r=fitting.e[,"e.logl"] - fitting.e[,"e.logl.static"] + priorLogOdds)
fitting.e <- cbind(fitting.e, e.lod.b=fitting.e[,"e.logl.static"] - fitting.e[,"e.logl.null"] + priorLogOdds)
fitting.e <- cbind(fitting.e, e.post.r=sapply(fitting.e[,"e.lod.r"], logistic))
fitting.e <- cbind(fitting.e, e.post.b=sapply(fitting.e[,"e.lod.b"], logistic))

# fitted <- cbind(mutants,fitting.w,fitting.m,fitting.e)
fitted <- cbind(mutants,fitting.w,fitting.e)
rownames(fitted) <- fitted$hgvs

#fix variant types in table
fitted$subject <- NULL
fitted$type[which(fitted$variant=="*")] <- "nonsense"
fitted$variant[fitted$type=="synonymous"] <- fitted$ancestral[fitted$type=="synonymous"]


###############################################################
# Create an interpolation model to deterimine non-linearity 
# correction factors for epistasis models
###############################################################

cat("Determining nonlinearity correction model...\n")

logl.cutoff <- -6
# delta.cutoff <- log(10)
# posteriorCutoff <- .95
# concs <- c(12,25,100,200)

filtered <- fitted[
	which(fitted[,"w.logl"] > logl.cutoff & fitted[,"e.logl"] > logl.cutoff),
]

interpolate <- function(x,y,binMids=seq(0,1.5,0.01),binWidth=0.1) {
	running <- sapply(binMids, function(mid) {
		median(y[which(abs(x-mid) < binWidth/2)],na.rm=TRUE)
	})

	polymat <- cbind(binMids,binMids^2,binMids^3,binMids^4)
	mainmodel <- lm(running~polymat)

	cutoff <- 1.2
	# tailMids <-  which(binMids > 1)
	# tailmodel <- lm(running[tailMids] ~ binMids[tailMids])
	tailIC <- (sum(coef(mainmodel)*(cutoff+0.05)^(0:ncol(polymat)))-sum(coef(mainmodel)*cutoff^(0:ncol(polymat))))/0.05
	tailParams <- c(
		intercept=sum(coef(mainmodel)*cutoff^(0:ncol(polymat)))-cutoff*tailIC,
		slope=tailIC
	)

	# modelfun <- function(.x) if (.x < 1) sum(coef(mainmodel)*.x^(0:ncol(polymat))) else sum(coef(tailmodel)*.x^(0:1))
	modelfun <- function(.x) if (.x < cutoff) sum(coef(mainmodel)*.x^(0:ncol(polymat))) else sum(tailParams*.x^(0:1))
	list(model=modelfun,data=cbind(x=binMids,median=running,model=sapply(binMids,modelfun)))
}

ipb <- interpolate(filtered$w.fitness,filtered$e.b)
ipr <- interpolate(filtered$w.fitness,filtered$e.r)
corrFuns <- list(b=ipb$model,r=ipr$model)
edmCorr <- function(wb,conc) corrFuns$b(wb) + corrFuns$r(wb)*conc

#export the correction function
save(corrFuns, edmCorr, file="results/edm_correction.Rdata")



########################################################################
# Run repeat fitting for genetic interactions with correction functions#
########################################################################

cat("Fitting corrected epsistasis model\n")
fitting.e <- do.call(rbind,pbmclapply(1:nrow(fitted), function(i) {
	expected <- function(conc) {
		if (fitted[i,"w.post"] > 0.5) {
			(fitted[i,"w.fitness"] + conc*fitted[i,"w.remediation"]) * 
			(fitted["p.Ala222Val","w.fitness"] + conc*fitted["p.Ala222Val","w.remediation"]) +
			corrFuns$b(fitted[i,"w.fitness"]) + corrFuns$r(fitted[i,"w.fitness"])*conc
		} else {
			mean(unlist(fitted[i,sprintf("w%d.score",concs)]),na.rm=TRUE) * 
			(fitted["p.Ala222Val","w.fitness"] + conc*fitted["p.Ala222Val","w.remediation"]) +
			corrFuns$b(fitted[i,"w.fitness"]) + corrFuns$r(fitted[i,"w.fitness"])*conc
		}
	}
	data <- data.frame(
		concentration=concs,
		score=unlist(fitted[i,sprintf("m%d.score",concs)]),
		se=unlist(fitted[i,sprintf("m%d.se",concs)])
	)
	if (sum(is.na(data$score)) > 2 || is.na(fitted[i,"w.fitness"])) {
		#if more than two datapoints are missing we can't do anything
		return(c(e.b=NA,e.r=NA,e.b.static=NA,e.logl=NA,e.logl.static=NA,e.logl.null=NA))
	} 
	calcLogl <- function(x) {
		b <- x[[1]]
		r <- if (length(x) > 1) x[[2]] else 0
		sum(apply(data,1,function(row) {
			pred <- b+row["concentration"]*r + expected(row["concentration"])
			dnorm(pred,row["score"],row["se"],log=TRUE) 
		}),na.rm=TRUE)
	}
	null.logl <- calcLogl(0)
	#default case in case of failure
	b <- r <- logl <- b.static <- logl.static <- NA
	try({
		start <- c(0,0)
		optim <- optim_nm(calcLogl,k=2,start=start,maximum=TRUE)
		stopifnot(!identical(optim$par,start)) 
		b <- optim$par[[1]]
		r <- optim$par[[2]]
		logl <- optim$function_value
	},silent=TRUE)
	try({
		optim <- optimize(calcLogl,c(-1,1),maximum=TRUE)
		b.static <- optim$maximum
		logl.static <- optim$objective
	},silent=TRUE)
	return(c(e.b=b,e.r=r,e.b.static=b.static,e.logl=logl,e.logl.static=logl.static,e.logl.null=null.logl))
},mc.cores=8))


#define assumed prior probability of epistasis
priorProb <- 0.01
priorLogOdds <- log(priorProb/(1-priorProb))

#normalize log likelihoods
loglmax <- max(fitting.e[,c("e.logl","e.logl.static","e.logl.null")],na.rm=TRUE)
fitting.e[,"e.logl"] <- fitting.e[,"e.logl"] - loglmax
fitting.e[,"e.logl.null"] <- fitting.e[,"e.logl.null"] - loglmax
fitting.e[,"e.logl.static"] <- fitting.e[,"e.logl.static"] - loglmax
#calculate log odds by adding prior to LLRs and convert to posterior probability
fitting.e <- cbind(fitting.e, e.lod.r=fitting.e[,"e.logl"] - fitting.e[,"e.logl.static"] + priorLogOdds)
fitting.e <- cbind(fitting.e, e.lod.b=fitting.e[,"e.logl.static"] - fitting.e[,"e.logl.null"] + priorLogOdds)
fitting.e <- cbind(fitting.e, e.post.r=sapply(fitting.e[,"e.lod.r"], logistic))
fitting.e <- cbind(fitting.e, e.post.b=sapply(fitting.e[,"e.lod.b"], logistic))

fitted <- cbind(fitted[,1:27],fitting.e)

dir.create("results",showWarnings=FALSE)
write.csv(fitted,"results/folate_response_model5.csv",row.names=FALSE)



################################################
# Graphically illustrate the interpolated nonlinearity correction
#################################################

cat("Drawing plots...\n")
modelRange <- seq(-0.1,2,0.01)

# plot(filtered[,c("w.fitness","e.b")])

dir.create("vis",showWarnings=FALSE)
pdf("vis/epsilon_normalization_model.pdf",8,4)
layout(cbind(1,2))
op <- par(mar=c(5,4,1,1)+.1)
plot(filtered[,c("w.fitness","e.b")],pch=20,col=yogitools::colAlpha("black",0.1))
abline(h=0,v=0:1,col="gray")
lines(ipb$data[,"x"],ipb$data[,"median"],col=2,lwd=2)
lines(modelRange,sapply(modelRange,ipb$model),col=3,lwd=2)

plot(filtered[,c("w.fitness","e.r")],pch=20,col=yogitools::colAlpha("black",0.1))
abline(h=0,v=0:1,col="gray")
lines(ipr$data[,"x"],ipr$data[,"median"],col=2,lwd=2)
lines(modelRange,sapply(modelRange,ipr$model),col=3,lwd=2)
par(op)
dev.off()


pdf("vis/edm_adjustment.pdf",5,5)
plot(function(x) sapply(x,function(.x)edmCorr(.x,12)),from=0,to=1.5,ylim=c(-.5,.5),ylab="Expected DM fitness adjustment",xlab="SM fitness")
plot(function(x) sapply(x,function(.x)edmCorr(.x,25)),from=0,to=1.5,add=TRUE,col=2)
plot(function(x) sapply(x,function(.x)edmCorr(.x,100)),from=0,to=1.5,add=TRUE,col=3)
plot(function(x) sapply(x,function(.x)edmCorr(.x,200)),from=0,to=1.5,add=TRUE,col=4)
legend("topright",paste(c(12.5,25,100,200),"ug/ml folate"),col=1:4,lty=1)
dev.off()


##############################################
# Print general statistics about the results
##############################################

ngood <- sum(fitted$w.logl > -10,na.rm=TRUE)
nresp <- with(fitted,sum(w.logl > -10 & w.post > 0.95, na.rm=TRUE))
nia <- with(fitted,sum(w.logl > -10 & e.logl.static > -10 & e.post.b > 0.95, na.rm=TRUE))
niar <- with(fitted,sum(w.logl > -10 & e.logl > -10 & e.post.r > 0.95,na.rm=TRUE))

sprintf("%d variants (=%.02f%%) received logL > -10.",ngood,100*ngood/nrow(fitted))
sprintf("%d variants (=%.02f%%) are folate-responsive",nresp,100*nresp/nrow(fitted))
sprintf("%d variants (=%.02f%%) have genetic interactions",nia,100*nia/nrow(fitted))
sprintf("%d variants (=%.02f%%) have folate-responsive interactions",niar,100*niar/nrow(fitted))


#####################################
# Visualize individual model fits
#####################################

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
	#draw datapoints with error bars
	for (cond in c("w","m")) {
		score <- unlist(fitted[i,sprintf("%s%d.score",cond,concs)])
		se <- unlist(fitted[i,sprintf("%s%d.se",cond,concs)])
		with(fitted[i,],segments(concs-2,score,concs+2,score,col=condCols[[cond]]))
		with(fitted,arrows(concs,score-se,concs,score+se,code=3,length=0.01,angle=90,col=condCols[[cond]]))
		abline(h=0:1,col=c("firebrick3","chartreuse3"),lty="dotted")
	}
	#select single mutant model
	if (!is.na(fitted[i,"w.post"]) && fitted[i,"w.post"] > 0.5) {
		sm <- function(conc) fitted[i,"w.fitness"]+fitted[i,"w.remediation"]*conc
	} else {
		meanf <- mean(unlist(fitted[i,sprintf("w%d.score",concs)]),na.rm=TRUE)
		sm <- function(conc) rep(meanf,length(conc))
	}
	#expected double-mutant fitness
	edm <- function(conc) {
		sm(conc) * (fitted["p.Ala222Val","w.fitness"] + conc*fitted["p.Ala222Val","w.remediation"])
	}
	#select double mutant model
	if (!is.na(fitted[i,"e.post.r"]) && fitted[i,"e.post.r"] > 0.5) {
		dm <- function(conc) edm(conc) + fitted[i,"e.b"] + fitted[i,"e.r"]*conc
	} else if (!is.na(fitted[i,"e.post.b"]) && fitted[i,"e.post.b"] > 0.5) {
		dm <- function(conc) edm(conc) + fitted[i,"e.b"]
	} else if (any(is.na(fitted[i,c("e.post.r","e.post.b")]))) {
		dm <- function(conc) rep(NA,length(conc))
	} else {
		dm <- edm
	}
	#draw curves
	curve(sm, from=0,to=200,col=condCols[[1]],add=TRUE)
	curve(edm,from=0,to=200,col="gold",add=TRUE)
	curve(dm,from=0,to=200,col=condCols[[2]],add=TRUE)
	return(invisible(NULL))
}


examples <- "p.Ala461Ser"

pdf("vis/example_trajectories.pdf",5,5)
# invisible(lapply(sample(nrow(fitted),9),drawDoubleTrajectory))
invisible(lapply(examples,drawDoubleTrajectory))
par(op)
dev.off()


###################################################
# Find a conversion function from logL to stderr
###################################################

#compare logL to average stderr
avstderr.w <- apply(fitted[,c("w12.se","w25.se","w100.se","w200.se")],1,median,na.rm=TRUE)
loglBins <- seq(-15,-1)
# bin.q95 <- sapply(loglBins, function(bin) quantile(avstderr.w[which(fitted$w.logl >= bin & fitted$w.logl < bin+1)],0.95,na.rm=TRUE))
bin.medians <- sapply(loglBins, function(bin) median(avstderr.w[which(fitted$w.logl >= bin & fitted$w.logl < bin+1)],na.rm=TRUE))
z <- lm(log(bin.medians) ~ loglBins)
co <- coefficients(z)
logl2stderr <- function(x) {
	exp(co[[1]] + co[[2]] * x)
}

################################################
# Draw logL histogram and conversion function
################################################

pdf("vis/loglVstderr.pdf",7,5)
layout(rbind(1,2),heights=c(1,2))
op <- par(mar=c(0,4,1,1)+.1)

x <- seq(-35,0,0.5)
y <- sapply(x,function(.x)sum(is.na(fitted$w.logl) | fitted$w.logl < .x,na.rm=TRUE))/nrow(fitted)
plot(x,y*100,
	type="l",axes=FALSE,lwd=2,col="steelblue3",
	xlab="",ylab="cdf (%)",
	ylim=c(0,100),xlim=c(-30,0)
)
grid(NULL,NULL)
abline(v=-10,lty="dotted")
axis(2)

op <- par(mar=c(5,4,0,1)+.1)
with(fitted,plot(
	w.logl,avstderr.w,
	log="y",pch=".",
	xlim=c(-30,0),
	xlab="log Likelihood",ylab="average stderr"
))
mtext("(across folate conditions)",side=2,line=2,cex=0.7)
grid(NULL,NULL)
abline(h=0.3,v=-10,lty="dotted")
curve(logl2stderr,from=-50,to=0,add=TRUE,col="gray",lty="dashed")
# lines(loglBins,bin.medians,col="green")
par(op)
invisible(dev.off())



#######################
#Export list of "compelling" folate responders
#######################

signif.idx <- which(with(fitted,w.post > postThr & w.logl > -10))
signif.response <- fitted[which(with(fitted,w.post > postThr & w.logl > -10)),]
signif.response <- signif.response[order(abs(signif.response$w.remediation),decreasing=TRUE),]

write.csv(signif.response,"results/folate_response_model_signif5.csv",row.names=FALSE)



cat("Execution complete\n")
