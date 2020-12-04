#####################################
# This script is written to process Amber MD output
# in NetCDF format. The data is not provided here
# as it comprises over 200GB and thus exceeds storage
# limits for github.
#####################################

options(stringsAsFactors=FALSE)

if (!("yogitools" %in% installed.packages()[,1])) {
	stop("This script requires the package 'yogitools'\n=> https://github.com/jweile/yogitools")
}

#Load NetCDF parser
library(ncdf4)
library(yogitools)
library(beeswarm)
library(mclust)


#function to parse the metadata file
parseMetadata <- function(mdFile) {
	#read all the lines in the file
	lines <- readLines(mdFile)
	#find the lines with section headers
	flaglines <- grep("^%FLAG",lines)
	#extract the section names
	sectionNames <- trimws(sub("^%FLAG ","",lines[flaglines]))
	#extract the section formatting control strings
	sectionFormats <- trimws(gsub("%FORMAT\\(|\\)","",lines[flaglines+1]))
	#determine section start and end line numbers
	sectionStarts <- flaglines+2
	sectionEnds <- c(flaglines[-1]-1,length(lines))

	#a helper function to parse a section block given its format
	#nf = number of fields; fl = field length
	parseBlock <- function(i,castFun=as.integer,nf=10,fl=8) {
		do.call(c,lapply(lines[sectionStarts[[i]]:sectionEnds[[i]]],function(line){
			na.omit(castFun(trimws(mapply(substr,line,0:(nf-1)*fl+1,1:nf*fl,SIMPLIFY=FALSE))))
		}))
	}

	#parse each section using the parseBlock() function while applying
	# the appropriate format as per the control string
	metaData <- setNames(lapply(1:length(flaglines), function(i) {
		switch(sectionFormats[[i]],
			`10I8`=parseBlock(i,as.integer,10,8),
			`1a80`=parseBlock(i,as.character,1,80),
			`20a4`=parseBlock(i,as.character,20,4),
			`3I8`=parseBlock(i,as.integer,3,8),
			`5E16.8`=parseBlock(i,as.numeric,5,16)
		)
	}),sectionNames)
}


#function to calculate time series of euclidian distances between atoms
#across frames from coordinate array
euDists <- function(coordinates,atoms1,atoms2) {
	#select the atom slices from the matrix and find their center of mass
	centerTrj1 <- if(length(atoms1) > 1) {
		apply(coordinates[,atoms1,],c(1,3),mean)
	} else coordinates[,atoms1,]
	centerTrj2 <- if(length(atoms2) > 2) {
		apply(coordinates[,atoms2,],c(1,3),mean)
	} else coordinates[,atoms2,]
	# and use apply on frame dimension
	apply(yogitools::zbind(centerTrj1,centerTrj2), 2, function(coord) {
		#distance = sqrt of sum of squared differences
		sqrt(sum((coord[,1]-coord[,2])^2))
	})
}

#function to calculate distance tables for a set of ncdf files
#atom group is a list of lists, containing pairs of atom index vectors
distances <- function(nc.files, atomGroups) {
	lapply(nc.files, function(nc.file) {
		#open NetCDF file connection
		ncCon <- nc_open(nc.file)
		#extract multidimensional coordinate array
		#dimension 1 : x-y-z; dimension 2: atom index; dimension 3: frames
		coordinates <- ncvar_get(ncCon,"coordinates")
		#extract time table
		timepoints <- ncvar_get(ncCon,"time")
		#close file connection
		nc_close(ncCon)

		#for each pair of atom sets
		dists <- do.call(cbind,lapply(atomGroups, function(grpPair) {
			#extract euclidian distance track for two example atoms
			euDists(coordinates, grpPair[[1]], grpPair[[2]])
		}))
		colnames(dists) <- paste0("dists",1:ncol(dists))
		
		#return table with timepoints and distances
		cbind(timepoints,dists)
	})
}


#function to determine the matrix index of a given atom
getAtomIndex <- function(metaData,resiPos,resiName,atomName,offset=35) {
	if (is.na(resiPos)) {
		resiPos <- which(metaData$RESIDUE_LABEL == resiName)+offset
		if (length(resiPos) != 1) {
			stop("No unique position for residue name!")
		} 
	} else if (metaData$RESIDUE_LABEL[[resiPos-offset]] != resiName){
		stop("Residue name does not match position!")
	} 
	rangeStart <- metaData$RESIDUE_POINTER[[resiPos-offset]]
	rangeEnd <- metaData$RESIDUE_POINTER[[resiPos-offset+1]]-1
	atomIdx <- which(metaData$ATOM_NAME[rangeStart:rangeEnd] == atomName)+rangeStart-1
	return(atomIdx)
}

#process a given directory with ncdf and metadata files
processDirectory <- function(parentdir,aa165="TRP",offset=35) {
	#select the metadata file
	# mdFile <- paste0(parentdir,"parm_strip")
	mdFile <- list.files(parentdir,pattern="^parm.+strip$",full.names=TRUE)
	#parse the metadata file
	cat("Parsing metadata...\n")
	metaData <- parseMetadata(mdFile)

	#atoms of interest 
	carbonBeta165 <- getAtomIndex(metaData,165,aa165,"CB",offset=offset)
	carbonC1FAD <- getAtomIndex(metaData,NA,"FAD","C1'",offset=offset)
	nitroN5FAD <- getAtomIndex(metaData,NA,"FAD","N5",offset=offset)

	#defining the top ring of the beta barrel
	barrelRim <- data.frame(
		resiName=c("TYR","ILE","THR","GLY","LEU","THR","VAL","PHE"),
		resiPos=c(321,256,227,196,156,129,93,64),
		atomName=rep("CA",8)
	)
	barrelRimIdx <- with(barrelRim,mapply(getAtomIndex,list(metaData),resiPos,resiName,atomName,offset))

	#define the list of atom groups to compare against each other
	atomGroups <- list(
		list(carbonBeta165, carbonC1FAD),
		list(nitroN5FAD, barrelRimIdx)
	)

	#find NetCDF files in parent directory
	nc.files <- list.files(parentdir,pattern="nc$",full.names=TRUE)
	#and calculate distance tables for C1' of FAD and CB of res165
	cat("Processing trajectories...\n")
	dtables <- distances(nc.files, atomGroups=atomGroups)

	return(dtables)
}


#select parent directory
glndist <- processDirectory("mthfr_simulations/mthfr/gln/","GLN")
trpdist <- processDirectory("mthfr_simulations/mthfr/model11/","TRP")
tyrdist <- processDirectory("mthfr_simulations/mthfr/tyr/","TYR")


glnfoldist <- processDirectory("mthfr_simulations/folate/gln_folate/","GLN",offset=36)
trpfoldist <- processDirectory("mthfr_simulations/folate/trp_folate/","TRP",offset=36)
tyrfoldist <- processDirectory("mthfr_simulations/folate/tyr_folate/","TYR",offset=36)




# max.time <- do.call(max,lapply(dtables,function(dt)dt[,"timepoints"]))
# distRange <- range(do.call(c,lapply(dtables,function(dt)dt[,"dists"])))

running <- function(x,y,w=10,n=100) {
	mids <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length.out=n)
	bins <- do.call(rbind,lapply(mids, function(m) {
		ys <- y[which(abs(x-m)<w)]
		c(median=median(ys,na.rm=TRUE),sd=sd(ys,na.rm=TRUE))
	}))
	return(cbind(binmid=mids,bins))
}
plotSmooth<- function(time,burial,col="steelblue3") {
		runfun <- running(time,burial)
		polygon(c(runfun[,1],rev(runfun[,1])),c(runfun[,2]+runfun[,3],rev(runfun[,2]-runfun[,3])),border=NA,col=colAlpha(col,0.1))
		lines(runfun[,1],runfun[,2],col=col,lwd=1)
}
plotTracks <- function(dtables, name, col="steelblue3",max.time=200,distRange=c(0,30),cutoff=9) {
	plot(NA,type="n",
		xlim=c(0,max.time),ylim=distRange,
		xlab="time (ns)",ylab="distance (A)"
	)
	invisible(lapply(dtables,function(dt) plotSmooth(dt[,"timepoints"]/1000,dt[,name],col)))
	abline(h=cutoff,lty="dashed")
	grid(NA,NULL)
}


#plot the distance over time
pdf("vis/dist_tracks.pdf",10,5)
layout(cbind(1:4,5:8),heights=c(.3,1,1,1))
op <- par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.5,"Loop - FAD",cex=2)
par(mar=c(5,4,0.1,3)+.1)
plotTracks(trpdist,"dists1",col="gray30")
mtext("WT",side=4,line=1)
plotTracks(glndist,"dists1",col="orange")
mtext("p.Trp165Gln",side=4,line=1)
plotTracks(tyrdist,"dists1",col="steelblue3")
mtext("p.Trp165Tyr",side=4,line=1)

par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.5,"FAD - Beta barrel",cex=2)
par(mar=c(5,4,0.1,3)+.1)
plotTracks(trpdist,"dists2",col="gray30",distRange=c(0,13),cutoff=7)
mtext("WT",side=4,line=1)
plotTracks(glndist,"dists2",col="orange",distRange=c(0,13),cutoff=7)
mtext("p.Trp165Gln",side=4,line=1)
plotTracks(tyrdist,"dists2",col="steelblue3",distRange=c(0,13),cutoff=7)
mtext("p.Trp165Tyr",side=4,line=1)
par(op)
invisible(dev.off())



#plot the distance over time
pdf("vis/dist_tracks_folate.pdf",10,5)
layout(cbind(1:4,5:8),heights=c(.3,1,1,1))
op <- par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.5,"Loop - FAD",cex=2)
par(mar=c(5,4,0.1,3)+.1)
plotTracks(trpfoldist,"dists1",col="gray30")
mtext("WT",side=4,line=1)
plotTracks(glnfoldist,"dists1",col="orange")
mtext("p.Trp165Gln",side=4,line=1)
plotTracks(tyrfoldist,"dists1",col="steelblue3")
mtext("p.Trp165Tyr",side=4,line=1)

par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.5,"FAD - Beta barrel",cex=2)
par(mar=c(5,4,0.1,3)+.1)
plotTracks(trpfoldist,"dists2",col="gray30",distRange=c(0,20),cutoff=7)
mtext("WT",side=4,line=1)
plotTracks(glnfoldist,"dists2",col="orange",distRange=c(0,20),cutoff=7)
mtext("p.Trp165Gln",side=4,line=1)
plotTracks(tyrfoldist,"dists2",col="steelblue3",distRange=c(0,20),cutoff=7)
mtext("p.Trp165Tyr",side=4,line=1)
par(op)
invisible(dev.off())


#calculate ratio of time spent beyond given distence
timeratios <- function(dtables, cutoffDist=9, name="dists1") {
	sapply(dtables, function(dt) {
		sum(dt[,name] > cutoffDist,na.rm=TRUE)/nrow(dt)
	})
}

timeTrpFAD <- list(
	nofol.trp=timeratios(trpdist,9,"dists1")*100,
	nofol.gln=timeratios(glndist,9,"dists1")*100,
	nofol.tyr=timeratios(tyrdist,9,"dists1")*100,
	fol.trp=timeratios(trpfoldist,9,"dists1")*100,
	fol.gln=timeratios(glnfoldist,9,"dists1")*100,
	fol.tyr=timeratios(tyrfoldist,9,"dists1")*100
)
timeFADbarrel <- list(
	nofol.trp=timeratios(trpdist,7,"dists2")*100,
	nofol.gln=timeratios(glndist,7,"dists2")*100,
	nofol.tyr=timeratios(tyrdist,7,"dists2")*100,
	fol.trp=timeratios(trpfoldist,7,"dists2")*100,
	fol.gln=timeratios(glnfoldist,7,"dists2")*100,
	fol.tyr=timeratios(tyrfoldist,7,"dists2")*100
)

pvalsTrpFAD <- do.call(c,lapply(2:6, function(i){
	labels <- paste0(names(timeTrpFAD)[[i]],"-",names(timeTrpFAD)[1:(i-1)])
	setNames(sapply(1:(i-1),function(j){
		p <- wilcox.test(timeTrpFAD[[i]],timeTrpFAD[[j]])$p.value
	}),labels)
}))

pvalsFADbarrel <- do.call(c,lapply(2:6, function(i){
	labels <- paste0(names(timeFADbarrel)[[i]],"-",names(timeFADbarrel)[1:(i-1)])
	setNames(sapply(1:(i-1),function(j){
		p <- wilcox.test(timeFADbarrel[[i]],timeFADbarrel[[j]])$p.value
	}),labels)
}))

pvalsTrpFAD < 0.05
pvalsFADbarrel < 0.05


pdf("vis/dist_timeShare.pdf",10,10)
layout(cbind(1:2))
op <- par(bty="n")
beeswarm(timeTrpFAD,
	col=c("gray30","orange","steelblue3"),pch=20,
	ylab="% time spent at distance > 9A",
	labels=c(
		"WT fol-","p.Trp165Gln fol-","p.Trp165Tyr fol-",
		"WT fol+","p.Trp165Gln fol+","p.Trp165Tyr fol+"
	),
	main=expression(165*C^alpha - "flavin")
)
grid(NA,NULL)
bxplot(timeTrpFAD,col="gray",add=TRUE)
beeswarm(timeFADbarrel,
	col=c("gray30","orange","steelblue3"),pch=20,
	ylab="% time spent at distance > 7A",
	labels=c(
		"WT fol-","p.Trp165Gln fol-","p.Trp165Tyr fol-",
		"WT fol+","p.Trp165Gln fol+","p.Trp165Tyr fol+"
	),
	main=expression("flavin" - "active site")
)
bxplot(timeFADbarrel,col="gray",add=TRUE)
grid(NA,NULL)
par(op)
invisible(dev.off())







################################
# Binding mode cluster analysis
################################


#define cross product
`%X%` <- function(a, b) {
    i1 <- c(2,3,1)
    i2 <- c(3,1,2)
    return (a[i1]*b[i2] - a[i2]*b[i1])
}
#L2 norm (Euclidian)
l2 <- function(x) {
	sqrt(sum(x^2))
}

#normalize vector unit length
normalize <- function(x) x/l2(x)

#transform rotation matrix to euler angles
mat2euler <- function(m) {
	phi <- acos((m[1,1]+m[2,2]+m[3,3]-1)/2)
	d <- (2*sin(phi))
	x <- (m[3,2]-m[2,3])/d
	y <- (m[1,3]-m[3,1])/d
	z <- (m[2,1]-m[1,2])/d
	c(phi=phi,x=x,y=x,z=z)
}

#transform rotation matrix to quaternion
mat2quat <- function(m) {
	qr <- sqrt(1+m[1,1]+m[2,2]+m[3,3])/2
	qi <- (m[3,2]-m[2,3])/(4*qr)
	qj <- (m[1,3]-m[3,1])/(4*qr)
	qk <- (m[2,1]-m[1,2])/(4*qr)
	c(qr=qr,qi=qi,qj=qj,qk=qk)
}


#function to extract geometry characterization from a simulation frame
#based on atomic coordinates
characterize <- function(fad.a, fad.b, fad.c, trp.a, trp.b, trp.c){

	fad.ab.norm <- normalize(fad.b-fad.a)
	fad.ac.norm <- normalize(fad.c-fad.a)
	fad.orthonorm <- normalize(fad.ab.norm %X% fad.ac.norm)
	fad.ac.prime <- normalize(fad.orthonorm %X% fad.ab.norm)

	A <- cbind(fad.ab.norm,fad.ac.prime,fad.orthonorm)

	trp.ab.norm <- normalize(trp.b-trp.a)
	trp.ac.norm <- normalize(trp.c-trp.a)
	trp.ortho <- trp.ab.norm %X% trp.ac.norm
	trp.ac.prime <- trp.ortho %X% trp.ab.norm

	B <- cbind(trp.ab.norm,trp.ac.prime,trp.ortho)

	#invert A
	Ainv <- solve(A)

	#rebase B in terms of A
	Bprime <- Ainv%*%B
	relative.location <- Ainv%*%(trp.a-fad.a)
	Bquat <- mat2quat(Bprime)

	return(c(relative.location,Bquat))
}

#function to determine the matrix index of a given atom
makeResiTable <- function(metaData,offset=35) {
	labels <- metaData$RESIDUE_LABEL[metaData$RESIDUE_LABEL != ""]
	atomPointer <- metaData$RESIDUE_POINTER
	pos <- seq_along(labels)+offset

	data.frame(label=labels,pos=pos,atomPointer=atomPointer)
}

makeAtomTable <- function(metaData,xyz) {
	aname <- metaData$ATOM_NAME[metaData$ATOM_NAME != ""]
	atoms <- data.frame(
		num=1:length(aname),atom=aname,
		x=xyz[,1],y=xyz[,2],z=xyz[,3],
		elem=substr(aname,1,1)
	)
}


#residues: label, pos, atomPointer
#atoms: x, y, z, 
write.pdb <- function(resis,atoms,outfile) {

	library("gdata")

	resiExpand <- with(resis,{
		resl <- c(
			sapply(2:length(atomPointer),function(i) {
				atomPointer[[i]]-atomPointer[[i-1]]
			}),
			nrow(atoms)+1-atomPointer[[length(atomPointer)]]
		)
		do.call(rbind,mapply(function(resl,label,pos) {
			data.frame(aa=rep(label,resl),seqpos=rep(pos,resl))
		},resl,label,pos,SIMPLIFY=FALSE))
	})

	atomLines <- data.frame(tag="ATOM",num=atoms$num,sep1="",
		atom.name=atoms$atom,
		aloc="",aa=resiExpand$aa,sep2="",chain="A",
		seq.pos=resiExpand$seqpos,
		ins="",sep3="",round(atoms[,c("x","y","z")],digits=3),occ=1.0,temp="",sep4="",
		segment="A",element=atoms$elem,charge=""
	)
	atomLines[atomLines$aa == "FAD","segment"] <- "B"

	#fixed-width definition of atom entries
	widths <- c(
		tag=6,num=5,sep1=1,
		atom.name=4,aloc=1,
		aa=3,sep2=1,chain=1,seq.pos=4,
		ins=1,sep3=3,x=8,y=8,z=8,
		occ=6,temp=6,sep4=6,segment=4,
		element=2,charge=2
	)
	# tag     num atom aa chn seqpos  x       y        z      occ  temp   segm    elem  charge
	# ATOM      1  N   THR A  89      -6.611  24.602   9.492  1.00 66.84           N  
	# ATOM      2  CA  THR A  89      -5.936  25.604   8.610  1.00 67.02           C  
	# ATOM      3  C   THR A  89      -5.630  26.953   9.294  1.00 66.96           C  
	# CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1   
	# ATOM      1 PA   FAD A  36      68.763  87.228  40.698     1            B   P   
	# ATOM      2 O1A  FAD A  36      69.680  88.290  41.199     1            B   O   


	#cryst1 entry for pdb files         1   "
	cryst1 <- function(x,y,z,phi,psi,xi) {
		sprintf("CRYST1    %.3f    %.3f    %.3f  %.2f  %.2f  %.2f P 1           1   ",x,y,z,phi,psi,xi)
	}

	crystLine <- cryst1(1,1,1,90,90,90)

	con <- file(outfile,open="w")
	writeLines(crystLine,con)
	write.fwf(atomLines, con, colnames=FALSE, width=widths, sep="",append=TRUE)
	writeLines("END",con)
	close(con)
	
}


#process WT simulations
parentdir <- "mthfr_simulations/mthfr/model11/"
#select the metadata file
mdFile <- paste0(parentdir,"parm_strip")
#parse the metadata file
cat("Parsing metadata...\n")
metaData <- parseMetadata(mdFile)

#atoms of interest 
fadN10 <- getAtomIndex(metaData,NA,"FAD","N10")
fadC9A <- getAtomIndex(metaData,NA,"FAD","C9A")
fadC10 <- getAtomIndex(metaData,NA,"FAD","C10")

trpCG <- getAtomIndex(metaData,165,"TRP","CG")
trpCD1 <- getAtomIndex(metaData,165,"TRP","CD1")
trpCD2 <- getAtomIndex(metaData,165,"TRP","CD2")

#define the list of atom groups to compare against each other
atomList <- c(fadN10,fadC9A,fadC10,trpCG,trpCD1,trpCD2)

#find NetCDF files in parent directory
nc.files <- list.files(parentdir,pattern="nc$",full.names=TRUE)
#and calculate distance tables for C1' of FAD and CB of res165
cat("Processing trajectories...\n")
# dtables <- distances(nc.files, atomGroups=atomGroups)

characteristics <- lapply(nc.files, function(nc.file) {
	#open NetCDF file connection
	ncCon <- nc_open(nc.file)
	#extract multidimensional coordinate array
	#dimension 1 : x-y-z; dimension 2: atom index; dimension 3: frames
	coordinates <- ncvar_get(ncCon,"coordinates")
	#extract time table
	timepoints <- ncvar_get(ncCon,"time")
	#close file connection
	nc_close(ncCon)

	nframes <- dim(coordinates)[[3]]
	characteristics <- do.call(rbind,lapply(1:nframes,function(fr) {
		atomCoords <- coordinates[,atomList,fr]
		do.call(characterize,lapply(1:length(atomList),function(i)atomCoords[,i]))
	}))

	return(characteristics)
})

allchar <- do.call(rbind,characteristics)
colnames(allchar) <- c("x","y","z","r","i","j","k")

######################
#Perform clustering

mc3 <- Mclust(allchar,G=3)
smallest <- which.min(table(mc3$classification))
smallest.idx <- which(mc3$classification == smallest)

subset <- allchar[smallest.idx,]
mcSub <- Mclust(subset,G=6)

combined.classes <- mc3$classification
combined.classes[smallest.idx] <- mcSub$classification+2

#Extract avatars for each cluster
combined.centroids <- cbind(mc3$parameters$mean[,1:2],mcSub$parameters$mean)
avatars <- apply(combined.centroids,2,function(cen) {
	dists <- apply(allchar,1,function(vec) sqrt(sum((vec-cen)^2)))
	which.min(dists)
})

cumuFrames <- c(0,cumsum(sapply(characteristics,nrow)))
avatarIdx <- do.call(rbind,lapply(avatars, function(x) {
	last <- max(which(cumuFrames < x))
	c(sim=last,frame=x-cumuFrames[[last]])
}))

lapply(1:length(avatars), function(i) {
	cat("Extracting avatar",i,"\n")
	nc.file <- nc.files[[avatarIdx[i,"sim"]]]
	ncCon <- nc_open(nc.file)
	xyz <- t(ncvar_get(ncCon,"coordinates")[,,avatarIdx[i,"frame"]])
	nc_close(ncCon)
	atoms <- makeAtomTable(metaData,xyz)
	resis <- makeResiTable(metaData)
	write.pdb(resis,atoms,paste0("results/cluster_avatar_",i,".pdb"))
})


#####################
# Visualize clusters
#####################
pc <- princomp(allchar,cor=TRUE)
pc2 <- princomp(allchar)
um <- umap::umap(allchar)

mycols <- c("black","firebrick3","darkolivegreen3","steelblue3","orange","purple3","gray50","magenta","cyan")
mycols <- sapply(mycols,colAlpha,0.2)

pdf("vis/clustering_pairs.pdf")
pairs(allchar,pch=".",col=mycols[combined.classes])
dev.off()

pdf("vis/clustering_pairs_small.pdf")
pairs(allchar[smallest.idx,],pch=".",col=mycols[mcSub$classification+2])
dev.off()

pdf("vis/clustering_pca_umap.pdf",15,5)
op <- par(mfrow=c(1,3))
plot(pc$scores[,1:2],pch=20,col=mycols[combined.classes])
plot(pc2$scores[,1:2],pch=20,col=mycols[combined.classes])
plot(um$layout,col=mycols[combined.classes],pch=20)
par(op)
dev.off()

########################
# Examine state transitions
############################
transitions <- do.call(c,lapply(2:length(combined.classes), function(i) {
	if (combined.classes[[i]] != combined.classes[[i-1]]) {
		paste0(combined.classes[[i-1]],"-",combined.classes[[i]])
	} else NULL
}))
trCounts <- table(transitions)

trTable <- matrix(0,nrow=8,ncol=8)
for (i in 1:length(trCounts)) {
	fromTo <- as.integer(strsplit(names(trCounts)[[i]],"-")[[1]])
	trTable[fromTo[[1]],fromTo[[2]]] <- trCounts[[i]]
}

probTable <- apply(trTable,1,function(x) x/sum(x))
nodeTable <- table(combined.classes)/length(combined.classes)

#########################
#draw transition graph

#function to calculate circle positions
circPos <- function(m,r,phis) {
	do.call(rbind,lapply(phis,function(phi) {
		c(
			x = m[[1]]+r*cos(phi),
			y = m[[2]]+r*sin(phi)
		)
	}))
}

#function to generate phi angles
genPhis <- function(n)seq(0,2,length.out=n+1)[-1]*pi

pdf("vis/transitionGraph.pdf")
ms <- circPos(m=c(0,0),r=3.5,phis=genPhis(8))

mycols <- c("black","firebrick3","darkolivegreen3","steelblue3","orange","purple3","gray50","magenta","cyan")
op <- par(mar=c(0,0,0,0))
plot(NA,type="n",xlim=c(-5,5),ylim=c(-5,5),axes=FALSE,xlab="",ylab="")

lapply(1:nrow(probTable), function(i) {
	lapply(1:ncol(probTable), function(j) {
		if (probTable[i,j] > 0) {
			segments(ms[i,1],ms[i,2],ms[j,1],ms[j,2],lwd=probTable[i,j]*10,col="gray")
		}
	})
})

lapply(1:length(nodeTable), function(i) {
	radius <- sqrt(10*nodeTable[[i]]/pi)
	xy <- circPos(ms[i,],radius,genPhis(50))
	polygon(xy[,1],xy[,2],col=mycols[[i]],border=NA)
	text(ms[i,1],ms[i,2],i,col="white")
})
par(op)
dev.off()

