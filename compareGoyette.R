#!/usr/bin/Rscript
options(stringsAsFactors=FALSE)

if (!("yogitools" %in% installed.packages()[,1])) {
	stop("This script requires the package 'yogitools'\n=> https://github.com/jweile/yogitools")
}

fitted <- read.csv("results/folate_response_model5.csv")
rownames(fitted) <- fitted$hgvs

goyette <- read.csv("reference_data/goyette_rozen.csv")

joint <- cbind(goyette,fitted[goyette$AA.change,c("w200.score","m200.score","w.fitness","e.b","e.r","e.post.b","e.post.r")])

pdf("goyette.pdf",5,3)
op <- par(mar=c(5,4,1,1))
plot(NA,type="n",xlim=c(0,5),ylim=c(0,1.5),xlab="Goyette & Rozen 2000",ylab="VE map at 200ug/ml folate")
with(joint,points(WT.rel,w200.score,pch=20))
with(joint,points(A222V.rel,m200.score,pch=20,col="gray50"))
with(joint,arrows(WT.rel,w200.score,A222V.rel,m200.score,length=.05))
grid(NULL,NULL)
abline(h=0:1,v=0:1,col=c("firebrick3","chartreuse3"),lty="dashed")
legend("right",c("WT","A222V"),col=c("black","gray50"),pch=20,bg="white")
par(op)
dev.off()

#calculate correlation
cor(yogitools::fin(with(joint,cbind(
	c(WT.rel,A222V.rel),c(w200.score,m200.score)
))),method="spearman")
