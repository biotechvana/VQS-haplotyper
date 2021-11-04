 
flnms <- list.files(flashDir)
snms <- sub("\\.fastq$","",flnms)
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
parts <- parts[,-c(3,4)]
colnames(parts) <- c("PatID","SmplID")

bxp.dt <- list(stats=matrix(0,nrow=5,ncol=length(flnms)),
               names=parts[,1])
res <- matrix(0,nrow=length(flnms),ncol=8)
colnames(res) <- c("TotReads","ShortRds","AvQV","SdQV",
                   "MedU33","MedU30","MedU20","MedU13")
rownames(res) <- paste(parts[,1],parts[,2],sep=".")
for(iii in 1:length(flnms))
{
  ###  Load fastq file
  sqq <- readFastq(dirPath=flashDir,patt=flnms[iii])

  ###  Sequences as a DNAStringSet
  seqs <- sread(sqq)
  ###  Sequence lengths
  lns <- width(seqs)
  ###  Quality score
  qs <- quality(sqq)

  ###  Codi ASCII-33
  phrsc <- as(qs,"matrix")

  ###  Eliminar seqs massa curtes
  fl <- lns > min.len
  n.tot <- length(lns)
  n.short <-  n.tot-sum(fl)
  seqs <- seqs[fl]
  lns <- lns[fl]
  phrsc <- phrsc[fl,]

  res[iii,1] <- n.tot
  res[iii,2] <- n.short
  
  pdf.flnm <- file.path(repDir,paste(snms[iii],"pdf",sep="."))
  pdf(pdf.flnm,paper="a4r",width=10,height=7)
  par(mfrow=c(2,2),mar = c(5,5,2,5))
  txt.flnm <- file.path(repDir,paste(snms[iii],"txt",sep="."))
  sink(txt.flnm)
  
  cat("\nFASTQ analysis: ",flnms[iii],"\n")
  cat("\nTotal number of reads: ",n.tot,"\n")
  msg <- paste("\nNumber of short reads (<",min.len,"bp) filtered out: ",
               sep="")
  cat(msg,n.short," (",round(n.short/n.tot*100,2),"%)\n",sep="")
  cat("\nRead lengths quantiles\n\n")

  ###  Estudi longituds
  print(quantile(lns,p=c(0.025,0.05,0.25,0.5,0.75,0.90,0.95,0.975)))

  frq <- table(lns)
  ln <- as.integer(names(frq))
  xmx <- quantile(lns,p=0.95)

  ###  Distribució de longituds amb quantils (per amplicons)
  plot(ln,as.vector(frq),type="h",xlab="Read length",ylab="Frequency",
       lwd=2,main="Read length distribution",xlim=c(0,max(ln)))
  abline(h=0)

  ###  Freqüència acumulada de longituds amb quantils
  yp <- rev(cumsum(rev(as.vector(frq))))
  xmx <- quantile(lns,p=0.98)
  plot(ln,yp,type="l",xlab="read length",ylab="Cumulative frequency",
       lwd=2,xlim=c(min(ln),xmx),ylim=c(0,max(yp)),
	   main="Read length distribution",)	   
  qln <- quantile(lns,p=c(0.5,0.25,0.1))
  abline(v=qln[1],lty=4,col="gray")
  abline(v=qln[2],lty=4,col="gray")
  abline(v=qln[3],lty=4,col="gray")
  text(qln[1],0.5*yp[1],paste("50% cov =",qln[1]),cex=0.6,font=2)
  text(qln[2],0.75*yp[1],paste("75% cov =",qln[2]),cex=0.6,font=2)
  text(qln[3],0.9*yp[1],paste("90% cov =",qln[3]),cex=0.6,font=2)
  text(min(ln),0,paste(yp[1],"reads"),cex=0.8,adj=c(0,0),font=2)

  ###  Inspecció de la qualitat dels reads en promig
  mean.qv <- apply(phrsc,1,mean,na.rm=TRUE)
  res[iii,3] <- round(mean(mean.qv),1)
  res[iii,4] <- round(sd(mean.qv),1)
  bxp.dt$stats[,iii] <- fivenum(mean.qv)
  
  cat("\nSummary of read mean QV\n\n")
  print(summary(mean.qv))
  Fn <- ecdf(mean.qv)
  o <- order(mean.qv)
  plot(mean.qv[o],1-Fn(mean.qv[o]),type="l",
       xlab="Mean read Phed score",
       ylab="Fraction of reads",
       main="Mean read QV")
  grid()
  abline(v=median(mean.qv),lty=5,col="blue")

  Pqv <- function(qv) 1-10^(-qv/10)
  p.err <- function(qv) 10^(-qv/10)

  pred.err <- apply(phrsc,1,function(qv) 
                sum(p.err(qv),na.rm=TRUE)/sum(!is.na(qv)))
	 
  cat("\nSummary of predicted consensus accuracy\n\n")
  print(signif(summary(1-pred.err),3))
  Fn <- ecdf(pred.err)
  o <- order(pred.err)

  hres <- hist(1-pred.err[o],breaks=100,col="lavender",freq=TRUE,main="",
               xlim=c(0.85,1),xlab="Predicted accuracy",ylab="Reads")
  abline(v=median(1-pred.err),lty=5,col="blue")
  ypos <- max(hres$counts)
  text(median(1-pred.err),ypos,adj=0.5,"median ",col="blue",cex=0.8,font=2)
  par(new=T)
  plot(1-pred.err[o],Fn(pred.err[o]),type="l",ylim=c(0,1),xlim=c(0.85,1),
       axes=F,xlab=NA,ylab=NA,col="blue",
       # xlab="Predicted accuracy",
       # ylab="Fraction of reads below accuracy",
       main="Consensus predicted accuracy")
  axis(side=4,col="blue",col.axis="blue",at=seq(0,1,0.1),las=2)
  grid()
  abline(h=seq(0,1,0.1),lty=3,col="gray")

  ###  Inspecció de la qualitat dels reads per la qualitat de les seves bases
  p.true <- rev(c(0.95,0.99,0.999,0.9995))
  phred.cut <- round(-10 * log10(1-p.true),1)
  cat("\nStatistics of number of bases per read below a QV\n")
  for(i in 1:length(phred.cut))
  { xtt <- paste("# bases below Phred ",phred.cut[i],
                 " (",p.true[i]*100,")",sep="")
    nb.low <- sapply(1:nrow(phrsc),function(k) 
                    sum(phrsc[k,!is.na(phrsc[k,])]<phred.cut[i]))
    tbl <- table(nb.low)
    plot(as.integer(names(tbl)),as.vector(tbl),type="h",xlab=xtt,ylab="# reads",
	       xlim=c(0,quantile(nb.low,p=0.95)))
 	par(new=T)
    frac <- cumsum(tbl)/sum(tbl)
    plot(as.integer(names(tbl)),frac,type="l",axes=F,xlab=NA,ylab=NA,col="blue",
	     xlim=c(0,quantile(nb.low,p=0.95)))
    axis(side=4,col="blue",col.axis="blue",at=seq(0,1,0.1),las=2)
	grid(ny=NA)
    abline(h=c(seq(0.5,1,0.1),seq(0.8,1,0.05)),lty=3,col="gray")	
    mtext(side=4,line=3,'Cumulative frequency',col="blue",cex=0.8)	

    cat("\nQuantiles for Phred ", phred.cut[i]," (",p.true[i]*100,")\n\n",sep="")
    print(quantile(nb.low,p=c(0.5,0.6,0.7,0.8,0.90,0.95,0.99)))
	res[iii,4+i] <- median(nb.low)
  }
  sink()		 
  dev.off()
}

rownames(bxp.dt$stats) <- c("low.w","low.h","median","high.h","high.w")
colnames(bxp.dt$stats) <- bxp.dt$names

txt.flnm <- file.path(repDir,paste("FLASH_",proj.nm,"_EDA.txt",sep=""))
sink(txt.flnm)
cat("\nExploratory data analysis of fastq files\n\n")
print(res)
cat("\nFivenums of mean read phred score by pool\n\n")
print(round(bxp.dt$stats,1))
sink()

pdf.flnm <- file.path(repDir,paste("FLASH_",proj.nm,"_EDA.pdf",sep=""))
pdf(pdf.flnm,paper="a4r",width=10,height=5)
par(mfrow=c(1,2))

barplot(res[,1],col="lavender",las=2,main="Reads by fastq file")

bxp(bxp.dt,pars=list(boxfilll="lavender"),border="navy",
    ylab="mean phred scores",las=2)
title(main="Read mean phred scores by pool")

par(mfrow=c(1,1))
library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")
mat <- res[,5:8]
ylm <- c(0,max(mat)*1.2)
barplot(t(mat),las=2,beside=TRUE,ylim=ylm,col=pal[1:4],
        main="median number of bases below QV by read")
legend("top",horiz=TRUE,fill=pal[1:4],legend=paste(p.true),cex=0.8)		
dev.off()
