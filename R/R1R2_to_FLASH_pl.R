
flnms <- list.files(runDir)
snms <- sub("\\.fastq$","",flnms)
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
parts <- parts[,-c(3,5)]
colnames(parts) <- c("PatID","SmplID","Read")

flnms <- file.path(runDir,flnms)
R1.flnms <- flnms[parts[,3]=="R1"]
R2.flnms <- flnms[parts[,3]=="R2"]
out.flnms <- paste(parts[parts[,3]=="R1",1],parts[parts[,3]=="R1",2],
                   "flash.fastq",sep="_")
out.flnms <- file.path(flashDir,out.flnms)

parts <- parts[parts[,3]=="R1",,drop=FALSE]
res <- matrix(0,nrow=length(out.flnms),ncol=2)
rownames(res) <- paste(parts[,1],parts[,2],sep="_")
colnames(res) <- c("Extended","NoExtd")

for(i in 1:length(out.flnms))
{ 
  ok <- file.exists(R1.flnms[i]) & file.exists(R2.flnms[i])
  if(!ok) next
  command <- paste(flash,flash.opts,R1.flnms[i],R2.flnms[i],collapse=" ")
  es <- system(command,intern=FALSE,wait=TRUE,show.output.on.console=FALSE,
              ignore.stdout=FALSE,invisible=TRUE)
  if(es!=0) next
  file.copy(from="out.extendedFrags.fastq",to=out.flnms[i],overwrite=TRUE)

  sqq <- readFastq(dirPath=".",patt="out.extendedFrags.fastq")
  res[i,1] <- length(sqq)
  sqq <- readFastq(dirPath=".",patt="out.notCombined_1.fastq")
  res[i,2] <- length(sqq)
}

df.res <- data.frame(res,Yield=round(res[,1]/(res[,1]+res[,2])*100,1))
txt.flnm <- file.path(repDir,"FLASH_report.txt")
sink(txt.flnm)
cat("\nExtending Illumina reads by FLASH\n")
cat("\nFLASH parameters:")
cat("\n    Minimum overlap:",min.ov)
cat("\n    Maximum overlap:",max.ov)
cat("\n        Error level:",err.lv,"\n\n")
print(df.res)
sink()

flash.res <- df.res
save(flash.res,file=file.path(repDir,"FLASH_table.RData"))

pdf.flnm <- file.path(repDir,"FLASH_barplot.pdf")
pdf(pdf.flnm,paper="a4",width=6,height=10)
par(mfrow=c(2,1))

library(RColorBrewer)
pal=brewer.pal(8,"Dark2")
M <- data.matrix(df.res[,1:2])
ymx <- max(M)*1.15
barplot(t(M),beside=TRUE,las=2,col=pal[1:2],ylim=c(0,ymx))
legend("top",horiz=TRUE,fill=pal[1:2],cex=0.8,
       legend=c("Extended","Not extended"))
title(main="FLASH results on paired-ends",line=2.5)
title(main="Yield in number of reads",cex.main=1,line=1)

bp <- barplot(df.res$Yield,col="Lavender",ylim=c(0,100),ylab="Yield",
              names.arg=rownames(df.res),las=2)	   
text(bp,10,paste(df.res$Yield,"%",sep=""),srt=90,col="navy")
title(main="Yield in percentage",cex.main=1,line=1)

dev.off()
