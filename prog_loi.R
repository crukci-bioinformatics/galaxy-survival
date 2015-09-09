# --- Command-line arguments (for use in Galaxy) -------------------------------

ca <- commandArgs(trailingOnly=TRUE)

infile = ca[1]
outfile = ca[2]

# ------------------------------------------------------------------------------

library(reshape)

xfs <- read.table(infile, sep="\t", header=T, as.is=TRUE)
xfs2 <- xfs[!is.na(xfs$HighIsGood),]

table.xfs <- table(xfs2$GeneName, xfs2$HighIsGood)
table.xfs <- data.frame(table.xfs[,1], table.xfs[,2])
colnames(table.xfs) <- c("False", "True")

both <- table.xfs[(table.xfs[,1]>0 & table.xfs[,2]>0),]
both.xfs <- row.names(both)
both.rep <- rep("Both", length(both.xfs))
both.list <- cbind(both.xfs, both.rep)

true <- table.xfs[(table.xfs[,1]==0 & table.xfs[,2]>0),]
true.xfs <- row.names(true)
true.rep <- rep("Good", length(true.xfs))
true.list <- cbind(true.xfs, true.rep)

false <- table.xfs[(table.xfs[,1]>0 & table.xfs[,2]==0),]
false.xfs<-row.names(false)
false.rep <- rep("Poor", length(false.xfs))
false.list <- cbind(false.xfs, false.rep)

all.list <- rbind(both.list, true.list, false.list)
colnames(all.list) <- c("GeneName", "Prognosis")

write.table(all.list, outfile, quote=FALSE, row.names=FALSE, sep="\t")
