# --- Command-line arguments (for use in Galaxy) -------------------------------

ca <- commandArgs(trailingOnly=TRUE)

datafile            <- ca[1]
infile              <- ca[2]
col                 <- as.numeric(ca[3])
if (is.na(col)) col <- 1
outfile_txt         <- ca[4]
outfile_html        <- ca[5]
outfile_path        <- ca[6]
pcutoff             <- as.numeric(ca[7])

# --- Hard-wired parameters ----------------------------------------------------

#datafile     <- 'GlnData.Rda'
#infile       <- 'genelist.txt'
#col          <- 1
#outfile_txt  <- 'results.txt'
#outfile_html <- 'results.html'
#outfile_path <- "results"
#pcutoff      <- 0.1
# ----------------------------------------------------

library(hgu133a.db)
library(party)
library(survival)

load(datafile)

dir.create(outfile_path)

cat("<html>\n", file=outfile_html)
cat("<body>\n", file=outfile_html, append=TRUE)
cat("<h3>Prostate Cancer Recursive Partitioning</h3>\n", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat("The following are the results of recursive partitioning for input genes", file=outfile_html, append=TRUE)
cat(" using a prostate cancer gene expression dataset from Glinsky <em>et al.</em>,", file=outfile_html, append=TRUE)
cat(" <em>J. Clin. Invest.</em> 113, 913-923 (2004) with accompanying", file=outfile_html, append=TRUE)
cat(" survival data.", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat(" This dataset comprises expression profiles for primary tumours from 79 patients.", file=outfile_html, append=TRUE)
cat(" The data were generated with Affymetrix Human Genome U133A microarrays. ", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat("Kaplan-Meier survival curves and recursive partitioning plots are",      file=outfile_html, append=TRUE)
cat(" displayed below for those genes which are predictive of recurrence",    file=outfile_html, append=TRUE)
cat(paste(" free survival with p &lt; ", pcutoff),                            file=outfile_html, append=TRUE)
cat(" (corrected for the testing of multiple cut-offs, but not genes). ",      file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat("Warning: P-values can vary considerably if a study is repeated from ",    file=outfile_html, append=TRUE)
cat("sampling alone. You may have a statistically significant result but ",    file=outfile_html, append=TRUE)
cat("is the effect size of importance? The p-value can help you answer the ",  file=outfile_html, append=TRUE)
cat("yes/no question but the question you often want to ask is the how much ", file=outfile_html, append=TRUE)
cat("question. See Halsey et al., Nat Methods 12(3):179-85, 2015.", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

geneList <- read.table(infile, sep="\t", as.is=TRUE)[,col]
geneList <- toupper(geneList)
geneList <- unique(geneList)
# geneList <-c("AURKA", "AURKB", "CDC20", "MELK", "PLK1", "CCNE1","CENPE", "CDC2", "RAD17","BRCA1", "HES6", "UBE2C", "CCNB1", "CCNB2")
probes   <- mget(geneList, hgu133aALIAS2PROBE, ifnotfound=NA)
AffyID   <- probes[[1]]

names    <- rep(geneList[1], length(probes[[1]]))

if(length(geneList)>1) {
  for(j in 2:length(geneList)) {
    AffyID <- c(AffyID, probes[[j]])
    names <- c(names, rep(geneList[j],length(probes[[j]])))
  }
}

probetable      <- cbind(names, AffyID)

uniqueProbeList <- na.omit(AffyID)
pvalList        <- NA
pvalList2       <- NA
splitList       <- NA
splitList2      <- NA
highIsGoodList  <- NA
symbolList      <- uniqueProbeList[1]
idlist          <- NA

if (length(uniqueProbeList) > 0) {
  for (j in 1:length(uniqueProbeList)) {
    # get gene symbol and affy ID#
    genenum       <- grep(uniqueProbeList[[j]], probetable[, 2], ignore.case=TRUE)
    symbol2       <- probetable[genenum, 1]
    symbol        <- symbol2[[1]]
    symbolList[j] <- symbol
    symbolProbe   <- paste(symbol, "-", uniqueProbeList[[j]])
    pngprefix     <- paste(symbol, "_", uniqueProbeList[[j]], sep="")
    rowID         <- match(uniqueProbeList[[j]], colnames(data.gln))
    idlist[j]     <- rowID
    rowID         <- na.omit(rowID)
    
    if(length(rowID)>0) {
      
      combined.data <- cbind(demo.gln, symbolProbe=data.gln[, rowID])
      surv.xfs      <- Surv((combined.data$DFI), combined.data$RELAPSE=="Y")
      mincriterion       <- 1 - pcutoff # mincriterion is 1-pvalue, this must exceeded to implement a split
      
      xlabel <- "Disease Free Interval (months)"
      ylabel <- "Probability of Relapse Free Survival"
      title  <- "Relapse Free Survival"
      
      ctree_xfs     <- ctree(surv.xfs~symbolProbe, data=combined.data, controls=ctree_control(mincriterion=mincriterion))
      pvalue      <- 1 - ctree_xfs@tree$criterion$maxcriterion
      newPval     <- signif(pvalue, digits = 2)
      pvalList[j] <- newPval
      ps3         <- NA
      highIsGood  <- NA
      
      if(newPval<pcutoff) {
        rp_pngfile <- paste(pngprefix, "rp.png", sep="_")
        png(filename = paste(outfile_path, rp_pngfile, sep="/"))
        plot(ctree(surv.xfs~symbolProbe, data=combined.data, control=ctree_control(mincriterion=mincriterion)), xlab=xlabel, ylab=ylabel, main=paste(title, symbolProbe))
        dev.off()
        
        ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
        ps  <- signif(ps2[1], digits = 3)
        
        if(length(ps2)==1) {
          combined.data$symbolProbe_cp <- combined.data$symbolProbe<=ps2[1]
          coxphRegressionModel         <- coxph(surv.xfs~combined.data$symbolProbe_cp)
          hazardRatio                  <- exp(coxphRegressionModel$coefficients)
          highIsGood                   <- hazardRatio > 1
          nt                           <- table(combined.data$symbolProbe_cp)
          symbolProbe.survfit.xfs      <- survfit(surv.xfs~combined.data$symbolProbe_cp)
          newPval2                     <- NA
          
          km_pngfile <- paste(pngprefix, "km.png", sep="_")
          png(filename = paste(outfile_path, km_pngfile, sep="/"))
          plot(symbolProbe.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(symbolProbe,", p=", newPval), col=c(2, 4))
          legend(0.1, 0.1, c(paste(symbol, ">", ps, "n=", nt[[1]]), paste(symbol, "<=", ps, "n=", nt[[2]])), col=c(2, 4), lty=1, lwd=1.5, bty="n")
          dev.off()
          
          cat("<img src=\"", file=outfile_html, append=TRUE)
          cat(km_pngfile, file=outfile_html, append=TRUE)
          cat("\">\n", file=outfile_html, append=TRUE)
        }
        
        if(length(ps2)==2) {
          if(ps2[1]>ps2[2]) {
            ps3                           <- round(ps2[2], digits=3)
            combined.data$symbolProbe_cp1 <- combined.data$symbolProbe<=ps2[1]
            combined.data$symbolProbe_cp2 <- combined.data$symbolProbe<=ps2[2]
            combined.data$symbolProbe_cp3 <- combined.data$symbolProbe_cp1+combined.data$symbolProbe_cp2
            nt                            <- table(combined.data$symbolProbe_cp3)
            symbolProbe.survfit.xfs       <- survfit(surv.xfs~combined.data$symbolProbe_cp3)
            pvalue2                       <- 1 - ctree_xfs@tree$left[[3]][[2]]
            newPval2                      <- signif(pvalue2, digits=2)
            
            if(newPval2>=pcutoff) {
              km_pngfile <- paste(pngprefix, "km.png", sep="_")
              png(filename = paste(outfile_path, km_pngfile, sep="/"))
              plot(symbolProbe.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(symbolProbe,", p=", newPval), col=c(2, 4, 5))
              legend(0.1, 0.15, c(paste(symbol, ">", ps,  "n=", nt[[1]]), paste(symbol, "<=", ps, "& >", ps3,  "n=", nt[[2]]), paste(symbol, "<=", ps3,  "n=", nt[[3]])), col=c(2, 4, 5), lty=1, lwd=1.5, bty="n")
              dev.off()
            } else {
              km_pngfile <- paste(pngprefix, "km.png", sep="_")
              png(filename = paste(outfile_path, km_pngfile, sep="/"))
              plot(symbolProbe.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(symbolProbe,", p=", newPval, "& p=", newPval2), col=c(2, 4, 5))
              legend(0.1, 0.15, c(paste(symbol, ">", ps,  "n=", nt[[1]]), paste(symbol, "<=", ps, "& >", ps3,  "n=", nt[[2]]), paste(symbol, "<=", ps3,  "n=", nt[[3]])), col=c(2, 4, 5), lty=1, lwd=1.5, bty="n")
              dev.off()
            }
            
            cat("<img src=\"", file=outfile_html, append=TRUE)
            cat(km_pngfile, file=outfile_html, append=TRUE)
            cat("\">\n", file=outfile_html, append=TRUE)
          } else {
            ps3                           <- round(ps2[2], digits=3)
            combined.data$symbolProbe_cp1 <- combined.data$symbolProbe<=ps2[1]
            combined.data$symbolProbe_cp2 <- combined.data$symbolProbe<=ps2[2]
            combined.data$symbolProbe_cp3 <- combined.data$symbolProbe_cp1+combined.data$symbolProbe_cp2
            nt                            <- table(combined.data$symbolProbe_cp3)
            symbolProbe.survfit.xfs       <- survfit(surv.xfs~combined.data$symbolProbe_cp3)
            pvalue2                       <- 1 - ctree_xfs@tree$right[[3]][[2]]
            newPval2                      <- signif(pvalue2, digits=2)
            
            if(newPval2>=pcutoff) {
              km_pngfile = paste(pngprefix, "km.png", sep="_")
              png(filename=paste(outfile_path, km_pngfile, sep="/"))
              plot(symbolProbe.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(symbolProbe,", p=", newPval), col=c(2, 4, 5))
              legend(0.1, 0.15, c(paste(symbol, ">", ps3,  "n=", nt[[1]]), paste(symbol, ">", ps, "& <=", ps3,  "n=", nt[[2]]), paste(symbol, "<=", ps,  "n=", nt[[3]])), col=c(2, 4, 5), lty=1, lwd=1.5, bty="n")
              dev.off()
            } else {
              km_pngfile <- paste(pngprefix, "km.png", sep="_")
              png(filename = paste(outfile_path, km_pngfile, sep="/"))
              plot(symbolProbe.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(symbolProbe,", p=", newPval, "& p=", newPval2), col=c(2, 4, 5))
              legend(0.1, 0.15, c(paste(symbol, ">", ps3,  "n=", nt[[1]]), paste(symbol, ">", ps, "& <=", ps3,  "n=", nt[[2]]), paste(symbol, "<=", ps,  "n=", nt[[3]])), col=c(2, 4, 5), lty=1, lwd=1.5, bty="n")
              dev.off()
            }
            
            cat("<img src=\"", file=outfile_html, append=TRUE)
            cat(km_pngfile, file=outfile_html, append=TRUE)
            cat("\">\n", file=outfile_html, append=TRUE)
          }
        }
        
        cat("<img src=\"", file=outfile_html, append=TRUE)
        cat(rp_pngfile, file=outfile_html, append=TRUE)
        cat("\">\n", file=outfile_html, append=TRUE)
        cat("<p>\n", file=outfile_html, append=TRUE)
      }
      
      if(newPval>=pcutoff) {
        ps       <- NA
        newPval2 <- NA
      }
      
      pvalList2[j]      <- newPval2
      splitList[j]      <- ps
      splitList2[j]     <- ps3
      highIsGoodList[j] <- highIsGood
    }
  }
  
  missinggene     <- setdiff(geneList, unique(symbolList))
  missinggenelist <- NULL
  
  if (length (missinggene)>0) {
    for (k in 1:length(missinggene)) {
      missinggenetemp<-c(missinggene[k], rep(NA, 6))
      missinggenetemp<-as.matrix(t(missinggenetemp))
      missinggenelist<-rbind(missinggenetemp, missinggenelist)
    }
  }
  
  GeneName   <- symbolList
  AffyID     <- uniqueProbeList
  PValue1    <- pvalList
  PValue2    <- pvalList2
  CutOff1    <- splitList
  CutOff2    <- splitList2
  HighIsGood <- highIsGoodList
  
  xfsList    <- cbind(GeneName, AffyID, CutOff1, CutOff2, PValue1, PValue2, HighIsGood)
  
  if (length (missinggene)>0) {
    xfsList  <- rbind(xfsList,missinggenelist)
  }
  
  write.table(xfsList, outfile_txt, quote=FALSE, row.names=FALSE, sep="\t")
  
} else {
  cat("<h3>No matches to Affymetrix probes found for query genes</h3><p>\n", file=outfile_html, append=TRUE)
}

cat("</body>\n", file=outfile_html, append=TRUE)
cat("</html>\n", file=outfile_html, append=TRUE)

