# --- Command-line arguments (for use in Galaxy) -------------------------------

ca <- commandArgs(trailingOnly=TRUE)

datafile            <- ca[1]
infile              <- ca[2]
col                 <- as.numeric(ca[3])
if (is.na(col)) col <- 1
outfile_txt         <- ca[4]
outfile_html        <- ca[5]
outfile_path        <- ca[6]
margins             <- ca[7]
groups              <- ca[8]
gleason             <- ca[9]
tstage              <- ca[10]

# --- Hard-wired parameters ----------------------------------------------------

#datafile     <- 'TaylorNew.Rdata'
#infile       <- 'genelist.txt'
#col          <- 1
#outfile_txt  <- 'results.txt'
#outfile_html <- 'results.html'
#outfile_path <- "."
#margins      <- 'negative'
#margins      <- 'positive'
#margins      <- 'all'
#groups       <- 'primary'
#groups       <- 'all'
#groups       <- 'metastasis'
#gleason      <- 'low'
#gleason      <- 'high'
#gleason      <- 'all'
#tstage       <- 'two'
#tstage       <- 'threefour'
#tstage       <- 'all'

# ------------------------------------------------------------------------------


load(datafile)

library(party)
library(survival)

dir.create(outfile_path)

cat("<html>\n", file=outfile_html)
cat("<body>\n", file=outfile_html, append=TRUE)
cat("<h3>Taylor Prostate Cancer Biochemical Recurrence</h3>\n", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat("The following are the results of recursive partitioning for input genes", file=outfile_html, append=TRUE)
cat(" using a prostate cancer gene expression dataset from Taylor <em>et al.</em>,", file=outfile_html, append=TRUE)
cat(" <em>Cancer Cell</em>  18, 11-22, (2010) with accompanying", file=outfile_html, append=TRUE)
cat(" biochemical recurrence data.", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat(" This dataset comprises expression profiles for prostate cancer from 140 patients.", file=outfile_html, append=TRUE)
cat(" Consisting of 131 primary tumours and 9 metastatic tumours,", file=outfile_html, append=TRUE)
cat(" 107 surgical margin negative and 33 surgical margin positive.", file=outfile_html, append=TRUE)
cat(" The expression data were generated with Affymetrix Human Genome Exon microarrays.", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat("Kaplan-Meier survival curves and recursive partitioning plots are", file=outfile_html, append=TRUE)
cat(" displayed below for those genes which are predictive of biochemical recurrence", file=outfile_html, append=TRUE)

cat(" with p &lt; 0.05", file=outfile_html, append=TRUE)
cat(" (corrected for the testing of multiple cut-offs).", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat("The following options were selected: Margins:", paste(margins), file=outfile_html, append=TRUE)
cat(", Groups:", paste(groups), file=outfile_html, append=TRUE)
cat(", Gleason Grade:", paste(gleason), file=outfile_html, append=TRUE)
cat(", T-Stage:", file=outfile_html, append=TRUE)

if (tstage == 'two') {
  cat(" T2", file=outfile_html, append=TRUE)

} else if (tstage == 'threefour') {
  cat(" T3/T4", file=outfile_html, append=TRUE)

} else {
  cat(" all", file=outfile_html, append=TRUE)
}

cat("<p>\n", file=outfile_html, append=TRUE)

#geneList<-c("HES6")
geneList     <- read.table(infile, sep="\t", as.is=TRUE)[, col]
geneList     <- toupper(geneList)
geneList     <- unique(geneList)
rowidtemp    <- NA
affyidtemp   <- NA
namestemp    <- NA
accessidtemp <- NA
rowid        <- NULL
affyid       <- NULL
accessid     <- NULL
names        <- NULL

if(length(geneList)>0) {
  for(i in 1:length(geneList)) {
    gene = paste("^", geneList[i], "$", sep="")
    affyidtemp   <- rownames(annot.bcr)[grep(gene, annot.bcr[, 7])]
    affyid       <- c(affyid, affyidtemp)
    namestemp    <- rep(geneList[i], length(affyidtemp))
    names        <- c(names, namestemp)
    rowidtemp    <- grep(gene, annot.bcr[, 7])
    rowid        <- c(rowid, rowidtemp)
    accessidtemp <- annot.bcr[rowidtemp, 2]
    accessid     <- c(accessid, accessidtemp)
  }
}

table          <- cbind(names,rowid,affyid)
pvalList       <- NA
pvalList2      <- NA
splitList      <- NA
splitList2     <- NA
highIsGoodList <- NA


if (length(affyid) > 0) {
  for (j in 1:length(affyid)) {
    #get gene symbol and affy ID#
    geneid    <- paste(names[j], "-", accessid[j])
    pngprefix <- paste(names[j], "_", accessid[j], sep="")
   
    if(length(rowid)>0) {
      #Extract data for PgR etc#
      combined.data <- cbind(demo.bcr, geneexp=as.numeric(data.bcr[,(rowid[j])]))
       #********
      if (margins =='negative') {
        combined.data <- combined.data[combined.data$SMS == "Negative",]
      } else if (margins =='positive') {
        combined.data <- combined.data[combined.data$SMS == "Positive",]              	
      } else { 
       	combined.data <- combined.data
      }

      combined.data$gleasonGrade2 <- as.numeric(combined.data$gleasonGrade)
      if (gleason =='low') {
        combined.data <- combined.data[combined.data$gleasonGrade2<4,]
      } else if (gleason =='high') {
        combined.data <- combined.data[combined.data$gleasonGrade2>3,]     
      } else {
        combined.data <- combined.data
      }

      combined.data$tstage  <- as.numeric(combined.data$pathologicalStage)
      if (tstage =='two') {
        combined.data       <- combined.data[combined.data$tstage<6,]
      } else if (tstage =='threefour') {
        combined.data       <- combined.data[combined.data$tstage>5,]     
      } else {
              combined.data <- combined.data
      }

      if (groups =='primary') {
        combined.data <- combined.data[combined.data$groups == "Primary",]
      } else if (groups=='metastasis') {
        combined.data <- combined.data[combined.data$groups == "Metastasis",]            
      } else {
        combined.data <- combined.data
      }

      surv.xfs <- Surv((combined.data$BCR_FreeTime/12), combined.data$BCR_Event=="BCR_Algorithm")
      xlabel   <- "Time to Biochemical Recurrence (years)"
      ylabel   <- "Probability of Freedom from Biochemical Recurrence"
      title    <- "Biochemical Recurrence"
      
      #Carry out RP  #
      ctree_xfs   <- ctree(surv.xfs~geneexp, data=combined.data)
      pvalue      <- 1 - ctree_xfs@tree$criterion$maxcriterion
      newPval     <- signif(pvalue, digits = 2)
      pvalList[j] <- newPval
      ps3         <- NA
      highIsGood  <- NA

      if(newPval<0.05) {
        rp_pngfile <- paste(pngprefix, "rp.png", sep="_")
        png(filename=paste(outfile_path, rp_pngfile, sep="/"))
        plot(ctree(surv.xfs~geneexp, data=combined.data), xlab=xlabel, ylab=ylabel, main=paste(title, geneid))
        dev.off()

        ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
        ps  <- signif(ps2[1], digits = 3)

        if(length(ps2)==1) {
          combined.data$geneexp_cp <- combined.data$geneexp<=ps2[1]
          nt                       <- table(combined.data$geneexp_cp)
          geneexp.survfit.xfs      <- survfit(surv.xfs~combined.data$geneexp_cp)
          newPval2                 <- NA

          km_pngfile <- paste(pngprefix, "km.png", sep="_")
          png(filename=paste(outfile_path, km_pngfile, sep="/"))
          plot(geneexp.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(geneid,", p=", newPval), col=c(2,4))
          legend(0.1, 0.1, c(paste(names[j], ">", ps, "n=", nt[[1]]), paste(names[j], "<=", ps, "n=", nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty="n")
          dev.off()

          cat("<img src=\"", file=outfile_html, append=TRUE)
          cat(km_pngfile, file=outfile_html, append=TRUE)
          cat("\">\n", file=outfile_html, append=TRUE)
        }

        if(length(ps2)==2) {
          if(ps2[1]>ps2[2]) {
            ps3                       <- round(ps2[2], digits=3)
            combined.data$geneexp_cp1 <- combined.data$geneexp<=ps2[1]
            combined.data$geneexp_cp2 <- combined.data$geneexp<=ps2[2]
            combined.data$geneexp_cp3 <- combined.data$geneexp_cp1+combined.data$geneexp_cp2
            nt                        <- table(combined.data$geneexp_cp3)
            geneexp.survfit.xfs       <- survfit(surv.xfs~combined.data$geneexp_cp3)
            pvalue2                   <- 1 - ctree_xfs@tree$left[[3]][[2]]
            newPval2                  <- signif(pvalue2, digits=2)

            if(newPval2>=0.05) {
              km_pngfile <- paste(pngprefix, "km.png", sep="_")
              png(filename=paste(outfile_path, km_pngfile, sep="/"))
              plot(geneexp.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(geneid,", p=", newPval), col=c(2, 4, 5))
              legend(0.1, 0.15, c(paste(names[j], ">", ps,  "n=", nt[[1]]), paste(names[j], "<=", ps, "& >", ps3,  "n=", nt[[2]]), paste(names[j], "<=", ps3,  "n=", nt[[3]])), col=c(2, 4, 5), lty=1, lwd=1.5, bty="n")
              dev.off()
            } else {
              km_pngfile <- paste(pngprefix, "km.png", sep="_")
              png(filename=paste(outfile_path, km_pngfile, sep="/"))
              plot(geneexp.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(geneid,", p=", newPval, "& p=", newPval2), col=c(2, 4, 5))
              legend(0.1, 0.15, c(paste(names[j], ">", ps,  "n=", nt[[1]]), paste(names[j], "<=", ps, "& >", ps3,  "n=", nt[[2]]), paste(names[j], "<=", ps3,  "n=", nt[[3]])), col=c(2, 4, 5), lty=1, lwd=1.5, bty="n")
              dev.off()
            }

            cat("<img src=\"", file=outfile_html, append=TRUE)
            cat(km_pngfile, file=outfile_html, append=TRUE)
            cat("\">\n", file=outfile_html, append=TRUE)
          } else {
            ps3                       <- round(ps2[2], digits=3)
            combined.data$geneexp_cp1 <- combined.data$geneexp<=ps2[1]
            combined.data$geneexp_cp2 <- combined.data$geneexp<=ps2[2]
            combined.data$geneexp_cp3 <- combined.data$geneexp_cp1+combined.data$geneexp_cp2
            nt                        <- table(combined.data$geneexp_cp3)
            geneexp.survfit.xfs       <- survfit(surv.xfs~combined.data$geneexp_cp3)
            pvalue2                   <- 1 - ctree_xfs@tree$right[[3]][[2]]
            newPval2                  <- signif(pvalue2, digits=2)

            if(newPval2>=0.05) {
              km_pngfile <- paste(pngprefix, "km.png", sep="_")
              png(filename=paste(outfile_path, km_pngfile, sep="/"))
              plot(geneexp.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(geneid,", p=", newPval), col=c(2,4,5))
              legend(0.1, 0.15, c(paste(names[j], ">", ps3,  "n=", nt[[1]]), paste(names[j], ">", ps, "& <=", ps3,  "n=", nt[[2]]), paste(names[j], "<=", ps,  "n=", nt[[3]])), col=c(2,4,5), lty=1, lwd=1.5, bty="n")
              dev.off()
            } else {
              km_pngfile <- paste(pngprefix, "km.png", sep="_")
              png(filename=paste(outfile_path, km_pngfile, sep="/"))
              plot(geneexp.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(geneid,", p=", newPval, "& p=", newPval2), col=c(2,4,5))
              legend(0.1, 0.15, c(paste(names[j], ">", ps3,  "n=", nt[[1]]), paste(names[j], ">", ps, "& <=", ps3,  "n=", nt[[2]]), paste(names[j], "<=", ps,  "n=", nt[[3]])), col=c(2,4,5), lty=1, lwd=1.5, bty="n")
              dev.off()
            }

            cat("<img src=\"", file=outfile_html, append=TRUE)
            cat(km_pngfile, file=outfile_html, append=TRUE)
            cat("\">\n", file=outfile_html, append=TRUE)
          }
        }

        combined.data$geneexp_cp <- combined.data$geneexp <= ps2[1]
        coxphRegressionModel     <- coxph(surv.xfs~combined.data$geneexp_cp)
        hazardRatio              <- exp(coxphRegressionModel$coefficients)
        highIsGood               <- hazardRatio > 1

        cat("<img src=\"", file=outfile_html, append=TRUE)
        cat(rp_pngfile, file=outfile_html, append=TRUE)
        cat("\">\n", file=outfile_html, append=TRUE)
        cat("<p>\n", file=outfile_html, append=TRUE)
      }

      if(newPval>=0.05) {
        ps       <- NA
        newPval2 <- NA
      }

      pvalList2[j]      <- newPval2
      splitList[j]      <- ps
      splitList2[j]     <- ps3
      highIsGoodList[j] <- highIsGood
    }
  }

  missinggene     <- setdiff(geneList, names)
  missinggenelist <- NULL

  if (length (missinggene)>0) {
    for (k in 1:length(missinggene)) {
      missinggenetemp <- c(missinggene[k], rep(NA, 6))
      missinggenetemp <- as.matrix(t(missinggenetemp))
      missinggenelist <- rbind(missinggenetemp, missinggenelist)
    }
  }

  Gene       <- names
  Accession  <- accessid
  PValue1    <- pvalList
  PValue2    <- pvalList2
  CutOff1    <- splitList
  CutOff2    <- splitList2
  HighIsGood <- highIsGoodList

  xfsList <- cbind(Gene, Accession, CutOff1, CutOff2, PValue1, PValue2, HighIsGood)
  if (length (missinggene)>0) {
    xfsList <- rbind(xfsList, missinggenelist)
  }

  write.table(xfsList, outfile_txt, quote=FALSE, row.names=FALSE, sep="\t")

} else {
  cat("<h3>No matches to probes found for query genes</h3><p>\n", file=outfile_html, append=TRUE)
}

cat("</body>\n", file=outfile_html, append=TRUE)
cat("</html>\n", file=outfile_html, append=TRUE)

