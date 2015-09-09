# --- Command-line arguments (for use in Galaxy) -------------------------------

ca <- commandArgs(trailingOnly=TRUE)

datafile = ca[1]
infile = ca[2]
col = as.numeric(ca[3])
if (is.na(col)) col = 1
outfile = ca[4]
outcome = ca[5]
signature = ca[6]

# --- Hard-wired parameters ----------------------------------------------------

#datafile = 'galaxy/tool-data/survival_analysis/Shedden3.RData'
#infile = 'genelist.txt'
#col = 1
#outfile = 'results.txt'
#outcome = 'rfs'
#outcome = 'ofs'
#signature = 'kras.pos'
#signature = 'kras.neg'
#signature = 'kras.indef'
#signature = 'all'

# ------------------------------------------------------------------------------

load(datafile)

library(hgu133a.db)
library(survival)

geneList <- read.table(infile, sep="\t", as.is=TRUE)[,col]
geneList <- toupper(geneList)
geneList <- unique(geneList)

probes <- mget(geneList, hgu133aALIAS2PROBE, ifnotfound=NA)

probelist <- probes[[1]]
names <- rep(geneList[1], length(probes[[1]]))

if (length(geneList) > 1)
{
  for (j in 2:length(geneList))
  {
    probelist <- c(probelist, probes[[j]])
    names <- c(names, rep(geneList[j], length(probes[[j]])))
  }
}

probetable <- cbind(names, probelist)

uniqueProbeList <- unique(na.omit(probelist))
pvalList <- NA
hrList <- NA
lblist <- NA
ublist <- NA
symbolList <- uniqueProbeList[1]

if (length(uniqueProbeList) > 0)
{

  for (j in 1:length(uniqueProbeList))
  {
    # get gene symbol and affy ID
    genenum <- grep(uniqueProbeList[[j]], probetable[,2], ignore.case=TRUE)
    symbol2 <- probetable[genenum, 1]
    symbol <- symbol2[[1]]
    symbolList[j] <- symbol
    symbolProbe <- paste(symbol, "-", uniqueProbeList[[j]])

    rowID <- match(uniqueProbeList[[j]], colnames(data.ordered))
    rowID <- na.omit(rowID)

    if(length(rowID)>0)
    {
      combined.data = cbind(demo.ordered, symbolProbe=as.numeric(data.ordered[,rowID]))

      # Extract data for RFS/DMFS, ER/PgR, etc
      if (signature == 'kras.pos')
      {
        combined.data <- combined.data[combined.data$kras == 1,]
      } else if (signature == 'kras.neg')
      {
        combined.data <- combined.data[combined.data$kras == 0,]
      } else if (signature == 'kras.indef')
      {
        combined.data <- combined.data[combined.data$kras == -1,]
      }
      #ASSUMED THAT ALL DOESN'T NEED A STATEMENT HERE

      # Create survival object.#
      if (outcome == 'rfs')
      {
        surv.xfs <- Surv((combined.data$time.rfs), combined.data$status.rfs == "Progressor")
      } else
      {
        surv.xfs <- Surv((combined.data$time.os), combined.data$status.os == "Dead")
      }
      x <- coxph(surv.xfs~combined.data$symbolProbe)
      se <- sqrt(diag(surv <-x$var))
      coef <- x$coefficients
      pval <- signif(1-pchisq((coef/se)^2, 1), 2)
      hr <- round((exp(coef)), digits=2)
      lb <- round(exp(coef-(1.96*se)), digits=2)
      ub <- round(exp(coef+(1.96*se)), digits=2)
      lblist[j] <- lb
      ublist[j] <- ub
      pvalList[j] <- pval
      hrList[j] <- hr
    }
  }

  GeneName <- symbolList
  AffyID <- uniqueProbeList
  PValue <- pvalList
  HR <- hrList
  LB <- lblist
  UB <- ublist

  results <- cbind(GeneName, AffyID, HR, LB, UB, PValue)

  write.table(results, outfile, quote=FALSE, row.names=FALSE, sep="\t")

} else
{
  print("No matches to Affymetrix probes found for query genes\n")
}

