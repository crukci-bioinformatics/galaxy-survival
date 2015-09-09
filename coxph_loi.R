# --- Command-line arguments (for use in Galaxy) -------------------------------

ca <- commandArgs(trailingOnly=TRUE)

datafile            <- ca[1]
infile              <- ca[2]
col                 <- as.numeric(ca[3])
if (is.na(col)) col <- 1
outfile             <- ca[4]
outcome             <- ca[5]
receptor            <- ca[6]
treatment           <- ca[7]

# --- Hard-wired parameters ----------------------------------------------------

#datafile  <- 'galaxy/tool-data/survival_analysis/LUMINAL.1.RData'
#datafile  <- 'LUMINAL.1.RData'
#infile    <- 'genelist.txt'
#col       <- 2
#outfile   <- 'results.txt'
#outcome   <- 'rfs'
#outcome   <- 'dmfs'
#receptor  <- 'er'
#receptor  <- 'pgr'
#receptor  <- 'all'
#treatment <- 'tamoxifen'
#treatment <- 'none'

# ------------------------------------------------------------------------------

load(datafile)

library(hgu133plus2.db)
library(survival)

geneList  <- read.table(infile, sep="\t", as.is=TRUE)[, col]
geneList  <- toupper(geneList)
geneList  <- unique(geneList)
probes    <- mget(geneList, hgu133plus2ALIAS2PROBE, ifnotfound=NA)
probelist <- probes[[1]]
names     <- rep(geneList[1], length(probes[[1]]))

if (length(geneList) > 1) {
  for (j in 2:length(geneList)) {
    probelist <- c(probelist, probes[[j]])
    names     <- c(names, rep(geneList[j], length(probes[[j]])))
  }
}

probetable      <- cbind(names, probelist)

uniqueProbeList <- unique(na.omit(probelist))
pvalList        <- NA
hrList          <- NA
lblist          <- NA
ublist          <- NA
symbolList      <- uniqueProbeList[1]

if (length(uniqueProbeList) > 0) {
  for (j in 1:length(uniqueProbeList)) {
    # get gene symbol and affy ID
    genenum       <- grep(uniqueProbeList[[j]], probetable[, 2], ignore.case=TRUE)
    symbol2       <- probetable[genenum, 1]
    symbol        <- symbol2[[1]]
    symbolList[j] <- symbol
    symbolProbe   <- paste(symbol, "-", uniqueProbeList[[j]])

    if (treatment == 'tamoxifen') {
      rowID <- match(uniqueProbeList[[j]], annot.tam[, 1])
    } else {
      rowID <- match(uniqueProbeList[[j]], annot.untreated[, 1])
    }
    rowID   <- na.omit(rowID)

    if(length(rowID)>0) {
      # Extract data for RFS/DMFS, ER/PgR, etc
      if (treatment == 'tamoxifen') {
        combined.data <- cbind(demo.tam, symbolProbe=as.numeric(data.tam[, rowID]))
      } else {
        combined.data <- cbind(demo.untreated, symbolProbe=as.numeric(data.untreated[, rowID]))
      }
      if (receptor == 'er') {
        combined.data <- combined.data[combined.data$er == 1, ]
      } else if (receptor == 'pgr') {
        combined.data <- combined.data[combined.data$pgr == 1, ]
      }
      if (outcome == 'rfs') {
        combined.data <- combined.data[!is.na(combined.data$e.rfs), ]
        combined.data <- combined.data[!is.na(combined.data$t.rfs), ]
      } else {
        combined.data <- combined.data[!is.na(combined.data$e.dmfs), ]
        combined.data <- combined.data[!is.na(combined.data$t.dmfs), ]
      }

      # Create survival object.#
      if (outcome == 'rfs') {
        surv.xfs  <- Surv((combined.data$t.rfs/365.25), combined.data$e.rfs == "1")
      } else {
        surv.xfs  <- Surv((combined.data$t.dmfs/365.25), combined.data$e.dmfs == "1")
      }
      x           <- coxph(surv.xfs~combined.data$symbolProbe)
      se          <- sqrt(diag(surv <-x$var))
      coef        <- x$coefficients
      pval        <- signif(1-pchisq((coef/se)^2, 1), 2)
      hr          <- round((exp(coef)), digits=2)
      lb          <- round(exp(coef-(1.96*se)), digits=2)
      ub          <- round(exp(coef+(1.96*se)), digits=2)
      lblist[j]   <- lb
      ublist[j]   <- ub
      pvalList[j] <- pval
      hrList[j]   <- hr
    }
  }

  GeneName <- symbolList
  AffyID   <- uniqueProbeList
  PValue   <- pvalList
  HR       <- as.character(hrList)
  LB       <- as.character(lblist)
  UB       <- as.character(ublist)
  results  <- cbind(GeneName, AffyID, HR, LB, UB, PValue)

  write.table(results, outfile, quote=FALSE, row.names=FALSE, sep="\t")

} else {
  print("<h3>No matches to Affymetrix probes found for query genes</h3>\n")
}

