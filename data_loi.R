# --- Command-line arguments (for use in Galaxy) -------------------------------

ca <- commandArgs(trailingOnly=TRUE)

datafile = ca[1]
infile = ca[2]
col = as.numeric(ca[3])
if (is.na(col)) col = 1
outfile = ca[4]
outcome = ca[5]
receptor = ca[6]
treatment = ca[7]

# --- Hard-wired parameters ----------------------------------------------------

#datafile = '~/galaxy-cri/tool-data/cri/survival_analysis/LUMINAL.RData'
#infile = 'genelist.txt'
#col = 2
#outfile = 'loi_data.txt'
#outcome = 'rfs'
#receptor = 'all'
#treatment = 'tamoxifen'

# ------------------------------------------------------------------------------


load(datafile)

library(hgu133plus2.db)
library(party)

geneList <- read.table(infile, sep="\t", as.is=TRUE)[,col]
geneList <- toupper(geneList)
geneList <- unique(geneList)
probes <- mget(geneList,hgu133plus2ALIAS2PROBE, ifnotfound=NA)
AffyID <- probes[[1]]

 if (treatment == 'tamoxifen')
      {
      	   combined.data<-demo.tam

      } else
      {
        combined.data <- demo.untreated
      }


names <- rep(geneList[1], length(probes[[1]]))

if(length(geneList)>1)
{
  for(j in 2:length(geneList))
  {
    AffyID <- c(AffyID, probes[[j]])
    names <- c(names,rep(geneList[j],length(probes[[j]])))
  }
}

probetable <- cbind(names, AffyID)

uniqueProbeList <- na.omit(AffyID)
symbolList <- uniqueProbeList[1]
idlist <- NA

if (length(uniqueProbeList) > 0)
{

  for (j in 1:length(uniqueProbeList))
  {
    #get gene symbol and affy ID#
    genenum <- grep(uniqueProbeList[[j]], probetable[,2], ignore.case=TRUE)
    symbol2 <- probetable[genenum, 1]
    symbol <- symbol2[[1]]
    symbolList[j] <- symbol
    symbolProbe <- paste(symbol, "-", uniqueProbeList[[j]])

    if (treatment == 'tamoxifen')
    {
      rowID <- match(uniqueProbeList[[j]], annot.tam[,1])
    } else
    {
      rowID <- match(uniqueProbeList[[j]], annot.untreated[,1])
    }

    idlist[j] <- rowID
    rowID <- na.omit(rowID)

    if(length(rowID)>0)
    {

      #Extract data for PgR etc#
      if (treatment == 'tamoxifen')
      {
      	combined.data[,ncol(combined.data)+1]=as.numeric(data.tam[,rowID])
        colnames(combined.data)[ncol(combined.data)]<-paste(symbolProbe)
       }
       else
       {
       	combined.data[,ncol(combined.data)+1]=as.numeric(data.untreated[,rowID])
        colnames(combined.data)[ncol(combined.data)]<-paste(symbolProbe)
       }
    }
  }
}
  if (receptor == 'er')
      {
        combined.data <- combined.data[combined.data$er == 1,]
      } else if (receptor == 'pgr')
      {
        combined.data <- combined.data[combined.data$pgr == 1,]
      } else if (receptor == 'erpgr')
      {
        combined.data <- combined.data[combined.data$er == 1,]
        combined.data <- combined.data[combined.data$pgr == 1,]
      }

      if (outcome == 'rfs')
      {
        combined.data <- combined.data[!is.na(combined.data$e.rfs),]
        combined.data <- combined.data[!is.na(combined.data$t.rfs),]
      } else if (outcome == 'dmfs')
      {
        combined.data <- combined.data[!is.na(combined.data$e.dmfs),]
        combined.data <- combined.data[!is.na(combined.data$t.dmfs),]
      } else
      {
        combined.data <- combined.data
      } 
write.table(combined.data, outfile, quote=FALSE, row.names=FALSE, sep="\t")
