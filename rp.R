# --------------------------------------------------------------------------------
# --- Command-line arguments (for use in Galaxy) 
# --------------------------------------------------------------------------------

ca <- commandArgs(trailingOnly=TRUE)

# Common arguments for all datasets
dataname = ca[1]
datafile = ca[2]
infile = ca[3]
col = as.numeric(ca[4])
if (is.na(col)) col = 1
outfile_txt = ca[5]
outfile_html = ca[6]
outfile_path = ca[7]

# All other arguments - if not in used set to "0"
outcome = ca[8]
receptor = ca[9]
treatment = ca[10]
er = ca[11]
age = ca[12]
grade = ca[13]
signature = ca[14]
margins = ca[15]
groups = ca[16]
gleason = ca[17]
tstage = ca[18]

# --------------------------------------------------------------------------------
# --- Methods
# --------------------------------------------------------------------------------

plot_km <- function(outfile_path, pngprefix, data, xlabel, ylabel, main, legend, col) {
    km_pngfile = paste(pngprefix, "km.png", sep="_")
    png(filename=paste(outfile_path, km_pngfile, sep="/"))
    plot(data, xlab=xlabel, ylab=ylabel, main=main, col=col)
    legend(0.1, 0.15, legend, col=col, lty=1, lwd=1.5, bty="n")
    dev.off()
}

plot_km_with_ps2_1 <- function(outfile_path, pngprefix, symbolProbe, symbol, combined.data, surv.xfs, xlabel, ylabel) {
    ctree_xfs <- ctree(surv.xfs~symbolProbe, data=combined.data)
        
    ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
    ps2_1 <- signif(ps2[1], digits=3)
    
    combined.data$symbolProbe_cp <- combined.data$symbolProbe<=ps2[1]
    nt <- table(combined.data$symbolProbe_cp)
    symbolProbe.survfit.xfs <- survfit(surv.xfs~combined.data$symbolProbe_cp)

    main = paste(symbolProbe,", p=", newPval)
    legend = c(paste(symbol, ">", ps2_1, "n=", nt[[1]]), paste(symbol, "<=", ps2_1 , "n=", nt[[2]]))
    km_pngfile <- plot_km(outfile_path, pngprefix, symbolProbe.survfit.xfs, xlabel, ylabel, main, legend, c(2,4))

    cat("<img src=\"", file=outfile_html, append=TRUE)
    cat(km_pngfile, file=outfile_html, append=TRUE)
    cat("\">\n", file=outfile_html, append=TRUE)
}

plot_km_with_ps2_2 <- function(outfile_path, pngprefix, symbolProbe, symbol, combined.data, surv.xfs, xlabel, ylabel) {
    ctree_xfs <- ctree(surv.xfs~symbolProbe, data=combined.data)
    
    pvalue <- 1 - ctree_xfs@tree$criterion$maxcriterion
    newPval <- signif(pvalue, digits = 2)
    
    pvalue2 <- 1 - ctree_xfs@tree$left[[3]][[2]]
    newPval2 <- signif(pvalue2, digits=2)
    
    ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
    ps2_1 <- signif(ps2[1], digits=3)
    ps2_2 <- round(ps2[2], digits=3)
    
    combined.data$symbolProbe_cp1 <- combined.data$symbolProbe<=ps2[1]
    combined.data$symbolProbe_cp2 <- combined.data$symbolProbe<=ps2[2]
    combined.data$symbolProbe_cp3 <- combined.data$symbolProbe_cp1+combined.data$symbolProbe_cp2
    nt <- table(combined.data$symbolProbe_cp3)
    symbolProbe.survfit.xfs <- survfit(surv.xfs~combined.data$symbolProbe_cp3)
    
    if (newPval2>=0.05) {
        main = paste(symbolProbe,", p=", newPval)
    } else {
        main = paste(symbolProbe,", p=", newPval, "& p=", newPval2)
    }
    if (ps2[1]>ps2[2]) {
        legend=c(paste(symbol, ">", ps2_1 ,  "n=", nt[[1]]), paste(symbol, "<=", ps2_1 , "& >", ps2_2,  "n=", nt[[2]]), paste(symbol, "<=", ps2_2,  "n=", nt[[3]]))
    } else {
        legend = c(paste(symbol, ">", ps2_2,  "n=", nt[[1]]), paste(symbol, ">", ps2_1 , "& <=", ps2_2,  "n=", nt[[2]]), paste(symbol, "<=", ps2_1 ,  "n=", nt[[3]]))
    }
    km_pngfile <- plot_km(outfile_path, pngprefix, symbolProbe.survfit.xfs, xlabel, ylabel, main, legend, c(2,4,5))
    
    cat("<img src=\"", file=outfile_html, append=TRUE)
    cat(km_pngfile, file=outfile_html, append=TRUE)
    cat("\">\n", file=outfile_html, append=TRUE)
}


# --------------------------------------------------------------------------------
# --- Loi Specific
# --------------------------------------------------------------------------------
load(datafile)

library(hgu133plus2.db)
library(party)
library(survival)

dir.create(outfile_path)

cat("<html>\n", file=outfile_html)
cat("<body>\n", file=outfile_html, append=TRUE)
cat("<title>\n", file=outfile_html, append=TRUE)
cat("Breast Cancer Survival Analysis\n", file=outfile_html, append=TRUE)
cat("</title>\n", file=outfile_html, append=TRUE)
cat("<h3>Breast Cancer Survival Analysis</h3>\n", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat("The following are the results of recursive partitioning for input genes", file=outfile_html, append=TRUE)
cat(" using a breast cancer gene expression dataset from Loi <em>et al.</em>,", file=outfile_html, append=TRUE)
cat(" <em>J. Clin. Oncol.</em> 25, 1239-1246 (2007) with accompanying", file=outfile_html, append=TRUE)
cat(" survival data.", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat(" This dataset comprises expression profiles for primary tumours from 414 patients", file=outfile_html, append=TRUE)
cat(" of which 277 were treated with Tamoxifen and 137 were untreated.", file=outfile_html, append=TRUE)
cat(" The data were generated with Affymetrix Human Genome U133A and U133B microarrays.", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

cat("Kaplan-Meier survival curves and recursive partitioning plots are", file=outfile_html, append=TRUE)
cat(" displayed below for those genes which are predictive of", file=outfile_html, append=TRUE)
if (outcome == 'rfs') {
    cat(" recurrence", file=outfile_html, append=TRUE)
} else {
    cat(" distant metastasis", file=outfile_html, append=TRUE)
}
cat(" free survival with p &lt; 0.05", file=outfile_html, append=TRUE)
cat(" (corrected for the testing of multiple cut-offs).", file=outfile_html, append=TRUE)
cat("<p>\n", file=outfile_html, append=TRUE)

geneList <- read.table(infile, sep="\t", as.is=TRUE)[,col]
geneList <- toupper(geneList)
geneList <- unique(geneList)
probes <- mget(geneList,hgu133plus2ALIAS2PROBE, ifnotfound=NA)
AffyID <- probes[[1]]

names <- rep(geneList[1], length(probes[[1]]))

if(length(geneList)>1) {
    for(j in 2:length(geneList)) {
        AffyID <- c(AffyID, probes[[j]])
        names <- c(names,rep(geneList[j],length(probes[[j]])))
    }
}

probetable <- cbind(names, AffyID)

uniqueProbeList <- na.omit(AffyID)
pvalList <- NA
pvalList2 <- NA
splitList <- NA
splitList2 <- NA
highIsGoodList <- NA
symbolList <- uniqueProbeList[1]
idlist <- NA

if (length(uniqueProbeList) > 0) {

    for (j in 1:length(uniqueProbeList)) {
        #get gene symbol and affy ID#
        genenum <- grep(uniqueProbeList[[j]], probetable[,2], ignore.case=TRUE)
        symbol2 <- probetable[genenum, 1]
        symbol <- symbol2[[1]]
        symbolList[j] <- symbol
        symbolProbe <- paste(symbol, "-", uniqueProbeList[[j]])

        pngprefix = paste(symbol, "_", uniqueProbeList[[j]], sep="")

        if (treatment == 'tamoxifen') {
            rowID <- match(uniqueProbeList[[j]], annot.tam[,1])
        } else {
            rowID <- match(uniqueProbeList[[j]], annot.untreated[,1])
        }

        idlist[j] <- rowID
        rowID <- na.omit(rowID)

        if(length(rowID)>0) {

            #Extract data for PgR etc#
            if (treatment == 'tamoxifen') {
                combined.data <- cbind(demo.tam, symbolProbe=as.numeric(data.tam[,rowID]))
            } else {
                combined.data <- cbind(demo.untreated, symbolProbe=as.numeric(data.untreated[,rowID]))
            }

            if (receptor == 'er') {
                combined.data <- combined.data[combined.data$er == 1,]
            } else if (receptor == 'pgr') {
                combined.data <- combined.data[combined.data$pgr == 1,]
            } else if (receptor == 'erpgr') {
                combined.data <- combined.data[combined.data$er == 1,]
                combined.data <- combined.data[combined.data$pgr == 1,]
            }

            if (outcome == 'rfs') {
                combined.data <- combined.data[!is.na(combined.data$e.rfs),]
                combined.data <- combined.data[!is.na(combined.data$t.rfs),]
            } else {
                combined.data <- combined.data[!is.na(combined.data$e.dmfs),]
                combined.data <- combined.data[!is.na(combined.data$t.dmfs),]
            }
            
            #Create survival object #
            if (outcome == 'rfs') {
                surv.xfs <- Surv((combined.data$t.rfs/365.25), combined.data$e.rfs=="1")
                xlabel = "Time to Recurrence (years)"
                ylabel = "Probability of Recurrence Free Survival"
                title = "Recurrence Free Survival"
            } else {
                surv.xfs <- Surv((combined.data$t.dmfs/365.25), combined.data$e.dmfs=="1")
                xlabel = "Time to Distant Metastasis (years)"
                ylabel = "Probability of Distant Metastasis Free Survival"
                title = "Distant Metastasis Free Survival"
            }

            #Carry out RP #
            ctree_xfs <- ctree(surv.xfs~symbolProbe, data=combined.data)
            
            pvalue <- 1 - ctree_xfs@tree$criterion$maxcriterion
            newPval <- signif(pvalue, digits = 2)

            newPval2 <- NA
            ps2_1 <- NA
            ps2_2 <- NA
            highIsGood <- NA

            if(newPval<0.05) {
                rp_pngfile = paste(pngprefix, "rp.png", sep="_")
                png(filename=paste(outfile_path, rp_pngfile, sep="/"))
                plot(ctree(surv.xfs~symbolProbe, data=combined.data), xlab=xlabel, ylab=ylabel, main=paste(title, symbolProbe))
                dev.off()

                ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
                ps2_1 <- signif(ps2[1], digits = 3)
                ps2_2 <- round(ps2[2], digits=3)

                if(length(ps2)==1) {
                    plot_km_with_ps2_1(outfile_path, pngprefix, symbolProbe, symbol, combined.data, surv.xfs, xlabel, ylabel)
                }

                if(length(ps2)==2) {
                    pvalue2 <- 1 - ctree_xfs@tree$left[[3]][[2]]
                    newPval2 <- signif(pvalue2, digits=2)
                    plot_km_with_ps2_2(outfile_path, pngprefix, symbolProbe, symbol, combined.data, surv.xfs, xlabel, ylabel)
                }

                combined.data$symbolProbe_cp <- combined.data$symbolProbe <= ps2[1]
                coxphRegressionModel <- coxph(surv.xfs~combined.data$symbolProbe_cp)
                hazardRatio <- exp(coxphRegressionModel$coefficients)
                highIsGood <- hazardRatio > 1

                cat("<img src=\"", file=outfile_html, append=TRUE)
                cat(rp_pngfile, file=outfile_html, append=TRUE)
                cat("\">\n", file=outfile_html, append=TRUE)
                cat("<p>\n", file=outfile_html, append=TRUE)
            }
            
            pvalList[j] <- newPval
            pvalList2[j] <- newPval2
            splitList[j] <- ps2_1
            splitList2[j] <- ps2_2
            highIsGoodList[j] <- highIsGood
        }
    }

    GeneName <- symbolList
    AffyID <- uniqueProbeList
    PValue1 <- pvalList
    PValue2 <- pvalList2
    CutOff1 <- splitList
    CutOff2 <- splitList2
    HighIsGood <- highIsGoodList

    xfsList <- cbind(GeneName, AffyID, CutOff1, CutOff2, PValue1, PValue2, HighIsGood)

    write.table(xfsList, outfile_txt, quote=FALSE, row.names=FALSE, sep="\t")

} else {
    cat("No matches to Affymetrix probes found for query genes<p>\n", file=outfile_html, append=TRUE)
}

cat("</body>\n", file=outfile_html, append=TRUE)
cat("</html>\n", file=outfile_html, append=TRUE)


