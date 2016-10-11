'Usage: makeReports.r [-h | --help] [-p </path/to/data>] [--counts <name>] [--cols <name>]

Example:
Rscript makeReports.r -p /path/to/project --counts CountName --cols ColumnsName

Options:
    -h --help           Show this
    -p <path>           Directory that contains Cols.dat and Counts.dat [default: ./]
    --counts <name>     Name of collected featureCounts data file [default: Counts.dat]
    --cols <name>       Name of Column file that describes Counts [default: Cols.dat]'-> doc

# Set up command line argument library
library('docopt')

# Retrieve Command Line Arguments
opts <- docopt(doc)

# change to the directory
setwd(opts$p)

# load the necessary libraries
library('DESeq2')
library('ReportingTools')
library('ggplot2')
library('regionReport')

# Load the datasets
countData <- read.csv(opts$counts,sep='\t',header=TRUE)
colData <- read.csv(opts$cols,sep='\t',header=TRUE)

# Make DESeq Data Frame
dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=colData,
                              design = ~ condition)

# Filter out rows that contain only 1 or 0
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Label variables in condition column
dds$condition <- factor(dds$condition, levels=c("Wtctrl","cngcHS","WtHS","cngcctrl"))

# Run DESeq2
dds <- DESeq(dds)

# Make simple PCA Plot
rld <- rlog(dds)
pdf('Simple_PCA.pdf')
plotPCA(rld,intgroup='condition')
dev.off()

# Make Detailed regionReport
dir.create('SimpleReport',showWarnings = FALSE,recursive = TRUE)
report <- DESeq2Report(dds, project = 'DESeq2 Report',
                        intgroup = c('condition'),
                        outdir = 'SimpleReport',
                        output = 'DESeq2Report',
                        theme = theme_bw(),
                        output_format='pdf_document',
                        device='pdf')


# Adjust DESeq Data Frame to include other features
ddsAdj <- dds
design(ddsAdj) <- formula(~ run + condition)
ddsAdj <- DESeq(ddsAdj)

# Make PCA Plot with Run Numbers
rldAdj <- rlog(ddsAdj)
pdf('Detailed_PCA.pdf')
plotPCA(rldAdj,intgroup=c("condition","run"))
dev.off()

# Make ReportingTools Report
RTReport <- HTMLReport(shortName='RNASeq_Analysis_with_DESeq2',
                     title='RNASeq analysis of differential expression with DESeq2',
                     reportDirectory="./ReportingToolsReport")
publish(ddsAdj,RTReport,pvalueCutoff=0.05,factor=ddsAdj$condition,
        reportDir="./ReportingToolsReport")
finish(RTReport)

# Make Detailed regionReport
dir.create('DetailedReport',showWarnings = FALSE,recursive = TRUE)
report <- DESeq2Report(ddsAdj, project = 'DESeq2 Report',
                        intgroup = c("condition","run"),
                        outdir = 'DetailedReport',
                        output = 'DESeq2Report',
                        theme = theme_bw(),
                        output_format='pdf_document',
                        device='pdf')

# unload the libraries
detach("package:DESeq2")
detach("package:ReportingTools")
detach("package:ggplot2")
detach("package:regionReport")
detach("package:docopt")

# Make file to notify when done
fileConn <- file("FINISHED.txt")
writeLines("We are finished", fileConn)
close(fileConn)
