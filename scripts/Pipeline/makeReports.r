# store the current directory
initial.dir<-getwd()

# change to the new directory
setwd(initial.dir)

# load the necessary libraries
library('DESeq2')
library('ReportingTools')
library('ggplot2')
library('regionReport')

# Load the datasets
countData <- read.csv('Counts.dat',sep='\t',header=TRUE)
colData <- read.csv('Cols.dat',sep='\t',header=TRUE)

# Make DESeq Data Frame
dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=colData,
                              design = ~ condition)

# Filter out rows that contain only 1 or 0
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Label variables in condition column
dds$condition <- factor(dds$condition, levels=c('cngctrl',
                                                'cngcHS',
                                                'Wtctrl',
                                                'WtHS'))

# Run DESeq2
dds <- DESeq(dds)

# Make simple PCA Plot
rld <- rlog(dds)
pdf('Simple_PCA.pdf')
plotPCA(rld,intgroup='condition')
dev.off()

# Adjust DESeq Data Frame to include run numbers
ddsAdj <- dds
design(ddsAdj) <- formula(~ run_number + condition)
ddsAdj <- DESeq(ddsAdj)

# Make PCA Plot with Run Numbers
rldAdj <- rlog(ddsAdj)
pdf('Detailed_PCA.pdf')
plotPCA(rldAdj,intgroup=c('condition','run_number'))
dev.off()

# Make ReportingTools Report
RTReport <- HTMLReport(shortName='RNASeq_Analysis_with_DESeq2',
                     title='RNASeq analysis of diff expression with DESeq2',
                     reportDirectory="./ReportingToolsReport")
publish(ddsAdj,RTReport,pvalueCutoff=0.1,factor=ddsAdj$condition,
        reportDir="./ReportingToolsReport")
finish(RTReport)

# Make regionReport
dir.create('regionReport',showWarnings = FALSE,recursive = TRUE)
report <- DESeq2Report(ddsAdj, project = 'DESeq2 Report',
                        intgroup = c('condition','run_number'),
                        outdir = 'regionReport',
                        output = 'DESeq2Report',
                        theme = theme_bw(),
                        output_format='pdf_document',
                        device='pdf')

# unload the libraries
detach("package:DESeq2")
detach("package:ReportingTools")
detach("package:ggplot2")
detach("package:regionReport")

# change back to the original directory
setwd(initial.dir)
