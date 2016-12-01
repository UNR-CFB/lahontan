#!/usr/bin/python3

'''Usage: makeReportr.py [-h | --help] [-j <jsonfile>] [-t <tofile>]
                                                                                      
Options:                                                                              
    -h --help           Show this screen                                              
    -j <jsonfile>           Optional name of JSON file to be read [default: Metadata.json]
    -t <tofile>         Optional name of R report script [default: makeReport.r]
'''                                                                                   

from docopt import docopt
import csv
import subprocess
import json
import makeCols
import os

def createRscript(jsontoRead,name):
    ''' Arguments:
            jsontoRead = str; name of JSON Metadata file
                [default: Metadata.json]
            name = str; name for edgeR script
        Returns:
            None

        Creates DESeq2 analysis script by filling in a template
    '''
    makeReportsTemplate = """'Usage: makeReports.r [-h | --help] [-p </path/to/data>] [--counts <name>] [--cols <name>]

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
countData <- read.csv(opts$counts,sep='\\t',header=TRUE)
colData <- read.csv(opts$cols,sep='\\t',header=TRUE)

# Make DESeq Data Frame
dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=colData,
                              design = ~ {1_mf})

# Filter out rows that contain only 1 or 0
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Label variables in condition column
dds${2_mf} <- factor(dds${3_mf}, levels=c({4_factorlist}))

# Run DESeq2
dds <- DESeq(dds)

# Make simple PCA Plot
rld <- rlog(dds)
pdf('Simple_PCA.pdf')
plotPCA(rld,intgroup='{5_mf}')
dev.off()

# Make Detailed regionReport
dir.create('SimpleReport',showWarnings = FALSE,recursive = TRUE)
report <- DESeq2Report(dds, project = 'DESeq2 Report',
                        intgroup = c('{5_mf}'),
                        outdir = 'SimpleReport',
                        output = 'DESeq2Report',
                        theme = theme_bw(),
                        output_format='pdf_document',
                        device='pdf')


# Adjust DESeq Data Frame to include other features
ddsAdj <- dds
design(ddsAdj) <- formula(~ {6_formulalistplus})
ddsAdj <- DESeq(ddsAdj)

# Make PCA Plot with Run Numbers
rldAdj <- rlog(ddsAdj)
pdf('Detailed_PCA.pdf')
plotPCA(rldAdj,intgroup=c({7_formulalistcomma}))
dev.off()

# Make ReportingTools Report
RTReport <- HTMLReport(shortName='RNASeq_Analysis_with_DESeq2',
                     title='RNASeq analysis of differential expression with DESeq2',
                     reportDirectory="./ReportingToolsReport")
publish(ddsAdj,RTReport,pvalueCutoff={8_pvaluecutoff},factor=ddsAdj${9_mf},
        reportDir="./ReportingToolsReport")
finish(RTReport)

# Make Detailed regionReport
dir.create('DetailedReport',showWarnings = FALSE,recursive = TRUE)
report <- DESeq2Report(ddsAdj, project = 'DESeq2 Report',
                        intgroup = c({10_formulalistcomma}),
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
#fileConn <- file("FINISHED.txt")
#writeLines("We are finished", fileConn)
#close(fileConn)
"""
    factorlist,formulalistplus,formulalistcomma,mainFeature = getContext(jsontoRead)
    Context = {
            "1_mf": mainFeature,
            "2_mf": mainFeature,
            "3_mf": mainFeature,
            "4_factorlist": factorlist,
            "5_mf": mainFeature,
            "6_formulalistplus": formulalistplus,
            "7_formulalistcomma": formulalistcomma,
            "8_pvaluecutoff": 0.05,
            "9_mf": mainFeature,
            "10_formulalistcomma": formulalistcomma,
            }
    with open(name,'w') as f:
        f.write(makeReportsTemplate.format(**Context))

def getContext(jsontoRead):
    ''' Arguments:
            jsontoRead = str; name of JSON Metadata file
                            [default: Metadata.json]
        Returns:
            factorlist,formulalistplus,formulalistcomma,mainFeature
            = tuple containing:
                factorlist = str; the list of features to be analyzed
                                by DESeq2
                formulalistplus = str; formula to be analyzed
                formulalistcomma = str; list of factors to be analyzed
                mainFeature = str; name of mainFeature

        Scrapes Metadata for context to be filled into DESeq2
        template
    '''
    Metadata = jsontoRead
    mainFeature = Metadata['MainFeature']

    #TODO Generalize 'Cols.dat' in case it's not named Cols.dat
    with open('Cols.dat','r') as f:
        colData = csv.reader(f,dialect='unix',delimiter='\t')
        x = []
        for row in colData:
            x.append(row)

    formulalist = x[0][:]
    formulalistp = list(reversed(formulalist))
    x[0].insert(0,'Samples')

    z = list(zip(*x))
    for group in z:
        if str(group[0]) == str(mainFeature):
            realgroup = list(group)[1:]
            uniquesetj = set(realgroup)
            uniquelist = list(uniquesetj)

    factorlist = ','.join('"{0}"'.format(a) for a in uniquelist) 
    formulalistcomma = ','.join('"{0}"'.format(b) for b in formulalist)
    formulalistplus = ' + '.join(formulalistp)
    return factorlist,formulalistplus,formulalistcomma,mainFeature

if __name__ == '__main__':
    arguments = docopt(__doc__,version='1.0')
    createRscript(makeCols.readJSON(arguments['-j']),arguments['-t'])
