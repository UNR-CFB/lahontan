#!/usr/bin/python3

'''Usage: makeBallgownScript.py [-h | --help] [-j <jsonfile>] [-t <tofile>]

Options:
    -h --help           Show this screen
    -j <jsonfile>       Optional name of JSON file to be read [default: Metadata.json]
    -t <tofile>         Optional name of R report script [default: runBallgown.r]
'''

################################################################
# Importations
################################################################

from docopt import docopt
import csv
import subprocess
import json
import makeCols
import os

################################################################

def createRBallgownScript(jsontoRead,name):
    ''' Arguments:
            jsontoRead = str; name of JSON Metadata file
                [default: Metadata.json]
            name = str; name for edgeR script
        Returns:
            None

        Creates Ballgown analysis script by filling in a template
    '''
    scriptTemplate = """'Usage: runBallgown.r [-h | --help] [-p </path/to/data>] [--counts <name>] [--cols <name>]

Example:
Rscript runBallgown.r -p /path/to/project --cols ColumnsName

Options:
    -h --help           Show this
    -p <path>           Directory that contains Cols.dat and Counts.dat [default: ./]
    --cols <name>       Name of phenotype file that describes experiment [default: Cols.dat]'-> doc

# Set up command line argument library
library('docopt')

# Retrieve command line arguments
opts <- docopt(doc)

# Change to directory
setwd(opts$p)

# Load the necessary libraries
library('ballgown')
library('RSkittleBrewer')
library('genefilter')
library('dplyr')
library('devtools')

# Read phenotype sample data
phenotypeData <- read.csv(opts$cols)

# Read in expression data
stringtieData <- ballgown(dataDir = "StringtieResults", samplePattern="sample", pData=phenotypeData)

# Filter low abundance genes
stringtieDataFiltered <- subset(stringtieData, "rowVars(texpr(stringtieData)) > 1", genomesubset=TRUE)

# DE by transcript
results_transcripts <- stattest(stringtieDataFiltered, feature='transcript', covariate={mainfeature},
         {secondaryfeatures}getFC=TRUE, meas='FPKM')

# DE by gene
results_genes <-  stattest(stringtieDataFiltered, feature='gene', covariate={mainfeature},
         {secondaryfeatures}getFC=TRUE, meas='FPKM')

# Add gene name
results_transcripts <- data.frame(geneNames=ballgown::geneNames(stringtieDataFiltered),
          geneIDs=ballgown::geneIDs(stringtieDataFiltered), results_transcripts)

# Sort results from smallest p-value
results_transcripts <- arrange(results_transcripts, pval)
results_genes <-  arrange(results_genes, pval)

# Write results to CSV
write.csv(results_transcripts, {transcriptresultsfile}, row.names=FALSE)
write.csv(results_genes, {generesultsfile}, row.names=FALSE)

# Filter for genes with q-val <0.05
subset(results_transcripts, results_transcripts$qval <=0.05)
subset(results_genes, results_genes$qval <=0.05)

# Plotting setup
#tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
#palette(tropical)

# Plotting gene abundance distribution
#fpkm <- texpr(stringtieData, meas='FPKM')
#fpkm <- log2(fpkm +1)
#boxplot(fpkm, col=as.numeric(phenotypeData$sex), las=2,ylab='log2(FPKM+1)')

# Plot individual transcripts
#ballgown::transcriptNames(stringtieData)[12]
#plot(fpkm[12,] ~ phenotypeData$sex, border=c(1,2),
#     main=paste(ballgown::geneNames(stringtieData)[12], ' : ',ballgown::transcriptNames(stringtieData)[12]),
#     pch=19, xlab="Sex", ylab='log2(FPKM+1)')
#points(fpkm[12,] ~ jitter(as.numeric(phenotypeData$sex)), col=as.numeric(phenotypeData$sex))

# Plot gene of transcript 1729
#plotTranscripts(ballgown::geneIDs(stringtieData)[1729], stringtieData,
#                main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

# Plot average expression
#plotMeans(ballgown::geneIDs(stringtieData)[203], stringtieDataFiltered, groupvar="sex", legend=FALSE)

# unload the libraries
detach('package:ballgown')
detach('package:RSkittleBrewer')
detach('package:genefilter')
detach('package:dplyr')
detach('package:devtools')

# Make file to notify when done
#fileConn <- file("FINISHED.txt")
#writeLines("We are finished", fileConn)
#close(fileConn)
"""
    mainFeature, secondaryFeatures = getContext(jsontoRead)
    Context = {
            "mainfeature": '"{}"'.format(mainFeature),
            "secondaryfeatures": secondaryFeatures,
            "transcriptresultsfile": '"ballgownTranscriptsResults.csv"',
            "generesultsfile":  '"ballgownGeneResults.csv"'
            }
    with open(name,'w') as f:
        f.write(scriptTemplate.format(**Context))

def getContext(jsontoRead):
    ''' Arguments:
            jsontoRead = str; name of JSON Metadata file
                            [default: Metadata.json]
        Returns:
            tuple containing:
                mainFeature = str; name of mainFeature

        Scrapes Metadata for context to be filled into ballgown
        template
    '''
    Metadata = jsontoRead
    mainFeature = Metadata['MainFeature']
    featureNames = Metadata['FeatureNames']

    featureNames.remove(mainFeature)
    if len(featureNames) > 0:
        secondaryFeatures = 'adjustvars=c("' + '","'.join(featureNames) + '"), '
    else:
        secondaryFeatures = ''

    return mainFeature, secondaryFeatures

################################################################
if __name__ == '__main__':
    arguments = docopt(__doc__,version='1.0')
    createRBallgownScript(makeCols.readJSON(arguments['-j']),arguments['-t'])
