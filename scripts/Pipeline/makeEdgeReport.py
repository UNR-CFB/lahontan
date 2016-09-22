#!/usr/bin/python3

'''Usage: makeEdgeReport.py [-h | --help] [-j <jsonfile>] [-t <tofile>]
                                                                                      
Options:                                                                              
    -h --help           Show this screen                                              
    -j <jsonfile>           Optional name of JSON file to be read [default: Metadata.json]
    -t <tofile>         Optional name of R report script [default: makeEdge.r]
'''                                                                                   

from docopt import docopt
from itertools import combinations,permutations
import csv
import subprocess
import json
import makeCols
import os

def createEdgeR(jsontoRead,name):
    makeReportsTemplate = """'Usage: makeEdge.r [-h | --help] [-p </path/to/data>] [--counts <name>] [--cols <name>]

Options:
    -h --help           Show this
    -p <path>           Directory that contains Cols.dat and Counts.dat [default: ./]
    --counts <name>     Name of collected featureCounts data file [default: Counts.dat]
    --cols <name>       Name of Column file that describes Counts [default: Cols.dat]

Example:
Rscript makeEdge.r -p /path/to/project --counts CountName --cols ColumnsName

(Or if in directory with "Counts.dat" and "Cols.dat", can just call with:
    Rscript makeEdge.r
)'-> doc

# Set up command line argument library
library('docopt')

# Retrieve Command Line Arguments
opts <- docopt(doc)

# change to the directory
setwd(opts$p)

# load the necessary libraries
library('edgeR')
library('ggplot2')
library('regionReport')

# Load the datasets
countData <- read.csv(opts$counts,sep='\\t',header=TRUE)
colData <- read.csv(opts$cols,sep='\\t',header=TRUE)

# Make initial data frame
y <- DGEList(counts=countData)

# Format colData for edgeR
Group <- factor(paste(colData${mainfeature}, sep="."))
cbind(colData, Group=Group)
design <- model.matrix(~ 0+Group, data=y$samples)
colnames(design) <- levels(Group)

# Make edgeR Data Frame
y <- DGEList(counts=countData, group=Group)
y <- calcNormFactors(y) 

# Filter out rows that contain only 1 or 0
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

# Calculating Dispersion
y <- estimateDisp(y, design)

################################################################
# Actual Analysis
################################################################

###########################
# Quasi-likelihood F-test #
###########################
# Fit for F-test
QLfit <- glmQLFit(y, design)
# Making contrast matrix for QL F-test
{Contrasts}
# QL F-test
anova <- glmQLFTest(QLfit, contrast=QLcon)
topTags(anova)

#############
# GLM Stuff #
#############
# Generalized Linear Model Fit
GLfit <- glmFit(y, design)

# Perform all possible comparisons
{GeneralLinearModel}

# Make regionReport
report <- edgeReport(dge = y, object = {object},
                       project='edgeR Analysis',
                       intgroup='group',
                       outdir='edgeR_Report',
                       ouput='index',
                       theme=theme_linedraw(),
                       output_format='pdf_document',
                       device='pdf')

# unload the libraries
detach("package:edgeR")
detach("package:ggplot2")
detach("package:regionReport")
detach("package:docopt")

# Make file to notify when done
fileConn <- file("FINISHED.txt")
writeLines("We are finished", fileConn)
close(fileConn)
"""
    Combination, GLM, RRN, mf = getContext(jsontoRead)
    Context = {
            "Contrasts": Combination,
            "GeneralLinearModel": GLM,
            "object": RRN,
            "mainfeature": mf
            }
    with open(name,'w') as f:
        f.write(makeReportsTemplate.format(**Context))

def getContext(jsontoRead):
    from pprint import pprint
    Metadata = jsontoRead
    mainFeature = Metadata['MainFeature']

    Contrasts = 'QLcon <- makeContrasts(\n\t' + ',\n\t'.join(["{} = {}".format('{}_vs_{}'.format(a,b),'{} - {}'.format(a,b)) for (a,b) in list(combinations(set([Metadata['Samples'][samp]['Features'][mainFeature] for samp in Metadata['Samples']]),2))]) + ', levels=design)'

    newlist = []
    feats = list(set([Metadata['Samples'][samp]['Features'][mainFeature] for samp in Metadata['Samples']]))
    numbers = [1] + [1] + [0 for z in range(len(feats)-2)]
    for g in list(set(permutations(numbers))):
        g = list(g)
        g[g.index(int(1))] = -1
        newlist.append(g)
    versus = ['{}vs{}'.format(feats[c.index(-1)],feats[c.index(1)]) for c in newlist]

    Mycontrasts = 'mycontrasts <- makeContrasts(\n\t{},\n\tlevels = design)'.format(',\n\t'.join(['{} = {}'.format(a,b) for a,b in list(zip(versus,['{} - {}'.format(feats[c.index(-1)],feats[c.index(1)]) for c in newlist]))]))
    Comparisons = '\n'.join(['lrt_{} <- glmLRT(GLfit, contrast=mycontrasts[,"{}"])'.format(a,a) for a in versus])
    Toptags = '\n'.join(['topTags(lrt_{})'.format(v) for v in versus])

    GLM = Mycontrasts + '\n\n' + Comparisons + '\n\n' + Toptags
    regionReportName = 'lrt_' + versus[0]

    return Contrasts, GLM, regionReportName, str(mainFeature)

if __name__ == '__main__':
    arguments = docopt(__doc__,version='1.0')
    #TODO Check for Cols.dat
    #subprocess.run(["python","makeCols.py"],check=True)
    createEdgeR(makeCols.readJSON(arguments['-j']),arguments['-t'])
