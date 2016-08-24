'''Usage: asdf.py [-h | --help] [-f <file>]
                                                                                      
Options:                                                                              
    -h --help           Show this screen                                              
    -f <file>           Optional name of JSON file to be read [default: Metadata.json]
'''                                                                                   
from docopt import docopt                                                             
import csv
import subprocess
import json                                                                           
import makeCols

def AAA(jsontoRead):
    makeReports = """
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

# Make regionReport
dir.create('regionReport',showWarnings = FALSE,recursive = TRUE)
report <- DESeq2Report(ddsAdj, project = 'DESeq2 Report',
                        intgroup = c({10_formulalistcomma}),
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
detach("package:docopt")
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
    with open('makeReportspy.r','w') as f:
        f.write(makeReports.format(**Context))

def getContext(jsontoRead):
    Metadata = jsontoRead
    mainFeature = Metadata['MainFeature']

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
            newgroup = list(group)[1:]
            uniquegroup = set(newgroup)
            uniquegroup2 = list(uniquegroup)

    factorlist = ','.join('"{0}"'.format(a) for a in uniquegroup2) 
    formulalistcomma = ','.join('"{0}"'.format(b) for b in formulalist)
    formulalistplus = ' + '.join(formulalistp)
    return factorlist,formulalistplus,formulalistcomma,mainFeature

if __name__ == '__main__':
    arguments = docopt(__doc__,version='1.0')
    #subprocess.run(["python","makeCols.py"],check=True)
    AAA(makeCols.readJSON(arguments['-f']))
