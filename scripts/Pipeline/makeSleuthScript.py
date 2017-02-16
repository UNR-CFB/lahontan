#!/usr/bin/python3

'''Usage: makeBallgownScript.py [-h | --help] [-j <jsonfile>] [-t <tofile>]

Options:
    -h --help           Show this screen
    -j <jsonfile>       Optional name of JSON file to be read [default: Metadata.json]
    -t <tofile>         Optional name of R report script [default: runSleuth.r]
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

def createRSleuthScript(jsontoRead,name):
    ''' Arguments:
            jsontoRead = str; name of JSON Metadata file
                [default: Metadata.json]
            name = str; name for edgeR script
        Returns:
            None

        Creates Ballgown analysis script by filling in a template
    '''
    scriptTemplate =''''Usage: runSleuth.r [-h | --help] [-p </path/to/data>] [--cols <name>] [--base <name>]

Example:
Rscript runSleuth.r -p /path/to/project --cols ColumnsName

Options:
    -h --help           Show this
    -p <path>           Directory that contains Cols file [default: ./]
    --cols <name>       Name of phenotype file that describes experiment [default: Cols.dat]
    --base <name>       Name of Base Directory that contains Kallisto ouput 
                        [default: KallistoResults]'-> doc

# Set up command line argument library
library('docopt')

# Retrieve command line arguments
opts <- docopt(doc)

# Change to directory
setwd(opts$p)

# Load the necessary libraries
library('devtools')
library('sleuth')

base_dir <- opts$p
sample_id <- dir(file.path(base_dir,opts$base))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, opts$base, id, "abundance.h5"))
dataFrame <- read.csv(opts$cols, header=TRUE, stringsAsFactors=FALSE)
dataFrame <- dplyr::mutate(dataFrame, path=kal_dirs)
dataFrame
model <- sleuth_prep(dataFrame, ~ {mainfeature})
model <- sleuth_fit(model)
model <- sleuth_fit(model, ~1, 'reduced')
model <- sleuth_lrt(model, 'reduced', 'full')
models(model)
#sleuth_live(model) # to get a shiny server online which can be used to analyze data
results_table <- sleuth_results(model, 'reduced:full', test_type='lrt')
write.csv(results_table, "sleuthKallistoAnalysis.csv", row.names=FALSE)

# unload the libraries
detach('package:docopt')
detach('package:sleuth')
detach('package:devtools')'''
    mainFeature = getContext(jsontoRead)
    Context = {
            "mainfeature": '{}'.format(mainFeature),
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
    return mainFeature

################################################################
if __name__ == '__main__':
    arguments = docopt(__doc__,version='1.0')
    createRSleuthScript(makeCols.readJSON(arguments['-j']),arguments['-t'])
