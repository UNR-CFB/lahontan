#!/usr/bin/python3

'''Usage: makeCols.py [-h | --help] [-j <jsonfile>] [-t <tofile>]

Options:
    -h --help           Show this screen
    -j <jsonfile>           Optional name of JSON file to be read [default: Metadata.json] 
    -t <tofile>         Optional name of Column file to write to [default: Cols.dat]
'''

from docopt import docopt
from pprint import pprint
from collections import OrderedDict
import json
import os

def readJSON(name):
    with open(name) as jsonFile:
        jsonData = json.load(jsonFile,object_pairs_hook=OrderedDict)
    return jsonData

def parseJSON(jsontoRead):
    projectName = jsontoRead['ProjectName']
    numberofFeatures = jsontoRead['NumberofFeatures']
    numberofSamples = jsontoRead['NumberofSamples']
    featureNames = jsontoRead['FeatureNames']
    return featureNames,projectName,numberofFeatures,numberofSamples

def parseSamples(jsontoRead):
    SampleList = []
    counter = 0
    for sample in jsontoRead['Samples']:
        SampleList.append(["aligned.{}.bam".format(sample)])
        for feature in jsontoRead['Samples'][sample]['Features']:
            SampleList[counter].append(jsontoRead['Samples'][sample]['Features'][feature])
        counter += 1
    return SampleList

def makeCols(jsontoRead,filetoWrite):
    featureNames,_,__,___ = parseJSON(jsontoRead)
    sampleList = parseSamples(jsontoRead)
    deseqCols = []
    deseqCols.append(featureNames)
    for sample in sampleList:
        deseqCols.append(sample)
    with open(filetoWrite,'w') as colFile:
        colFile.writelines('\t'.join(row) + '\n' for row in deseqCols)

if __name__ == '__main__':
    arguments = docopt(__doc__,version='1.0')
    exec(open(os.environ["Input_Field"]).read())
    os.chdir(Postprocessing)
    makeCols(readJSON(arguments['-j']),arguments['-t'])
