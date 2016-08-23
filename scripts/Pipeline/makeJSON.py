'''Usage: makeJSON.py [-h | --help] [-f <file>]

Options:
    -h --help           Show this screen
    -f <file>           Optional name of JSON file to be saved to [default: Metadata.json] 
'''
from docopt import docopt
import json

def getProjectInfo():
    projectName = str(input("What is the name of the project? "))
    numberofSamples = int(input("How many samples are there? "))
    numberofFeatures = int(input("How many features does each sample have? "))
    return projectName,numberofSamples,numberofFeatures

def getSampleInfo():
    projectName,numberofSamples,numberofFeatures = getProjectInfo()

    SampleData = {}
    featureNames = []
    for feature in range(numberofFeatures):
        featureNames.append(str(input('What is the name of feature #%d? '%(feature+1))))

    for sample in range(numberofSamples):
        SampleData['sample_%.2d'%(sample+1)] = {}
        SampleData['sample_%.2d'%(sample+1)]['Features'] = {}
        for feature in featureNames:
            feat = str(input('What is the %s feature for sample #%d? '%(feature,sample+1)))
            SampleData['sample_%.2d'%(sample+1)]['Features'][feature] = feat
    return projectName,numberofSamples,numberofFeatures,SampleData

def makeJSON():
    MetaData = {}

    projectName,numberofSamples,numberofFeatures,sampleData = getSampleInfo()

    MetaData['ProjectName'] = projectName
    MetaData['NumberofSamples'] = numberofSamples
    MetaData['NumberofFeatures'] = numberofFeatures
    MetaData['Samples'] = sampleData

    return MetaData

def writeJSON(name):
    MetaDict = makeJSON()

    with open(name,'w') as File:
        json.dump(MetaDict, File, sort_keys=True,indent=4)

if __name__ == '__main__':
    arguments = docopt(__doc__,version='1.0')
    writeJSON(arguments['-f'])
