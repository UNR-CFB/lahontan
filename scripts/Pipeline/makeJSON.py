#!/usr/bin/python3

'''Usage:
    runPipe mj [-h | --help] [-j <jsonfile> | --jsonfile <jsonfile>]

Options:
    -h, --help
        Show this screen
    -j <jsonfile>, --jsonfile <jsonfile>
        Optional name of JSON file to be saved to [default: Metadata.json]'''

################################################################
# Importations
################################################################

from docopt import docopt
import json
import os

################################################################

def getProjectInfo():
    ''' Arguments:
            None
        Returns:
            projectName, numberofSamples, numberofFeatures
            = tuple containing:
                projectName = str; Name of the Project
                numberofSamples = int; number of samples
                numberofFeatures = int; numberofFeatures

        Asks for certain experiment characteristics
    '''
    projectName = str(input("What is the name of the project? "))
    numberofSamples = int(input("How many samples are there? "))
    numberofFeatures = int(input("How many features does each sample have? "))
    return projectName,numberofSamples,numberofFeatures

def getSampleInfo():
    ''' Arguments:
            None
        Returns:
            projectName,numberofSamples,numberofFeatures,
            featureNames,mainFeature,SampleData = tuple containing:
                projectName = str; Name of project
                numberofSamples = int; number of samples
                numberofFeatures = int; number of features
                featureNames = list; name of features
                mainFeature = str; name of main feature
                SampleData = dict; dict containing all sample names
                                and respective features

        Asks for specific sample characteristics
    '''
    projectName,numberofSamples,numberofFeatures = getProjectInfo()

    SampleData = {}
    featureNames = []
    for feature in range(numberofFeatures):
        featureNames.append(str(input('What is the name of feature #%d? '%(feature+1))))
    while True:
        mainFeature = str(input("What is the name of the main feature? "))
        if mainFeature not in featureNames:
            print('Need a valid feature name.')
        else:
            break

    for sample in range(numberofSamples):
        SampleData['sample_%.2d'%(sample+1)] = {}
        SampleData['sample_%.2d'%(sample+1)]['Features'] = {}
        for feature in featureNames:
            feat = str(input('What is the %s feature for sample #%d? '%(feature,sample+1)))
            SampleData['sample_%.2d'%(sample+1)]['Features'][feature] = feat
    return projectName,numberofSamples,numberofFeatures,featureNames,mainFeature,SampleData

def makeJSON():
    ''' Arguments:
            None
        Returns:
            MetaData = dict; dict of all Experimental Data

        Creates dictionary of metadata
    '''
    MetaData = {}

    projectName,numberofSamples,numberofFeatures,featureNames,mainFeature,sampleData = getSampleInfo()

    MetaData['ProjectName'] = projectName
    MetaData['NumberofSamples'] = numberofSamples
    MetaData['NumberofFeatures'] = numberofFeatures
    MetaData['FeatureNames'] = sorted(featureNames)
    MetaData['MainFeature'] = mainFeature
    MetaData['Samples'] = sampleData

    return MetaData

def writeJSON(name):
    ''' Arguments:
            name = str; name of JSON file to be written to
        Returns:
            None

        Converts Metadata dictionary to a JSON file
    '''
    MetaDict = makeJSON()
    with open(name,'w') as File:
        json.dump(MetaDict, File, sort_keys=True,indent=4)
    print('Done making {}'.format(name))

################################################################
if __name__ == '__main__':
    arguments = docopt(__doc__,version='1.0')
    writeJSON(arguments['--jsonfile'])
