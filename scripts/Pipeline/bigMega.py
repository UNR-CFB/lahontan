#!/usr/bin/python3

'''Usage: Mega.py [-h | --help] [-j <jsonFile> | --jsonfile <jsonFile>] [--noconfirm] [-c <placeToClean> | --clean <placeToClean>] [-s <sampleName> | --sampleclean <sampleName>] <pathtoInput>

Options:                                                                                           
    -h --help               Show this screen and exit
    -j <jsonFile>, --jsonfile <jsonFile>            Ignores JSON file creation and uses specified
                                                    path to JSON
    --noconfirm             Ignore all user prompts except JSON file creation
    -c <placeToClean>, --clean <placeToClean>       Cleans <placeToClean>: Possible places include:
                                                    Reference, Data, Postprocessing, All
    -s <sampleName>, --sampleclean <sampleName>     Similar to -c,--clean; but instead just cleans a
                                                    single sample directory <sampleName>
'''

################################################################
# Importations
################################################################

from docopt import docopt
import Mega
import os
import subprocess
import shutil

################################################################
# Utilities
################################################################

def exportVariablesforClass(pathtoInput):
    ''' Arguments:
            pathtoInput = string; path to INPUT file
        Returns:
            Variables
    '''
    Project, Reference, Original = Mega.sourceInput(pathtoInput)
    Gtf, Cdna, Genome = Mega.getReferenceVariables(Reference)
    Basename = Mega.getBasename(Genome)
    Fastq = Mega.getFastq(Original)
    Procs = Mega.countCPU()
    Variables = {
            "Projectpath": Project,
            "ogReference": Reference,
            "ogOriginal": Original,
            "Gtf": Gtf,
            "Cdna": Cdna,
            "Genome": Genome,
            "Basename": Basename,
            "Fastq": Fastq,
            "Procs": Procs
            }
    return Variables

def checkJSON(jsonFile):
    import json
    with open(jsonFile) as JF:
        try:
            json.load(JF)
        except ValueError as e:
            print("Invalid JSON: %s"%(e))
            print("Check your JSON, exiting now")
            raise SystemExit

################################################################
# Defining Experiment Class
################################################################

class Experiment:
    ''' Experiment is a Paired-End Sequencing Experiment.
        This class provides tools to analyze data
    '''

    def __init__(self, inputPath):
        variables = exportVariablesforClass(inputPath)
        self.Project = variables["Projectpath"]
        self.ogReference = variables["ogReference"]
        self.ogOriginal = variables["ogOriginal"]
        self.Reference = str(variables["Projectpath"] + '/Reference')
        self.Original = str(variables["Projectpath"] + '/Original')
        self.Data = str(variables["Projectpath"] + '/Data')
        self.Postprocessing = str(variables["Projectpath"] + '/Postprocessing')
        self.Gtf = str(variables["Gtf"])
        self.Cdna = str(variables["Cdna"])
        self.Genome = str(variables["Genome"])
        self.Basename = str(variables["Basename"])
        self.Fastq = str(variables["Fastq"])
        self.Procs = int(variables["Procs"])

    def __repr__(self):
        return 'Experiment(%r)'%(self.Project)

    ###############################################################
    # Utilities
    ###############################################################

    def getNumberofSamples(self):
        if Mega.getNumberofFiles(self.ogOriginal)%2 != 0:
            print('There are not an even number of files in {}'.format(self.ogOriginal))
        else:
            numSamp = Mega.getNumberofFiles(self.ogOriginal)/2
        return int(numSamp)

    def makeStructure(self):
        Mega.createStructure(self.Project, self.ogOriginal)

    def makeSyms(self):
        Mega.createSymLinks(self.Project, self.ogOriginal, self.ogReference)
        
    def qcRef(self):
        Mega.qcReference(self.Reference, self.Genome)

    def ppRef(self):
        Mega.preProcessingReference(self.Reference, self.Cdna, self.Gtf, self.Genome, self.Basename)

    def deployPipe(self):
        print("Pipeline is running...")
        Mega.runPipe(self.Data, self.Fastq, self.Procs, self.Reference,\
                self.Genome, self.Basename, self.Project, self.Gtf)

    def findPipeFinish(self):
        Mega.findFinish(self.Project, self.ogOriginal)

    def createJsonMetadata(self):
        Mega.createMetaData(self.Postprocessing)

    def createNiceCounts(self):
        Mega.makeNiceCounts(self.Postprocessing, self.Data)

    def getPipeTime(self):
        Mega.makeTotalTime(self.Postprocessing, self.Data)

    def createRCounts(self):
        os.chdir(self.Postprocessing)
        Mega.createCountFile(self.Postprocessing)

    def createRCols(self):
        os.chdir(self.Postprocessing)
        Mega.createColumnFile(self.Postprocessing)

    def createRProgram(self):
        Mega.createRScript(self.Postprocessing)
        
    def runRProgram(self):
        Mega.makeRreports(self.Postprocessing)

    def findRFinish(self):
        Mega.notifyEnding(self.Postprocessing)

    ################################################################
    # Cleaning Functions
    ################################################################

    def clean(self, thingToClean, sampleName=False):
        ''' Arguments:
                thingToClean = string:
                                Reference or Data or Postprocessing or All or Sample
                Note: If thingToClean=Sample, require second argument:
                sampleName = string; name of Sample to be cleaned
            Returns:
                None; cleans directories
        '''
        assert type(thingToClean) == str, '{} is not a valid argument'.format(thingToClean)
        if thingToClean == 'Reference':
            arg = '-r'
        elif thingToClean == 'Data':
            arg = '-d'
        elif thingToClean == 'Postprocessing':
            arg = '-p'
        elif thingToClean == 'All':
            arg = '-a'
        elif thingToClean == 'Sample':
            if sampleName == False:
                print('Need a valid sample name')
                raise SystemExit
            arg = sampleName
        else:
            print('Need a valid argument: Reference, Data, Postprocessing, All, or Sample')
            raise SystemExit
        if thingToClean == 'Sample':
            subprocess.run(['clean.sh','-s',arg, self.Genome, self.Cdna, self.Gtf, self.Reference, self.Data, self.Postprocessing],check=True)
        else:
            subprocess.run(['clean.sh',arg, self.Genome, self.Cdna, self.Gtf, self.Reference, self.Data, self.Postprocessing],check=True)

    def nukeProject(self):
        '''
        Removes entire Project Directory Structure
        '''
        while True:
            answer = input('Are you sure you want to remove entire Project?(y,n) ')
            if answer == 'y':
                shutil.rmtree(self.Project)
                break
            elif answer == 'n':
                break
            else:
                print('Please answer y or n')

    ################################################################
    # Useful Functions
    ################################################################

    def runStage1(self):
        ''' Note: does not need to be run more than once

        Prepare for Pipeline
        '''
        self.makeStructure()
        self.makeSyms()

    def runStage2(self):
        '''
        Prepare Reference Data
        '''
        self.qcRef()
        self.ppRef()

    def runStage3(self):
        ''' Note: In future deployPipe will be split up into
            multiple functions

        Run Pipeline
        '''
        self.deployPipe()
        self.findPipeFinish()

    def runStage4(self):
        ''' Note: Creating a metadata file should only need
            to be done once

        Prepare for DESeq2 Analysis
        '''
        self.createJsonMetadata()
        self.createNiceCounts()
        self.getPipeTime()
        self.createRCounts()
        self.createRCols()
        self.createRProgram()

    def runStage5(self):
        '''
        Run R pipeline
        '''
        self.runRProgram()
        self.findRFinish()

    def runAll(self):
        '''
        Runs Stage 1 through Stage 4
        '''
        self.runStage1()
        self.runStage2()
        self.runStage3()
        self.runStage4()
        self.runStage5()

################################################################
# Handling Command line arguments
################################################################

if __name__ == '__main__':
    arguments = docopt(__doc__, version='Version 1.1\nAuthor: Alberto')

    ######################
    ## Area for Testing
    #PROJ = Experiment(arguments['<pathtoInput>'])
    #PROJ.nukeProject()
    #raise SystemExit
    #####################

    # Handling Cleaning Arguments
    possibleCleanArguments = ['All','Reference','Data','Postprocessing']
    if arguments['--clean'] != None:
        assert arguments['--clean'] in possibleCleanArguments, 'Invalid Cleaning Argument: Run bigMega.py -h for available arguments'
        PROJ = Experiment(arguments['<pathtoInput>'])
        if arguments['--clean'] != 'All':
            if arguments['--clean'] == 'Data':
                if os.path.isdir(PROJ.Data) == False:
                    print('{} does not exist'.format(PROJ.Data))
                    raise SystemExit
            if arguments['--clean'] == 'Reference':
                if os.path.isdir(PROJ.Reference) == False:
                    print('{} does not exist'.format(PROJ.Reference))
                    raise SystemExit
            if arguments['--clean'] == 'Postprocessing':
                if os.path.isdir(PROJ.Postprocessing) == False:
                    print('{} does not exist'.format(PROJ.Postprocessing))
                    raise SystemExit
        PROJ.clean(arguments['--clean'])
        raise SystemExit
    elif arguments['--sampleclean'] != None:
        PROJ = Experiment(arguments['<pathtoInput>'])
        if os.path.isdir(str(PROJ.Data + '/' + arguments['--sampleclean'])) == False:
            print('{} does not exist'.format(str(PROJ.Data + '/' + arguments['--sampleclean'])))
            raise SystemExit
        PROJ.clean('Sample',arguments['--sampleclean'])
        raise SystemExit

    global noconfirm
    noconfirm = arguments['--noconfirm']
    global JSFI
    JSFI = arguments['--jsonfile']

    if JSFI != None:
        checkJSON(JSFI)

    Mega.JSFI = JSFI
    Mega.noconfirm = noconfirm

    PROJ = Experiment(arguments['<pathtoInput>'])
    PROJ.runAll()
