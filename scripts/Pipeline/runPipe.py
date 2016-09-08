#!/usr/bin/python3

'''Usage: runPipe.py [-h | --help] [-j <jsonFile> | --jsonfile <jsonFile>] [--noconfirm] [-c <placeToClean> | --clean <placeToClean>] [-s <sampleName> | --sampleclean <sampleName>] [--NUKE] [-r <runlogPath> | --runtime <runlogPath>] <pathtoInput>

Options:                                                                                           
    -h --help                                       Show this screen and exit
    -j <jsonFile>, --jsonfile <jsonFile>            Ignores JSON file creation and uses specified
                                                    path to JSON
    --noconfirm                                     Ignore all user prompts except JSON file creation
    -c <placeToClean>, --clean <placeToClean>       Cleans <placeToClean>: Possible places include:
                                                    Reference, Data, Postprocessing, All
    -s <sampleName>, --sampleclean <sampleName>     Similar to -c,--clean; but instead just cleans a
                                                    single sample directory <sampleName>
    --NUKE                                          Removes entire project Directory
    -r <runlogPath>, --runtime <runlogPath>         Optional directory path for Runtime Log file to 
                                                    be created in [default: $Project]
'''

################################################################
# Importations
################################################################

from docopt import docopt
from timeit import default_timer as timer
import time
import pipeUtils
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
    Project, Reference, Original = pipeUtils.sourceInput(pathtoInput)
    Gtf, Cdna, Genome = pipeUtils.getReferenceVariables(Reference)
    Basename = pipeUtils.getBasename(Genome)
    Fastq = pipeUtils.getFastq(Original)
    Procs = pipeUtils.countCPU()
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

def funTime(function):
    def wrapper(*args):
        t1 = timer()
        stuff = function(*args)
        t2 = timer()
        with open(RUNTIMELOG, 'a') as R:
            R.write('{} function took {:.3f} seconds\n'.format(function.__name__, t2 - t1))
        return stuff
    return wrapper

def makeTimeFile(logPath):
    with open(logPath, 'w') as R:
        R.write('runPipe.py Runtime File\n-----------------------------------------\n\n')

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
        if pipeUtils.getNumberofFiles(self.ogOriginal)%2 != 0:
            print('There are not an even number of files in {}'.format(self.ogOriginal))
        else:
            numSamp = pipeUtils.getNumberofFiles(self.ogOriginal)/2
        return int(numSamp)

    @funTime
    def makeStructure(self):
        print("Creating Structure...")
        pipeUtils.createStructure(self.Project, self.ogOriginal)
        makeTimeFile(RUNTIMELOG)

    @funTime
    def makeSyms(self):
        pipeUtils.createSymLinks(self.Project, self.ogOriginal, self.ogReference)
        
    @funTime
    def qcRef(self):
        pipeUtils.qcReference(self.Reference, self.Genome)

    @funTime
    def ppRef(self):
        print("Preprocessing Data...")
        pipeUtils.preProcessingReference(self.Reference, self.Cdna, self.Gtf, self.Genome, self.Basename)

    @funTime
    def deployPipe(self):
        print("Pipeline is running...")
        pipeUtils.runPipe(self.Data, self.Fastq, self.Procs, self.Reference,\
                self.Genome, self.Basename, self.Project, self.Gtf)

    @funTime
    def findPipeFinish(self):
        pipeUtils.findFinish(self.Project, self.ogOriginal)

    @funTime
    def createJsonMetadata(self):
        print('Making Metadata file; require user input...')
        pipeUtils.createMetaData(self.Postprocessing)
        print('Waiting for Pipeline to finish...')

    @funTime
    def createNiceCounts(self):
        pipeUtils.makeNiceCounts(self.Postprocessing, self.Data)

    @funTime
    def getPipeTime(self):
        pipeUtils.makeTotalTime(self.Postprocessing, self.Data)

    @funTime
    def createRCounts(self):
        os.chdir(self.Postprocessing)
        pipeUtils.createCountFile(self.Postprocessing)

    @funTime
    def createRCols(self):
        os.chdir(self.Postprocessing)
        pipeUtils.createColumnFile(self.Postprocessing)

    @funTime
    def createRProgram(self):
        pipeUtils.createRScript(self.Postprocessing)
        
    @funTime
    def runRProgram(self):
        print("R program is running...")
        pipeUtils.makeRreports(self.Postprocessing)

    @funTime
    def findRFinish(self):
        pipeUtils.notifyEnding(self.Postprocessing)

    def checkX11(self):
        while True:
            answer = input('Are you on a local X session or is X11 forwarding enabled? You wont be able to run DESeq without it...(y,n) '
            if answer == 'y':
                print('Running...')
            elif answer == 'n':
                print('Exiting now')
                raise SystemExit
            else:
                print('Please answer y or n')
            

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
        if noconfirm == True:
            shutil.rmtree(self.Project)
        else:
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

    @funTime
    def runStage1(self):
        ''' Note: does not need to be run more than once

        Prepare for Pipeline
        '''
        if noconfirm == False:
            self.checkX11()
        self.makeStructure()
        self.makeSyms()

    @funTime
    def runStage2(self):
        '''
        Prepare Reference Data
        '''
        self.qcRef()
        self.ppRef()

    @funTime
    def runStage3(self):
        ''' Note: In future deployPipe will be split up into
            multiple functions

        Run Pipeline
        '''
        self.deployPipe()
        self.createJsonMetadata()
        self.findPipeFinish()

    @funTime
    def runStage4(self):
        ''' Note: Creating a metadata file should only need
            to be done once

        Prepare for DESeq2 Analysis
        '''
        print('Preparing for DESeq2...')
        self.createNiceCounts()
        self.getPipeTime()
        self.createRCounts()
        self.createRCols()
        self.createRProgram()

    @funTime
    def runStage5(self):
        '''
        Run R pipeline
        '''
        self.runRProgram()
        self.findRFinish()

    @funTime
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
    t1 = timer()

    global noconfirm
    noconfirm = arguments['--noconfirm']
    global JSFI
    JSFI = arguments['--jsonfile']


    # Handling --NUKE argument
    if arguments['--NUKE'] == True:
        if noconfirm == True:
            PROJ = Experiment(arguments['<pathtoInput>'])
            PROJ.nukeProject()
            raise SystemExit
        else:
            PROJ = Experiment(arguments['<pathtoInput>'])
            while True:
                answer = input('Are you sure you wish to remove {}?(y/n) '.format(PROJ.Project))
                if answer == 'y':
                    PROJ.nukeProject()
                    raise SystemExit
                elif answer == 'n':
                    raise SystemExit
                else:
                    print('Please answer y or n')

    # Handling Cleaning Arguments
    possibleCleanArguments = ['All','Reference','Data','Postprocessing']
    if arguments['--clean'] != None:
        assert arguments['--clean'] in possibleCleanArguments, 'Invalid Cleaning Argument: Run runPipe.py -h for available arguments'
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

    if JSFI != None:
        checkJSON(JSFI)

    pipeUtils.JSFI = JSFI
    pipeUtils.noconfirm = noconfirm

    PROJ = Experiment(arguments['<pathtoInput>'])

    global RUNTIMELOG
    if arguments['--runtime'] == '$Project':
        RUNTIMELOG = str(PROJ.Project + '/Runtime.log')
    else:
        if os.path.isdir(arguments['--runtime']) == False:
            print('Need a valid directory to save Runtime file to')
            raise SystemExit
        else:
            RUNTIMELOG = str(arguments['--runtime'] + '/Runtime.log')


    # Call for actually doing stuff
    ############################################################
    # 
    PROJ.runAll()
    #
    ############################################################
    #

    t2 = timer()

    timeused = str(time.strftime('%H:%M:%S', time.gmtime(t2-t1)))
    print('Total time elapsed: {}'.format(timeused))
