#!/usr/bin/python3

'''Usage: Mega.py [-h | --help] [-j <jsonFile> | --jsonfile <jsonFile>] [--noconfirm] <pathtoInput>

Options:                                                                                           
    -h --help               Show this screen and exit
    -j <jsonFile>, --jsonfile <jsonFile>            Ignores JSON file creation and uses specified
                                                    path to JSON
    --noconfirm             Ignore all user prompts except JSON file creation
'''
#TODO Make a cleanup function or fix the other one

################################################################
# Importations
################################################################

from docopt import docopt
import Mega
import os

################################################################
# Defining Experiment Class
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
    # Useful Functions
    ################################################################

    def runStage1(self):
        ''' Note: does not need to be run more than once

        Prepare for Pipeline
        '''
        self.makeStructure()
        self.makeSyms()
        self.qcRef()
        self.ppRef()

    def runStage2(self):
        ''' Note: In future deployPipe will be split up into
            multiple functions

        Run Pipeline
        '''
        self.deployPipe()
        self.findPipeFinish()

    def runStage3(self):
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

    def runStage4(self):
        '''
        Run R pipeline
        '''
        self.runRProgram()
        self.findRFinish()

    def runStage5(self):
        print(':)')

    def runAll(self):
        self.runStage1()
        self.runStage2()
        self.runStage3()
        self.runStage4()
        self.runStage5()


################################################################

if __name__ == '__main__':
    arguments = docopt(__doc__, version='Version 1.1\nAuthor: Alberto')

    global noconfirm
    noconfirm = arguments['--noconfirm']
    global JSFI
    JSFI = arguments['--jsonfile']

    Mega.JSFI = JSFI
    Mega.noconfirm = noconfirm

    PROJ = Experiment(arguments['<pathtoInput>'])
    PROJ.runAll()
