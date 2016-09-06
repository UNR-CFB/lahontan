#!/usr/bin/python3

'''Usage: Mega.py [-h | --help] [-j <jsonFile> | --jsonfile <jsonFile>] [--noconfirm] <pathtoInput>

Options:
    -h --help               Show this screen and exit
    -j <jsonFile>, --jsonfile <jsonFile>            If already available, path to json file 
    --noconfirm             Ignore all user prompts except JSON file creation
'''

################################################################
# Importations
################################################################

from docopt import docopt
from twilio.rest import TwilioRestClient
import subprocess
import os
import multiprocessing
import shutil
import sys
import time
import makeJSON
import makeCols
import makeReportr

################################################################
# Useful Utilities
################################################################

def textMe(toPhoneNumber,message):
    ''' Arguments:
            toPhoneNumber = string; phone number that you wish message to be sent to
            message = string; message that you wish to send   
        Returns:
            None; Texts toPhoneNumber the message
    '''
    accountSID = 'AC17039a9f9b17f2ae2d188ca54db201ac'
    authToken = 'da2a8dbeca04139270ed22e4d0224b5f'
    twilioNumber = '+17756004137'
    twilioCli = TwilioRestClient(accountSID, authToken)
    message = twilioCli.messages.create(body=str(message), from_=twilioNumber, to=toPhoneNumber)

    Notification = 'This message was sent to {}:\n{}'.format(toPhoneNumber,message.body)
    print(Notification)

def getNumberofFiles(path):
    ''' Arguments:
            path = string of the path to a directory where you
                    want to count files
        Returns:
            counter = the number of files in path
    '''
    if os.path.isdir(path) == False:
        print("Path for getNumberofFiles is not a directory:\n{}".format(path))
        raise SystemExit

    num = len(os.listdir(path))
    return num

################################################################
# Getting INPUT variables
################################################################

def sourceInput(pathtoInput):
    ''' Arguments:
            pathtoInput = a string of the path to your input file
        Returns:
            Project, Reference, Original paths
        Note: INPUT file should contain:
                "Project" path
                "Reference" path
                "Original" path
    '''
    import imp

    assert type(pathtoInput) == str, "Path to INPUT not a string"
    if os.path.exists(pathtoInput) == False:
        print("Path is not a file:\n{}".format(pathtoInput))
        raise SystemExit

    with open(pathtoInput) as Input:
        global inputVars
        inputVars = imp.load_source('inputVars', '', Input)

    return inputVars.Project, inputVars.Reference, inputVars.Original

def getReferenceVariables(referencePath):
    ''' Arguments:
            referencePath = string; path to Reference files
        Returns:
            gtf, cdna, genome filenames
    '''
    if os.path.isdir(referencePath) == False:
        print("Path is not a directory:\n{}".format(path))
        raise SystemExit

    ls = os.listdir(referencePath)
    if len(ls) != 3:
        print("There are an incorrect number of files in:\n{}".format(referencePath))
        raise SystemExit

    for filename in ls:
        if filename.split(".")[-1] == 'gtf':
            gtf = str(filename)
        elif "cdna" in filename:
            cdna = str(filename)
        else:
            genome = str(filename)
    return gtf,cdna,genome

def getBasename(genomeName):
    ''' Arguments:
            genomeName = string; name of Genome file
        Returns:
            basename
    '''
    basename = str(genomeName.split(".")[0])
    return basename

def getFastq(originalPath):
    ''' Arguments:
            originalPath = string; path to ogOriginal
        Returns:
            fastq
    '''
    if os.path.isdir(originalPath) == False:
        print("Path is not a directory:\n{}".format(originalPath))
        raise SystemExit

    allFastq = []
    for filename in os.listdir(originalPath):
        allFastq.append(filename.split(".")[-2])

    uniqueFastq = set(allFastq)
    if len(uniqueFastq) != 1:
        raise NameError('File Endings in ogOriginal are not consistent')
    else:
        fastq = list(uniqueFastq)[0]

    if fastq != 'fastq' and fastq != 'fq':
        print('Unknown fastq extension: {}\nSetting default to "fastq"')
        fastq = 'fastq'

    return fastq

def countCPU():
    ''' Arguments:
            None
        Returns:
            PROCS (number of processors)
    '''
    procs = int(multiprocessing.cpu_count())
    return procs

def exportVariables(pathtoInput):
    ''' Arguments:
            pathtoInput = string; path to INPUT file
        Returns:
            Variables
        Note: INPUT file must contain Project, Reference, and Original paths
    '''
    if os.path.isdir(pathtoInput) == False:
        print("Path is not a directory:\n{}".format(pathtoInput))
        raise SystemExit

    Project, Reference, Original = sourceInput(pathtoInput)
    Gtf, Cdna, Genome = getReferenceVariables(Reference)
    Basename = getBasename(Genome)
    Fastq = getFastq(Original)
    Procs = countCPU()
    Variables = [Project,Reference,Original,
                Gtf,Cdna,Genome,
                Basename,Fastq,Procs]
    return Variables

################################################################
# Making Structure
################################################################

def createStructure(projectPath, originalPath):
    ''' Arguments:
            projectPath = string of the path to the Project Directory
            originalPath = string of the path to the Original Directory
        Returns:
            None; creates directory structure
    '''
    if os.path.isdir(originalPath) == False:
        print("Path is not a directory:\n{}".format(originalPath))
        raise SystemExit

    if getNumberofFiles(originalPath)%2 != 0:
        print('There are not an even number of files in {}'.format(originalPath))
        raise SystemExit
    else:
        numSamp = str(getNumberofFiles(originalPath)/2)

    subprocess.run(["makeStructure.sh", projectPath, numSamp],check=True)

    if noconfirm == False:
        if os.path.isdir(projectPath + '/Postprocessing') == False:
            print("Path is not a directory:\n{}".format(projectPath + '/Postprocessing'))
            raise SystemExit
        if os.path.isdir(projectPath + '/Reference') == False:
            print("Path is not a directory:\n{}".format(projectPath + '/Reference'))
            raise SystemExit
        if os.path.isdir(projectPath + '/Data') == False:
            print("Path is not a directory:\n{}".format(projectPath + '/Data'))
            raise SystemExit

################################################################
# Make Pointers to Data
################################################################

def createSymLinks(projectPath,originalPath,referencePath):
    ''' Arguments:
            projectPath = string of the path to the Project Directory
            originalPath = string of the path to the Original Directory
            referencePath = string of the path to the Reference Directory
        Returns:
            None; creates Sym links
    '''
    if noconfirm == False:
        if os.path.isdir(projectPath) == False:
            print("Path is not a directory:\n{}".format(projectPath))
            raise SystemExit
        if os.path.isdir(originalPath) == False:
            print("Path is not a directory:\n{}".format(originalPath))
            raise SystemExit
        if os.path.isdir(referencePath) == False:
            print("Path is not a directory:\n{}".format(referencePath))
            raise SystemExit

    if JSFI == None:
        subprocess.run(["makeSyms.sh",projectPath,originalPath,referencePath,'false'],check=True)
    else:
        subprocess.run(["makeSyms.sh",projectPath,originalPath,referencePath,JSFI],check=True)

    if noconfirm == False:
        if getNumberofFiles(originalPath)%2 != 0:
            print('There are not an even number of files in {}!'.format(originalPath))
        else:
            numSamp = int(getNumberofFiles(originalPath)/2)
        if len(os.listdir(projectPath + '/Reference')) != 3:
            print("There are not an appropriate amount of Symlinks in {}".format(projectPath + '/Reference'))
            raise SystemExit
        if len(os.listdir(projectPath + '/Data')) != numSamp:
            print("There are not an appropriate number of Samples in {}".format(projectPath + '/Data'))
            raise SystemExit
        if len(os.listdir(projectPath + '/Original')) != numSamp*2:
            print("There are not an appropriate number of samples in {}".format(projectPath + '/Original'))
            raise SystemExit


################################################################
# Quality Control of Reference Data
################################################################

def qcReference(referencePath,genomeName):
    ''' Arguments:
            referencePath = string of the path to the Reference Directory
            genomeName = string; name of the genome file
        Returns:
            None; outputs Reference_Report.txt to STDOUT
    '''
    if os.path.isdir(referencePath) == False:
        print("Path is not a directory:\n{}".format(referencePath))
        raise SystemExit
    if genomeName not in os.listdir(referencePath):
        print("{} not in {}".format(genomeName, referencePath))
        raise SystemExit

    subprocess.run(["QCofRef.sh",referencePath,genomeName],check=True)

    if noconfirm == False:
        os.chdir(referencePath)
        with open('Reference_Report.txt', 'r') as Report:
            print(Report.read())
        while True:
            answer = str(input("Would you like to review your reference data?(y,n) "))
            if answer == 'y':
                print('See Reference_Report.txt to view initial Diagnostics;\
                        exiting now')
                raise SystemExit
            elif answer == 'n':
                break
            else:
                print('Please answer y or n')

################################################################
# Pre-Processing of Reference Data once QC complete
################################################################

def preProcessingReference(referencePath,cdnaName,gtfName,genomeName,baseName):
    ''' Arguments:
            referencePath = string of the path to the Reference Directory
            cdnaName = string; name of Cdna file
            gtfName = string; name of Gtf file
            genomeName = string; name of the genome file
            baseName = string; base name created
        Returns:
            None; outputs Reference_Report.txt to STDOUT
    '''
    if noconfirm == False:
        if os.path.isdir(referencePath) == False:
            print("Path is not a directory:\n{}".format(referencePath))
            raise SystemExit
        if genomeName not in os.listdir(referencePath):
            print("{} not in {}".format(genomeName, referencePath))
            raise SystemExit
        if cdnaName not in os.listdir(referencePath):
            print("{} not in {}".format(cdnaName, referencePath))
            raise SystemExit
        if gtfName not in os.listdir(referencePath):
            print("{} not in {}".format(gtfName, referencePath))
            raise SystemExit

    os.chdir(referencePath)
    with open('Preprocessing.log','w') as PPlog:
        subprocess.run(["PPofRef.sh",referencePath,cdnaName,gtfName,genomeName,baseName],stdout=PPlog,stderr=subprocess.STDOUT,check=True)

    if noconfirm == False:
        if 'Reference_Report.txt' not in os.listdir(referencePath):
            print("{} not in {}".format('Reference_Report.txt', referencePath))
            raise SystemExit
        if 'splice_sites.txt' not in os.listdir(referencePath):
            print("{} not in {}".format('splice_sites.txt', referencePath))
            raise SystemExit
        if 'known_exons.txt' not in os.listdir(referencePath):
            print("{} not in {}".format('known_exons.txt', referencePath))
            raise SystemExit

################################################################
# Run Pipeline
################################################################
#Split up timeit3 into components
#       Definitely into the QC and actual alignment stuff

def runPipe(dataPath,fastqVar,numProcs,referencePath,genomeName,baseName,projectPath,gtfName):
    ''' Arguments:
            dataPath = string; path to Data directory
            fastqVar = string; fq vs fastq
            numProcs = integer; number of processors
            referencePath = string of the path to the Reference Directory
            genomeName = string; name of the genome file
            baseName = string; base name created
            projectPath = string; path to Project directory
            Gtf = string; name of Gtf file
        Returns:
            None; Runs Pipeline
    '''
    subprocess.run(["runall.sh",dataPath,fastqVar,str(numProcs),referencePath,genomeName,baseName,projectPath,gtfName],check=True)

def findFinish(projectPath, originalPath):
    ''' Arguments:
            projectPath = string; path to Project Directory
            originalPath = string; path to Original Directory
        Returns:
    '''
    if getNumberofFiles(originalPath)%2 != 0:
        print('There are not an even number of files in {}!'.format(originalPath))
    else:
        numSamp = int(getNumberofFiles(originalPath)/2)

    if os.path.isdir(projectPath + '/runPipeNotify') == False:
        print("Path is not a directory:\n{}".format(projectPath + '/runPipeNotify'))
        raise SystemExit
    else:
        os.chdir(projectPath + '/runPipeNotify')

    while True:
        if len([name for name in os.listdir('.') if os.path.isfile(name)]) == numSamp:
            #textMe('+17756227884','This is sent from your python script; The first part of the \
            #        Pipeline has finished. Please return to computer')
            print('At this point textMe would send you a text')
            os.chdir(projectPath)
            shutil.rmtree(projectPath + '/runPipeNotify')
            break
        else:
            time.sleep(30)

################################################################
# Making JSON file
################################################################

def createMetaData(postProcessingPath,jsonName='Metadata.json'):
    ''' Arguments:
            postProcessingPath = string; path to the Postprocessing Directory
        Returns:
            None; Creates Metadata.json
    '''
    if JSFI == None:
        if os.path.isdir(postProcessingPath) == False:
            print("Path is not a directory:\n{}".format(postProcessingPath))
            raise SystemExit
        os.chdir(postProcessingPath)
        makeJSON.writeJSON(jsonName)
        while True:
            answer = str(input("Would you like to review {}?(y,n) ".format(jsonName)))
            if answer == 'y':
                os.chdir(postProcessingPath)
                subprocess.run(['vim',jsonName],check=True)
            elif answer == 'n':
                break
            else:
                print('Please answer y or n')

################################################################
# Evaluating Results of Pipeline
################################################################

def makeNiceCounts(postProcessingPath,dataPath):
    ''' Arguments:
            postProcessingPath = string; path to the Postprocessing Directory
            dataPath = string; path to Data Directory
        Returns:
            None; Creates NiceCounts.dat
    '''
    if os.path.isdir(postProcessingPath) == False:
        print("Path is not a directory:\n{}".format(postProcessingPath))
        raise SystemExit
    if os.path.isdir(dataPath) == False:
        print("Path is not a directory:\n{}".format(dataPath))
        raise SystemExit

    subprocess.run(["formatFeatures.sh",postProcessingPath,dataPath],check=True)
    
def makeTotalTime(postProcessingPath,dataPath):
    ''' Arguments:
            postProcessingPath = string; path to the Postprocessing Directory
            dataPath = string; path to Data Directory
        Returns:
            None; Creates totalTime.dat
    '''
    if os.path.isdir(postProcessingPath) == False:
        print("Path is not a directory:\n{}".format(postProcessingPath))
        raise SystemExit
    if os.path.isdir(dataPath) == False:
        print("Path is not a directory:\n{}".format(dataPath))
        raise SystemExit

    subprocess.run(["getTotalTime.sh",postProcessingPath,dataPath],check=True)

################################################################
# Preparing for DESeq2
################################################################

def createColumnFile(postProcessingPath,jsonName='Metadata.json',columnName='Cols.dat'):
    ''' Arguments:
            postProcessingPath = string; path to the Postprocessing Directory
            jsonName = string; name of JSON file
            columnName = string; name of Column file
        Returns:
            None; Creates Cols.dat
    '''
    if os.path.isdir(postProcessingPath) == False:
        print("Path is not a directory:\n{}".format(postProcessingPath))
        raise SystemExit

    os.chdir(postProcessingPath)
    makeCols.makeCols(makeCols.readJSON(jsonName),columnName)

def createCountFile(postProcessingPath,countName='Counts.dat'):
    ''' Arguments:
            postProcessingPath = string; path to the Postprocessing Directory
            countName = string; name of Counts file
        Returns:
            None; Creates Counts.dat
    '''
    if os.path.isdir(postProcessingPath) == False:
        print("Path is not a directory:\n{}".format(postProcessingPath))
        raise SystemExit

    os.chdir(postProcessingPath)
    subprocess.run(["makeCounts.sh",postProcessingPath,countName],check=True)

def createRScript(postProcessingPath,jsonName='Metadata.json',rName='makeReport.r'):
    ''' Arguments:
            postProcessingPath = string; path to the Postprocessing Directory
            dataPath = string; path to Data Directory
        Returns:
            None; Creates makeReportr.dat
    '''
    if os.path.isdir(postProcessingPath) == False:
        print("Path is not a directory:\n{}".format(postProcessingPath))
        raise SystemExit

    os.chdir(postProcessingPath)
    makeReportr.createRscript(makeCols.readJSON(jsonName),rName)

################################################################
# Running DESeq2
################################################################

def makeRreports(postProcessingPath):
    ''' Arguments:
            postProcessingPath = string; path to Postprocessing Directory
        Returns:
            None; Creates R reports
    '''
    if noconfirm == False:
        if os.path.isdir(postProcessingPath) == False:
            print("Path is not a directory:\n{}".format(postProcessingPath))
            raise SystemExit
        if 'totalTime.dat' not in os.listdir(postProcessingPath):
            print("{} is not in {}".format('totalTime.dat',postProcessingPath))
            raise SystemExit
        if 'Cols.dat' not in os.listdir(postProcessingPath):
            print("{} is not in {}".format('Cols.dat',postProcessingPath))
            raise SystemExit
        if 'Counts.dat' not in os.listdir(postProcessingPath):
            print("{} is not in {}".format('Counts.dat',postProcessingPath))
            raise SystemExit
        if 'makeReport.r' not in os.listdir(postProcessingPath):
            print("{} is not in {}".format('makeReport.r',postProcessingPath))
            raise SystemExit

    subprocess.run(["runDESeq.sh",postProcessingPath],check=True)

def notifyEnding(postProcessingPath):
    ''' Arguments:
            postProcessingPath = string; path to Postprocessing Directory
        Returns:
            None
    '''
    if os.path.isdir(postProcessingPath) == False:
        print("Path is not a directory:\n{}".format(postProcessingPath))
        raise SystemExit

    while True:
        if 'FINISHED.txt' in os.listdir(postProcessingPath):
            #textMe('+17756227884','This is sent from your python script; The first part of the \
            #        Pipeline has finished. Please return to computer')
            print('DESeq2 reports are finished')
            os.remove(str(postProcessingPath + "/FINISHED.txt"))
            break
        else:
            time.sleep(30)


#################################################################
#def exportVariablesforClass(pathtoInput):
#    ''' Arguments:
#        Returns:
#    '''
#    Project, Reference, Original = sourceInput(pathtoInput)
#    Gtf, Cdna, Genome = getReferenceVariables(Reference)
#    Basename = getBasename(Genome)
#    Fastq = getFastq(Original)
#    Procs = countCPU()
#    Variables = {
#            "Projectpath": Project,
#            "ogReference": Reference,
#            "ogOriginal": Original,
#            "Gtf": Gtf,
#            "Cdna": Cdna,
#            "Genome": Genome,
#            "Basename": Basename,
#            "Fastq": Fastq,
#            "Procs": Procs
#            }
#    return Variables
#
#class Experiment:
#    def __init__(self, inputPath):
#        variables = exportVariablesforClass(inputPath)
#        self.Project = variables["Projectpath"]
#        self.ogReference = variables["ogReference"]
#        self.ogOriginal = variables["ogOriginal"]
#        self.Reference = str(variables["Projectpath"] + '/Reference')
#        self.Original = str(variables["Projectpath"] + '/Original')
#        self.Data = str(variables["Projectpath"] + '/Data')
#        self.Postprocessing = str(variables["Projectpath"] + '/Reference')
#        self.Gtf = str(variables["Gtf"])
#        self.Cdna = str(variables["Cdna"])
#        self.Genome = str(variables["Genome"])
#        self.Basename = str(variables["Basename"])
#        self.Fastq = str(variables["Fastq"])
#        self.Procs = int(variables["Procs"])
#    def __repr__(self):
#        return 'Experiment(%r)'%(self.Project)
#    ###############################################################
#    # Utilities
#    ###############################################################
#    def getNumberofSamples(self):
#        if getNumberofFiles(self.ogOriginal)%2 != 0:
#            print('There are not an even number of files in {}'.format(self.ogOriginal))
#        else:
#            numSamp = getNumberofFiles(self.ogOriginal)/2
#        return int(numSamp)
#    def makeStructure(self):
#        createStructure(self.Project, self.ogOriginal)
#    def makeSyms(self):
#        createSymLinks(self.Project, self.ogOriginal, self.ogReference)
#    def qcRef(self):
#        qcReference(self.Reference, self.Genome)
#    def ppRef(self):
#        preProcessingReference(self.Reference, self.Cdna, self.Gtf, self.Genome, self.Basename)
#    def deployPipe(self):
#        runPipe(self.Data, self.Fastq, self.Procs, self.Reference,\
#                self.Genome, self.Basename, self.Project, self.Gtf)
#    def findPipeFinish(self):
#        findFinish(self.Project, self.ogOriginal)
#    def createJsonMetadata(self):
#        createMetaData(self.Postprocessing)
#    def createNiceCounts(self):
#        makeNiceCounts(self.Postprocessing, self.Data)
#    def getPipeTime(self):
#        makeTotalTime(self.Postprocessing, self.Data)
#    def createRCounts(self):
#        createCountFile(self.Postprocessing)
#    def createRCols(self):
#        createColumnFile(self.Postprocessing)
#    def createRProgram(self):
#        createRScript(self.Postprocessing)
#    def runRProgram(self):
#        makeRreports(self.Postprocessing)
################################################################
def Main(pathtoInput):
    ''' Arguments:
            pathtoInput = string; path to INPUT file
        Returns:
            :)
    '''
    Project,Reference,Original,Gtf,Cdna,Genome,Basename,Fastq,Procs = exportVariables(pathtoInput)
    createStructure(Project,Original)
    createSymLinks(Project,Original,Reference)
    qcReference(Project + '/Reference', Genome)
    preProcessingReference(Project + '/Reference', Cdna, Gtf, Genome, Basename)
    runPipe(Project + '/Data', Fastq, Procs, Project + '/Reference', Genome, Basename,Project,Gtf)
    createMetaData(Project + '/Postprocessing')
    findFinish(Project, Original)
    makeNiceCounts(Project + '/Postprocessing', Project + '/Data')
    makeTotalTime(Project + '/Postprocessing', Project + '/Data')
    createColumnFile(Project + '/Postprocessing')
    createCountFile(Project + '/Postprocessing')
    createRScript(Project + '/Postprocessing')
    # TODO check for X11 forwarding
    makeRreports(Project + '/Postprocessing')
    notifyEnding(Project + '/Postprocessing')

if __name__ == '__main__':
    arguments = docopt(__doc__,version='Version 1.1\nAuthor: Alberto')
    print(arguments)

    global noconfirm
    noconfirm = arguments['--noconfirm']
    global JSFI
    JSFI = arguments['--jsonfile']

    Main(arguments['<pathtoInput>'])
