#!/usr/bin/python3

'''Usage: Mega.py [-h | --help] <pathtoInput>

Options:
    -h --help               Show this screen and exit
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
    counter = 0
    for filename in os.listdir(path):
        counter += 1
    return counter

################################################################
# Getting INPUT variables
################################################################
IdealInput="""
Project="/home/alberton/Project"
Reference="/path/to/3/reference/files"
Original="/path/to/actual/data"
"""

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
    assert type(pathtoInput) == str, "Path to INPUT not a string"
    exec(open(pathtoInput).read())
    return project, reference, original

def getReferenceVariables(referencePath):
    ''' Arguments:
        Returns:
            Gtf, Cdna, Genome filenames
    '''
    for filename in os.listdir(referencePath):
        if filename.split(".")[-1] == 'gtf':
            gtf = filename
        elif "cdna" in filename:
            cdna = filename
        else:
            genome = filename
    return gtf,cdna,genome

def getBasename(genomeName):
    ''' Arguments:
        Returns:
    '''
    basename = genomeName.split(".")[0]
    return basename

def getFastq(originalPath):
    ''' Arguments:
        Returns:
            fastq
    '''
    allFastq = []
    for filename in os.listdir(originalPath):
        allFastq.append(filename.split(".")[-2])
    uniqueFastq = set(allFastq)
    if len(uniqueFastq) != 1:
        raise NameError('File Endings in $Original are not consistent')
    else:
        fastq = list(uniqueFastq)[0]
    return fastq

def countCPU():
    ''' Arguments:
            None
        Returns:
            PROCS
    '''
    procs = multiprocessing.cpu_count()
    return procs

def exportVariables(pathtoInput):
    ''' Arguments:
        Returns:
    '''
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
    if getNumberofFiles(originalPath)%2 != 0:
        print('There are not an even number of files in {}'.format(originalPath))
    else:
        numSamp = str(getNumberofFiles(originalPath)/2)
    subprocess.call(["makeStructure.sh",projectPath,numSamp])

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
    subprocess.call(["makeSyms.sh",projectPath,originalPath,referencePath])

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
    subprocess.call(["QCofRef.sh",referencePath,genomeName])

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
    subprocess.call(["PPofRef.sh",referencePath,cdnaName,gtfName,genomeName,baseName])

################################################################
# Run Pipeline
################################################################
#Split up timeit3 into components
#       Definitely into the QC and actual alignment stuff

def runPipe(dataPath,fastqVar,numProcs,referencePath,genomeName,baseName,projectPath):
    ''' Arguments:
            dataPath = string; path to Data directory
            fastqVar = string; fq vs fastq
            numProcs = integer; number of processors
            referencePath = string of the path to the Reference Directory
            genomeName = string; name of the genome file
            baseName = string; base name created
            projectPath = string; path to Project directory
        Returns:
            None; Runs Pipeline
    '''
    subprocess.call(["runall.sh",dataPath,fastqVar,str(numProcs),referencePath,genomeName,baseName,projectPath])

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
    os.chdir(postProcessingPath)
    makeJSON.writeJSON(jsonName)

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
    subprocess.call(["formatFeatures.sh",postProcessingPath,dataPath])
    
def makeTotalTime(postProcessingPath,dataPath):
    ''' Arguments:
            postProcessingPath = string; path to the Postprocessing Directory
            dataPath = string; path to Data Directory
        Returns:
            None; Creates totalTime.dat
    '''
    subprocess.call(["getTotalTime.sh",postProcessingPath,dataPath])

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
    os.chdir(postProcessingPath)
    makeCols.makeCols(makeCols.readJSON(jsonName),columnName)

def createCountFile(postProcessingPath,countName='Counts.dat'):
    ''' Arguments:
            postProcessingPath = string; path to the Postprocessing Directory
            countName = string; name of Counts file
        Returns:
            None; Creates Counts.dat
    '''
    subprocess.call(["makeCounts.sh",postProcessingPath,countName])

def createRScript(postProcessingPath,jsonName='Metadata.json',rName='makeReport.r'):
    ''' Arguments:
            postProcessingPath = string; path to the Postprocessing Directory
            dataPath = string; path to Data Directory
        Returns:
            None; Creates makeReportr.dat
    '''
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
    subprocess.call(["runDESeq.sh",postProcessingPath])

def notifyEnding(postProcessingPath):
    ''' Arguments:
        Returns:
    '''
    if getNumberofFiles(originalPath)%2 != 0:
        print('There are not an even number of files in {}!'.format(originalPath))
    else:
        numSamp = int(getNumberofFiles(originalPath)/2)

    os.chdir(projectPath + '/runPipeNotify')
    while True:
        if 'FINISHED.txt' in os.listdir(postProcessingPath):
            #textMe('+17756227884','This is sent from your python script; The first part of the \
            #        Pipeline has finished. Please return to computer')
            print('At this point textMe would send you a text')
            os.chdir(postProcessingPath)
            os.remove("FINISHED.txt")
            break
        else:
            time.sleep(30)


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
    runPipe(Project + '/Data', Fastq, Procs, Project + '/Reference', Genome, Basename,Project)
    createMetaData(Project + '/Postprocessing')
    findFinish(Project, Original)
    makeNiceCounts(Project + '/Postprocessing', Project + '/Data')
    makeTotalTime(Project + '/Postprocessing', Project + '/Data')
    createColumnFile(Project + '/Postprocessing')
    createCountFile(Project + '/Postprocessing')
    createRScript(Project + '/Postprocessing')
    makeRreports(Project + '/Postprocessing')
    print(':)')

if __name__ == '__main__':
    arguments = docopt(__doc__,version='Version 1.1\nAuthor: Alberto')
    print('In development :)')
    #Main(arguments['<pathtoInput>'])
