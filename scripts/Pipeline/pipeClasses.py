#!/usr/bin/python3

'''Usage: runPipe.py [-h | --help] [-j <jsonFile> | --jsonfile <jsonFile>]
                    [--noconfirm] [-c <placeToClean> | --clean <placeToClean>]
                    [-s <sampleName> | --sampleclean <sampleName>]
                    [--NUKE] [-r <runlogPath> | --runtime <runlogPath>] <pathtoInput>

Options:                                                                                           
    -h --help                                  :    Show this screen and exit
    -j <jsonFile>, --jsonfile <jsonFile>       :    Ignores JSON file creation and uses specified
                                                    path to JSON
    --noconfirm                                :    Ignore all user prompts except JSON file creation
    -c <placeToClean>, --clean <placeToClean>  :    Cleans <placeToClean>: Possible places include:
                                                    Reference, Data, Postprocessing, All
    -s <sampleName>, --sampleclean <sampleName>:    Similar to -c,--clean; but instead just cleans a
                                                    single sample directory <sampleName>
    --NUKE                                     :    Removes entire project Directory
    -r <runlogPath>, --runtime <runlogPath>    :    Optional directory path for Runtime Log file to 
                                                    be created in [default: $Project]
'''

################################################################
# Importations
################################################################

from docopt import docopt
from timeit import default_timer as timer
import multiprocessing
import time
import pipeUtils
import os
import subprocess
import shutil
import glob

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

#def funTimeN(logName):
#    def Timer(function):
#        def Wrapper(*args,**kwargs):
#            with open(logName, 'a') as R:
#                R.write('{} function started\n'.format(function.__name__))
#            t1 = timer()
#            stuff = function(*args)
#            t2 = timer()
#            with open(logName, 'a') as R:
#                R.write('{} function took {:.3f} seconds\n'.format(function.__name__, t2 - t1))
#                R.write('{} function finished\n'.format(function.__name__))
#            return stuff
#        return Wrapper
#    return Timer


def makeTimeFile(logPath):
    with open(logPath, 'w') as R:
        R.write('runPipe.py Runtime File\n-----------------------------------------\n\n')

def unwrap_self_runSample(arg, **kwarg):
    return Experiment.runSample(*arg, **kwarg)

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
        self.inputPath = str(inputPath)

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
    def createSampleClasses(self):
        Samples = [Experiment.Sample(n, self.inputPath) for n in range(1,self.getNumberofSamples() + 1)]
        return Samples

    @funTime
    def runSample(self, sample):
        sample.runParts()

    @funTime
    def GO(self):
        self.makeNotifyFolder()
        Samples = self.createSampleClasses()
        with multiprocessing.Pool(self.Procs) as p:
            p.map(unwrap_self_runSample, zip([self]*len(Samples), Samples))
 
    @funTime
    def deployPipe(self):
        print("Pipeline is running...")

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
            answer = input('Are you on a local X session or is X11 forwarding enabled? You wont be able to run DESeq without it...(y,n) ')
            if answer == 'y':
                print('Running...')
                break
            elif answer == 'n':
                print('Exiting now')
                raise SystemExit
            else:
                print('Please answer y or n')
            
    def makeNotifyFolder(self):
        if os.path.isdir(self.Project + '/runPipeNotify') == False:
            os.mkdir(self.Project + '/runPipeNotify')
        else:
            shutil.rmtree(self.Project + '/runPipeNotify')
            os.mkdir(self.Project + '/runPipeNotify')

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
    ################################################################

    ################################################################
    # Sample Class
    ################################################################

    class Sample:
        '''For a specific sample: run analyses'''

        def __init__(self,sampleNumber,inputFile):
            Experiment.__init__(self,inputFile)
            assert sampleNumber <= Experiment.getNumberofSamples(self)
            assert sampleNumber > 0
            self.sampleName = 'sample_{}'.format(str(sampleNumber).zfill(2))
            self.samplePath = '{}/{}'.format(self.Data, self.sampleName)
            self.Read1 = self.getReadNames()[0]
            self.Read2 = self.getReadNames()[1]
            self.logPath = '{}/Runtime.{}.log'.format(self.samplePath, self.sampleName)
        def __repr__(self):
            return 'Sample(%r)'%(self.sampleName)

        ########################################################
        # Utilities
        ########################################################

        def formatCommand(self,Command):
            ''' Arguments:
                    Command = string; a command that you wish to be funTimed.
                Returns:
                    correctCommand = string; the command to be subprocessed
                Example:
                    Command = r'fastqc -t 48 -o /home/alberton/Version1.1_test2/Data/sample_01/fastqc.sample_01 *.fastq.gz'
                    returns:
                        '{ time fastqc -t 48 -o /home/alberton/Version1.1_test2/Data/sample_01/fastqc.sample_01 *.fastq.gz; } >> /home/alberton/Version1.1_test2/Data/sample_01/Runtime.sample_01.log 2>&1'
            '''
            correctCommand = r'{{ time {0}; }} >> {1} 2>&1'.format(Command, self.logPath)
            return correctCommand

        def getHostname(self):
            hostFile = '{}/host.txt'.format(self.samplePath)
            command = ['hostname']
            with open(hostFile,'w') as host:
                subprocess.run(command, stdout=host,
                                        stderr=subprocess.STDOUT,
                                        check=True)

        def getReadNames(self):
            files = sorted([os.path.basename(thing) for thing in glob.glob(self.samplePath + '/*.gz')])
            return files

        def writeToLog(self, message):
            with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'a') as LOG:
                LOG.write(message)

        def writeFunctionHeader(self, function):
            self.writeToLog('\n{} started\n\n'.format(function))

        def writeFunctionTail(self, function):
            self.writeToLog('\n{} done\n\n'.format(function))

        def initializeLog(self):
            os.chdir(self.samplePath)
            with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'w') as LOG:
                LOG.write('\t\tRuntime Log for {}\n'.format(sampleName))
                LOG.write('----------------------------------------\n\n')

        ########################################################
        # Quality Control
        ########################################################

        def runFastqc(self,runNumber):
            fastqcFolder = '{}/fastqc{}.{}'.format(self.samplePath,
                                                int(runNumber),
                                                self.sampleName)
            if not os.path.exists(fastqcFolder):
                os.makedirs(fastqcFolder)
            command1 = r'fastqc -t {0} -o {1} {2} {3}'.format(self.Procs,
                                                                fastqcFolder,
                                                                self.Read1,
                                                                self.Read2)
            command2 = r'unzip \*.zip'.format(fastqcFolder)
            goodCommand1 = self.formatCommand(command1)
            goodCommand2 = self.formatCommand(command2)

            # Executing Commands
            os.chdir(self.samplePath)
            subprocess.run(goodCommand1,
                                shell=True,
                                check=True)
            os.chdir(fastqcFolder)
            subprocess.run(goodCommand2,
                                shell=True,
                                check=True)

        def getPhred(self):
            fastqcFolder = '{}/fastqc1.{}'.format(self.samplePath, self.sampleName)
            os.chdir(glob.glob('{}/*/'.format(fastqcFolder))[0])
            command = r"grep -E 'Encoding' fastqc_data.txt | sed 's/[^0-9.]*//g'"
            proc = subprocess.Popen(command,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT)
            illumVN = float(proc.communicate()[0])
            phred = 33
            if illumVN >= 1.3 and illumVN < 1.8:
                phred = 64
            else:
                phred = 33
            return phred

        def findOverrepSeq(self,runNumber):
            fastqcFolder = '{}/fastqc{}.{}'.format(self.samplePath, runNumber, self.sampleName)
            os.chdir(fastqcFolder)
            availableDirs = glob.glob('{}/*/'.format(fastqcFolder))
            for directory in availableDirs:
                os.chdir(directory)
                command = r"awk '/Overrepresented sequences/,/>>END_MODULE/' fastqc_data.txt | tail -n +2 | head -n -1 > Overrepresented_sequences.txt"
                subprocess.run(command,
                                shell=True,
                                check=True)
                with open('Overrepresented_sequences.txt','r') as Over:
                    G = Over.readlines()
                if len(G) != 0:
                    thing = 'some'
                else:
                    thing = 'no'
                with open(self.logPath,'a') as LOG:
                    LOG.write('There are {} overrepresented sequence in {}'.format(thing,
                                str('{}/{}'.format(fastqcFolder, directory))))

        def runQCheck(self, runNumber):
            self.runFastqc(runNumber)
            self.findOverrepSeq(runNumber)

        def runTrimmomatic(self):
            os.chdir(self.samplePath)
            Phred = self.getPhred()
            Reads = self.getReadNames()

            # Making Command
            command = r'java -jar $RNASEQDIR/Trimmomatic/trimmomatic-0.35.jar PE -threads {procs} -phred{phred} {Read1} {Read2} read1.P.trim.{fastq}.gz read1.U.trim.{fastq}.gz read2.P.trim.{fastq}.gz read2.U.trim.{fastq}.gz ILLUMINACLIP:$RNASEQDIR/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25'
            Context = {
                    "procs": self.Procs,
                    "phred": Phred,
                    "Read1": self.Read1,
                    "Read2": self.Read2,
                    "fastq": self.Fastq
                    }
            commandWithContext = command.format(**Context)
            goodCommand = self.formatCommand(commandWithContext)

            # Executing
            os.chdir(self.samplePath)
            subprocess.run(goodCommand,
                                shell=True,
                                check=True)

        def smartQC(self):
            os.chdir(self.samplePath)
            run = 1
            while True:
                STOP = False
                self.runFastqc(run)
                self.findOverrepSeq(run)
                while True:
                    print('Pipeline paused. You should review sample {}'.format(self.sampleName))
                    answer1 = input("Are you ready to run Trimmomatic on {}?(y,n) ".format(self.sampleName))
                    if answer1 == 'y':
                        self.runTrimmomatic()
                        run += 1
                        break
                    elif answer1 == 'n':
                        answer2 = input("Would you like to proceed to Stage 2?(y,n) ")
                        if answer2 == 'y':
                            break
                        elif answer2 == 'n':
                            answer3 = input("Would you like to run QC again?(y,n) ")
                        else:
                            print('Please answer y or n')
                    else:
                        print('Please answer y or n')
                if STOP == True:
                    break

        def runPart1(self):
            os.chdir(self.samplePath)
            with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'w') as LOG:
                LOG.write('\t\tRuntime Log Part 1 for {}\n'.format(self.sampleName))
                LOG.write('----------------------------------------\n\n')
            self.runQCheck(1)
            self.runTrimmomatic()
            self.runQCheck(2)
            self.writeFunctionTail('runPart1')

        ########################################################
        # Pipeline
        ########################################################

        def runSeqtk(self):
            self.writeFunctionHeader('runSeqtk')
            # Making Command
            command1 = r'seqtk sample -s100 read1.P.trim.{0}.gz 10000 | seqtk seq -A - > sampled.read1.fa'.format(
                                                        self.Fastq)
            command2 = r'seqtk sample -s100 read2.P.trim.{0}.gz 10000 | seqtk seq -A - > sampled.read2.fa'.format(
                                                        self.Fastq)
            goodCommand1 = self.formatCommand(command1)
            goodCommand2 = self.formatCommand(command2)
            # Executing
            os.chdir(self.samplePath)
            subprocess.run(goodCommand1,
                                shell=True,
                                check=True)
            subprocess.run(goodCommand2,
                                shell=True,
                                check=True)
            self.writeFunctionTail('runSeqtk')

        def runBlastn(self):
            self.writeFunctionHeader('runBlastn')
            # Making Command
            db = self.Reference + '/' + self.Basename + '.cdna.all'
            command1 = r'blastn -query sampled.read1.fa -db {0} -out sampled.read1_vscdna.out -task blastn-short -outfmt "6 std sstrand" -max_target_seqs 1 -num_threads {1}'.format(db, self.Procs)
            command2 = r'blastn -query sampled.read2.fa -db {0} -out sampled.read2_vscdna.out -task blastn-short -outfmt "6 std sstrand" -max_target_seqs 1 -num_threads {1}'.format(db, self.Procs)
            goodCommand1 = self.formatCommand(command1)
            goodCommand2 = self.formatCommand(command2)
            # Executing
            os.chdir(self.samplePath)
            subprocess.run(goodCommand1,
                                shell=True,
                                check=True)
            subprocess.run(goodCommand2,
                                shell=True,
                                check=True)
            self.writeFunctionTail('runBlastn')

        def findStranded(self):
            # Running stranded_classifier.py
            # Making Command
            with open('{}/Runtime.{}.log'.format(self.samplePath,
                                                self.sampleName), 'a') as LOG:
                LOG.write('findStranded started\n\n')
            command1 = r'stranded_classifier.py -1 sampled.read1_vscdna.out -2 sampled.read2_vscdna.out'
            goodCommand1 = self.formatCommand(command1)
            # Executing
            os.chdir(self.samplePath)
            subprocess.run(goodCommand1,
                                shell=True,
                                check=True)
            with open('{}/Runtime.{}.log'.format(self.samplePath,
                                                self.sampleName), 'a') as LOG:
                LOG.write('\nfindStranded done\n\n')
            # Scraping stranded_classifier.py output
            command2 = r"awk '/findStranded started/,/findStranded done/' Runtime.{}.log | grep -q False && stranded=0 || stranded=1 && echo $stranded".format(self.sampleName)
            # Executing
            os.chdir(self.samplePath)
            proc = subprocess.Popen(command2,
                                        shell=True,
                                        universal_newlines=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT)
            strandedBool = int(proc.communicate()[0])
            if strandedBool == 1:
                return True
            elif strandedBool == 0:
                return False
            else:
                print('There was an error with findStranded')
                raise SystemExit

        def runHisat(self):
            self.writeFunctionHeader('runHisat')
            # hisat2 manual line 553
            # Making Command
            command = r"""hisat2 -k 5 -p {numProcs} --dta --phred{phred} --known-splicesite-infile {ref}/splice_sites.txt -x {ref}/{basename} -1 read1.P.trim.{fastq}.gz -2 read2.P.trim.{fastq}.gz -S aligned.{sample}.sam"""
            Phred = self.getPhred()
            context = {
                    "numProcs": self.Procs,
                    "phred": str(Phred),
                    "ref": self.Reference,
                    "basename": self.Basename,
                    "fastq": self.Fastq,
                    "sample": self.sampleName
                    }
            
            goodCommand = self.formatCommand(command.format(**context)) 
            # Executing
            os.chdir(self.samplePath)
            subprocess.run(goodCommand,
                                shell=True,
                                check=True)
            self.writeFunctionTail('runHisat')

        def runCompression(self):
            self.writeFunctionHeader('runCompression')
            # Making Command
            command = r"samtools view -bT {ref}/{genome} -@{procs} aligned.{sample}.sam -o aligned.{sample}.bam"
            context = {
                    "ref": self.Reference,
                    "genome": self.Genome,
                    "procs": self.Procs,
                    "sample": self.sampleName,
                    }
            goodCommand = self.formatCommand(command.format(**context))
            # Executing
            os.chdir(self.samplePath)
            subprocess.run(goodCommand,
                                shell=True,
                                check=True)
            os.chdir(self.samplePath)
            os.remove('aligned.{}.sam'.format(self.sampleName))
            self.writeFunctionTail('runCompression')

        def runFeatureCounts(self):
            self.writeFunctionHeader('runFeatureCounts')
            # Making Command
            command = r"featureCounts -T {procs} -p -C --primary --ignoreDup -t exon -g gene_id -a {ref}/{gtf} -o aligned.{sample}.counts aligned.{sample}.bam"
            context = {
                    "procs": self.Procs,
                    "ref": self.Reference,
                    "gtf": self.Gtf,
                    "sample": self.sampleName,
                    }
            goodCommand = self.formatCommand(command.format(**context))
            # Executing
            os.chdir(self.samplePath)
            subprocess.run(goodCommand,
                                shell=True,
                                check=True)
            self.writeFunctionTail('runFeatureCounts')

        def getNiceColumns(self):
            self.writeFunctionHeader('getNiceColumns')
            # Making Command
            command = r'''tail -n +2 aligned.{sample}.counts | awk '{{printf ("%5s\t%s\t%s\n", $1, $6, $7)}}' > aligned.{sample}.counts.three'''
            context = {"sample": self.sampleName}
            goodCommand = self.formatCommand(command.format(**context))
            # Executing
            os.chdir(self.samplePath)
            subprocess.run(goodCommand,
                                shell=True,
                                check=True)
            self.writeFunctionTail('getNiceColumns')

        def getAlignedColumn(self):
            self.writeFunctionHeader('getAlignedColumn')
            # Making Command
            command = r'''tail -n +2 aligned.{sample}.counts | awk '{{print $7}}' > aligned.{sample}.counts.one'''
            context = {"sample": self.sampleName}
            goodCommand = self.formatCommand(command.format(**context))
            # Executing
            os.chdir(self.samplePath)
            subprocess.run(goodCommand,
                                shell=True,
                                check=True)
            self.writeFunctionTail('getAlignedColumn')

        def runPart2(self):
            '''
            Should be relatively smooth, meaning should not need
            user input
            '''
            os.chdir(self.samplePath)
            with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'a') as LOG:
                LOG.write('\t\tRuntime Log Part 2 for {}\n'.format(self.sampleName))
                LOG.write('----------------------------------------\n\n')
            self.runSeqtk()
            self.runBlastn()
            self.runHisat()
            self.runCompression()
            self.runFeatureCounts()
            self.getNiceColumns()
            self.getAlignedColumn()
            self.writeFunctionTail('runPart2')
            with open(self.Project + '/runPipeNotify/{}'.format('done'+self.sampleName), 'w') as N:
                N.write('{} is done'.format(self.samplePath))

        def runParts(self):
            self.runPart1()
            self.runPart2()

    ################################################################
    ################################################################

    ################################################################
    # Stage Run Functions
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
        self.GO()
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
    # End
    ################################################################

################################################################
# Handling Command line arguments
################################################################

if __name__ == '__main__':
    #global RUNTIMELOG
    #RUNTIMELOG='/home/alberton/testlog/Runtime.log'
    #E = Experiment('/home/alberton/idealINPUT')
    #E.runAll()
    ##with Pool(48) as p:
    ##    print(p.map(fa, [g for g in range(10)]))
    ##E.dostuff()
    ##pool = Pool(48)
    ##read = pool.map(fa, [1,2,3,4,5])
    ##pool.close()
    ##print(read)

    ##print('In development')
    ##raise SystemExit

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
    #PROJ.runAll()
    PROJ.runStage3()
    PROJ.runStage4()
    PROJ.runStage5()
    #
    ############################################################
    #

    t2 = timer()

    timeused = str(time.strftime('%H:%M:%S', time.gmtime(t2-t1)))
    print('Total time elapsed: {}'.format(timeused))
