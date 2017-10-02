#!/usr/bin/python3

################################################################
# Importations
################################################################

from docopt import docopt
from timeit import default_timer as timer
from math import ceil as ceiling
from functools import partial
import multiprocessing
import time
import os
import subprocess
import shutil
import glob
import json
import re
import makeBallgownScript
import makeSleuthScript
import optPath
import makeJSON
import makeCols
import makeReportr
import makeEdgeReport

################################################################
# Utilities
################################################################

def testFileExistence(filePath):
    if not os.path.exists(filePath):
        raise SystemExit("File does not exist: {}".format(filePath))

def testDirExistence(dirPath):
    if not os.path.isdir(dirPath):
        raise SystemExit("Directory does not exist: {}".format(dirPath))

def expandPath(path):
    """ Arguments:
            path : str; any path
        Returns:
            expandedPath : str; path with expanded variables and path
    """
    return os.path.abspath(os.path.expandvars(path))

def getFastq(originalPath):
    ''' Arguments:
            originalPath = string; path to ogOriginal
        Returns:
            fastq = str; either 'fastq' or 'fq'

        Scrapes Original Data to figure out whether
        fastq notation or fq notation is used
        ex: thingsample1.read1.fastq.gz
            fastq = 'fastq'
            thingsample1.read1.fq.gz
            fastq = 'fq'
    '''
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

def getReferenceVariables(referencePath):
    ''' Arguments:
            referencePath = string; path to Reference files
        Returns:
            gtf, cdna, genome filenames

        Scrapes Reference directory to determine name of GTF, cdna, and
        genome files
    '''
    ls = os.listdir(referencePath)
    if len(ls) != 3:
        raise SystemExit(
                "There are an incorrect number of files in:\n" +
                "{}\n".format(referencePath) + "Need only: a gtf file with " +
                "'gtf' in filename, a cdna with 'cdna' in filename', and a " +
                "genome file")

    for filename in ls:
        if "gtf" in filename.split("."):
            gtf = str(filename)
        elif "cdna" in filename:
            cdna = str(filename)
        else:
            genome = str(filename)
    return gtf,cdna,genome

def readInit(referencePath):
    ''' Arguments:
            referencePath = string; path to Reference files
        Returns:
            gtf, cdna, genome filenames

        Scrapes .init file inside of Reference directory to determine name of
        GTF, cdna, and genome files
    '''
    Init = os.path.join(referencePath, '.init')
    with open(Init,'r') as f:
        Stuff = f.readlines()
    Gtf = Stuff[0].rstrip('\n')
    Cdna = Stuff[1].rstrip('\n')
    Genome = Stuff[2].rstrip('\n')
    return Gtf, Cdna, Genome

def exportVariablesforClass(globalArgs):
    ''' Arguments:
            globalArgs = dict; contains manifest combined with CLI options
        Returns:
            Variables = dict; Experiment parameters

        Gathers instance variables for Experiment Class
    '''
    Project = expandPath(globalArgs["locations"]["Project"])
    Reference = expandPath(globalArgs["locations"]["Reference"])
    Original = expandPath(globalArgs["locations"]["Original"])
    testDirExistence(Reference)
    testDirExistence(Original)
    if IS_REFERENCE_PREPARED:
        Gtf, Cdna, Genome = readInit(Reference)
    elif os.path.exists(os.path.join(Reference, '.init')):
        Gtf, Cdna, Genome = readInit(Reference)
    else:
        Gtf, Cdna, Genome = getReferenceVariables(Reference)
        with open(os.path.join(Reference, '.init'),'w') as f:
            f.write('\n'.join([Gtf, Cdna, Genome]))
    Basename = str(Genome.split(".")[0])
    Fastq = getFastq(Original)
    maxCPU = globalArgs['--maxcpu']
    Procs = os.cpu_count() if maxCPU == None else maxCPU
    NumSamples = int(len(os.listdir(Original))/2)
    if len(os.listdir(Original))%2 != 0:
        raise SystemExit("There are not an even number of files in Original " +
                "directory meaning that there may be paired-end files missing")
    Variables = {
            "Projectpath": Project,
            "ogReference": Reference,
            "ogOriginal": Original,
            "Gtf": Gtf,
            "Cdna": Cdna,
            "Genome": Genome,
            "Basename": Basename,
            "Fastq": Fastq,
            "Procs": Procs,
            "NumSamples": NumSamples
            }
    return Variables

def checkJSON(jsonFile, behavior='default'):
    ''' Arguments:
            jsonFile = string; path to JSON file
            *behavior = 'default' or other str; non default also
                        checks whether or not jsonFile argument
                        was given
        Returns:
            None

        Tests if JSON file has correct syntax
    '''
    if behavior != 'default':
        assert jsonFile != None, 'Require a JSON file as an argument'
    with open(jsonFile) as JF:
        try:
            json.load(JF)
        except ValueError as e:
            raise SystemExit("Need --json argument with valid JSON file: {}".format(e))

def funTime(function):
    ''' Arguments:
            function = str; name of function
        Returns:
            None

        Wrapper to time a function and write into RUNTIMELOG
    '''
    def wrapper(*args):
        t1 = timer()
        stuff = function(*args)
        t2 = timer()
        with open(RUNTIMELOG, 'a') as R:
            R.write('{} function took {:.3f} seconds\n'.format(
                                    function.__name__, t2 - t1))
        return stuff
    return wrapper

def makeTimeFile(logPath):
    ''' Arguments:
            logPath = string; path to log file
        Returns:
            None

        Initializes log file
    '''
    with open(logPath, 'w') as R:
        R.write('runPipe Runtime File\n-----------------------------------------\n\n')

def unwrap_self_runSample_fc(arg, **kwarg):
    ''' Magic for multiprocessing
        i.e. Multiprocessing needs to occur in top level functions.
             -> multiprocessing won't work if called from a class method'''
    return FCountsExperiment.runSample(*arg, **kwarg)

def unwrap_self_runSample_st(arg, **kwarg):
    ''' Magic for multiprocessing
        i.e. Multiprocessing needs to occur in top level functions.
             -> multiprocessing won't work if called from a class method'''
    return StringtieExperiment.runSample(*arg, **kwarg)

def unwrap_self_runSample_ka(arg, **kwarg):
    ''' Magic for multiprocessing
        i.e. Multiprocessing needs to occur in top level functions.
             -> multiprocessing won't work if called from a class method'''
    return KallistoExperiment.runSample(*arg, **kwarg)

################################################################
# Defining Experiment Class
################################################################

class Experiment:
    ''' Experiment is a Paired-End RNA Sequencing Experiment.
        This class provides tools to analyze data
    '''

    def __init__(self, globalArgs):
        ''' Arguments:
                globalArgs = dict; contains manifest combined with CLI options
            Returns:
                None

            Initializes variables that get scraped from
            exportVariablesforClass(inputPath)
        '''
        variables = exportVariablesforClass(globalArgs)
        # 3 Location Variables
        self.Project = variables["Projectpath"]
        self.ogReference = variables["ogReference"]
        self.ogOriginal = variables["ogOriginal"]
        # 4 Project Directories
        self.Reference = os.path.join(variables["Projectpath"], 'Reference')
        self.Original = os.path.join(variables["Projectpath"], 'Original')
        self.Data = os.path.join(variables["Projectpath"], 'Data')
        self.Postprocessing = os.path.join(variables["Projectpath"], 'Postprocessing')
        # 3 Reference Files
        self.Gtf = str(variables["Gtf"])
        self.Cdna = str(variables["Cdna"])
        self.Genome = str(variables["Genome"])
        # Miscellaneous
        self.Basename = str(variables["Basename"])
        self.Fastq = str(variables["Fastq"])
        self.Procs = int(variables["Procs"])
        self.inputPath = str(expandPath(globalArgs['<input>']))
        self.Blacklist = str(expandPath(globalArgs['--use-blacklist']))
        self.Numsamples = int(variables["NumSamples"])
        # Global Args
        self.GlobalArgs = globalArgs

    def __repr__(self):
        ''' Arguments:
                None
            Returns:
                A string representation of class
                ex: >print(Experiment)
                        Experiment('/path/to/Project')

            Creates a string representation of class
        '''
        return 'Experiment(%r)'%(self.Project)

    ###############################################################
    # Utilities
    ###############################################################

    def checkMan(self, default, manifest):
        """ Arguments:
             default : 
            manifest : 
            Returns:
                None
        """
        return default if manifest == None else manifest

    def checkManBool(self, default, manifest):
        """ Arguments:
             default : 
            manifest : 
            Returns:
                None
        """
        return default if manifest == None else ""

    def isStructurePrepared(self):
        ''' Arguments:
                None
            Returns:
                A boolean {
                        True : if .init file exists and contains
                                an 'S' meaning the Project Structure
                                has been created
                        False : otherwise}

            Checks to see if directory structure has been made by checking
            if there is an 'S' in '.init' which is in the Project folder.
            The 'S' is added after createSymLinks has completed
        '''
        if os.path.exists(os.path.join(self.Project,'.init')):
            with open(os.path.join(self.Project, '.init'), 'r') as f:
                G = f.read()
            if 'S' in G:
                return True
        return False

    def isReferencePrepared(self):
        ''' Arguments:
                None
            Returns:
                A boolean {
                        True : if .init file exists and contains
                                an 'P' meaning Reference Data has been
                                processed
                        False : otherwise}

            Checks to see if reference data has been processed by checking
            if there is a 'P' in '.init' which is in the Project folder.
            The 'P' is added after preProcessingReference has
            completed
        '''

        if os.path.exists(os.path.join(self.Project, '.init')):
            with open(os.path.join(self.Project, '.init'), 'r') as f:
                G = f.read()
            print(G)
            if 'P' in G:
                return True
        return False

    def findPipeFinish(self):
        ''' Arguments:
                None
            Returns:
                None

            Figure out when self.GO()
            function has finished analysis on all samples
        '''
        doneGlob = os.path.join(self.Data, 'sample*/.done')
        notifyFolder = os.path.join(self.Project, 'runPipeNotify')
        while True:
            if os.path.isdir(notifyFolder):
                if (len(os.listdir(notifyFolder)) == self.Numsamples or 
                    len(glob.glob(doneGlob)) == self.Numsamples):
                    shutil.rmtree(notifyFolder)
                    break
                else:
                    pass
            elif len(glob.glob(doneGlob)) == self.Numsamples:
                break
            else:
                pass
            time.sleep(2)

    def is3Finished(self):
        ''' Arguments:
                None
            Returns:
                boolean = True if finished, False if not finished

            Figures out when self.GO()
            function has finished analysis on all samples
        '''
        doneGlob = os.path.join(self.Data, 'sample*/.done')
        notifyFolder = os.path.join(self.Project, 'runPipeNotify')
        finished = False
        if os.path.isdir(notifyFolder):
            if (len(os.listdir(notifyFolder)) == self.Numsamples or 
                len(glob.glob(doneGlob)) == self.Numsamples):
                finished = True
            else:
                finished = False
        elif len(glob.glob(doneGlob)) == self.Numsamples:
            finished = True
        else:
            finished = False
        return finished

    def getOptimal(self, cluster, jsonPath=None):
        ''' Arguments:
                cluster = list; a list of integers that describes number of
                            CPUs on each machine in a cluster
                *jsonPath = boolean; if there already exists a json with
                            desired path
            Returns:
                dict; optimal path for slurm batch file

            Figures out best way to run your samples
        '''
        Arguments = {
                "--maxcpu": str(self.Procs),
                "<numsamples>": str(self.Numsamples),
                "<cluster>": ','.join(str(node) for node in cluster),
                "--customize": False,
                "--tofile": "GeneratedOptimalPath.dat"
                }
        if jsonPath:
            with open(jsonPath) as JF:
                jsonData = json.load(JF,object_pairs_hook=makeCols.OrderedDict)
            return jsonData
        else:
            optPath.main(Arguments)
            with open(Arguments['--tofile'], "r") as JF:
                jsonData = json.load(JF,object_pairs_hook=makeCols.OrderedDict)
            return jsonData

    @funTime
    def createJsonMetadata(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls makeJSON.writeJSON() to make JSON file if not already
            provided
        '''
        print('Making Metadata file; require user input...')
        if os.path.exists(os.path.join(self.Postprocessing, 'Metadata.json')):
            pass
        elif JSFI == None:
            testDirExistence(self.Postprocessing)
            os.chdir(self.Postprocessing)
            makeJSON.writeJSON('Metadata.json')
        print('Waiting for Pipeline to finish...')

    def checkX11(self):
        ''' Arguments:
                None
            Returns:
                None

            Asks user if X server is available. Required for R stuff.
        '''
        while True:
            answer = input('Are you on a local X session or is X11 forwarding enabled? You' + 
                            ' wont be able to run DESeq without it...(y,n) ')
            if answer == 'y':
                print('Running...')
                break
            elif answer == 'n':
                print('Exiting now')
                raise SystemExit
            else:
                print('Please answer y or n')

    def makeNotifyFolder(self):
        ''' Arguments:
                None
            Returns:
                None

            Initializes notification folder that is used to determine when
            Pipeline has finished running
            'runPipeNotify/' in self.Project
        '''
        runPipeNotify = os.path.join(self.Project, 'runPipeNotify')
        if not os.path.isdir(runPipeNotify):
            os.mkdir(runPipeNotify)
        else:
            shutil.rmtree(runPipeNotify)
            os.mkdir(runPipeNotify)

    def makeNotifyFolder2(self):
        ''' Arguments:
                None
            Returns:
                None

            Initializes notification folder that is used to determine when
            Pipeline has finished running
            'runPipeNotify/' in self.Project
        '''
        runPipeNotify = os.path.join(self.Project, 'runPipeNotify')
        if not os.path.isdir(runPipeNotify):
            try:
                os.mkdir(runPipeNotify)
            except:
                pass

    def gatherAllSampleOverrep(self, runNumber):
        ''' Arguments:
                runNumber = int; the fastqc trial number that you want
                                to examine
            Returns:
                None

            Collects all Sample Overrepresented Sequences for all samples
            for a specific Fastqc trial. File will be in Postprocessing
        '''
        with open('{}/Run{}-OverrepSeq.txt'.format(self.Postprocessing, runNumber),'w') as R:
            R.write('Overrepresented Sequences for fastqc number {}\n'.format(runNumber))
            R.write('-------------------------------------------------------------\n\n')
        samples = glob.glob(self.Data + '/sa*')
        for sample in samples:
            OverSeq = glob.glob(sample + '/Sample*Run{}*'.format(runNumber))
            for ov in OverSeq:
                header = '-'.join(ov.split('/')[-3:-1])
                with open('{}/Run{}-OverrepSeq.txt'.format(self.Postprocessing,
                                                            runNumber),'a') as R:
                    R.write('\n' + header + '\n')
                    with open(ov,'r') as O:
                        contents = O.read()
                    R.write(contents + '\n')

    def checkEach(self):
        ''' Arguments:
                None
            Returns:
                None

            Checks if there exists a .done file in each sample folder.
            .done file is created at the completion of each respective
            sample
        '''
        Check = [os.path.exists(sample+'/.done')
                    for sample in glob.glob(os.path.join(self.Data,'sample*'))]
        if False in Check:
            return False
        else:
            return True

    def redirectSTDERR(self,command,logfile):
        ''' Arguments:
                command = string; a command that you wish to control stderr
                logfile = string; file for stderr and stdout to be saved to
            Returns:
                correctCommand = string; the command to be subprocessed
        '''
        correctCommand = r'{{ time -p {0}; }} >> {1} 2>&1'.format(command, logfile)
        return correctCommand

    def quickTrimStats(self):
        ''' Arguments:
                None
            Returns:
                None

            Scrapes trimmomatic stats from Runtime logs into one file
        '''
        command = """grep 'Input Read Pairs: ' {runtimelogglob} > {trimlog}"""
        Context = {
                "runtimelogglob": os.path.join(self.Data,
                                    'sample_*/Runtime.sample_*.log'),
                "trimlog": os.path.join(self.Postprocessing, 'trim.summary.log')
                }
        goodCommand = command.format(**Context)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True,
                            executable="/bin/bash")

    ###############################################################
    # Stage 1 and 2 Functions
    ###############################################################
    # @funTime is just a wrapper that records time of function
    @funTime
    def makeStructure(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls makeStructure.sh to make Directory Structure
        '''
        if not self.isStructurePrepared():
            print("Creating Structure...")
            command = r'''makeStructure.sh {} {}'''.format(self.Project,
                    self.Numsamples)
            subprocess.run(command,
                shell=True,
                check=True,
                executable="/bin/bash")
            makeTimeFile(RUNTIMELOG)

    @funTime
    def makeSyms(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls makeSyms.sh to make Symbolic Links for
            Original Data, Reference Data; and Metadata if available
        '''
        if not self.isStructurePrepared():
            command = r'''makeSyms.sh {} {} {} {}'''.format(self.Project,
                    self.ogOriginal, self.ogReference,
                    'false' if JSFI == None else JSFI)
            subprocess.run(command,
                shell=True,
                check=True,
                executable="/bin/bash")
            with open(os.path.join(self.Project,'.init'), 'w') as F:
                F.write('S')

    @funTime
    def qcRef(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls QCofRef to make
            Reference_Report.txt in Reference Folder
        '''
        if not self.isReferencePrepared():
            print("Running Quality Control Check on Reference Data...")
            command = r'''QCofRef.sh {} {}'''.format(self.Reference,
                    self.Genome)
            subprocess.run(command,
                shell=True,
                check=True,
                executable="/bin/bash")
            if not NOCONFIRM:
                with open(os.path.join(self.Reference, 'Reference_Report.txt'), 
                        'r') as Report:
                    print(Report.read())
                print('\nSee Reference_Report.txt to view initial Diagnostics')

    @funTime
    def ppRef(self):
        ''' Arguments:
                None
            Returns:
                None

            Pre-Process Reference Data using blast, hisat2, and samtools
        '''
        if not self.isReferencePrepared():
            print("Preprocessing Reference Data...")
            ppLog = os.path.join(self.Reference, 'Preprocessing.log')
            Context = {
                    "cdna": self.Cdna,
                    "basename": self.Basename,
                    "gtf": self.Gtf,
                    "genome": self.Genome,
                    "cpu": self.Procs
                    }
            makeBlastdb = """time -p makeblastdb -in {cdna} -dbtype nucl -out {basename}.cdna.all""".format(**Context)
            extractSpliceSites = """time -p extract_splice_sites.py {gtf} > splice_sites.txt""".format(**Context)
            extractExons = """time -p extract_exons.py {gtf} > known_exons.txt""".format(**Context)
            hisatBuild = """time -p hisat2-build -p {cpu} --ss splice_sites.txt --exon known_exons.txt {genome} {basename}""".format(**Context)
            samtoolsFaidx = """time -p samtools faidx {genome}""".format(**Context)
            os.chdir(self.Reference)
            with open(ppLog, 'w') as PPlog:
                PPlog.write('\n{}\n{}'.format(makeBlastdb,'='*50))
                subprocess.run(makeBlastdb,
                                    shell=True,
                                    check=True,
                                    executable="/bin/bash",
                                    stdout=PPlog,
                                    stderr=subprocess.STDOUT)
                PPlog.write('\n{}\n{}'.format(extractSpliceSites,'='*50))
                subprocess.run(extractSpliceSites,
                                    shell=True,
                                    check=True,
                                    executable="/bin/bash",
                                    stdout=PPlog,
                                    stderr=subprocess.STDOUT)
                PPlog.write('\n{}\n{}'.format(extractExons,'='*50))
                subprocess.run(extractExons,
                                    shell=True,
                                    check=True,
                                    executable="/bin/bash",
                                    stdout=PPlog,
                                    stderr=subprocess.STDOUT)
                PPlog.write('\n{}\n{}'.format(hisatBuild,'='*50))
                subprocess.run(hisatBuild,
                                    shell=True,
                                    check=True,
                                    executable="/bin/bash",
                                    stdout=PPlog,
                                    stderr=subprocess.STDOUT)
                PPlog.write('\n{}\n{}'.format(samtoolsFaidx,'='*50))
                subprocess.run(samtoolsFaidx,
                                    shell=True,
                                    check=True,
                                    executable="/bin/bash",
                                    stdout=PPlog,
                                    stderr=subprocess.STDOUT)
            with open(os.path.join(self.Project, '.init'), 'a') as F:
                F.write('P')

    ################################################################
    # Cleaning Functions
    ################################################################

    def clean(self, thingToClean, sampleName=None):
        ''' Arguments:
                thingToClean = string:
                                Reference or Data or Postprocessing or All or Sample
                Note: If thingToClean=Sample, require second argument:
                *sampleName = string; name of Sample to be cleaned
            Returns:
                None; cleans directories

            Cleaning function for Experiment.
            See runPipe --help for help
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
            if sampleName == None:
                raise SystemExit('Need a valid sample name')
            arg = sampleName
        else:
            raise SystemExit('Need a valid argument: ' +
                            'Reference, Data, Postprocessing, All, or Sample')
        if thingToClean == 'Sample':
            subprocess.run(['clean.sh','-s',arg, self.Genome, self.Cdna, self.Gtf, self.Reference,
                            self.Data, self.Postprocessing],check=True)
        else:
            subprocess.run(['clean.sh',arg, self.Genome, self.Cdna, self.Gtf, self.Reference,
                            self.Data, self.Postprocessing],check=True)

    ################################################################
    # Stage Run Functions
    ################################################################
    @funTime
    def runStage1(self):
        self.makeStructure()
        self.makeSyms()

    @funTime
    def runAll(self):
        '''
        Runs Stage 1 through Stage 5
        '''
        self.runStage1()
        self.runStage2()
        self.runStage3()
        self.runStage4()
        self.runStage5()

    ################################################################
    # End of Experiment class
    ################################################################

#@
class FCountsExperiment(Experiment):
    '''
        Inherit from general experiment class. Can
        apply featureCounts specific methods
    '''

    def __init__(self,globalArgs):
        Experiment.__init__(self,globalArgs)

    def __repr__(self):
        ''' Arguments:
                None
            Returns:
                A string representation of class
                ex: >print(Experiment)
                        Experiment('/path/to/Project')

            Creates a string representation of class
        '''
        return 'FCountsExperiment(%r)'%(self.Project)

    ###############################################################
    # Utilities
    ###############################################################

    def runSample(self, sample):
        ''' Arguments:
                sample =  class instance; a sample to execute
            Returns:
                None

            Executes a sample by calling its runParts() method
        '''
        if not sample.isFinished():
            sample.runParts()

    def GO(self, subject=0):
        ''' Arguments:
                subject = int; if is 0, then run all available samples;
                            if any other number, execute it individually
            Returns:
                None

            Execute Stage 3 for either a single sample or all samples
        '''
        if subject == 0:
            self.makeNotifyFolder()
            # Initialize all samples
            Samples = self.createAllSampleClasses()
            with multiprocessing.Pool(self.Procs) as p:
                p.map(unwrap_self_runSample_fc, zip([self]*len(Samples), Samples))
        else:
            name = 'sample_{:02g}'.format(subject)
            if os.path.exists(os.path.join(self.Data, name)):
                # Initialize only one sample
                experimentSample = self.createSampleClassNumber(subject)
                self.runSample(experimentSample)

    def makeBatch(self, cluster, jsonFile=None):
        ''' Arguments:
                cluster = list; list of integers that describe number
                            of CPUs in each node of your cluster
                *jsonFile = boolean; if optimal path already exists
                            in a file
            Returns:
                None

            Creates Pipeline batch script to be used with slurm
        '''
        batch =  """#!/bin/bash
#SBATCH --nodes={NODES}
#SBATCH --time=400
#SBATCH --cpus-per-task={CPT}
#SBATCH --ntasks={NTASKS}
#SBATCH --job-name="featureCounts Pipeline"
#SBATCH --export=PATH,RNASEQDIR,HOME,LC_ALL="en_US.UTF-8"

inputFile='{INPUT}'
jsonFile='{JSON}'

# Stage 1
{STAGE1}

wait

# Stage 2
{STAGE2}

wait

# Stage 3
{STAGE3}
wait

# Stage 4
{STAGE4}

wait

# Stage 5
{STAGE5}

wait
scontrol show job $SLURM_JOB_ID
wait
"""
        if JSFI == None:
            raise SystemExit('Need to specify a Metadata file with "--jsonfile"\nCan use "runPipe mj" to create a metadata file')
        numSamps = self.Numsamples
        if not IS_REFERENCE_PREPARED:
            ref = ''
        else:
            ref = ' --use-reference'
        command1 = 'srun -N1 -c1 -n1 runPipe fcounts --noconfirm{} --jsonfile "${{jsonFile}}" --execute 1 "${{inputFile}}"'.format(ref)
        command2 = 'srun -N1 -c{1} -n1 runPipe fcounts --noconfirm{0} --jsonfile "${{jsonFile}}" --maxcpu {1} --execute 2 "${{inputFile}}"'.format(ref,max(cluster))
        bestPath = self.getOptimal(cluster,jsonPath=jsonFile)
        command3,counter,sampleNum = '',1,1
        for path in sorted(bestPath):
            Pstep = bestPath[path]['Procs']
            Sstep = bestPath[path]['Samps']
            for S in range(sampleNum,sampleNum + Sstep):
                com = 'srun -N1 -c{0} -n1 --exclusive runPipe fcounts --noconfirm{2} --use-blacklist {3} --jsonfile "${{jsonFile}}" --maxcpu {0} -e 3 -r {1} "${{inputFile}}" &\n'.format(Pstep,S,ref,self.Blacklist)
                command3 += com
            if counter != len(bestPath):
                command3 += 'wait\n'
            counter += 1
            sampleNum += Sstep
        command4 = 'srun -N1 -c1 -n1 runPipe fcounts --noconfirm{0} --jsonfile "${{jsonFile}}" --execute 4 "${{inputFile}}"'.format(ref)
        command5a = 'srun -N1 -c1 -n1 --exclusive runPipe fcounts --noconfirm{0} --jsonfile "${{jsonFile}}" --execute 5 --edger "${{inputFile}}" &'.format(ref)
        command5b = 'srun -N1 -c1 -n1 --exclusive runPipe fcounts --noconfirm{0} --jsonfile "${{jsonFile}}" --execute 5 --deseq "${{inputFile}}" &'.format(ref)
        command5 = command5a + '\n' + command5b
        Context = {
                "NODES": len(cluster),
                "CPT": bestPath['Step 1']['Procs'],
                "NTASKS": bestPath['Step 1']['Samps'],
                "INPUT": self.inputPath,
                "JSON": JSFI,
                "STAGE1": command1,
                "STAGE2": command2,
                "STAGE3": command3,
                "STAGE4": command4,
                "STAGE5": command5,
                }
        batchScript = batch.format(**Context)
        with open('pipeBatch','w') as f:
            f.write(batchScript)

    def makeINTC(self):
        ''' Arguments:
                None
            Returns:
                None

            Scrapes GTF and writes ID,Name,BioType,and Chr
            for each gene into an INTC file
        '''
        geneID = r'gene_id [\S]*;'
        geneName = r'gene_name [\S]*;'
        geneType = r'gene_biotype [\S]*;'
        GI = re.compile(geneID)
        GN = re.compile(geneName)
        GT = re.compile(geneType)
        with open(os.path.join(self.Reference,self.Gtf),'r') as F:
            All = F.readlines()[5:]
        updatedGTF = []
        for line in All:
            try:
                ID = re.search(GI,line).group(0).split(' ')[1].split('"')[1]
            except AttributeError: # line doesn't have a field for name, id, or type
                ID = 'N/A'
            try:
                Name = re.search(GN,line).group(0).split(' ')[1].split('"')[1]
            except AttributeError: # line doesn't have a field for name, id, or type
                Name = 'N/A'
            try:
                BioType = re.search(GT,line).group(0).split(' ')[1].split('"')[1]
            except AttributeError: # line doesn't have a field for name, id, or type
                BioType = 'N/A'
            Chr = line.split('\t')[0]
            updatedGTF.append('\t'.join([ID,Name,BioType,Chr]))
        with open(os.path.join(self.Postprocessing,'INTC'),'w') as F:
            F.write('\n'.join(sorted(set(updatedGTF),key=lambda x: x.split('\t')[0])))

    def makeGoodCounts(self):
        ''' Arguments:
                None
            Returns:
                None

            Scrapes GTF and writes ID,Name,BioType,and Chr
            for each gene into an INTC file
        '''
        self.makeINTC()
        INTC = os.path.join(self.Postprocessing,'INTC')
        if os.path.exists(INTC):
            with open(INTC,'r') as I:
                leftSide = I.readlines()
            oldCounts = os.path.join(self.Postprocessing,'NiceCounts.dat')
            with open(oldCounts,'r') as F:
                unsortedData = F.readlines()
            Header = unsortedData[0]
            sortedData = sorted(unsortedData[1:],key=lambda x: x.split('\t')[0])
            newHeader = ['\t'.join(['Geneid','gene_name','gene_biotype','Chr'] +
                                   Header.split('\t'))]
            newHeader[0] = newHeader[0].strip()
            checkStuff = newHeader + ['\t'.join([a.strip(),b.strip()])
                                        for a,b in zip(leftSide,sortedData)]
            GoodStuff = []
            for line in checkStuff:
                splitLine = line.split('\t')
                if splitLine[0] != splitLine[4]:
                    print('GoodCounts.dat is off at line {}'.format(checkStuff.index(line)))
                else:
                    GoodStuff.append('\t'.join(splitLine[:4]+splitLine[5:]))
            newCounts = os.path.join(self.Postprocessing,'GoodCounts.dat')
            with open(newCounts,'w') as F:
                F.write('\n'.join(GoodStuff))
        else:
            print('{} does not exist'.format(INTC))

    @funTime
    def createNiceCounts(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls formatFeatures.sh to make the composite counts file
            with length column; Not used for my R analysis
        '''
        os.chdir(self.Postprocessing)
        command = r'''formatFeatures.sh {} {}'''.format(self.Postprocessing,
                self.Data)
        subprocess.run(command,
            shell=True,
            check=True,
            executable="/bin/bash")

    def createRCounts(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls makeCounts.sh to make the Count file(Counts.dat)
            that works in my R analysis
        '''
        os.chdir(self.Postprocessing)
        command = r'''makeCounts.sh {} {}'''.format(self.Postprocessing,
                'Counts.dat')
        subprocess.run(command,
            shell=True,
            check=True,
            executable="/bin/bash")

    def createRCols(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls makeCols.makeCols to make column file(Cols.dat)
            that describes Experiment. Uses Metadata JSON file for this
        '''
        os.chdir(self.Postprocessing)
        makeCols.makeCols(makeCols.readJSON('Metadata.json'), 'Cols.dat')

    def createRProgram(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls makeReportr.createRscript and makeEdgeReport.createEdgeR()
            to make DESeq2 R analysis program and edgeR analysis program
        '''
        makeReportr.createRscript(makeCols.readJSON('Metadata.json'),
                'makeReport.r')
        makeEdgeReport.createEdgeR(makeCols.readJSON('Metadata.json'),
                'makeEdge.r')

    @funTime
    def runRProgram(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls runDESeq() and runEdgeR()
            to create DESeq2 and edgeR reports
        '''
        print("R program is running...")
        self.runDESeq()
        self.runEdgeR()

    @funTime
    def runDESeq(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs DESeq2
        '''
        deseqCommand = r'''{ time Rscript "makeReport.r"; } > makeEdgeTime.log 2>&1'''
        subprocess.run(deseqCommand,
            shell=True,
            check=True,
            executable="/bin/bash")

    @funTime
    def runEdgeR(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs edgeR
        '''
        edgeCommand = r'''{ time Rscript "makeEdge.r"; } > makeEdgeTime.log 2>&1'''
        subprocess.run(edgeCommand,
            shell=True,
            check=True,
            executable="/bin/bash")

    def createAllSampleClasses(self):
        ''' Arguments:
                None
            Returns:
                experimentSamples = list; list contains Sample Classes

            Creates all Sample Classes and returns them in a list
        '''
        numSamples = self.Numsamples
        # Don't want to oversubscribe computer with multiprocessing
        sampCpuMax = self.Procs//numSamples if self.Procs//numSamples != 0 else 1
        experimentSamples = [FCountsSample(n, self.GlobalArgs)
                            for n in range(1,numSamples + 1)]
        return experimentSamples

    def createSampleClassNumber(self,subject):
        ''' Arguments:
                subject = int; the sample to create, default is to create all sample classes
            Returns:
                experimentSample = sample instance

            Creates specified sample class
        '''
        assert (type(subject) == int 
                and subject > 0 
                and subject <= self.Numsamples
                ), 'Need an int argument as a subject to create'
        experimentSample = FCountsSample(subject, self.GlobalArgs)
        return experimentSample

    def fullFCStats(self):
        ''' Arguments:
                None
            Returns:
                None

            Scrapes featureCounts summary logs into one file
        '''
        command = """paste {fcsummaryglob} | cut -sf {cols} > {summarylog}"""
        Context = {
                "fcsummaryglob": os.path.join(self.Data,
                                    'sample_*/aligned.sample_*.counts.summary'),
                "cols": ','.join(str(i) for i in [1]+[2*sampleNum for sampleNum in range(1, self.Numsamples+1)]),
                "summarylog": os.path.join(self.Postprocessing, 'fc.summary.log')
                }
        goodCommand = command.format(**Context)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True,
                            executable="/bin/bash")

    def quickFCStats(self):
        ''' Arguments:
                None
            Returns:
                None

            Scrapes quick featureCounts stats into one file
        '''
        command = """grep 'fragments :' {runtimelogglob} > {quicklog}"""
        fixFormat = """paste <(awk -F '| |' '{{print $1}}' {quicklog} | awk -F '/' '{{print $NF}}') <(awk -F ' ' '{{for (i=2; i<=NF; i++) printf $i " "; print $NF}}' {quicklog}) > {quicklog}"""
        Context = {
                "runtimelogglob": os.path.join(self.Data,
                                    'sample_*/Runtime.sample_*.log'),
                "quicklog": os.path.join(self.Postprocessing, 'fc.quickstats.log')
                }
        goodCommand = command.format(**Context)
        goodFixFormat = fixFormat.format(**Context)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True,
                            executable="/bin/bash")
        #subprocess.run(goodFixFormat,
        #                    shell=True,
        #                    check=True,
        #                    executable="/bin/bash")

    def fullH2Stats(self):
        ''' Arguments:
                None
            Returns:
                None

            Scrapes hisat2 stats from Runtime logs into one file
        '''
        command = """grep -A 14 'reads; of these:' {runtimelogglob} > {h2log}"""
        Context = {
                "runtimelogglob": os.path.join(self.Data,
                                    'sample_*/Runtime.sample_*.log'),
                "h2log": os.path.join(self.Postprocessing, 'h2.summary.log')
                }
        goodCommand = command.format(**Context)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True,
                            executable="/bin/bash")

    ################################################################
    # Stage Run Functions
    ################################################################

    def runStage2(self):
        if not IS_REFERENCE_PREPARED:
            self.qcRef()
            self.ppRef()
        else:
            with open(os.path.join(self.Project, '.init'), 'a') as F:
                F.write('P')

    def runStage3(self):
        print("Pipeline is running...")
        self.makeNotifyFolder()
        self.GO()
        self.findPipeFinish()

    def executeSample(self, number):
        self.makeNotifyFolder2() #TODO Fix makeNotifyFolder1 vs 2
        print("Pipeline is running for sample_{} on {}...".format(number,os.uname()[1]))
        self.GO(int(number))

    def runStage4(self):
        while True:
            if self.is3Finished() or self.checkEach():
                time.sleep(1)
                print('Preparing for R analysis...')
                if os.path.isdir(os.path.join(self.Project,'runPipeNotify')):
                    shutil.rmtree(os.path.join(self.Project,'runPipeNotify'))
                self.gatherAllSampleOverrep(1)
                self.gatherAllSampleOverrep(2)
                self.fullFCStats()
                self.quickFCStats()
                self.fullH2Stats()
                self.quickTrimStats()
                self.createJsonMetadata()
                self.createNiceCounts()
                self.createRCounts()
                self.createRCols()
                self.createRProgram()
                self.makeGoodCounts()
                break
            else:
                time.sleep(1)

    def runStage5(self):
        self.runRProgram()

#@@
class StringtieExperiment(Experiment):
    '''
        Inherit from general experiment class. Can
        apply Stringtie specific methods
    '''

    def __init__(self,globalArgs):
        Experiment.__init__(self,globalArgs)

    def __repr__(self):
        ''' Arguments:
                None
            Returns:
                A string representation of class
                ex: >print(Experiment)
                        Experiment('/path/to/Project')

            Creates a string representation of class
        '''
        return 'StringtieExperiment(%r)'%(self.Project)

    ###############################################################
    # Utilities
    ###############################################################

    def runSample(self, sample, stPhase):
        ''' Arguments:
                sample = class instance; a sample to execute
                stPhase = str; phase of Stringtie to run
            Returns:
                None

            Executes a sample by calling its runParts() method
        '''
        if not sample.isFinished():
            sample.runParts(stPhase)

    def GO(self, phases, subject=0):
        ''' Arguments:
                phases = list; list of phases given on command line
                subject = int; if is 0, then run all available samples;
                            if any other number, execute it individually
                            Note: phase 'b' must be ran with subject=0
            Returns:
                None

            Execute Stage 3 for either a single sample or all samples
        '''
        for phase in phases:
            nextPhase = self.runStringtiePhase(phases)
            if nextPhase == 'DONE':
                break
            else:
                if subject == 0:
                    self.makeNotifyFolder()
                    if nextPhase == 'b':
                        self.stringtiePart2b()
                    else:
                        Samples = self.createAllSampleClasses()
                        with multiprocessing.Pool(self.Procs) as p:
                            # Have to use partial in order to get optional argument into
                            # multiprocessing
                            p.map(partial(unwrap_self_runSample_st, stPhase=nextPhase),
                                    zip([self]*len(Samples), Samples))
                else:
                    if nextPhase == 'b':
                        raise SystemExit('Cannot specify a single sample with --stringtie b')
                    else:
                        name = 'sample_{:02g}'.format(subject)
                        if os.path.exists(os.path.join(self.Data,name)):
                            experimentSample = self.createSampleClassNumber(subject)
                            self.runSample(experimentSample,stPhase=nextPhase)

    def makeStringtieBatch(self, cluster, jsonFile=None):
        ''' Arguments:
                cluster = list; list of integers that describe number
                            of CPUs in each node of your cluster
                *jsonFile = boolean; if optimal path already exists
                            in a file
            Returns:
                None

            Creates Pipeline batch script to be used with slurm
        '''
        batch =  """#!/bin/bash
#SBATCH --nodes={NODES}
#SBATCH --time=400
#SBATCH --cpus-per-task={CPT}
#SBATCH --ntasks={NTASKS}
#SBATCH --job-name="Stringtie Pipeline"
#SBATCH --export=PATH,RNASEQDIR,HOME,LC_ALL="en_US.UTF-8"

inputFile='{INPUT}'
jsonFile='{JSON}'

# Stage 1
{STAGE1}
wait

# Stage 2
{STAGE2}
wait

# Stage 3
{STAGE3a}
wait

{STAGE3b}
wait

{STAGE3c}
wait

# Stage 4
{STAGE4}
kait

# Stage 5
{STAGE5}
wait

scontrol show job $SLURM_JOB_ID
wait
"""
        if JSFI == None:
            raise SystemExit('Need to specify a Metadata file with "--jsonfile"\nCan use "runPipe mj" to create a metadata file')
        numSamps = self.Numsamples
        if not IS_REFERENCE_PREPARED:
            ref = ''
        else:
            ref = ' --use-reference'
        command1 = 'srun -N1 -c1 -n1 runPipe string --noconfirm{} --jsonfile "${{jsonFile}}" --execute 1 "${{inputFile}}"'.format(ref)
        command2 = 'srun -N1 -c{1} -n1 runPipe string --noconfirm{0} --jsonfile "${{jsonFile}}" --maxcpu {1} --execute 2 "${{inputFile}}"'.format(ref,max(cluster))
        bestPath = self.getOptimal(cluster,jsonPath=jsonFile)
        def getStage3(phase):
            # For Phase a and c
            command3,counter,sampleNum = '',1,1
            for path in sorted(bestPath):
                Pstep = bestPath[path]['Procs']
                Sstep = bestPath[path]['Samps']
                for S in range(sampleNum,sampleNum + Sstep):
                    com = 'srun -N1 -c{0} -n1 --exclusive runPipe string --noconfirm{2} --use-blacklist {4} --jsonfile "${{jsonFile}}" --maxcpu {0} -e 3 -r {1} --phase {3} "${{inputFile}}" &\n'.format(Pstep,S,ref,phase,self.Blacklist)
                    command3 += com
                if counter != len(bestPath):
                    command3 += 'wait\n'
                counter += 1
                sampleNum += Sstep
            return command3
        command3b = 'srun -N1 -c{1} -n1 --exclusive runPipe string --noconfirm{0} --use-blacklist {2} --maxcpu {1} --jsonfile "${{jsonFile}}" --execute 3 --phase b "${{inputFile}}"'.format(ref, max(cluster), self.Blacklist)
        command4 = 'srun -N1 -c1 -n1 runPipe string --noconfirm{0} --jsonfile "${{jsonFile}}" --execute 4 "${{inputFile}}"'.format(ref)
        command5 = 'srun -N1 -c{1} -n1 --exclusive runPipe string --noconfirm{0} --maxcpu {1} --jsonfile "${{jsonFile}}" --execute 5 "${{inputFile}}"'.format(ref, max(cluster))
        Context = {
                "NODES": len(cluster),
                "CPT": bestPath['Step 1']['Procs'],
                "NTASKS": bestPath['Step 1']['Samps'],
                "INPUT": self.inputPath,
                "JSON": JSFI,
                "STAGE1": command1,
                "STAGE2": command2,
                "STAGE3a": getStage3('a'),
                "STAGE3b": command3b,
                "STAGE3c": getStage3('c'),
                "STAGE4": command4,
                "STAGE5": command5
                }
        batchScript = batch.format(**Context)
        with open('pipeBatch','w') as f:
            f.write(batchScript)

    def createAllSampleClasses(self):
        ''' Arguments:
                None
            Returns:
                experimentSamples = list; list contains Sample Classes

            Creates all Sample Classes and returns them in a list
        '''
        numSamples = self.Numsamples
        # Don't want to oversubscribe computer with multiprocessing
        sampCpuMax = self.Procs//numSamples if self.Procs//numSamples != 0 else 1
        experimentSamples = [StringtieSample(n, self.GlobalArgs)
                            for n in range(1,numSamples + 1)]
        return experimentSamples

    def createSampleClassNumber(self,subject):
        ''' Arguments:
                subject = int; the sample to create, default is to create all sample classes
            Returns:
                experimentSample = sample instance

            Creates specified sample class
        '''
        assert (type(subject) == int 
                and subject > 0 
                and subject <= self.Numsamples
                ), 'Need an int argument as a subject to create'
        experimentSample = StringtieSample(subject, self.GlobalArgs)
        return experimentSample

    def fullH2Stats(self):
        ''' Arguments:
                None
            Returns:
                None

            Scrapes hisat2 stats from Runtime logs into one file
        '''
        command = """grep -A 14 'reads; of these:' {runtimelogglob} > {h2log}"""
        Context = {
                "runtimelogglob": os.path.join(self.Data,
                                    'sample_*/Runtime.sample_*.log'),
                "h2log": os.path.join(self.Postprocessing, 'h2.summary.log')
                }
        goodCommand = command.format(**Context)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True,
                            executable="/bin/bash")

    ############################################################
    # Stringtie Utilities
    ############################################################

    def makeStringtieMergelist(self):
        ''' Arguments:
                None
            Returns:
                None

            Creates StringtieMerge directory and mergelist.txt
        '''
        stMergeDir = self.Postprocessing + '/StringtieMerge'
        if not os.path.isdir(stMergeDir):
            try:
                os.mkdir(stMergeDir)
            except:
                pass
        sampleNames = glob.glob(os.path.join(self.Data,'sample*'))
        mergeList = '\n'.join([os.path.join(sample,os.path.basename(sample)+'.st.gtf')
                                for sample in sampleNames])
        with open(os.path.join(stMergeDir,'mergelist.txt'),'w') as MergeList:
            MergeList.write(mergeList)

    def stringtieMerge(self):
        ''' Arguments:
                None
            Returns:
                None
        '''
        self.makeStringtieMergelist()
        stMergeDir = self.Postprocessing + '/StringtieMerge'
        mergeList = os.path.join(stMergeDir,'mergelist.txt')
        logFile = os.path.join(stMergeDir,'StringtieRuntime.log')
        # Making Command
        #command = r"stringtie --merge -p {procs} -G {ref}/{gtf} -o {mergedir}/{projectname}.stmerged.gtf {mergelist}"
        command = r"stringtie --merge -p {procs} -o {mergedir}/{projectname}.stmerged.gtf {mergelist}"
        context = {
                "procs": self.Procs,
                "ref": self.Reference,
                "gtf": self.Gtf,
                "mergedir": stMergeDir,
                "projectname": os.path.basename(self.Project),
                "mergelist": mergeList
                }
        goodCommand = self.redirectSTDERR(command.format(**context), logFile)
        # Executing
        os.chdir(stMergeDir)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True)

    def compareTranscripts(self):
        ''' Arguments:
                None
            Returns:
                None
        '''
        projectName = os.path.basename(self.Project)
        stMergeDir = self.Postprocessing + '/StringtieMerge'
        stMergedFile = os.path.join(stMergeDir,"{}.stmerged.gtf".format(projectName))
        logFile = os.path.join(stMergeDir,'GFFCompareRuntime.log')
        # Making Command
        command = r"gffcompare -r {ref}/{gtf} -G -o {projectname}.merged {stmerged}"
        context = {
                "ref": self.Reference,
                "gtf": self.Gtf,
                "projectname": projectName,
                "stmerged": stMergedFile
                }
        goodCommand = self.redirectSTDERR(command.format(**context), logFile)
        # Executing
        os.chdir(stMergeDir)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True)

    def allSampleLookup(self, fileTemplate):
        ''' Arguments:
                fileTemplate = str; with {} inside which represents
                                sampleName Returns:
                None
        '''
        boolTable = []
        for sample in glob.glob(os.path.join(self.Data,'*')):
            if os.path.exists(os.path.join(sample,fileTemplate.format(os.path.basename(sample)))):
                boolTable.append(True)
            else:
                boolTable.append(False)
        return all(boolTable)

    def checkStringtiePhases(self):
        ''' Arguments:
                None
            Returns:
                Status = str; string of phases completed(a,b,c)
        '''
        projectName = os.path.basename(self.Project)
        stMergeDir = self.Postprocessing + '/StringtieMerge'
        stMergedFile = os.path.join(stMergeDir,"{}.stmerged.gtf".format(projectName))

        Status = ' '
        if self.allSampleLookup("{}.st.gtf"):
            Status += 'a'
        if os.path.exists(stMergedFile):
            Status += 'b'
        if self.allSampleLookup("{}.good.st.gtf"):
            Status += 'c'
        return Status

    def nextStringtiePhase(self):
        ''' Arguments:
                None
            Returns:
                nextPhase = str; string of next phase that hasn't
                been completed yet i.e. "a", "b", or "c"
        '''
        Status = self.checkStringtiePhases()
        nextPhase = ''
        if Status[-1] == ' ':
            nextPhase = 'a'
        elif Status[-1] == 'a':
            nextPhase = 'b'
        elif Status[-1] == 'b':
            nextPhase = 'c'
        elif Status[-1] == 'c':
            nextPhase = 'DONE'
        return nextPhase

    def runStringtiePhase(self, stringtiePhases, behavior='default'):
        ''' Arguments:
                None
            Returns:
                None

            Compares stringtie phases needing to be ran with phases
            given on command line
        '''
        runPhase = ''
        nextPhase = self.nextStringtiePhase()
        if nextPhase == 'DONE':
            runPhase = nextPhase
        elif len(stringtiePhases) == 1:
            if stringtiePhases == nextPhase:
                runPhase = nextPhase
            else:
                if behavior == 'default':
                    raise SystemExit('Cannot run --stringtie {}\nPlease run --stringtie {} first for all samples'.format(stringtiePhases,nextPhase))
        elif nextPhase in stringtiePhases:
            runPhase = nextPhase
        else:
            runPhase = 'DONE'
        return runPhase

    @funTime
    def stringtiePart2b(self):
        ''' Arguments:
                None
            Returns:
                None

            To be run in between sample stringtie part 2a and 2c.
            Requires results from phase a for each sample.
        '''
        self.makeStringtieMergelist()
        self.stringtieMerge()
        self.compareTranscripts()

    ########################################################
    # Gathering Data (for stringtie)
    ########################################################

    def organizeStringtieOutput(self):
        ''' Arguments:
                None
            Returns:
                None

            Gathers stringtie results into Postprocessing
        '''
        resultsDirectory = os.path.join(self.Postprocessing, 'StringtieResults')
        if not os.path.isdir(resultsDirectory):
            try:
                os.mkdir(resultsDirectory)
            except:
                pass
        for sample in glob.glob(os.path.join(self.Data,'sample*')):
            sampleResults = os.path.join(resultsDirectory, os.path.basename(sample))
            if not os.path.isdir(sampleResults):
                try:
                    os.mkdir(sampleResults)
                except:
                    pass
            for ctab in glob.glob(os.path.join(sample, '*ctab')):
                try:
                    os.symlink(ctab, os.path.join(sampleResults, os.path.basename(ctab)))
                except FileExistsError:
                    print('Symbolic link failed since' + 
                            ' {} already exists in StringtieResults'.format(ctab))
            goodSampleGtf = os.path.join(sample, '{}.good.st.gtf'.format(
                                                os.path.basename(sample)))
            if os.path.exists(goodSampleGtf):
                try:
                    os.symlink(goodSampleGtf,
                                os.path.join(sampleResults, os.path.basename(goodSampleGtf)))
                except FileExistsError:
                    print('Symbolic link failed since' + 
                            ' {} already exists in StringtieResults'.format(ctab))

    def createBallgownCols(self, jsonName='Metadata.json', columnName='Cols.dat'):
        ''' Arguments:
                None
            Returns:
                None

            Creates table for ballgown analysis
        '''
        jsonFile = os.path.join(self.Postprocessing, jsonName)
        colFile = os.path.join(self.Postprocessing, columnName)
        makeCols.makeStringtieCols(makeCols.readJSON(jsonFile), colFile)

    def createBallgownScript(self, jsonName='Metadata.json', programName='runBallgown.r'):
        ''' Arguments:
                None
            Returns:
                None

            Creates Ballgown R script
        '''
        jsonFile = os.path.join(self.Postprocessing, jsonName)
        programFile = os.path.join(self.Postprocessing, programName)
        makeBallgownScript.createRBallgownScript(makeCols.readJSON(jsonFile),
                                                 programFile)

    @funTime
    def runBallgownAnalysis(self):
        ''' Arguments:
                None
            Returns:
                None

            Executes Ballgown R script
        '''
        os.chdir(self.Postprocessing)
        ballgownCommand = r'''{ time -p Rscript "runBallgown.r"; } > runBallgownTime.log 2>&1'''
        subprocess.run(ballgownCommand,
                            shell=True,
                            check=True)

    ################################################################
    # Stage Run Functions
    ################################################################

    def runStage2(self):
        if not IS_REFERENCE_PREPARED:
            self.qcRef()
            self.ppRef()
        else:
            with open(os.path.join(self.Project, '.init'), 'a') as F:
                F.write('P')

    def runStage3(self, phases):
        print("Pipeline is running...")
        self.makeNotifyFolder()
        self.GO(phases)
        stringtieStatus = self.runStringtiePhase(phases, behavior='non-default')
        if stringtieStatus == 'DONE' or stringtieStatus == 'c':
            self.findPipeFinish()

    def executeSample(self, number, phases):
        ''' Arguments:
                number = int; sample number
                phases = list; list of Stringtie phases to be executed
                               given on command line
            Returns:
                None

        Run Stage 3 for Sample #number
        '''
        self.makeNotifyFolder2()
        print("Pipeline is running for sample_{} on {}...".format(number,os.uname()[1]))
        self.GO(phases, int(number))

    def runStage4(self):
        while True:
            if self.is3Finished() or self.checkEach():
                time.sleep(1)
                print('Preparing for R analysis...')
                if os.path.isdir(os.path.join(self.Project,'runPipeNotify')):
                    shutil.rmtree(os.path.join(self.Project,'runPipeNotify'))
                self.gatherAllSampleOverrep(1)
                self.gatherAllSampleOverrep(2)
                self.fullH2Stats()
                self.quickTrimStats()
                self.createJsonMetadata()
                self.organizeStringtieOutput()
                self.createBallgownCols()
                self.createBallgownScript()
                break
            else:
                time.sleep(1)

    def runStage5(self):
        self.runBallgownAnalysis()

#@@@
class KallistoExperiment(Experiment):
    '''
        Inherit from general experiment class. Can
        apply Kallisto specific methods
    '''

    def __init__(self,globalArgs):
        Experiment.__init__(self,globalArgs)

    def __repr__(self):
        ''' Arguments:
                None
            Returns:
                A string representation of class
                ex: >print(Experiment)
                        Experiment('/path/to/Project')

            Creates a string representation of class
        '''
        return 'StringtieExperiment(%r)'%(self.Project)

    ###############################################################
    # Utilities
    ###############################################################

    def runSample(self, sample):
        ''' Arguments:
                sample =  class instance; a sample to execute
            Returns:
                None

            Executes a sample by calling its runParts() method
        '''
        if not sample.isFinished():
            sample.runParts()

    def GO(self, subject=0):
        ''' Arguments:
                subject = int; if is 0, then run all available samples;
                            if any other number, execute it individually
            Returns:
                None

            Execute Stage 3 for either a single sample or all samples
        '''
        if subject == 0:
            self.makeNotifyFolder()
            Samples = self.createAllSampleClasses()
            with multiprocessing.Pool(self.Procs) as p:
                p.map(unwrap_self_runSample_ka,
                        zip([self]*len(Samples), Samples))
        else:
            name = 'sample_{:02g}'.format(subject)
            if os.path.exists(self.Data + '/' + name):
                experimentSample = self.createSampleClassNumber(subject)
                self.runSample(experimentSample)

    def makeKallistoBatch(self, cluster, jsonFile=None):
        ''' Arguments:
                cluster = list; list of integers that describe number
                            of CPUs in each node of your cluster
                *jsonFile = boolean; if optimal path already exists
                            in a file
            Returns:
                None

            Creates Pipeline batch script to be used with slurm
        '''
        batch =  """#!/bin/bash
#SBATCH --nodes={NODES}
#SBATCH --time=400
#SBATCH --cpus-per-task={CPT}
#SBATCH --ntasks={NTASKS}
#SBATCH --job-name="Kallisto Pipeline"
#SBATCH --export=PATH,RNASEQDIR,HOME,LC_ALL="en_US.UTF-8"

inputFile='{INPUT}'
jsonFile='{JSON}'

# Stage 1
{STAGE1}
wait

# Stage 2
{STAGE2}
wait

# Stage 3
{STAGE3}
wait

# Stage 4
{STAGE4}
wait

# Stage 5
{STAGE5}
wait

scontrol show job $SLURM_JOB_ID
wait
"""
        if JSFI == None:
            raise SystemExit('Need to specify a Metadata file with "--jsonfile"\nCan use "runPipe mj" to create a metadata file')
        numSamps = self.Numsamples
        if not IS_REFERENCE_PREPARED:
            ref = ''
        else:
            ref = ' --use-reference'
        command1 = 'srun -N1 -c1 -n1 runPipe kall --noconfirm{} --jsonfile "${{jsonFile}}" --execute 1 "${{inputFile}}"'.format(ref)
        command2 = 'srun -N1 -c{1} -n1 runPipe kall --noconfirm{0} --jsonfile "${{jsonFile}}" --maxcpu {1} --execute 2 "${{inputFile}}"'.format(ref,max(cluster))
        bestPath = self.getOptimal(cluster,jsonPath=jsonFile)
        def getStage3():
            command3,counter,sampleNum = '',1,1
            for path in sorted(bestPath):
                Pstep = bestPath[path]['Procs']
                Sstep = bestPath[path]['Samps']
                for S in range(sampleNum,sampleNum + Sstep):
                    com = 'srun -N1 -c{0} -n1 --exclusive runPipe kall --noconfirm{2} --use-blacklist {3} --jsonfile "${{jsonFile}}" --maxcpu {0} -e 3 -r {1} "${{inputFile}}" &\n'.format(Pstep,S,ref,self.Blacklist)
                    command3 += com
                if counter != len(bestPath):
                    command3 += 'wait\n'
                counter += 1
                sampleNum += Sstep
            return command3
        command4 = 'srun -N1 -c1 -n1 runPipe kall --noconfirm{0} --jsonfile "${{jsonFile}}" --execute 4 "${{inputFile}}"'.format(ref)
        command5 = 'srun -N1 -c{1} -n1 --exclusive runPipe kall --noconfirm{0} --maxcpu {1} --jsonfile "${{jsonFile}}" --execute 5 "${{inputFile}}"'.format(ref, max(cluster))
        Context = {
                "NODES": len(cluster),
                "CPT": bestPath['Step 1']['Procs'],
                "NTASKS": bestPath['Step 1']['Samps'],
                "INPUT": self.inputPath,
                "JSON": JSFI,
                "STAGE1": command1,
                "STAGE2": command2,
                "STAGE3": getStage3(),
                "STAGE4": command4,
                "STAGE5": command5
                }
        batchScript = batch.format(**Context)
        with open('pipeBatch','w') as f:
            f.write(batchScript)

    def createAllSampleClasses(self):
        ''' Arguments:
                None
            Returns:
                Samples = list; list contains Sample Classes

            Creates all Sample Classes and returns them in a list
        '''
        numSamples = self.Numsamples
        # Don't want to oversubscribe computer with multiprocessing
        sampCpuMax = self.Procs//numSamples if self.Procs//numSamples != 0 else 1
        experimentSamples = [KallistoSample(n, self.GlobalArgs)
                            for n in range(1,numSamples + 1)]
        return experimentSamples

    def createSampleClassNumber(self,subject):
        ''' Arguments:
                subject = int; the sample to create, default is to create all sample classes
            Returns:
                experimentSample = sample instance

            Creates specified sample class
        '''
        assert (type(subject) == int 
                and subject > 0 
                and subject <= self.Numsamples
                ), 'Need an int argument as a subject to create'
        experimentSample = KallistoSample(subject, self.GlobalArgs)
        return experimentSample

    ################################################################
    # Kallisto Utilities
    ################################################################

    def needToBuildKaliIndex(self):
        ''' Arguments:
                None
            Returns:
                None

            KaliIndexBuilt is a file that gets written right after 
            Kallisto index gets built
        '''
        if os.path.exists(os.path.join(self.Reference, 'KaliIndexBuilt')):
            return False
        else:
            return True

    @funTime
    def buildKallistoIndex(self):
        ''' Arguments:
                None
            Returns:
                None

            Build Kallisto Index and save it in Reference directory
        '''
        # ? What is behavior of kallisto index, will it build failed index? Assume no
        if self.needToBuildKaliIndex():
            logFile = os.path.join(self.Reference, 'KallistoRuntime.log')
            # Making Command
            command = r"kallisto index -i {basename}.kali.cdna.fa.idx {cdna}"
            context = {
                    "cdna": self.Cdna,
                    "basename": self.Basename
                    }
            goodCommand = self.redirectSTDERR(command.format(**context), logFile)
            # Executing
            print('Building kallisto index...')
            os.chdir(self.Reference)
            subprocess.run(goodCommand,
                                shell=True,
                                check=True)
            with open(os.path.join(self.Reference, 'KaliIndexBuilt'),'w') as F:
                F.write('True')

    def organizeKallistoOutput(self):
        ''' Arguments:
                None
            Returns:
                None

            Gather results of kallisto analysis for each sample into
            Postprocessing
        '''
        resultsDirectory = os.path.join(self.Postprocessing, 'KallistoResults')
        if not os.path.isdir(resultsDirectory):
            try:
                os.mkdir(resultsDirectory)
            except:
                pass
        for sample in glob.glob(os.path.join(self.Data,'sample*')):
            sampleResults = os.path.join(resultsDirectory, os.path.basename(sample))
            kallistoOutput = os.path.join(sample, "KallistoOutput.{}".format(os.path.basename(sample)))
            linkCommand = r"ln -sr {} {}".format(kallistoOutput, sampleResults)
            if not os.path.isdir(sampleResults):
                try:
                    subprocess.run(linkCommand,
                                        shell=True,
                                        check=True)
                except:
                    pass

    def createSleuthCols(self, jsonName='Metadata.json', columnName='Cols.dat'):
        ''' Arguments:
                *jsonName = str; name of Metadata file if it exists
                *columnName = str; name of Sleuth table to create
            Returns:
                None

            Creates table for Sleuth R analysis
        '''
        jsonFile = os.path.join(self.Postprocessing, jsonName)
        colFile = os.path.join(self.Postprocessing, columnName)
        makeCols.makeKallistoCols(makeCols.readJSON(jsonFile), colFile)

    def createSleuthScript(self, jsonName='Metadata.json', programName='runSleuth.r'):
        ''' Arguments:
                *jsonName = str; name of Metadata file if it exists
                *programName = str; name of Sleuth R script to create
            Returns:
                None

            Creates Sleuth R script
        '''
        jsonFile = os.path.join(self.Postprocessing, jsonName)
        programFile = os.path.join(self.Postprocessing, programName)
        makeSleuthScript.createRSleuthScript(makeCols.readJSON(jsonFile), programFile) 

    @funTime
    def runSleuthAnalysis(self):
        ''' Arguments:
                None
            Returns:
                None

            Execute Sleuth script
        '''
        os.chdir(self.Postprocessing)
        sleuthCommand = r'''{ time -p Rscript "runSleuth.r"; } > runSleuthTime.log 2>&1'''
        subprocess.run(sleuthCommand,
                            shell=True,
                            check=True,
                            executable="/bin/bash")

    ################################################################
    # Stage Run Functions
    ################################################################

    def runStage2(self):
    #TODO Check if statement
        if self.needToBuildKaliIndex():
            self.buildKallistoIndex()
        if not IS_REFERENCE_PREPARED:
            self.qcRef()
            self.ppRef()
        else:
            with open(os.path.join(self.Project, '.init'), 'a') as F:
                F.write('P')

    def runStage3(self):
        print("Pipeline is running...")
        self.makeNotifyFolder()
        self.GO()
        self.findPipeFinish()

    def executeSample(self, number):
        ''' Arguments:
                number = int; sample number
            Returns:
                None

        Run Stage 3 for Sample #number
        '''
        self.makeNotifyFolder2()
        print("Pipeline is running for sample_{} on {}...".format(number,os.uname()[1]))
        self.GO(int(number))

    def runStage4(self):
        while True:
            if self.is3Finished() or self.checkEach():
                time.sleep(1)
                print('Preparing for R analysis...')
                if os.path.isdir(os.path.join(self.Project,'runPipeNotify')):
                    shutil.rmtree(os.path.join(self.Project,'runPipeNotify'))
                self.gatherAllSampleOverrep(1)
                self.gatherAllSampleOverrep(2)
                self.quickTrimStats()
                self.createJsonMetadata()
                self.organizeKallistoOutput()
                self.createSleuthCols()
                self.createSleuthScript()
                break
            else:
                time.sleep(1)

    def runStage5(self):
        self.runSleuthAnalysis()


################################################################
# Defining Sample Class
################################################################
class Sample:
    '''
    For a specific sample: run analyses
    Inherits from Experiment Parent Class
    '''

    def __init__(self,sampleNumber,globalArgs):
        ''' Arguments:
                sampleNumber = int; sample number used for naming
                globalArgs = dict; contains options from cli and manifest
            Returns:
                None

            Initializes class variables from Parent class and also
            initializes some sample specific variables
        '''
        Experiment.__init__(self,globalArgs)
        assert sampleNumber <= Experiment.Numsamples and sampleNumber > 0
        self.sampleNumber = sampleNumber
        self.sampleName = 'sample_{}'.format(str(sampleNumber).zfill(2))
        self.samplePath = '{}/{}'.format(self.Data, self.sampleName)
        self.Read1 = self.getReadNames()[0]
        self.Read2 = self.getReadNames()[1]
        self.logPath = '{}/Runtime.{}.log'.format(self.samplePath, self.sampleName)
        if globalArgs['--maxcpu'] != None:
            self.Procs = globalArgs['--maxcpu']

    ########################################################
    # Utilities
    ########################################################

    def formatCommand(self,Command):
        ''' Arguments:
                Command = string; a command that you wish to be funTimed.
            Returns:
                correctCommand = string; the command to be subprocessed
            Example:
                Command = r'fastqc -t 48 \
                -o /home/alberton/Version1.1_test2/Data/sample_01/fastqc.sample_01 \
                *.fastq.gz'
                returns:
                '{ time fastqc -t 48 \
                -o /home/alberton/Version1.1_test2/Data/sample_01/fastqc.sample_01 *.fastq.gz; } \
                >> /home/alberton/Version1.1_test2/Data/sample_01/Runtime.sample_01.log 2>&1'
        '''
        correctCommand = r'{{ time -p {0}; }} >> {1} 2>&1'.format(Command, self.logPath)
        return correctCommand

    def getHostname(self):
        ''' Arguments:
                None
            Returns:
                None

            Writes to host.txt the hostname
        '''
        hostFile = '{}/host.txt'.format(self.samplePath)
        command = ['hostname']
        with open(hostFile,'w') as host:
            subprocess.run(command, stdout=host,
                                    stderr=subprocess.STDOUT,
                                    check=True)

    def getReadNames(self):
        ''' Arguments:
                None
            Returns:
                files = list; contains read1 and read2 original
                        data for sample

            Finds the name of the original data names and returns
            them in a list
        '''
        files = sorted([os.path.basename(thing) for thing in glob.glob(self.samplePath + '/*.gz')])
        return files

    def writeToLog(self, message):
        ''' Arguments:
                message = str; a line to be written to specific sample
                            Runtime.log
            Returns:
                None

            Writes to sample log a string(message)
        '''
        with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'a') as LOG:
            LOG.write(message)

    def writeFunctionHeader(self, function):
        ''' Arguments:
                function = str; name of function
            Returns:
                None

            Calls self.writeToLog with name of function
            Used to identify specific function output and times
        '''
        self.writeToLog('\n{} started\n\n'.format(function))

    def writeFunctionTail(self, function):
        ''' Arguments:
                function = str; name of function
            Returns:
                None

            Calls self.writeToLog with name of function
            Used to identify specific function output and times
        '''
        self.writeToLog('\n{} done\n\n'.format(function))

    def writeFunctionCommand(self, command):
        ''' Arguments:
                command= str; command used for execution
            Returns:
                None

            Calls self.writeToLog with name of function
            Used to identify specific function output and times
        '''
        self.writeToLog('\nCommand Used:\n{}\n'.format(command))

    def initializeLog(self):
        ''' Arguments:
                None
            Returns:
                None

            Initializes sample log
        '''
        os.chdir(self.samplePath)
        with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'w') as LOG:
            LOG.write('\t\tRuntime Log for {}\n'.format(sampleName))
            LOG.write('----------------------------------------\n\n')

    def isFinished(self):
        ''' Arguments:
                None
            Returns:
                bool: True = sample has finished running
                      False = sample not finished

            Checks to see if '.done' exists in sample directory
        '''
        if os.path.exists(os.path.join(self.samplePath, '.done')):
            return True
        return False

    ########################################################
    # Quality Control
    ########################################################

    def runFastqc(self,runNumber,read1,read2):
        ''' Arguments:
                runNumber = int; trial number for QC
                read1 = str; first read to run fastQC on
                read2 = str; second read to run fastQC on
            Returns:
                None

            Runs fastQC on read1 and read2
        '''
        fastqcFolder = '{}/fastqc{}.{}'.format(self.samplePath,
                                            int(runNumber),
                                            self.sampleName)
        if not os.path.exists(fastqcFolder):
            os.makedirs(fastqcFolder)
        localArgs = self.GlobalArgs['runQCheck']
        command1 = r'''{fastqc} -t {-t} -o {-o} {other} {read1} {read2}'''
        c1Context = {
                "fastqc": self.checkMan("fastqc", localArgs['fastqc']),
                "-t": self.checkMan(self.Procs, localArgs['-t']),
                "-o": self.checkMan(fastqcFolder, localArgs['-o']),
                "other": self.checkMan('', localArgs['other']),
                "read1": self.checkMan(read1, localArgs['read1']),
                "read2": self.checkMan(read2, localArgs['read2']),
                }
        command1.format(**c1Context)
        command2 = r'unzip \*.zip'.format(fastqcFolder)
        goodCommand1 = self.formatCommand(command1)
        goodCommand2 = self.formatCommand(command2)
        # Executing Commands
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand1)
        subprocess.run(goodCommand1,
                            shell=True,
                            check=True,
                            executable='/bin/bash')
        os.chdir(fastqcFolder)
        self.writeFunctionCommand(goodCommand2)
        subprocess.run(goodCommand2,
                            shell=True,
                            check=True,
                            executable='/bin/bash')

    def getPhred(self):
        ''' Arguments:
                None
            Returns:
                phred = int; Phred Version Number

            Scrapes Phred Version Number from fastqc output
        '''
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
        ''' Arguments:
                runNumber = int; trial number for QC
            Returns:
                None

            Scrapes Overrepresented Sequences from Fastqc output
            and writes to a file in sample fastqc directory
        '''
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

    def initOverrepLog(self, runNumber):
        ''' Arguments:
                runNumber = int; trial number for QC
            Returns:
                None

            Initializes Log for Overepresented Sequences
        '''
        with open('{}/Sample{}Run{}-OverrepSeq.txt'.format(self.samplePath,
                                                            self.sampleNumber,
                                                            runNumber), 'w') as Log:
            Log.write('Overrepresented Sequences for {} and fastqc number {}\n'.format(
                                                    self.sampleName,runNumber))
            Log.write('-------------------------------------------------------------\n\n')

    def gatherOverrep(self, runNumber):
        ''' Arguments:
                runNumber = int; trial number for QC
            Returns:
                None

            Collects Overrepresented sequences together into a specific sample file
        '''
        QCs = glob.glob(self.samplePath + '/fastqc*')
        goodQC = [qc for qc in QCs if 'fastqc{}'.format(runNumber) in qc]
        if len(goodQC) != 1:
            with open('{}/Sample{}Run{}-OverrepSeq.txt'.format(self.samplePath,
                                                                self.sampleNumber,
                                                                runNumber), 'a') as Log:
                Log.write('Error: Unable to gather Overrepresented Sequences for {}\n'.format(
                            self.sampleName))
        else:
            unzipped = glob.glob(goodQC[0] + '/*/')
            for unzip in unzipped:
                Overrep = glob.glob(unzip + '/Overrep*')
                header = '-'.join(unzip.split('/')[-3:-1])
                with open('{}/Sample{}Run{}-OverrepSeq.txt'.format(self.samplePath,
                                                                    self.sampleNumber,
                                                                    runNumber), 'a') as Log:
                    Log.write('\n' + header + '\n')
                    with open(Overrep[0],'r') as O:
                        contents = O.read()
                    Log.write(contents + '\n')

    def runQCheck(self, runNumber, read1, read2):
        ''' Arguments:
                runNumber = int; trial number for QC
                read1 = str; first read to run fastQC on
                read2 = str; second read to run fastQC on
            Returns:
                None

            Packages self.runFastqc and self.findOverrepSeq together
        '''
        self.runFastqc(runNumber, read1, read2)
        self.findOverrepSeq(runNumber)
        self.gatherOverrep(runNumber)

    def runTrimmomatic(self, read1, read2):
        ''' Arguments:
                read1 = str; first read to run Trimmomatic on
                read2 = str; second read to run Trimmomatic on
            Returns:
                None

            Runs Trimmomatic on read1 and read2
        '''
        os.chdir(self.samplePath)
        Phred = self.getPhred()
        Reads = self.getReadNames()
        localArgs = self.GlobalArgs['runTrimmomatic']
        # Making Command
        #command = r'java -jar $RNASEQDIR/Trimmomatic/trimmomatic-0.36.jar PE -threads {procs} 
        # -phred{phred} {Read1} {Read2} read1.P.trim.{fastq}.gz read1.U.trim.{fastq}.gz 
        # read2.P.trim.{fastq}.gz read2.U.trim.{fastq}.gz 
        # ILLUMINACLIP:{black}:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:35'

        command = (r'''{java} -jar {-jar} PE -threads {-threads}'''+
            ''' -phred{-phred} {other} {read1} {read2} {read1pout} {read1uout}'''+
            ''' {read2pout} {read2uout} ILLUMINACLIP:{blacklist}:{ILLUMINACLIP}'''+
            ''' LEADING:{LEADING} TRAILING:{TRAILING}'''+
            ''' SLIDINGWINDOW:{SLIDINGWINDOW} MINLEN:{MINLEN}''')
        Context = {
                "java": self.checkMan("java", localArgs['java']),
                "-jar": self.checkMan(
                    "$RNASEQDIR/Trimmomatic/trimmomatic-0.36.jar",
                    localArgs['-jar']),
                "-threads": self.checkMan(self.Procs, localArgs['-threads']),
                "-phred": self.checkMan(Phred, localArgs['-phred']),
                "other": self.checkMan("", localArgs['other']),
                "read1": self.checkMan(read1, localArgs['java']),
                "read2": self.checkMan(read2, localArgs['java']),
                "read1pout": self.checkMan(
                    "read1.P.trim.{}.gz".format(self.Fastq),
                    localArgs['read1pout']),
                "read1uout": self.checkMan(
                    "read1.U.trim.{}.gz".format(self.Fastq),
                    localArgs['read1uout']),
                "read2pout": self.checkMan(
                    "read2.P.trim.{}.gz".format(self.Fastq),
                    localArgs['read2pout']),
                "read2uout": self.checkMan(
                    "read2.U.trim.{}.gz".format(self.Fastq),
                    localArgs['read2uout']),
                "blacklist": self.checkMan(self.Blacklist,
                    localArgs['blacklist']),
                "ILLUMINACLIP": self.checkMan("2:30:10",
                    localArgs['ILLUMINACLIP']),
                "LEADING": self.checkMan("5", localArgs['LEADING']),
                "TRAILING": self.checkMan("5", localArgs['TRAILING']),
                "SLIDINGWINDOW": self.checkMan("4:5",
                    localArgs['SLIDINGWINDOW']),
                "MINLEN": self.checkMan("35", localArgs['MINLEN']),
                }
        commandWithContext = command.format(**Context)
        goodCommand = self.formatCommand(commandWithContext)
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True,
                            executable='/bin/bash')

    def runPart1(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs Fastqc once, then Trimmomatic, and then fastqc again
        '''
        os.chdir(self.samplePath)
        with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'w') as LOG:
            LOG.write('\t\tRuntime Log Part 1 for {}\n'.format(self.sampleName))
            LOG.write('----------------------------------------\n\n')
        if self.GlobalArgs['fcounts/main']['runQCheck']:
            self.runQCheck(1, self.Read1, self.Read2)
        if self.GlobalArgs['fcounts/main']['runTrimmomatic']:
            self.runTrimmomatic(self.Read1, self.Read2)
        if (self.GlobalArgs['fcounts/main']['runQCheck'] and 
            self.GlobalArgs['fcounts/main']['runTrimmomatic']):
            self.runQCheck(2, 'read1.P.trim.{}.gz'.format(self.Fastq),
                'read2.P.trim.{}.gz'.format(self.Fastq))
        self.writeFunctionTail('runPart1')

    ########################################################
    # Determine strandicity
    ########################################################

    def runSeqtk(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs seqtk on reads from trimmomatic
        '''
        self.writeFunctionHeader('runSeqtk')
        localArgs = self.GlobalArgs['runSeqtk']
        # Making Command
        #command1 = r'seqtk sample -s100 read1.P.trim.{0}.gz 10000 |'+
        #' seqtk seq -A - > sampled.read1.fa'.format(self.Fastq)
        #command2 = r'seqtk sample -s100 read2.P.trim.{0}.gz 10000 |'+
        #' seqtk seq -A - > sampled.read2.fa'.format(self.Fastq)
        command1 = ("""{seqtk} sample -s{-s} {read1} {samples} |"""+
                """ seqtk seq {-A} {other} - > {read1out}""")
        command2 = ("""{seqtk} sample -s{-s} {read2} {samples} |"""+
                """ seqtk seq {-A} {other} - > {read2out}""")
        Context = {
                "seqtk": self.checkMan("seqtk", localArgs['seqtk']),
                "-s": self.checkMan("100", localArgs['-s']),
                "read1": self.checkMan("read1.P.trim.{}.gz".format(self.Fastq),
                    localArgs['read1']),
                "read2": self.checkMan("read2.P.trim.{}.gz".format(self.Fastq),
                    localArgs['read2']),
                "samples": self.checkMan("10000", localArgs['samples']),
                "-A": self.checkManBool("-A", localArgs['-A']),
                "read1out": self.checkMan("sampled.read1.fa",
                    localArgs['read1out']),
                "read2out": self.checkMan("sampled.read2.fa",
                    localArgs['read2out']),
                }
        goodCommand1 = self.formatCommand(command1.format(**Context))
        goodCommand2 = self.formatCommand(command2.format(**Context))
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand1)
        subprocess.run(goodCommand1,
                            shell=True,
                            check=True,
                            executable='/bin/bash')
        self.writeFunctionCommand(goodCommand2)
        subprocess.run(goodCommand2,
                            shell=True,
                            check=True,
                            executable='/bin/bash')
        self.writeFunctionTail('runSeqtk')

    def runBlastn(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs Blastn on seqtk output
        '''
        self.writeFunctionHeader('runBlastn')
        # Making Command
        db = os.path.join(self.Reference,'{}.cdna.all'.format(self.Basename))
        localArgs = self.GlobalArgs['runBlastn']
        #command1 = r'blastn -query sampled.read1.fa -db {0}
        #-out sampled.read1_vscdna.out -task blastn-short
        #-outfmt "6 std sstrand" -max_target_seqs 1 
        #-num_threads {1}'.format(db, self.Procs)
        #
        #command2 = r'blastn -query sampled.read2.fa -db {0}
        #-out sampled.read2_vscdna.out -task blastn-short
        #-outfmt "6 std sstrand" -max_target_seqs 1
        #-num_threads {1}'.format(db, self.Procs)
        command1 = ('''{blastn} -query {read1} -db {-db} -out {out1}'''+
                ''' -task {-task} -outfmt {-outfmt}'''+
                ''' -max_target_seqs {-max_target_seqs}'''+
                ''' -num_threads {-num_threads} {other}''')
        command2 = ('''{blastn} -query {read2} -db {-db} -out {out2}'''+
                ''' -task {-task} -outfmt {-outfmt}'''+
                ''' -max_target_seqs {-max_target_seqs}'''+
                ''' -num_threads {-num_threads} {other}''')
        Context = {
                "blastn": self.checkMan("blastn", localArgs['blastn']),
                "read1": self.checkMan("sampled.read1.fa", localArgs['read1']),
                "read2": self.checkMan("sampled.read2.fa", localArgs['read2']),
                "-db": self.checkMan(db, localArgs['-db']),
                "out1": self.checkMan("sampled.read1_vscdna.out",
                    localArgs['out1']),
                "out2": self.checkMan("sampled.read2_vscdna.out",
                    localArgs['out2']),
                "-task": self.checkMan("blastn-short", localArgs['-task']),
                "-outfmt": self.checkMan('"6 std sstrand"',
                    localArgs['-outfmt']),
                "-max_target_seqs": self.checkMan("1",
                    localArgs['-max_target_seqs']),
                "-num_threads": self.checkMan(self.Procs,
                    localArgs['-num_threads']),
                "other": self.checkMan("", localArgs['other']),
                }
        goodCommand1 = self.formatCommand(command1.format(**Context))
        goodCommand2 = self.formatCommand(command2.format(**Context))
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand1)
        subprocess.run(goodCommand1,
                            shell=True,
                            check=True,
                            executable='/bin/bash')
        self.writeFunctionCommand(goodCommand2)
        subprocess.run(goodCommand2,
                            shell=True,
                            check=True,
                            executable='/bin/bash')
        self.writeFunctionTail('runBlastn')

    def findStranded(self):
        ''' Arguments:
                None
            Returns:
                Boolean:{
                            True: if data is stranded
                            False: if data is not stranded
                        }
            Runs stranded_classifier.py and then scrapes output
            to determine if data is stranded or not
        '''
        # Running stranded_classifier.py
        # Making Command
        command1 = r'stranded_classifier.py -1 sampled.read1_vscdna.out -2 sampled.read2_vscdna.out'
        goodCommand1 = self.formatCommand(command1)
        self.writeFunctionCommand(goodCommand1)
        with open('{}/Runtime.{}.log'.format(self.samplePath,
                                            self.sampleName), 'a') as LOG:
            LOG.write('findStranded started\n\n')
        # Executing
        os.chdir(self.samplePath)
        subprocess.run(goodCommand1,
                            shell=True,
                            check=True)
        with open('{}/Runtime.{}.log'.format(self.samplePath,
                                            self.sampleName), 'a') as LOG:
            LOG.write('\nfindStranded done\n\n')
        # Scraping stranded_classifier.py output
        command2 = r"awk '/read1 plus percent/,/findStranded done/' Runtime.{}.log | grep -q True && stranded=1 || stranded=0 && echo $stranded".format(self.sampleName)
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(command2)
        proc = subprocess.Popen(command2,
                                    shell=True,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT)
        strandedBool = int(proc.communicate()[0])
        if strandedBool == 1:
            Data = []
            with open('Runtime.{}.log'.format(self.sampleName),'r') as F:
                copy = False
                for line in F:
                    if 'read1 plus percent' in line.strip():
                        copy = True
                        Data.append(line)
                    elif line.strip() == 'findStranded done':
                        copy = False
                    elif copy:
                        Data.append(line)
            for line in Data:
                if line.startswith('read1 plus percent'):
                    R1 = float(re.findall(r'\d+\.\d+', line)[0])
                if line.startswith('read2 plus percent'):
                    R2 = float(re.findall(r'\d+\.\d+', line)[0])
            if R1 > R2:
                self.writeToLog('! Going to use "--rna-strandness FR" for hisat and "-s 1" for FC or --fr-stranded for kallisto quant')
                return 1
            else:
                self.writeToLog('! Going to use "--rna-strandness RF" for hisat and "-s 2" for FC or --rf-stranded for kallisto quant')
                return 2
        elif strandedBool == 0:
            self.writeToLog('! Reads not stranded')
            return 0
        else:
            print('There was an error with findStranded')

    ########################################################
    # End of Sample class
    ########################################################

class FCountsSample(Sample):

    def __init__(self,sampleNumber,globalArgs):
        Sample.__init__(self, sampleNumber, globalArgs)

    def __repr__(self):
        ''' Arguments:
                None
            Returns:
                String representation print
                ex: >print(Sample)
                        Sample(sample_01)

            Returns string representation of class
        '''
        return 'FCountsSample(%r)'%(self.sampleName)

    ########################################################
    # Pipeline
    ########################################################

    def runHisat(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs hisat2 on data
        '''
        self.writeFunctionHeader('runHisat')
        strandedVar = self.findStranded()
        if strandedVar == 0:
            FR = ''
        elif strandedVar == 1:
            FR = ' --rna-strandness FR'
        else:
            FR = ' --rna-strandness RF'
        Phred = str(self.getPhred())
        localArgs = self.GlobalArgs['fcounts/runHisat']
        # Making Command
        #command = r"""hisat2 -k 5 -p {numProcs}{FRoRF} --dta 
        # --phred{phred} --known-splicesite-infile {ref}/splice_sites.txt
        # -x {ref}/{basename} -1 read1.P.trim.{fastq}.gz 
        # -2 read2.P.trim.{fastq}.gz -S aligned.{sample}.sam"""
        command = ("""{hisat2} -k {-k} -p {-p} {--rna-strandedness}"""+
            """ {--dta} --phred{phred} {other}"""+
            """ --known-slicesite-infile {--known-splicesite-infile}"""+
            """ -x {-x} -1 {-1} -2 {-2} -S {-S}""")
        Context = {
                "hisat2": self.checkMan("hisat2", localArgs['hisat2']),
                "-k": self.checkMan("5", localArgs['-k']),
                "-p": self.checkMan(self.Procs, localArgs['-p']),
                "--rna-strandedness": (FR 
                    if localArgs['--rna-strandedness'] == None
                    else "--rna-strandedness "+localArgs['--rna-strandedness']),
                "--dta": self.checkManBool("--dta", localArgs['--dta']),
                "phred": self.checkMan(Phred, localArgs['phred']),
                "other": self.checkMan("", localArgs['other']),
                "--known-splicesite-infile": self.checkMan(
                    "{}/splice_sites.txt".format(self.Reference),
                    localArgs['--known-splicesite-infile']),
                "-x": self.checkMan(
                    "{}/{}".format(self.Reference, self.Basename),
                    localArgs['-x']),
                "-1": self.checkMan(
                    "read1.P.trim.{}.gz".format(self.Fastq),
                    localArgs['-1']),
                "-2": self.checkMan(
                    "read2.P.trim.{}.gz".format(self.Fastq),
                    localArgs['-2']),
                "-S": self.checkMan(
                    "aligned.{}.sam".format(self.sampleName),
                    localArgs['-S']),
                }
        goodCommand = self.formatCommand(command.format(**Context))
        # Executing
        self.writeFunctionCommand(goodCommand)
        os.chdir(self.samplePath)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True,
                            executable=True)
        self.writeFunctionTail('runHisat')

    def runCompression(self):
        ''' Arguments:
                None
            Returns:
                None

            Compresses output of hisat2 with samtools
        '''
        self.writeFunctionHeader('runCompression')
        localArgs = self.GlobalArgs['fcounts/runCompression']
        # Making Command
        #command = r"samtools view -bT {ref}/{genome} -@{procs}
        # aligned.{sample}.sam -o aligned.{sample}.bam"
        command = ("""{samtools} view {-b} -T {-T} -@{-@} {other}"""+
            """ {in} -o {-o}""")
        Context = {
                "samtools"  : self.checkMan("samtools",
                    localArgs['samtools']),
                "-b"        : self.checkManBool("-b",
                    localArgs['-b']),
                "-T"        : self.checkMan(
                    "{}/{}".format(self.Reference, self.Genome),
                    localArgs['-T']),
                "-@"        : self.checkMan(self.Procs,
                    localArgs['-@']),
                "other"     : self.checkMan("",
                    localArgs['other']),
                "in"        : self.checkMan(
                    "aligned.{}.sam".format(self.sampleName),
                    localArgs['in']),
                "-o"        : self.checkMan(
                    "aligned.{}.bam".format(self.sampleName),
                    localArgs['-o']),
                }
        goodCommand = self.formatCommand(command.format(**Context))
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True,
                            executable='/bin/bash')
        os.chdir(self.samplePath)
        os.remove(self.checkMan("aligned.{}.sam".format(self.sampleName),
                    localArgs['in']))
        self.writeFunctionTail('runCompression')

    ########################################################
    # Feature Counts
    ########################################################

    def runFeatureCounts(self):
        ''' Arguments:
                None
            Returns:
                None

            Collects counts from data with featureCounts
        '''
        self.writeFunctionHeader('runFeatureCounts')
        strandedVar = str(self.findStranded())
        localArgs = self.GlobalArgs['fcounts/runFeatureCounts']
        # Making Command
        #command = r"featureCounts -T {procs} -s {stranded} -p -C --primary
        # --ignoreDup -t exon -g gene_id -a {ref}/{gtf} 
        # -o aligned.{sample}.counts aligned.{sample}.bam"
        command = ("""{featureCounts} -T {-T}"""+
            """ -s {-s}{-p}{-C}{--primary}{--ignoreDup} -t {-t} -g {-g}"""+
            """ -a {-a} {other} -o {-o} {in}""")
        Context = {
                "featureCounts" : self.checkMan("featureCounts",
                                    localArgs['featureCounts']),
                "-T"            : self.checkMan(self.Procs,
                                    localArgs['-T']),
                "-s"            : self.checkMan(strandedVar,
                                    localArgs['-s']),
                "-p"            : self.checkManBool(" -p",
                                    localArgs['-p']),
                "-C"            : self.checkManBool(" -C",
                                    localArgs['-C']),
                "--primary"     : self.checkManBool(" --primary",
                                    localArgs['--primary']),
                "--ignoreDup"   : self.checkManBool(" --ignoreDup",
                                    localArgs['--ignoreDup']),
                "-t"            : self.checkMan("exon",
                                    localArgs['-t']),
                "-g"            : self.checkMan("gene_id",
                                    localArgs['-g']),
                "-a"            : self.checkMan(
                                    os.path.join(self.Reference, self.Gtf),
                                    localArgs['-a']),
                "-o"            : self.checkMan(
                                    "aligned.{}.counts".format(self.sampleName),
                                    localArgs['-o']),
                "in"            : self.checkMan(
                                    "aligned.{}.bam".format(self.sampleName),
                                    localArgs['in']),
                "other"         : self.checkMan("",
                                    localArgs['other']),
                }
        goodCommand = self.formatCommand(command.format(**Context))
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand)
        subprocess.run(goodCommand,
            shell=True,
            check=True,
            executable='/bin/bash')
        self.writeFunctionTail('runFeatureCounts')

    ########################################################
    # Gathering Data (for non-stringtie)
    ########################################################

    def getNiceColumns(self):
        ''' Arguments:
                None
            Returns:
                None

            Scrapes featureCounts output for gene id, length of gene, and count
        '''
        self.writeFunctionHeader('getNiceColumns')
        # Making Command
        command = r'''tail -n +2 aligned.{sample}.counts | awk '{{printf ("%5s\t%s\t%s\n", $1, $6, $7)}}' > aligned.{sample}.counts.three'''
        context = {"sample": self.sampleName}
        goodCommand = self.formatCommand(command.format(**context))
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True)
        self.writeFunctionTail('getNiceColumns')

    def getAlignedColumn(self):
        ''' Arguments:
                None
            Returns:
                None

            Scrapes featureCounts output to get count column
        '''
        self.writeFunctionHeader('getAlignedColumn')
        # Making Command
        command = r'''tail -n +2 aligned.{sample}.counts | awk '{{print $7}}' > aligned.{sample}.counts.one'''
        context = {"sample": self.sampleName}
        goodCommand = self.formatCommand(command.format(**context))
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True)
        self.writeFunctionTail('getAlignedColumn')

    ########################################################
    # Run Sample
    ########################################################

    def runPart2(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs Pipeline sequentially
            Note: Does not include Quality control steps: Fastqc and trimmomatic
        '''
        os.chdir(self.samplePath)
        with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'a') as LOG:
            LOG.write('\t\tRuntime Log Part 2 for {}\n'.format(self.sampleName))
            LOG.write('----------------------------------------\n\n')
        if self.GlobalArgs['fcounts/main']['runQCheck']:
            self.runSeqtk()
        if self.GlobalArgs['fcounts/main']['runBlastn']:
            self.runBlastn()
        if self.GlobalArgs['fcounts/main']['runHisat']:
            self.runHisat()
        if self.GlobalArgs['fcounts/main']['runCompression']:
            self.runCompression()
        if self.GlobalArgs['fcounts/main']['runFeatureCounts']:
            self.runFeatureCounts()
        if self.GlobalArgs['fcounts/main']['getNiceColumns']:
            self.getNiceColumns()
        if self.GlobalArgs['fcounts/main']['getAlignedColumn']:
            self.getAlignedColumn()
        self.writeFunctionTail('runPart2')
        with open(self.Project + '/runPipeNotify/{}'.format('done'+self.sampleName), 'w') as N:
            N.write('{} is done'.format(self.samplePath))
        with open(self.samplePath + '/.done', 'w') as N:
            N.write('{} is done'.format(self.samplePath))

    def runParts(self):
        ''' Arguments:
                None
            Returns:
                None

            Run QC and Pipeline part sequentially
        '''
        self.runPart1()
        self.runPart2()

    ########################################################
    # End of FCountsSample class
    ########################################################
class StringtieSample(Sample):

    def __init__(self,sampleNumber,globalArgs):
        Sample.__init__(self, sampleNumber, globalArgs)

    def __repr__(self):
        ''' Arguments:
                None
            Returns:
                String representation print
                ex: >print(Sample)
                        Sample(sample_01)

            Returns string representation of class
        '''
        return 'StringtieSample(%r)'%(self.sampleName)

    ########################################################
    # Pipeline
    ########################################################

    def runHisat(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs hisat2 on data
        '''
        self.writeFunctionHeader('runHisat')
        strandedVar = self.findStranded()
        if strandedVar == 0:
            FR = ''
        elif strandedVar == 1:
            FR = ' --rna-strandness FR'
        else:
            FR = ' --rna-strandness RF'
        # Making Command
        command = r"""hisat2 -k 5 -p {numProcs}{FRoRF} --dta --phred{phred} --known-splicesite-infile {ref}/splice_sites.txt -x {ref}/{basename} -1 read1.P.trim.{fastq}.gz -2 read2.P.trim.{fastq}.gz -S aligned.{sample}.sam"""
        Phred = self.getPhred()
        context = {
                "numProcs": self.Procs,
                "phred": str(Phred),
                "ref": self.Reference,
                "basename": self.Basename,
                "fastq": self.Fastq,
                "sample": self.sampleName,
                "FRoRF": FR
                }
        goodCommand = self.formatCommand(command.format(**context))
        # Executing
        self.writeFunctionCommand(goodCommand)
        os.chdir(self.samplePath)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True)
        self.writeFunctionTail('runHisat')

    def runCompression(self):
        ''' Arguments:
                None
            Returns:
                None

            Compresses output of hisat2 with samtools
        '''
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
        self.writeFunctionCommand(goodCommand)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True)
        os.chdir(self.samplePath)
        os.remove('aligned.{}.sam'.format(self.sampleName))
        self.writeFunctionTail('runCompression')

    ########################################################
    # Stringtie
    ########################################################

    def assembleTranscripts(self):
        ''' Arguments:
                None
            Returns:
                None

            Assemble gene transcripts using Reference GTF
        '''
        self.writeFunctionHeader('assembleTranscripts')
        # Making Command
        #command = r"stringtie -p {procs} -G {ref}/{gtf} -o {sample}.st.gtf -l {sample} aligned.{sample}.bam"
        command = r"stringtie -p {procs} -o {sample}.st.gtf -l {sample} aligned.{sample}.bam"
        context = {
                "procs": self.Procs,
                "ref": self.Reference,
                "gtf": self.Gtf,
                "sample": self.sampleName,
                }
        goodCommand = self.formatCommand(command.format(**context))
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True)
        self.writeFunctionTail('assembleTranscripts')

    def estimateTranscriptAbundances(self):
        ''' Arguments:
                None
            Returns:
                None

            Using the stringtie merged GTF file, estimate the
            transcript abundances and create a table counts
        '''
        self.writeFunctionHeader('estimateTranscriptAbundances')
        projectName = os.path.basename(self.Project)
        stMergeDir = self.Postprocessing + '/StringtieMerge'
        stMergedFile = os.path.join(stMergeDir,"{}.stmerged.gtf".format(projectName))
        # Making Command
        command = r"stringtie -e -B -p {procs} -G {stmerged} -o {sample}.good.st.gtf aligned.{sample}.bam"
        context = {
                "procs": self.Procs,
                "stmerged": stMergedFile,
                "sample": self.sampleName
                }
        goodCommand = self.formatCommand(command.format(**context))
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True)
        self.writeFunctionTail('estimateTranscriptAbundances')

    def runAltCompression(self):
        ''' Arguments:
                None
            Returns:
                None

            Compresses output of hisat2 with samtools
        '''
        self.writeFunctionHeader('runCompression')
        # Making Command
        command = r"samtools sort -@ {procs} -o aligned.{sample}.bam aligned.{sample}.sam"
        context = {
                "procs": self.Procs,
                "sample": self.sampleName,
                }
        goodCommand = self.formatCommand(command.format(**context))
        # Executing
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True)
        os.chdir(self.samplePath)
        os.remove('aligned.{}.sam'.format(self.sampleName))
        self.writeFunctionTail('runCompression')

    ########################################################
    # Run Sample
    ########################################################

    def stringtiePart2a(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs Pipeline sequentially
        '''
        self.writeFunctionHeader('stringtiePart2a')
        os.chdir(self.samplePath)
        with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'a') as LOG:
            LOG.write('\t\tRuntime Log Part 2 for {}\n'.format(self.sampleName))
            LOG.write('----------------------------------------\n')
            LOG.write('Part 2a\n')
            LOG.write('----------------------------------------\n\n')
        self.runSeqtk()
        self.runBlastn()
        self.runHisat()
        self.runAltCompression()
        self.assembleTranscripts()
        self.writeFunctionTail('stringtiePart2a')

    def stringtiePart2c(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs Pipeline sequentially
        '''
        self.writeFunctionHeader('stringtiePart2c')
        os.chdir(self.samplePath)
        with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'a') as LOG:
            LOG.write('Part 2c\n')
            LOG.write('----------------------------------------\n\n')
        self.estimateTranscriptAbundances()
        self.writeFunctionTail('stringtiePart2c')
        with open(self.Project + '/runPipeNotify/{}'.format('done'+self.sampleName), 'w') as N:
            N.write('{} is done'.format(self.samplePath))
        with open(self.samplePath + '/.done', 'w') as N:
            N.write('{} is done'.format(self.samplePath))

    def runParts(self,stringtiePhase=None):
        ''' Arguments:
                None
            Returns:
                None

            Run QC and Pipeline part sequentially
        '''
        if stringtiePhase == 'a':
            self.runPart1()
            self.stringtiePart2a()
        elif stringtiePhase == 'c':
            self.stringtiePart2c()
        else:
            raise ValueError('stringtiePhase must be specified for stringtie sample calculations')

    ########################################################
    # End of StringtieSample class
    ########################################################

class KallistoSample(Sample):

    def __init__(self,sampleNumber,globalArgs):
        Sample.__init__(self, sampleNumber, globalArgs)

    def __repr__(self):
        ''' Arguments:
                None
            Returns:
                String representation print
                ex: >print(Sample)
                        Sample(sample_01)

            Returns string representation of class
        '''
        return 'KallistoSample(%r)'%(self.sampleName)

    ########################################################
    # Kallisto
    ########################################################

    def runKallisto(self):
        ''' Arguments:
                None
            Returns:
                None

        '''
        self.writeFunctionHeader('runKallisto')
        kallistoOutputDir = os.path.join(self.samplePath, 'KallistoOutput.{}'.format(
                                                            self.sampleName))
        strandedVar = self.findStranded()
        if strandedVar == 0:
            FR = ''
        elif strandedVar == 1:
            FR = ' --fr-stranded'
        else:
            FR = ' --rf-stranded'
        # Making Command
        command = r"kallisto quant -i {transcriptindex} -o {outputdir} --threads {procs}{FRoRF} -b 100 <( zcat read1.P.trim.{fastq}.gz ) <( zcat read2.P.trim.{fastq}.gz )" 
        context = {
                "transcriptindex": os.path.join(self.Reference, '{}.kali.cdna.fa.idx'.format(
                                                                            self.Basename)),
                "outputdir": kallistoOutputDir,
                "fastq": self.Fastq,
                "FRoRF": FR,
                "procs": self.Procs
                }
        goodCommand = self.formatCommand(command.format(**context))
        # Executing
        if not os.path.isdir(kallistoOutputDir):
            try:
                os.mkdir(kallistoOutputDir)
            except:
                pass
        os.chdir(self.samplePath)
        self.writeFunctionCommand(goodCommand)
        subprocess.run(goodCommand,
                            shell=True,
                            check=True,
                            executable="/bin/bash")
        self.writeFunctionTail('runKallisto')

    ########################################################
    # Run Sample
    ########################################################

    def runPart2Kallisto(self):
        ''' Arguments:
                None
            Returns:
                None

            Runs Pipeline sequentially
            Note: Does not include Quality control steps: Fastqc and trimmomatic
        '''
        os.chdir(self.samplePath)
        with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'a') as LOG:
            LOG.write('\t\tRuntime Log Part 2 for {}\n'.format(self.sampleName))
            LOG.write('----------------------------------------\n\n')
        self.runSeqtk()
        self.runBlastn()
        self.runKallisto()
        self.writeFunctionTail('runPart2')
        with open(self.Project + '/runPipeNotify/{}'.format('done'+self.sampleName), 'w') as N:
            N.write('{} is done'.format(self.samplePath))
        with open(self.samplePath + '/.done', 'w') as N:
            N.write('{} is done'.format(self.samplePath))

    def runParts(self):
        ''' Arguments:
                None
            Returns:
                None

            Run QC and Pipeline part sequentially
        '''
        self.runPart1()
        self.runPart2Kallisto()

    ########################################################
    # End of KallistoSample class
    ########################################################

#################################################################
if __name__ == '__main__':
    raise SystemExit('Please use runPipe\nSee runPipe --help')
