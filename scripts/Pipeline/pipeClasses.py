#!/usr/bin/python3

'''Usage: pipeClasses.py [-h | --help] [-j <jsonFile> | --jsonfile <jsonFile>]
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
from math import ceil as ceiling
import multiprocessing
import time
import pipeUtils
import os
import subprocess
import shutil
import glob
import json

################################################################
# Utilities
################################################################

def exportVariablesforClass(pathtoInput,maxCPU=None):
    ''' Arguments:
            pathtoInput = string; path to INPUT file
        Returns:
            Variables

        Gathers instance variables for Experiment Class
    '''
    Project, Reference, Original = pipeUtils.sourceInput(pathtoInput)
    Gtf, Cdna, Genome = pipeUtils.getReferenceVariables(Reference)
    Basename = pipeUtils.getBasename(Genome)
    Fastq = pipeUtils.getFastq(Original)
    Procs = pipeUtils.countCPU(maxCPU)
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

def checkJSON(jsonFile, behavior='default'):
    ''' Arguments:
            jsonFile = string; path to JSON file
        Returns:
            None

        Tests if JSON file has correct syntax
    '''
    if behavior == 'default':
        with open(jsonFile) as JF:
            try:
                json.load(JF)
            except ValueError as e:
                raise SystemExit("Invalid JSON: %s"%(e))
    else:
        if jsonFile == None:
            raise SystemExit("Need --json argument")
        else:
            with open(jsonFile) as JF:
                try:
                    json.load(JF)
                except ValueError:
                    raise SystemExit("Need --json argument with valid JSON file")

def funTime(function):
    ''' Arguments:
            None
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
        R.write('runPipe.py Runtime File\n-----------------------------------------\n\n')

def unwrap_self_runSample(arg, **kwarg):
    ''' Magic for multiprocessing '''
    return Experiment.runSample(*arg, **kwarg)

################################################################
# Defining Experiment Class
################################################################

class Experiment:
    ''' Experiment is a Paired-End RNA Sequencing Experiment.
        This class provides tools to analyze data
    '''

    def __init__(self, inputPath, maxCPU=None):
        ''' Arguments:
                inputPath = string; path to INPUT file
            Returns:
                None

            Initializes variables that get scraped from
            exportVariablesforClass(inputPath)
        '''
        variables = exportVariablesforClass(inputPath,maxCPU)
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

    def getNumberofSamples(self):
        ''' Arguments:
                None
            Returns:
                numSamp = int; The number of samples in experiment

            Finds number of samples in Experiment by counting number
            of files in the Original Data folder and dividing by 2
            because there are 2 reads to a sample
        '''
        if pipeUtils.getNumberofFiles(self.ogOriginal)%2 != 0:
            print('There are not an even number of files in {}'.format(
                                                    self.ogOriginal))
        else:
            numSamp = pipeUtils.getNumberofFiles(self.ogOriginal)/2
        return int(numSamp)

    def checkStructure(self):
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
            The 'S' is added after pipeUtils.createSymLinks has completed
        '''
        if os.path.exists(self.Project + '/.init'):
            with open(self.Project + '/.init', 'r') as f:
                G = f.read()
            if 'S' in G:
                return True
        return False

    def checkReference(self):
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
            The 'P' is added after pipeUtils.preProcessingReference has completed
        '''

        if os.path.exists(self.Project + '/.init'):
            with open(self.Project + '/.init', 'r') as f:
                G = f.read()
            if 'P' in G:
                return True
        return False

    # @funTime is just a wrapper that times function
    @funTime
    def makeStructure(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.createStructure to make Directory Structure
        '''
        if noconfirm == True:
            print("Creating Structure...")
            pipeUtils.createStructure(self.Project,
                                        self.ogOriginal)
            makeTimeFile(RUNTIMELOG)
        else:
            if self.checkStructure() == False:
                print("Creating Structure...")
                pipeUtils.createStructure(self.Project,
                                            self.ogOriginal)
                makeTimeFile(RUNTIMELOG)

    @funTime
    def makeSyms(self,exists=False):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.createSymLinks to make Symbolic Links for
            Original Data, Reference Data; and Metadata if available
        '''
        if noconfirm == True:
            pipeUtils.createSymLinks(self.Project,
                                        self.ogOriginal,
                                        self.ogReference,
                                        exists)
        else:
            if self.checkStructure() == False:
                pipeUtils.createSymLinks(self.Project,
                                            self.ogOriginal,
                                            self.ogReference,
                                            exists)

    @funTime
    def qcRef(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.qcReference to make Reference_Report.txt in Reference
            Folder
        '''
        if self.checkReference() == False:
            pipeUtils.qcReference(self.Reference,
                                    self.Genome)

    @funTime
    def ppRef(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.preProcessingReference to Pre-Process Reference Data
        '''
        if self.checkReference() == False:
            print("Preprocessing Data...")
            pipeUtils.preProcessingReference(self.Reference,
                                                self.Cdna,
                                                self.Gtf,
                                                self.Genome,
                                                self.Basename)

    def createSampleClasses(self,subject='all'):
        ''' Arguments:
                None
            Returns:
                Samples = list; list contains Sample Classes

            Creates all Sample Classes necessary and returns them in a list
            Number of samples is based on self.getNumberofSamples() function
        '''
        if self.Procs < os.cpu_count():
            if subject == 'all':
                Samples = [Experiment.Sample(n, self.inputPath,self.Procs) for n
                        in range(1,self.getNumberofSamples() + 1)]
                return Samples
            elif type(subject) == int:
                Sample = Experiment.Sample(subject, self.inputPath, self.Procs)
                return Sample
        else:
            if subject == 'all':
                Samples = [Experiment.Sample(n, self.inputPath) for n
                        in range(1,self.getNumberofSamples() + 1)]
                return Samples
            elif type(subject) == int:
                Sample = Experiment.Sample(subject, self.inputPath)
                return Sample

    def runSample(self, sample):
        ''' Arguments:
                sample = class; runs pipeline on sample
            Returns:
                None

            Function calls Sample class function runParts()
            Function used in multiprocessing as map function on
                self.createSampleClasses() list
        '''
        sample.runParts()

    @funTime
    def GO(self,subject='all'):
        ''' Arguments:
                None
            Returns:
                None

            Uses python multiprocessing to distribute sample analysis to
            CPUs
        '''
        if subject == 'all':
            self.makeNotifyFolder()
            Samples = self.createSampleClasses()
            with multiprocessing.Pool(self.Procs) as p:
                p.map(unwrap_self_runSample, zip([self]*len(Samples), Samples))
        elif type(subject) == int:
            name = 'sample_{:02g}'.format(subject)
            if os.path.exists(self.Data + '/' + name):
                Sample = self.createSampleClasses(subject)
                self.runSample(Sample)
        else:
            raise SystemExit('There is a problem with executing sample')

    def getOptimal(self, cluster):
        ''' Arguments:
                cluster = list; a list of integers that describes number of
                            CPUs on each machine in a cluster
                numSamples = int; the number of samples that you need to run
            Returns:
                dict; what you need to run

            Figures out best way to run your samples
        '''
        numSamples = self.getNumberofSamples()
        def atatime(cpuPerSample,limits=cluster):
            return sum([x//cpuPerSample for x in limits])
        def iterationsTilFinish(cpuPerSample, numSamp=numSamples, limits=cluster):
            return ceiling(numSamp/atatime(cpuPerSample,limits))
        def timePerSample(cpuPerSample):
            return 2221.85*cpuPerSample**-0.831626
        def totalTime(cpuPerSample, numSamp=numSamples, limits=cluster):
            ''' If your cpuPerSample is held constant for all iterations required'''
            return iterationsTilFinish(cpuPerSample,numSamp,limits)*timePerSample(cpuPerSample)
        def oneIteration(cpuPerSample, numSamp=numSamples, limits=cluster):
            return [numSamp-atatime(cpuPerSample,limits),[cpuPerSample,timePerSample(cpuPerSample)]]
        def chokePoints(limits=cluster):
            chokes,final = [],[]
            for amount in [(n,atatime(n,limits)) for n in range(1,max(limits)+1)]:
                if amount[1] not in chokes:
                    chokes.append(amount[1])
                    final.append(amount)
                else:
                    final.pop(-1)
                    final.append(amount)
            return [a[0] for a in final]
        def smart(numSamp=numSamples,limits=cluster):
            Paths,Finished,Best,counter = [],[1],9E9,0
            available = chokePoints(limits)
            while True:
                if counter == 0:
                    for M in available:
                        ITER = oneIteration(M,numSamp,limits)
                        if ITER[0] <= 0:
                            test = [sum(i) for i in zip(*ITER[1:])][1]
                            if test < Best:
                                Best = test
                                del Finished[0]
                                Finished.append(ITER[1:])
                        else: Paths.append(ITER)
                    counter += 1
                else:
                    for path in Paths:
                        for n in available:
                            ITER = path[:]
                            ITER[0] = ITER[0] - atatime(n)
                            ITER.append(oneIteration(n,numSamp,limits)[1])
                            if ITER[0] <= 0:
                                test = [sum(i) for i in zip(*ITER[1:])][1]
                                if test < Best:
                                    Best = test
                                    del Finished[0]
                                    Finished.append(ITER[1:])
                            else:
                                Paths.append(ITER)
                        Paths.remove(path)
                    counter += 1
                if len(Paths) == 0:
                    break
            return Best,Finished
        #def smartMultiprocessing(M,numSamp=numSamples,limits=cluster):
        #    Paths,Finished,Best,counter = [],[1],9E9,0
        #    available = chokePoints(limits)
        #    while True:
        #        if counter == 0:
        #            ITER = oneIteration(M,numSamp,limits)
        #            if ITER[0] <= 0:
        #                test = [sum(i) for i in zip(*ITER[1:])][1]
        #                if test < Best:
        #                    Best = test
        #                    del Finished[0]
        #                    Finished.append(ITER[1:])
        #            else:
        #                Paths.append(ITER)
        #            counter += 1
        #        else:
        #            for path in Paths:
        #                for n in available:
        #                    ITER = path[:]
        #                    ITER[0] = ITER[0] - atatime(n)
        #                    ITER.append(oneIteration(n,numSamp,limits)[1])
        #                    if ITER[0] <= 0:
        #                        test = [sum(i) for i in zip(*ITER[1:])][1]
        #                        if test < Best:
        #                            Best = test
        #                            del Finished[0]
        #                            Finished.append(ITER[1:])
        #                    else:
        #                        Paths.append(ITER)
        #                Paths.remove(path)
        #            #counter += 1
        #        if len(Paths) == 0:
        #            break
        #    return Best,Finished
        def scrapeResults(Best,numSamps=numSamples,limits=cluster):
            resultDict,counter,totalSamps = {},0,numSamps
            while True:
                S = atatime(Best[1][0][counter][0],limits)
                if S > totalSamps:
                    S = totalSamps
                else:
                    totalSamps -= S
                resultDict['Step {}'.format(counter+1)] = {}
                resultDict['Step {}'.format(counter+1)]['Procs'] = Best[1][0][counter][0]
                resultDict['Step {}'.format(counter+1)]['Samps'] = S
                counter += 1
                if counter+1 > len(Best[1][0]):
                    break
            return resultDict
        #with multiprocessing.Pool(self.Procs) as p:
        #    available = chokePoints(cluster)
        #    Results = p.map(smartMultiprocessing, available)
        #    bestResult = list(sorted(Results)[0])
        #    return scrapeResults(bestResult)
        return scrapeResults(smart())

    def makeBatchBiox(self):
        ''' Arguments:
                None
            Returns:
                None

            Creates Pipeline batch script to be used with slurm
        '''
        batch =  """#!/bin/bash
#SBATCH --nodes={NODES}
#SBATCH --time=120
#SBATCH --cpus-per-task={CPT}
#SBATCH --ntasks={NTASKS}
#SBATCH --job-name="Pipeline"
#SBATCH --export=PATH,RNASEQDIR

inputFile='{INPUT}'
jsonFile='{JSON}'

# Stage 1 and 2
{STAGE12}

wait

# Stage 3
{STAGE3}
wait

# Stage 4
{STAGE4}

wait

# Stage 5
{STAGE5}

scontrol show job $SLURM_JOB_ID
wait
"""
        numSamps = self.getNumberofSamples()
        command12 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 1,2 "${inputFile}"'
        bestPath = self.getOptimal([48,32])
        command3,counter,sampleNum = '',1,1
        for path in bestPath:
            Pstep = bestPath[path]['Procs']
            Sstep = bestPath[path]['Samps']
            for S in range(sampleNum,sampleNum + Sstep):
                com = 'srun -N1 -c{0} -n1 --exclusive runPipe.py --noconfirm --maxcpu {0} -e 3 -r {1} "${{inputFile}}" &\n'.format(Pstep,S)
                command3 += com
            if counter != len(bestPath):
                command3 += 'wait\n'
            counter += 1
            sampleNum += Sstep
        command4 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 4 "${inputFile}"'
        command5a = 'srun -N1 -c1 -n1 --exclusive runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 5 --edger "${inputFile}" &'
        command5b = 'srun -N1 -c1 -n1 --exclusive runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 5 --deseq "${inputFile}" &'
        command5 = command5a + '\n' + command5b
        Context = {
                "NODES": 2,
                "CPT": 16,
                "NTASKS": 5,
                "INPUT": self.inputPath,
                "JSON": JSFI,
                "STAGE12": command12,
                "STAGE3": command3,
                "STAGE4": command4,
                "STAGE5": command5,
                }
        batchScript = batch.format(**Context)
        with open('pipeBatch','w') as f:
            f.write(batchScript)

    def makeBatch(self, cluster):
        ''' Arguments:
                cluster = list; list of integers that describe number
                            of CPUs in each node of your cluster
            Returns:
                None

            Creates Pipeline batch script to be used with slurm
        '''
        batch =  """#!/bin/bash
#SBATCH --nodes={NODES}
#SBATCH --time=120
#SBATCH --cpus-per-task={CPT}
#SBATCH --ntasks={NTASKS}
#SBATCH --job-name="Pipeline"
#SBATCH --export=PATH,RNASEQDIR

inputFile='{INPUT}'
jsonFile='{JSON}'

# Stage 1 and 2
{STAGE12}

wait

# Stage 3
{STAGE3}
wait

# Stage 4
{STAGE4}

wait

# Stage 5
{STAGE5}

scontrol show job $SLURM_JOB_ID
wait
"""
        numSamps = self.getNumberofSamples()
        command12 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 1,2 "${inputFile}"'
        bestPath = self.getOptimal(cluster)
        command3,counter,sampleNum = '',1,1
        for path in bestPath:
            Pstep = bestPath[path]['Procs']
            Sstep = bestPath[path]['Samps']
            for S in range(sampleNum,sampleNum + Sstep):
                com = 'srun -N1 -c{0} -n1 --exclusive runPipe.py --noconfirm --maxcpu {0} -e 3 -r {1} "${{inputFile}}" &\n'.format(Pstep,S)
                command3 += com
            if counter != len(bestPath):
                command3 += 'wait\n'
            counter += 1
            sampleNum += Sstep
        command4 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 4 "${inputFile}"'
        command5a = 'srun -N1 -c1 -n1 --exclusive runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 5 --edger "${inputFile}" &'
        command5b = 'srun -N1 -c1 -n1 --exclusive runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 5 --deseq "${inputFile}" &'
        command5 = command5a + '\n' + command5b
        Context = {
                "NODES": len(cluster),
                "CPT": bestPath['Step 1']['Procs'],
                "NTASKS": bestPath['Step 1']['Samps'],
                "INPUT": self.inputPath,
                "JSON": JSFI,
                "STAGE12": command12,
                "STAGE3": command3,
                "STAGE4": command4,
                "STAGE5": command5,
                }
        batchScript = batch.format(**Context)
        with open('pipeBatch','w') as f:
            f.write(batchScript)

#    def makeBatchScript(self,cpuLimit=None,behavior='default'):
#        ''' Arguments:
#                None
#            Returns:
#                None
#
#            Creates Pipeline batch script to be used with slurm
#        '''
#        batch = """#!/bin/bash
##SBATCH --nodes=2
##SBATCH --time=120
##SBATCH --cpus-per-task={CPT}
##SBATCH --ntasks={NT}
##SBATCH --job-name="Pipeline"
##SBATCH --export=PATH,RNASEQDIR
#
#inputFile='{INPUT}'
#{COMMANDS}
#
#wait"""
#        numSamps = self.getNumberofSamples()
#        if behavior == 'default':
#            if cpuLimit == None:
#                com = 'srun -N1 -c48 -n1 --exclusive runPipe.py --noconfirm -e 3 -r {} "${{inputFile}}" &'
#                commands = '\n'.join([com.format(num) for num in range(1,numSamps+1)])
#                cpt = 48
#                nt = 2
#            else:
#                com = 'srun -N1 -c{0} -n1 --exclusive runPipe.py --noconfirm --maxcpu {0} -e 3 -r {1} "${{inputFile}}" &'
#                commands = '\n'.join([com.format(self.Procs,num) for num in range(1,numSamps+1)])
#                cpt = self.Procs
#                nt = 80//cpt
#        elif behavior == 'biox':
#            cpt = 16
#            nt = 5
#            com = 'srun -N1 -c16 -n1 --exclusive runPipe.py --noconfirm --maxcpu 16 -e 3 -r {0} "${{inputFile}}" &'
#            commands = '\n'.join([com.format(num) for num in range(1,numSamps+1)])
#        Context = {
#                "NUMSAMPLES": numSamps,
#                "INPUT": self.inputPath,
#                "COMMANDS": commands,
#                "CPT": cpt,
#                "NT": nt
#                }
#        batchScript = batch.format(**Context)
#        with open('pipeBatch','w') as f:
#            f.write(batchScript)
#
#    def makeRScript(self,cpuLimit=None):
#        ''' Arguments:
#                None
#            Returns:
#                None
#
#            Creates R batch script to be used with slurm
#        '''
#        batch = """#!/bin/bash
##SBATCH --nodes=1
##SBATCH --time=130
##SBATCH --cpus-per-task=1
##SBATCH --ntasks=2
##SBATCH --job-name="R-Analysis"
##SBATCH --export=PATH,RNASEQDIR
#
#inputFile='{INPUT}'
#jsonFile='{JSON}'
#{COMMAND1}
#{COMMAND2}
#{COMMAND3}
#
#wait"""
#        if cpuLimit == None:
#            command1 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 4 "${inputFile}"'
#            command2 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 5 --edger "${inputFile}" &'
#            command3 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 5 --deseq "${inputFile}" &'
#        else:
#            command1 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --maxcpu {} --jsonfile "${{jsonFile}}" --execute 4 "${{inputFile}}"'.format(self.Procs)
#            command2 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --maxcpu {} --jsonfile "${{jsonFile}}" --execute 5 --edger "${{inputFile}}" &'.format(self.Procs)
#            command3 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --maxcpu {} --jsonfile "${{jsonFile}}" --execute 5 --deseq "${{inputFile}}" &'.format(self.Procs)
#
#        Context = {
#                "INPUT": self.inputPath,
#                "JSON": JSFI,
#                "COMMAND1": command1,
#                "COMMAND2": command2,
#                "COMMAND3": command3,
#                }
#        batchScript = batch.format(**Context)
#        with open('rBatch','w') as f:
#            f.write(batchScript)
#
#    def makeMasterScript(self,cpuLimit=None):
#        ''' Arguments:
#                None
#            Returns:
#                None
#
#            Creates R batch script to be used with slurm
#        '''
#        batch = """#!/bin/bash
##SBATCH --nodes=1
##SBATCH --time=10
##SBATCH --cpus-per-task=1
##SBATCH --ntasks=1
##SBATCH --job-name="RNA-Seq"
##SBATCH --export=PATH,RNASEQDIR
#
#inputFile='{INPUT}'
#jsonFile='{JSON}'
#{COMMAND1}
#
#step1ID=$SLURM_JOB_ID
#sbatch -d afterok:$step1ID "{PWD}/pipeBatch"
#
#step2ID=$(( step1ID++ ))
#sbatch -d afterok:$step2ID "{PWD}/rBatch"
#
#wait"""
#        if cpuLimit == None:
#            command1 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --jsonfile "${jsonFile}" --execute 1,2 "${inputFile}"'
#        else:
#            command1 = 'srun -N1 -c1 -n1 runPipe.py --noconfirm --maxcpu {} --jsonfile "${{jsonFile}}" --execute 1,2 "${{inputFile}}"'.format(self.Procs)
#        cd = os.getcwd()
#        Context = {
#                "INPUT": self.inputPath,
#                "JSON": JSFI,
#                "COMMAND1": command1,
#                "PWD": cd
#                }
#        batchScript = batch.format(**Context)
#        with open('MasterBatch','w') as f:
#            f.write(batchScript)
#
#    def makeSlurm(self,cpuLimit,batch='default'):
#        checkJSON(JSFI,behavior='makeSlurm')
#        self.makeRScript(cpuLimit=None)
#        self.makeBatchScript(cpuLimit,batch)
#        self.makeMasterScript(cpuLimit=None)

    def findPipeFinish(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.findFinish() to figure out when self.GO()
            function has finished analysis on all samples
        '''
        pipeUtils.findFinish(self.Project, self.ogOriginal)

    def is3Finished(self):
        ''' Arguments:
                None
            Returns:
                boolean = True if finished, False if not finished

            Calls pipeUtils.findFinish() to figure out when self.GO()
            function has finished analysis on all samples
        '''
        return pipeUtils.findFinish(self.Project, self.ogOriginal, behavior='non-default')

    @funTime
    def createJsonMetadata(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.createMetaData() to make JSON file if not already
            provided
        '''
        print('Making Metadata file; require user input...')
        pipeUtils.createMetaData(self.Postprocessing)
        print('Waiting for Pipeline to finish...')

    @funTime
    def createNiceCounts(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.makeNiceCounts() to make the composite counts file
            with length column; Not used for my R analysis
        '''
        pipeUtils.makeNiceCounts(self.Postprocessing, self.Data)

    @funTime
    def getPipeTime(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.makeTotalTime() to scrape Sample times from
            respective sample Runtime.log files
        '''
        pipeUtils.makeTotalTime(self.Postprocessing, self.Data)

    def createRCounts(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.createCountFile() to make the Count file(Counts.dat)
            that works in my R analysis
        '''
        os.chdir(self.Postprocessing)
        pipeUtils.createCountFile(self.Postprocessing)

    def createRCols(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.createColumnFile() to make column file(Cols.dat)
            that describes Experiment. Uses Metadata JSON file for this
        '''
        os.chdir(self.Postprocessing)
        pipeUtils.createColumnFile(self.Postprocessing)

    @funTime
    def createRProgram(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.createRScript() and pipeUtils.createEdgeRScript()
            to make DESeq2 R analysis program and edgeR analysis program
        '''
        pipeUtils.createRScript(self.Postprocessing)
        pipeUtils.createEdgeRScript(self.Postprocessing)

    @funTime
    def runRProgram(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.makeEdgeRreport() and pipeUtils.makeRreports()
            to create DESeq2 and edgeR reports
        '''
        print("R program is running...")
        pipeUtils.makeEdgeRreport(self.Postprocessing)
        pipeUtils.makeRreports(self.Postprocessing)

    @funTime
    def runDESeq(self):
        pipeUtils.makeRreports(self.Postprocessing)
    @funTime
    def runEdgeR(self):
        pipeUtils.makeEdgeRreport(self.Postprocessing)

    @funTime
    def findRFinish(self):
        ''' Arguments:
                None
            Returns:
                None

            Calls pipeUtils.notifyEnding() to figure out when self.runRProgram()
            has finished making reports
        '''
        pipeUtils.notifyEnding(self.Postprocessing)

    def checkX11(self):
        ''' Arguments:
                None
            Returns:
                None

            Asks user if X server is available. Required for R stuff.
        '''
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
        ''' Arguments:
                None
            Returns:
                None

            Initializes notification folder that is used to determine when
            Pipeline has finished running
            'runPipeNotify/' in self.Project
        '''
        if os.path.isdir(self.Project + '/runPipeNotify') == False:
            os.mkdir(self.Project + '/runPipeNotify')
        else:
            shutil.rmtree(self.Project + '/runPipeNotify')
            os.mkdir(self.Project + '/runPipeNotify')

    def makeNotifyFolder2(self):
        ''' Arguments:
                None
            Returns:
                None

            Initializes notification folder that is used to determine when
            Pipeline has finished running
            'runPipeNotify/' in self.Project
        '''
        if not os.path.isdir(self.Project + '/runPipeNotify'):
            try:
                os.mkdir(self.Project + '/runPipeNotify')
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
        Check = [os.path.exists(sample+'/.done') for sample in glob.glob(self.Data + '/sample*')]
        if False in Check:
            return False
        else:
            return True

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

            Cleaning function for Experiment.
            See runPipe.py --help for help
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
        ''' Arguments:
                None
            Returns:
                None

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
        '''
        For a specific sample: run analyses
        Inherits from Experiment Parent Class
        '''

        def __init__(self,sampleNumber,inputFile, maxCPU=None):
            ''' Arguments:
                    sampleNumber = int; sample number used for naming
                    inputFile = str; path to inputFile
                Returns:
                    None

                Initializes class variables from Parent class and also
                initializes some sample specific variables
            '''
            Experiment.__init__(self,inputFile)
            assert sampleNumber <= Experiment.getNumberofSamples(self)
            assert sampleNumber > 0
            self.sampleNumber = sampleNumber
            self.sampleName = 'sample_{}'.format(str(sampleNumber).zfill(2))
            self.samplePath = '{}/{}'.format(self.Data, self.sampleName)
            self.Read1 = self.getReadNames()[0]
            self.Read2 = self.getReadNames()[1]
            self.logPath = '{}/Runtime.{}.log'.format(self.samplePath, self.sampleName)
            if maxCPU != None:
                self.Procs = maxCPU

        def __repr__(self):
            ''' Arguments:
                    None
                Returns:
                    String representation print
                    ex: >print(Sample)
                            Sample(sample_01)

                Returns string representation of class
            '''
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
            command1 = r'fastqc -t {0} -o {1} {2} {3}'.format(self.Procs,
                                                                fastqcFolder,
                                                                read1,
                                                                read2)
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
                    Log.write('Error: Unable to gather Overrepresented Sequences for {}\n'.format(self.sampleName))
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
            # Making Command
            command = r'java -jar $RNASEQDIR/Trimmomatic/trimmomatic-0.36.jar PE -threads {procs} -phred{phred} {Read1} {Read2} read1.P.trim.{fastq}.gz read1.U.trim.{fastq}.gz read2.P.trim.{fastq}.gz read2.U.trim.{fastq}.gz ILLUMINACLIP:$RNASEQDIR/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25'
            Context = {
                    "procs": self.Procs,
                    "phred": Phred,
                    "Read1": read1,
                    "Read2": read2,
                    "fastq": self.Fastq
                    }
            commandWithContext = command.format(**Context)
            goodCommand = self.formatCommand(commandWithContext)
            # Executing
            os.chdir(self.samplePath)
            subprocess.run(goodCommand,
                                shell=True,
                                check=True)

        #def smartQC(self):
        #    os.chdir(self.samplePath)
        #    run = 1
        #    while True:
        #        STOP = False
        #        self.runFastqc(run)
        #        self.findOverrepSeq(run)
        #        while True:
        #            print('Pipeline paused. You should review sample {}'.format(self.sampleName))
        #            answer1 = input("Are you ready to run Trimmomatic on {}?(y,n) ".format(self.sampleName))
        #            if answer1 == 'y':
        #                self.runTrimmomatic()
        #                run += 1
        #                break
        #            elif answer1 == 'n':
        #                answer2 = input("Would you like to proceed to Stage 2?(y,n) ")
        #                if answer2 == 'y':
        #                    break
        #                elif answer2 == 'n':
        #                    answer3 = input("Would you like to run QC again?(y,n) ")
        #                else:
        #                    print('Please answer y or n')
        #            else:
        #                print('Please answer y or n')
        #        if STOP == True:
        #            break

        def runPart1(self):
            ''' Arguments:
                    None
                Returns:
                    None

                Runs Fastqc once, then Trimmomatic, and then fastqc again
            '''
            #TODO Fix this to desired behavior
            os.chdir(self.samplePath)
            with open('{}/Runtime.{}.log'.format(self.samplePath, self.sampleName), 'w') as LOG:
                LOG.write('\t\tRuntime Log Part 1 for {}\n'.format(self.sampleName))
                LOG.write('----------------------------------------\n\n')
            self.runQCheck(1, self.Read1, self.Read2)
            self.runTrimmomatic(self.Read1, self.Read2)
            self.runQCheck(2, 'read1.P.trim.{}.gz'.format(self.Fastq),
                                'read2.P.trim.{}.gz'.format(self.Fastq))
            self.writeFunctionTail('runPart1')

        ########################################################
        # Pipeline
        ########################################################

        def runSeqtk(self):
            ''' Arguments:
                    None
                Returns:
                    None

                Runs seqtk on reads from trimmomatic
            '''
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
            ''' Arguments:
                    None
                Returns:
                    None

                Runs Blastn on seqtk output
            '''
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

        def runHisat(self):
            ''' Arguments:
                    None
                Returns:
                    None

                Runs hisat2 on data
            '''
            self.writeFunctionHeader('runHisat')
            strandedVar = self.findStranded()
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
            subprocess.run(goodCommand,
                                shell=True,
                                check=True)
            os.chdir(self.samplePath)
            os.remove('aligned.{}.sam'.format(self.sampleName))
            self.writeFunctionTail('runCompression')

        def runFeatureCounts(self):
            ''' Arguments:
                    None
                Returns:
                    None

                Collects counts from data with featureCounts
            '''
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
            subprocess.run(goodCommand,
                                shell=True,
                                check=True)
            self.writeFunctionTail('getAlignedColumn')

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
        # End of Sample class
        ########################################################

    ################################################################
    # Stage Run Functions
    ################################################################

    @funTime
    def runStage1(self,exists=False):
        ''' Note: does not need to be run more than once

        Prepare for Pipeline
        '''
        if noconfirm == False:
            self.checkX11()
        self.makeStructure()
        self.makeSyms(exists)

    @funTime
    def runStage2(self,exists=False):
        '''
        Prepare Reference Data
        '''
        if not exists:
            self.qcRef()
            self.ppRef()

    @funTime
    def runStage3(self):
        '''
        Run Pipeline
        '''
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

    @funTime
    def runStage4(self):
        ''' Note: Creating a metadata file should only need
            to be done once

        Prepare for DESeq2 Analysis
        '''
        while True:
            if self.is3Finished() or self.checkEach():
                time.sleep(10)
                print('Preparing for DESeq2...')
                shutil.rmtree(self.Project + '/runPipeNotify')
                self.gatherAllSampleOverrep(1)
                self.gatherAllSampleOverrep(2)
                self.createJsonMetadata()
                self.createNiceCounts()
                self.getPipeTime()
                self.createRCounts()
                self.createRCols()
                self.createRProgram()
                break
            else:
                time.sleep(30)

    @funTime
    def runStage5(self):
        '''
        Run R Analysis
        '''
        self.runRProgram()
        self.findRFinish()

    @funTime
    def runAll(self,exists=False):
        '''
        Runs Stage 1 through Stage 5
        '''
        self.runStage1(exists)
        self.runStage2(exists)
        self.runStage3()
        self.runStage4()
        self.runStage5()

    ################################################################
    # End of Experiment class
    ################################################################

#################################################################
if __name__ == '__main__':
    print('Please use runPipe.py\nSee runPipe --help')
    raise SystemExit
    #main()
