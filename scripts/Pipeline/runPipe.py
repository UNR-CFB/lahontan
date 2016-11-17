#!/usr/bin/python3

'''Usage:
    runPipe.py [options] <pathtoInputFile>

Options:
    -h, --help
        Show this screen and exit
    -j <jsonFile>, --jsonfile <jsonFile>
        Ignores JSON file creation and uses specified path
        to JSON
    -c <placeToClean>, --clean <placeToClean>
        Cleans <placeToClean>; Possible places include:
            Reference, Data, Postprocessing, All
    -s <sampleName>, --sampleclean <sampleName>
        Similar to --clean; but instead just cleans a
        single sample directory, <sampleName>
    -e <stage>, --execute <stage>
        Comma-separated list of stages to be executed.
        Possible stages include:
            1: Creating Project Structure
            2: Preparing Reference Data
            3: Running actual Pipeline
            4: Preparing for R analysis
            5: Running R analysis
            A: (1,2,3,4,5); A=all i.e. runs entire pipeline
        [default: A]
    -r <integer>, --runsample <integer>
        Runs Stage 3 of the pipeline on the sample specified
        by the integer
    --maxcpu <CPUs>
        Limits number of CPUs used by Pipeline. Without
        argument, default is to use all available CPUs
    --reference-qc <pathtoReferences>
        Runs Quality Control check on Reference files
    --reference-pp <pathtoReferences>
        Pre processes Reference data
    --makebatch <cluster>
        Makes batch file to be used with slurm. The argument
        it takes is a comma-separated list of CPUs on each
        node in your cluster
    --makebatchbiox
        Modifies behavior of --makebatch. Makes batch file
        with best behavior for our cluster(compute-1,compute-2)
    --noconfirm
        Ignore all user prompts except JSON file creation
    --NUKE
        Removes entire project Directory
    --edger
        Runs edgeR analysis only. Default is to run both 
    --deseq
        Runs DESeq2 analysis only. Default is to run both

Examples:
    runPipe.py /path/to/INPUT
        This will run the pipeline using the input variables
        from /path/to/INPUT with all the default behavior

    runPipe.py --execute 1,2 /path/to/INPUT
        This will run the first and second stage of the
        pipeline using the input variables from /path/to/INPUT

    runPipe.py --runsample 1,2 /path/to/INPUT
        This will run third stage on the first and second
        sample. Note that it will also run first and second
        stage if needed

    runPipe.py --jsonfile /path/to/Metadata /path/to/INPUT
        This will run the pipeline using the input variables
        from /path/to/INPUTfile and will use the JSON located
        at /path/to/Metadata for DESeq2 and edgeR analysis

    runPipe.py --clean Data /path/to/INPUT
        This will wipe all of the pipeline's actions on all
        samples within the respective Project Directory
        leaving only the symlinks that were created

    runPipe.py --NUKE /path/to/INPUT
        This will remove entire Project directory denoted by
        the INPUT variable "Project".
        Note: Will leave original Reference and Data folders
        as they were originally made
'''

################################################################
# Importations
################################################################

from docopt import docopt
from timeit import default_timer as timer
import time
import pipeClasses
import os

################################################################
# Handling Command line arguments
################################################################

def main(arguments):
    ''' Arguments:                   
            None                     
        Returns:                     
            None                     
                                     
        Handles arguments from docopt
        Creates Experiment Class     
        Runs Experiment methods      
        Tip: search for @@@          
    '''                              
    t1 = timer()

    global noconfirm
    noconfirm = arguments['--noconfirm']
    global JSFI
    JSFI = arguments['--jsonfile']

    pipeClasses.JSFI = JSFI
    pipeClasses.noconfirm = noconfirm
    pipeClasses.pipeUtils.JSFI = JSFI
    pipeClasses.pipeUtils.noconfirm = noconfirm

    # Handling --NUKE argument
    if arguments['--NUKE'] == True:
        if noconfirm == True:
            PROJ = pipeClasses.Experiment(arguments['<pathtoInputFile>'])
            PROJ.nukeProject()
            raise SystemExit
        else:
            PROJ = pipeClasses.Experiment(arguments['<pathtoInputFile>'])
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
        PROJ = pipeClasses.Experiment(arguments['<pathtoInputFile>'])
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
        PROJ = pipeClasses.Experiment(arguments['<pathtoInputFile>'])
        if os.path.isdir(str(PROJ.Data + '/' + arguments['--sampleclean'])) == False:
            print('{} does not exist'.format(str(PROJ.Data + '/' + arguments['--sampleclean'])))
            raise SystemExit
        PROJ.clean('Sample',arguments['--sampleclean'])
        raise SystemExit

    if JSFI != None:
        pipeClasses.checkJSON(JSFI)

    def makeBatch(ExperimentClass):
        ''' Arguments:
                ExperimentClass = class; experiment to run analysis on
            Returns:
                None

            If --makebatch argument given, then make batch script to be
            used with slurm, and then exit
        '''
        if arguments['--makebatchbiox']:
            ExperimentClass.makeBatchBiox()
            raise SystemExit('Batch file successfully created:\n\t{0}/pipeBatch'.format(os.getcwd()))
        elif arguments['--makebatch'] != None:
            preNodes,Nodes = arguments['--makebatch'].split(','),[]
            for node in preNodes:
                if not node.isdigit() or int(node) <= 0:
                    raise SystemExit('Not a valid argument to --makebatch: {}'.format(
                                            arguments['--makebatch']))
                Nodes.append(int(node))
            ExperimentClass.makeBatch(cluster=Nodes)
            raise SystemExit('Batch file successfully created:\n\t{0}/pipeBatch'.format(os.getcwd()))

    def checkMaxCPU():
        if arguments['--maxcpu'] == None:
            return None
        else:
            if arguments['--maxcpu'].isdigit():
                Max = int(arguments['--maxcpu'])
                if Max <= os.cpu_count() and Max > 0:
                    return Max
                else:
                    raise SystemExit('--maxcpu greater than available CPUs')
            else:
                raise SystemExit('Invalid value to --maxcpu: {}'.format(
                            arguments['--maxcpu']))


    PROJ = pipeClasses.Experiment(arguments['<pathtoInputFile>'], maxCPU=checkMaxCPU())

    global RUNTIMELOG
    RUNTIMELOG = str(PROJ.Project + '/Runtime.log')

    pipeClasses.RUNTIMELOG = RUNTIMELOG

    def runR(ExperimentClass):
        ''' Arguments:
                ExperimentClass = class; experiment to run analysis on
            Returns:
                None

            Runs R analysis, if:
                --edger  :   Then only edgeR analysis done
                --deseq  :   Then only DESeq2 analysis done
                neither argument is provided then both are done
        '''
        if arguments['--edger'] or arguments['--deseq']:
            if arguments['--edger']:
                ExperimentClass.runEdgeR()
            if arguments['--deseq']:
                ExperimentClass.runDESeq()
        else:
            ExperimentClass.runStage5()

    def checkdashr():
        if arguments['--runsample'] == None:
            return None
        else:
            samples = arguments['--runsample'].split(',')
            possibleSamples = [str(a+1) for a in range(PROJ.getNumberofSamples())]
            if len(samples) == 0:
                print('Argument not valid: {}'.format(
                        arguments['--runsample']))
                print('Possible arguments for --runsample: {}'.format(
                        ','.join(possibleSamples)))
                raise SystemExit
            if len(set(samples)) != len(samples):
                print('Argument not valid: {}'.format(
                        arguments['--runsample']))
                print('Possible arguments for --runsample: {}'.format(
                        ','.join(possibleSamples)))
                raise SystemExit('No repeats allowed')
            for sample in samples:
                if sample not in possibleSamples:
                    print('Argument not valid: {}'.format(
                            arguments['--runsample']))
                    print('Possible arguments for --runsample: {}'.format(
                            ','.join(possibleSamples)))
                    raise SystemExit('Argument not possible')
            return samples

    def checkdashe():
        executionStages = arguments['--execute'].split(',')
        possibleStages = ["1","2","3","4","5","A"]
        if len(executionStages) == 0:
            print('Argument not valid: {}'.format(
                    arguments['--execute']))
            print('Possible arguments for --execute: {}'.format(
                    ','.join(possibleStages)))
            raise SystemExit
        if len(set(executionStages)) != len(executionStages):
            print('Argument not valid: {}'.format(
                    arguments['--execute']))
            print('Possible arguments for --execute: {}'.format(
                    ','.join(possibleStages)))
            raise SystemExit('No repeats allowed')
        for stage in executionStages:
            if stage not in possibleStages:
                print('Argument not valid: {}'.format(
                        arguments['--execute']))
                print('Possible arguments for --execute: {}'.format(
                        ','.join(possibleStages)))
                raise SystemExit
        return executionStages

    def executeProject(ExperimentClass):
        ''' Arguments:
                ExperimentClass = class; experiment to run analysis on
            Returns:
                None

            Runs Experiment methods based on --execute,--runsample or
            just entire pipeline if not specified
        '''
        makeBatch(ExperimentClass)
        executionSamples = checkdashr()
        executionStages = checkdashe()

        if executionSamples == None:
            if 'A' in executionStages:
                ExperimentClass.runAll()
            else:
                for stageNumber in executionStages:
                    if str(stageNumber) == '5':
                        runR(ExperimentClass)
                    else:
                        exec('{}.runStage{}()'.format('ExperimentClass',stageNumber))
        else:
            if not os.path.exists(ExperimentClass.Data):
                ExperimentClass.runStage1()
            if not os.path.exists(ExperimentClass.Reference + '/{}'.format(
                                                        'Reference_Report.txt')):
                ExperimentClass.runStage2()
            for sample in executionSamples:
                ExperimentClass.executeSample(int(sample))
            if ExperimentClass.is3Finished():
                if 'A' in executionStages:
                    ExperimentClass.runStage4()
                    runR(ExperimentClass)
                else:
                    if '4' in executionStages:
                        ExperimentClass.runStage4()
                    if '5' in executionStages:
                        runR(ExperimentClass)
            #print('Finished Running specified samples\n' + 
            #'If you wish to run R analysis, you will need to run all samples\n' +
            #'If you have run all samples, you can run R analysis by running:\n' +
            #'runPipe.py --execute 4,5 /path/to/INPUT/file\n' +
            #'along with any other arguments you wish')

    ############################################################
    # Call for actually doing stuff
    ############################################################
    # @@@
    executeProject(PROJ)
    ############################################################

    t2 = timer()

    timeused = str(time.strftime('%H:%M:%S', time.gmtime(t2-t1)))
    print('Total time elapsed: {}'.format(timeused))

def side(arguments):
    ''' Arguments:                   
            None                     
        Returns:                     
            None                     

        Runs Alternate Analysis
    '''                              
    if arguments['--reference-qc'] != None:
        pipeClasses.pipeUtils.qcReference(arguments['--reference-qc'],pipeClasses.pipeUtils.getReferenceVariables(arguments['--reference-qc'])[2])
        raise SystemExit
    if arguments['--reference-pp'] != None:
        Gtf,Cdna,Genome = pipeClasses.pipeUtils.getReferenceVariables(arguments['--reference-pp'])
        Basename = pipeClasses.pipeUtils.getBasename(Genome)
        pipeClasses.pipeUtils.preProcessingReference(arguments['--reference-pp'],Cdna,Gtf,Genome,Basename)
        raise SystemExit

def MainVsSide(arguments):
    if arguments['--reference-pp'] != None:
        return 'side'
    if arguments['--reference-qc'] != None:
        return 'side'
    return 'main'

################################################################
if __name__ == '__main__':
    argument = docopt(__doc__, version='Version 0.99\nAuthor: Alberto')
    if MainVsSide(argument) == 'main':
        main(argument)
    else:
        side(argument)
