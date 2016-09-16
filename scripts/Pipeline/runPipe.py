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
import time
import pipeClasses
import os

################################################################
# Handling Command line arguments
################################################################

def main():
    arguments = docopt(__doc__, version='Version 1.0\nAuthor: Alberto')
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
            PROJ = pipeClasses.Experiment(arguments['<pathtoInput>'])
            PROJ.nukeProject()
            raise SystemExit
        else:
            PROJ = pipeClasses.Experiment(arguments['<pathtoInput>'])
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
        PROJ = pipeClasses.Experiment(arguments['<pathtoInput>'])
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
        PROJ = pipeClasses.Experiment(arguments['<pathtoInput>'])
        if os.path.isdir(str(PROJ.Data + '/' + arguments['--sampleclean'])) == False:
            print('{} does not exist'.format(str(PROJ.Data + '/' + arguments['--sampleclean'])))
            raise SystemExit
        PROJ.clean('Sample',arguments['--sampleclean'])
        raise SystemExit

    if JSFI != None:
        pipeClasses.checkJSON(JSFI)

    PROJ = pipeClasses.Experiment(arguments['<pathtoInput>'])

    global RUNTIMELOG
    if arguments['--runtime'] == '$Project':
        RUNTIMELOG = str(PROJ.Project + '/Runtime.log')
    else:
        if os.path.isdir(arguments['--runtime']) == False:
            print('Need a valid directory to save Runtime file to')
            raise SystemExit
        else:
            RUNTIMELOG = str(arguments['--runtime'] + '/Runtime.log')

    pipeClasses.RUNTIMELOG = RUNTIMELOG

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

################################################################

if __name__ == '__main__':
    main()
