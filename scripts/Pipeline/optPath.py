#!/usr/bin/python3

'''Usage:
    lahontan fo [options] <numsamples> <cluster>

Options:
    -h, --help
        Show this screen and exit
    -t <filename, --tofile <filename>
        Name of file to be saved to
        [default: ./OptimalPath.dat]
    --maxcpu <CPUs>
        Limits number of CPUs used in calculating optimal path
        Default is to use all available CPUs
    -c, --customize
        If the calculation is too slow, can use this option
        to specify the path yourself

Arguments:
    numSamples  = int; number of samples in experiment
    cluster     = comma-separated list of ints; list of CPUs per
                       machine in cluster

Used to create optimal path of CPU usage for pipeline

Default is to write optimal path to ./OptimalPath.dat as JSON
format'''

################################################################
# Importations
################################################################

import multiprocessing
import sys
import json
from pprint import pprint
from functools import partial

################################################################

def ceiling(a,b):
    return -(-a//b)
def atatime(cpuPerSample,limits):
    return sum([x//cpuPerSample for x in limits])
def iterationsTilFinish(cpuPerSample, numSamp, limits):
    return ceiling(numSamp,atatime(cpuPerSample,limits))
def timePerSample(cpuPerSample):
    return 2221.85*cpuPerSample**-0.831626
def totalTime(cpuPerSample, numSamp, limits):
    ''' If your cpuPerSample is held constant for all iterations required'''
    return iterationsTilFinish(cpuPerSample,numSamp,limits)*timePerSample(cpuPerSample)
def oneIteration(cpuPerSample, numSamp, limits):
    return [numSamp-atatime(cpuPerSample,limits),[cpuPerSample,timePerSample(cpuPerSample)]]
def chokePoints(limits):
    chokes,final = [],[]
    for amount in [(n,atatime(n,limits)) for n in range(1,max(limits)+1)]:
        if amount[1] not in chokes:
            chokes.append(amount[1])
            final.append(amount)
        else:
            final.pop(-1)
            final.append(amount)
    #return [2,3,4,5,6,8,16,48]
    return [a[0] for a in final[1:]]
def allIteration(numSamp, limits):
    Total = [[numSamp-atatime(cpuPerSample,limits),[cpuPerSample,timePerSample(cpuPerSample)]] for cpuPerSample in chokePoints(limits)]
    Final = {}
    for test in Total:
        Final[str(test[1][0])] = test[1]
    return Final
def smartMultiprocessing(M,numSamp,limits):
    Paths,Finished,Best,counter = [],[1],9E9,0
    available = chokePoints(limits)
    calcIter = allIteration(numSamp, limits)
    while True:
        if counter == 0:
            ITER = oneIteration(M,numSamp,limits)
            if ITER[0] <= 0:
                test = [sum(i) for i in zip(*ITER[1:])][1]
                if test < Best:
                    Best = test
                    del Finished[0]
                    Finished.append(ITER[1:])
            else:
                Paths.append(ITER)
            counter += 1
        else:
            #if counter == 4:
            #    break
            for path in Paths:
                for n in available:
                    ITER = path[:]
                    ITER[0] = ITER[0] - atatime(n,limits)
                    ITER.append(calcIter[str(n)][:])
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

def scrapeResults(Best,numSamps,limits):
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

def ceiling(a,b):
    return -(-a//b)

def checkMaxCPU(arguments):
    # Returning a maximum CPU value if given
    if arguments['--maxcpu']:
        if arguments['--maxcpu'].isdigit():
            Max = int(arguments['--maxcpu'])
            if Max <= multiprocessing.cpu_count() and Max > 0:
                return Max
            else:
                raise SystemExit('--maxcpu greater than available CPUs')
        else:
            raise SystemExit('Invalid value to --maxcpu: {}'.format(
                        arguments['--maxcpu']))
    else:
        return multiprocessing.cpu_count()

def customizePath(numSamples, limits, filename):
    optimalPath = {}
    totalAvailableCPU = sum(limits)
    print("Type 0 to exit")
    print("Total available CPU: {}\n".format(totalAvailableCPU))
    remainingSamples = numSamples
    step = 1
    while True:
        optimalPath['Step {}'.format(step)] = {}
        print('Step {}:\nRemaining Samples: {}\n'.format(step, remainingSamples))
        s = int(input('Please specify a number of samples to run for Step {}: '.format(step)))
        if s == 0:
            raise SystemExit('Exiting...')
        elif s > remainingSamples:
            print('Given samples is more than remaining\n')
            continue
        p = int(input('Please specify the number of CPU each sample should use for Step {}: '.format(step)))
        if p == 0:
            raise SystemExit('Exiting...')
        elif p*s > totalAvailableCPU:
            print('Too many processors given: {}*{} > Total Available CPU'.format(p,s))
            continue
        elif p > min(limits):
            print('Cannot use {} processors on smallest node'.format(p))
            continue
        print('\n')
        optimalPath['Step {}'.format(step)]['Procs'] = p
        optimalPath['Step {}'.format(step)]['Samps'] = s
        step += 1
        remainingSamples -= s
        if remainingSamples == 0:
            print('Customized Path:')
            pprint(optimalPath)
            print('Saving to {}'.format(filename))
            with open(filename,'w') as JF:
                json.dump(optimalPath, JF, sort_keys=True, indent=4)
            break
    print('Optimal path saved')
    
def main(arguments):
    numSamples = int(arguments['<numsamples>'])
    cluster = [int(machine) for machine in arguments['<cluster>'].split(',')]
    procs = checkMaxCPU(arguments)
    if not arguments['--customize']:
        #pprint(optimize(numSamples, cluster, procs))
        available = chokePoints(cluster)
        with multiprocessing.Pool(procs) as p:
            Results = p.map(partial(smartMultiprocessing,numSamp=numSamples,limits=cluster), 
                    available)
        sortedResults = list(sorted(Results))
        bestResult = sortedResults[0]
        prettyResult = scrapeResults(bestResult, numSamps=numSamples, limits=cluster)
        print('Optimal Path:')
        pprint(prettyResult)
        print('Saving to {}'.format(arguments['--tofile']))
        with open(arguments['--tofile'],'w') as JF:
            json.dump(prettyResult, JF, sort_keys=True, indent=4)
        print('Optimal path saved')
    else:
        customizePath(numSamples, cluster, arguments['--tofile'])

if __name__ == '__main__':
    argument = docopt(__doc__)
    main(argument)
