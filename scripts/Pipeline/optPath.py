#!/usr/bin/python3

'''
Used to create optimal path of CPU usage for pipeline

Takes three arguments in order:
    numSamples = int; number of samples in experiment
    procs = int; number of processors to use to calculate
                optimal path
    cluster = comma-separated list of ints; list of CPUs per
                machine in cluster

Writes optimal path to ./OptimalPath.dat as JSON format
'''

################################################################
# Importations
################################################################

from math import ceil as ceiling
import multiprocessing
import sys
import json
from pprint import pprint

################################################################

numSamples = int(sys.argv[1])
procs = int(sys.argv[2])
cluster = [int(a) for a in sys.argv[3].split(',')]
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
    return [2,3,4,5,6,8,16,48]
    #return [a[0] for a in final[1:]]
def allIteration(numSamp=numSamples, limits=cluster):
    Total = [[numSamp-atatime(cpuPerSample,limits),[cpuPerSample,timePerSample(cpuPerSample)]] for cpuPerSample in chokePoints(limits)]
    Final = {}
    for test in Total:
        Final[str(test[1][0])] = test[1]
    return Final
def smartMultiprocessing2(M,numSamp=numSamples,limits=cluster):
    Paths,Finished,Best,counter = [],[1],9E9,0
    available = chokePoints(limits)
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
def smartMultiprocessing(M,numSamp=numSamples,limits=cluster):
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
            if counter == 8:
                break
            for path in Paths:
                for n in available:
                    ITER = path[:]
                    ITER[0] = ITER[0] - atatime(n)
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
def scrapeResults(Best,numSamps=numSamples,limits=cluster):
    resultDict,counter,totalSamps = {},0,numSamps
    while True:
        print(Best)
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
    with open('OptimalPath.dat','w') as JF:
        json.dump(resultDict, JF, sort_keys=True, indent=4)
    return resultDict
with multiprocessing.Pool(procs) as p:
    available = chokePoints(cluster)
    Results = p.map(smartMultiprocessing, available)
sortedResults = list(sorted(Results))
bestResult = sortedResults[0]
print(scrapeResults(bestResult))
