#!/usr/bin/env python

# jobs are indexed staring from 0

import os
import json
import numpy as np
import subprocess
from collections import OrderedDict

###############################################################################

### Modify these values and the generateJobs() function ###

# Set DRY = False to actually submit jobs
DRY = False

maxJobs = 500

sweepParamDict = OrderedDict([
    ('relativeN',(0.5,1.5)),
    ('tropicFractionI0',(0,1)),
    ('relativeR0',(0.75,1.25)),
    ('relativeTurnover',(0.5,2)),
    ('seasonalAmplitude',(0,0.25)),
])
steps = 10

NUM_RUNS = 20
totalJobs = NUM_RUNS*len(sweepParamDict)*steps
N_PER_JOB = totalJobs/maxJobs

totalN = 45000000
birthRate = 0.000091
Reff = 1.0
R0 = 1.8
nu = 0.2
beta = R0*nu
initialTraitA = -4.0
smithConversion = 0.07
muPhenotype = 0.0001
betweenDemePro = 0.001
deltaT = 0.1
meanStep = 0.6
sdStep = 0.3

def generateJobs(sweepParamDict, steps):

    initialPrS = Reff/R0

    initialPrI = birthRate/beta*(R0-1.0)
    totalInitialIs = int(initialPrI * totalN)


    initialPrR = (1.0 - initialPrS - initialPrI)/(1.0-smithConversion*abs(initialTraitA))
    initialIs = [totalInitialIs/2, totalInitialIs/2, 0,0,0,0]

    meanStep = 0.6
    sdStep = 0.3
    
    jobNum = 0

    for paramName in sweepParamDict:
        # set default sweep parameters
        relativeN = 1.0
        tropicFractionI0 = 1.0
        relativeR0 = 1.0
        relativeTurnover = 1.0
        seasonalAmplitude = 0.1
        sweepID = paramName

        for i in np.linspace(sweepParamDict[paramName][0],sweepParamDict[paramName][1],steps):
            if paramName=='relativeN':
                relativeN = i
            if paramName=='tropicFractionI0':
                tropicFractionI0 = i
            if paramName=='relativeR0':
                relativeR0 = i
            if paramName=='relativeTurnover':
                relativeTurnover = i
            if paramName=='seasonalAmplitude':
                seasonalAmplitude = i
                
            birthRates = [birthRate*relativeTurnover, birthRate*relativeTurnover,
                            birthRate,birthRate,
                            birthRate,birthRate]
                            
            deathRates = birthRates
            demeAmplitudes = [0,0,seasonalAmplitude,seasonalAmplitude,seasonalAmplitude,seasonalAmplitude]
            
            demeBaselines = [relativeR0, relativeR0, 1,1,1,1]
            
            #start tropics at equilibrium
            tropicsR0 = relativeR0 * R0
            tropicsBeta = tropicsR0 * (nu + birthRates[0])
        
            initialPrS = Reff/(tropicsR0)

            initialPrI = birthRates[2]/tropicsBeta*((tropicsR0)-1.0)
            totalInitialIs = int(initialPrI * totalN)

            initialPrR = (1.0 - initialPrS - initialPrI)/(1.0-smithConversion*abs(initialTraitA))
            
            initialNs = [int(relativeN*(totalN/(2.0+relativeN))/2),
                        int(relativeN*(totalN/(2.0+relativeN))/2),
                        int(totalN/(2.0 + relativeN)/2), 
                        int(totalN/(2.0 + relativeN)/2), 
                        int(totalN/(2.0 + relativeN)/2), 
                        int(totalN/(2.0 + relativeN)/2)]
                        
            initialIs = [int(totalInitialIs*tropicFractionI0/2), 
                        int(totalInitialIs*tropicFractionI0/2), 
                        int((totalInitialIs - totalInitialIs*tropicFractionI0)/2), 
                        int((totalInitialIs - totalInitialIs*tropicFractionI0)/2),
                        int((totalInitialIs - totalInitialIs*tropicFractionI0)/2),
                        int((totalInitialIs - totalInitialIs*tropicFractionI0)/2)]

            for runNum in range(NUM_RUNS):
                # This is the name that SLURM uses to identify the job
                jobName = jobNum

                # This is the subdirectory inside 'results' used to run the job
                jobSubdir = '{}'.format(jobName)

                # This dictionary is added to constant_parameters.json and written to
                # the parameters file

                paramDict = OrderedDict([
                    ('smithConversion',smithConversion),
                    ('meanStep',meanStep),
                    ('sdStep',sdStep),
                    ('initialIs',initialIs),
                    ('initialPrR',initialPrR),
                    ('totalN',totalN),
                    ('birthRate',birthRates),
                    ('deathRate',deathRates),
                    ('muPhenotype',muPhenotype),
                    ('seasonalAmplitude', seasonalAmplitude),
                    ('demeAmplitudes', demeAmplitudes),
                    ('demeBaselines', demeBaselines),
                    ('initialNs',initialNs),
                    ('relativeN', relativeN),
                    ('tropicFractionI0', tropicFractionI0),
                    ('relativeR0', relativeR0),
                    ('relativeTurnover', relativeTurnover),
                    ('betweenDemePro',betweenDemePro),
                    ('deltaT',deltaT),
                    ('sweepID',sweepID),
                ])
                yield (jobName, jobSubdir, paramDict)
                jobNum += 1
        
        


###############################################################################

def writeParameters(resultsDir, jobSubdir, paramDict):
    jobDir = os.path.join(resultsDir, jobSubdir)
    os.makedirs(jobDir)
    paramsFilename = os.path.join(jobDir, 'parameters.json')
    with open(paramsFilename, 'w') as paramsFile:
        json.dump(paramDict, paramsFile, indent=2)
        paramsFile.write('\n')

def submitJob(rootDir, sweepDir, resultsDir):
    # Construct SLURM command
    sbatchFilename = os.path.join(sweepDir, 'job.sbatch')
    submitCommand = [
        'sbatch',
        '-J{0}'.format('Antigen'),
        '-D{0}'.format(resultsDir),
        '--array=0-{0}'.format(maxJobs-1),
        sbatchFilename
    ]
    
    # Print command to terminal
    process = subprocess.Popen(['echo'] + submitCommand)
    process.wait()
    
    # Construct environment variables
    env = dict(os.environ)
    env['ANTIGEN_ROOT'] = rootDir
    env['N_PER_JOB'] = str(N_PER_JOB)
    env['SWEEP_DIR'] = sweepDir
    
    # Actually run command
    if not DRY:
        process = subprocess.Popen(submitCommand, env=env)
        process.wait()

def loadUncommentedJsonString(filename):
    lines = list()
    with open(filename) as jsonFile:
        for line in jsonFile:
            lines.append(line.rstrip('\n').split('//')[0])
    return '\n'.join(lines)

if __name__ == '__main__':
    sweepDir = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    rootDir = os.path.abspath(os.path.join(sweepDir, '..'))
    resultsDir = os.path.join(sweepDir, 'results')
    
    cParamsFilename = os.path.join(sweepDir, 'constant_parameters.json')
    cParamsJson = loadUncommentedJsonString(cParamsFilename)
    cParamsDict = json.loads(cParamsJson, object_pairs_hook=OrderedDict)
    
    for jobName, jobSubdir, paramDict in generateJobs(sweepParamDict, steps):
        allParamsDict = OrderedDict()
        allParamsDict.update(cParamsDict)
        allParamsDict.update(paramDict)
        writeParameters(resultsDir, jobSubdir, allParamsDict)
        
    submitJob(rootDir, sweepDir, resultsDir)
