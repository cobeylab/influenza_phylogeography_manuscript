#!/usr/bin/env python

# jobs are indexed staring from 0

import csv
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
MUPHENOTYPES = np.array([5e-5,1e-4,2e-4])
MEANSTEPS = np.linspace(0.1,0.9,9)
SDSTEPS = np.linspace(0.1,0.6,6)

NUM_RUNS = 10
totalJobs = NUM_RUNS*len(MUPHENOTYPES)*len(MEANSTEPS)*len(SDSTEPS)
if totalJobs < maxJobs:
    maxJobs = totalJobs
N_PER_JOB = totalJobs/maxJobs+1

totalN = 45000000
birthRate = 0.000091
Reff = 1.0
R0 = 1.8
nu = 0.2
beta = R0 * (nu+birthRate)
initialTraitA = -4.0
smithConversion = 0.07
muPhenotype = 0.0001
betweenDemePro = 0.001
deltaT = 0.1
startAtEquilibriumInfected = [False,False,False]
startAtEquilibriumImmune = [False,False,False]

relativeN = 1.0
tropicFractionI0 = 1.0
relativeR0 = 1.0
relativeTurnover = 1.0
seasonalAmplitude = 0.1

def generateJobs():    
    jobNum = 0

    for muPhenotype in MUPHENOTYPES:
        for meanStep in MEANSTEPS:
            for sdStep in SDSTEPS:
                birthRates = [birthRate,birthRate*relativeTurnover,birthRate]
                deathRates = birthRates
                demeAmplitudes = [seasonalAmplitude,0,seasonalAmplitude]
            
                demeBaselines = [1,relativeR0,1]
            
                #start tropics at equilibrium
                tropicsR0 = relativeR0 * R0
                tropicsBeta = tropicsR0 * (nu + birthRates[2])
        
                initialPrS = Reff/(tropicsR0)

                initialPrI = birthRates[2]/tropicsBeta*((tropicsR0)-1.0)
                totalInitialIs = int(initialPrI * totalN)

                initialPrR = (1.0 - initialPrS - initialPrI)/(1.0-smithConversion*abs(initialTraitA))
                initialIs = [0,totalInitialIs,0]
                initialNs = [int(totalN/(2.0 + relativeN)), int(relativeN*(totalN/(2.0+relativeN))),int(totalN/(2.0 + relativeN))]
                initialIs = [int((totalInitialIs - totalInitialIs*tropicFractionI0)/2), int(totalInitialIs*tropicFractionI0), int((totalInitialIs - totalInitialIs*tropicFractionI0)/2)]


        
                for runNum in range(NUM_RUNS):
                    # This is the name that SLURM uses to identify the job
                    jobName = jobNum

                    # This is the subdirectory inside 'results' used to run the job
                    jobSubdir = '{}'.format(jobName)

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
                        ('startAtEquilibriumImmune',startAtEquilibriumImmune),
                        ('startAtEquilibriumInfected',startAtEquilibriumInfected),
                        ('R0',R0),
                        ('initialTraitA',initialTraitA),
                        ('beta',beta),
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
    
    for jobName, jobSubdir, paramDict in generateJobs():
        allParamsDict = OrderedDict()
        allParamsDict.update(cParamsDict)
        allParamsDict.update(paramDict)
        writeParameters(resultsDir, jobSubdir, allParamsDict)
        
    submitJob(rootDir, sweepDir, resultsDir)
