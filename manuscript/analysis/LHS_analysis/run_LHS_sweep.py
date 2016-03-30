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
LHSsize = 500

NUM_RUNS = 20
totalJobs = NUM_RUNS*LHSsize
if totalJobs < maxJobs:
    maxJobs = totalJobs + 1
N_PER_JOB = totalJobs/maxJobs + 1

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

def generateJobs():

    
    jobNum = 0

    thisDir = os.path.abspath(os.path.dirname(__file__))
    LHSparams = csv.DictReader(open(os.path.join(thisDir,'LHSparams.csv')))    
    for row in LHSparams:
                
        relativeN = float(row['relativeN'])
        tropicFractionI0 = float(row['tropicFractionI0'])
        relativeR0 = float(row['relativeR0'])
        relativeTurnover = float(row['relativeTurnover'])
        seasonalAmplitude = float(row['seasonalAmplitude'])
        demeBaselines = [1,relativeR0,1]
        
        initialPrS = Reff/(R0*relativeR0)

        initialPrI = birthRate/beta*((R0*relativeR0)-1.0)
        totalInitialIs = int(initialPrI * totalN)


        initialPrI = totalInitialIs/totalN ##birthRate/beta*(R0-1.0)

        initialPrR = (1.0 - initialPrS - initialPrI)/(1.0-smithConversion*abs(initialTraitA))
        initialIs = [0,totalInitialIs,0]
        initialNs = [int(totalN/(2.0 + relativeN)), int(relativeN*(totalN/(2.0+relativeN))),int(totalN/(2.0 + relativeN))]
        initialIs = [int((totalInitialIs - totalInitialIs*tropicFractionI0)/2), int(totalInitialIs*tropicFractionI0), int((totalInitialIs - totalInitialIs*tropicFractionI0)/2)]

    
        birthRates = [birthRate,birthRate*relativeTurnover,birthRate]
        deathRates = birthRates
        demeAmplitudes = [seasonalAmplitude,0,seasonalAmplitude]

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
