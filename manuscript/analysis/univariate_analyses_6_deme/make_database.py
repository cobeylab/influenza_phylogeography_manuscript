#!/usr/bin/env python

import os
import sys
import csv
import sqlite3
import json
import numpy as np
import subprocess
from collections import OrderedDict

csvFiles = ['meanFluxRate','antigenicLagAG1','trunkProportions','meanAnnualIncidence','meanTMRCA']

PARAMS = ['totalN', 'initialPrR', 'muPhenotype', 'meanStep', 'sdStep', 'relativeN', 'relativeTurnover', 'tropicFractionI0', 'relativeR0', 'seasonalAmplitude', 'betweenDemePro', 'initialTraitA','beta','smithConversion','sweepID']
PARAM_TYPES = {
    'totalN' : 'INTEGER',
    'initialPrR' : 'REAL',
    'muPhenotype' : 'REAL',
    'meanStep' : 'REAL',
    'sdStep' : 'REAL',
    'relativeN' : 'REAL',
    'relativeTurnover' : 'REAL',
    'tropicFractionI0' : 'REAL',
    'relativeR0' : 'REAL',
    'seasonalAmplitude' : 'REAL',
    'betweenDemePro' : 'REAL',
    'initialTraitA' : 'REAL',
    'beta' : 'REAL',
    'smithConversion' : 'REAL',
    'sweepID' : 'TEXT',
}

def insert(db, table, values):
    return db.execute(
        '''INSERT INTO {0} VALUES ({1})'''.format(
            table,
            ','.join(['?'] * len(values))
        ),
        values
    )
    
# load mathematica outputs from csv files
def loadCsvData(fileName, jobDir, runId, comboDb):
    csvFile = fileName + '.csv'
    with open(os.path.join(jobDir,csvFile)) as f:
        cw = list(csv.reader(f, delimiter=','))
        insert(comboDb,fileName, ([runId] + [float(x) for y in cw for x in y]))

def getData(fileName, jobDir, runId, comboDb):
    csvFile = fileName + '.csv'
    with open(os.path.join(jobDir,csvFile)) as f:
        cw = list(csv.reader(f, delimiter=','))
        dat = [float(x) for y in cw for x in y]
        if(fileName == 'trunkProportions'):
            return([dat[0],dat[1]])
        elif(fileName == 'antigenicLagAG1'):
            return([dat[1], dat[4]])
        else:
            return([dat[0]])

# Load mathematica outputs from SQLite
def loadMathSqlite(fileName, jobDir, runId, mathDb, comboDb):
    mathOut = mathDb.execute("SELECT * FROM {0}".format(fileName))
    cw = list(mathOut)
    insert(comboDb,fileName, ([runId] + [float(x) for y in cw for x in y]))

def loadJob(comboDb, resultsDir, runId, jobDir):
    # Load dictionary of parameter values from parameters.json file
    paramsJson = loadUncommentedJsonString(os.path.join(jobDir, 'parameters.json'))
    paramsDict = json.loads(paramsJson, object_pairs_hook=OrderedDict)
    paramVals = paramsDict
    
    insert(comboDb, 'parameters', [runId] + [paramVals[paramName] for paramName in PARAMS])
    
    try:
        with open(os.path.join(jobDir, 'out.summary')) as summaryFile:
            cw = csv.reader(summaryFile, delimiter='\t')
            cw.next()
            for key, value in cw:
                insert(comboDb,'summary',[runId, key, value])
    except:
        pass
        
    try:
        extinct = 0
        excessDiversity = 0
        fluLike = 0
        if os.path.isfile(os.path.join(jobDir, 'out.extinct')):
            extinct = 1
        if os.path.isfile(os.path.join(jobDir, 'out.tmrcaLimit')):
            excessDiversity = 1
        if os.path.isfile(os.path.join(jobDir, 'out.branches')):
            fluLike = 1
        insert(comboDb, 'status',[runId, extinct, excessDiversity,fluLike])
    except:
        pass    
    
    try:
        if ((not os.path.isfile(os.path.join(jobDir, 'out.extinct'))) and (not os.path.isfile(os.path.join(jobDir, 'out.tmrcaLimit'))) and os.path.isfile(os.path.join(jobDir, 'out.branches'))):
            for fileName in csvFiles:
                loadCsvData(fileName,jobDir,runId, comboDb)
            insert(comboDb, 'pooled_results', [runId] + [paramVals[paramName] for paramName in PARAMS] + [extinct, excessDiversity, fluLike] + [x for y in [getData(fileName,jobDir,runId, comboDb) for fileName in csvFiles] for x in y])
        else:
            insert(comboDb, 'pooled_results', [runId] + [paramVals[paramName] for paramName in PARAMS] + [extinct, excessDiversity, fluLike] + [None, None, None, None, None, None, None])
    except:
        pass
        
    comboDb.commit()


def loadJobs(comboDb, sweepDir, resultsDir):
    # Extract parameters names
    subDirs = os.listdir(resultsDir)
    runId = 0
    for subDir in subDirs:
        jobDir = os.path.join(resultsDir, subDir)
        runId = subDir
        if os.path.isdir(jobDir):
            loadJob(comboDb, resultsDir, runId, jobDir)

def loadUncommentedJsonString(filename):
    lines = list()
    with open(filename) as jsonFile:
        for line in jsonFile:
            lines.append(line.rstrip('\n').split('//')[0])
    return '\n'.join(lines)

if __name__ == "__main__":
    sweepDir = os.path.abspath(os.path.dirname(__file__))
    resultsDir = os.path.join(sweepDir, 'results')
    
    comboDb = sqlite3.connect(os.path.join(sweepDir, 'results.sqlite'))
    
    # Index of runs with parameters
    comboDb.execute(
        "CREATE TABLE parameters (runId INTEGER, {0})".format(
            ', '.join([paramName + ' ' + PARAM_TYPES[paramName] for paramName in PARAMS])
        )
    )
    
    # Database from summary file
    comboDb.execute("CREATE TABLE summary (runId INTEGER, key TEXT, value TEXT);")
    
    # Database from results
    comboDb.execute("CREATE TABLE status (runId INTEGER, extinct INTEGER, excessDiversity INTEGER, fluLike INTEGER);")
    comboDb.execute("CREATE TABLE meanTMRCA (runId INTEGER, meanTMRCA REAL);")
    comboDb.execute("CREATE TABLE meanFluxRate (runId INTEGER, meanFluxRate REAL);")
    comboDb.execute("CREATE TABLE antigenicLagAg1 (runId INTEGER, lowerTropics1 REAL, tropics1 REAL, upperTropics1 REAL, lowerTropics2 REAL, tropics2 REAL, upperTropics2 REAL, lowerNorth1 REAL, north1 REAL, upperNorth1 REAL, lowerNorth2 REAL, north2 REAL, upperNorth2 REAL, lowerSouth1 REAL, south1 REAL, upperSouth1 REAL, lowerSouth2 REAL, south2 REAL, upperSouth2 REAL);")
    comboDb.execute("CREATE TABLE trunkProportions (runId INTEGER, tropics1 REAL, tropics2 REAL, north1 REAL, north2 REAL, south1 REAL, south2 REAL);" )
    comboDb.execute("CREATE TABLE meanAnnualIncidence (runId INTEGER, meanAnnualIncidence REAL);")
    comboDb.execute(
        "CREATE TABLE pooled_results (runId INTEGER, {0}, extinct INTEGER, excessDiversity INTEGER, fluLike INTEGER, meanFluxRate REAL, tropics1AgLag REAL, tropics2AgLag REAL, tropics1TrunkPro REAL, tropics2TrunkPro REAL, meanAnnualIncidence REAL, meanTMRCA REAL);".format(
            ', '.join([paramName + ' ' + PARAM_TYPES[paramName] for paramName in PARAMS])
        )
    )
    
    comboDb.commit()
    
    loadJobs(comboDb, sweepDir, resultsDir)

    comboDb.close()
    