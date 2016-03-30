#!/usr/bin/env python

import os
import sys
import csv
import sqlite3
import json
import numpy as np
import subprocess
from collections import OrderedDict

def summarize(resultsDb, parName1, parName2, parName3, parName4, parName5, parValue1, parValue2, parValue3, parValue4, parValue5):
	runIds = [row[0] for row in resultsDb.execute(
		"SELECT runId FROM parameters WHERE {0} = ? AND {1} = ? AND {2} = ? AND {3} = ? AND {4} = ?".format(parName1, parName2, parName3, parName4, parName5),
		[parValue1, parValue2, parValue3, parValue4, parValue5]
	)]
	
	for runId in runIds:
		try:
			northAgLag = float(resultsDb.execute('SELECT north FROM antigenicLagAg1 WHERE runId = ?', [runId]).next()[0])
			tropicsAgLag= float(resultsDb.execute('SELECT tropics FROM antigenicLagAg1 WHERE runId = ?', [runId]).next()[0])
			southAgLag = float(resultsDb.execute('SELECT south FROM antigenicLagAg1 WHERE runId = ?', [runId]).next()[0])
			northTrunkPro =	float(resultsDb.execute('SELECT north FROM trunkProportions WHERE runId = ?', [runId]).next()[0])
			tropicsTrunkPro = float(resultsDb.execute('SELECT tropics FROM trunkProportions WHERE runId = ?', [runId]).next()[0])
			southTrunkPro = float(resultsDb.execute('SELECT south FROM trunkProportions WHERE runId = ?', [runId]).next()[0])
			
			resultsDb.execute(
				'INSERT INTO pooled_results VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
				[runId, parValue1, parValue2, parValue3, parValue4, parValue5, northAgLag, tropicsAgLag, southAgLag, northTrunkPro, tropicsTrunkPro, southTrunkPro]
			)
		except:
			pass
	

	
if __name__ == '__main__':
	os.chdir(os.path.dirname(__file__))
	
	parName1 = "relativeN"
	parName2 = "tropicFractionI0"
	parName3 = "relativeR0"
	parName4 = "relativeTurnover"
	parName5 = "seasonalAmplitude"
	
	resultsDb = sqlite3.connect('results.sqlite')
	
	parameterSets = [row for row in resultsDb.execute(
		"SELECT DISTINCT {0}, {1}, {2}, {3}, {4} FROM parameters".format(parName1,parName2,parName3,parName4,parName5)
	)]
	
	# Create table containing various metrics for each rate, sd pair
	resultsDb.execute(
		"CREATE TABLE pooled_results ({0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11})".format('runId',parName1,parName2,parName3,parName4,parName5,"northAgLag","tropicsAgLag","southAgLag","northTrunkPro","tropicsTrunkPro","southTrunkPro")
	)
	
	for parValue1, parValue2, parValue3, parValue4, parValue5 in parameterSets:
		summarize(resultsDb, parName1, parName2, parName3, parName4, parName5, parValue1, parValue2, parValue3, parValue4, parValue5)
	
	resultsDb.commit()
	
