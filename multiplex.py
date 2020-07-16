#!/usr/bin/env python

# multiplexJSON.py #################
#
# A "multiplexer" to combine many JSON files (like for different SRs) and outputs
# a single JSON with the "best" SR written out. Default figure of merit is expected CLs
# but is configurable.
#
# By: Larry Lee - Jan 2018

import json
import argparse
import glob
import os,sys
import pandas as pd

def main(args,inputDataList):


	# Print out the settings
	for setting in dir(args):
		if not setting[0]=="_":
			pass

	databaseList = []

	databaseListTheoryUp = []
	databaseListTheoryDown = []
	databaseListUpperLimit = []

	for (fID,data) in inputDataList:
		if not data:
			continue
		df = pd.DataFrame(data, columns=data[0].keys())
		df['fID'] = fID
		databaseList.append(df)


		if "Nominal" in fID and not args.ignoreTheory:

			try:
				with open(filename.replace("Nominal","Up")) as data_file:
					data = json.load(data_file)
					df = pd.DataFrame(data, columns=data[0].keys())
					df['fID'] = filename
					databaseListTheoryUp.append(df)
			except:
				pass
			try:
				with open(filename.replace("Nominal","Down")) as data_file:
					pass
					data = json.load(data_file)
					df = pd.DataFrame(data, columns=data[0].keys())
					df['fID'] = filename
					databaseListTheoryDown.append(df)
			except:
				print(">>> Can't find %s"%filename.replace("Nominal","Down"))

		if "Nominal" in fID and not args.ignoreUL:

			try:
				with open(filename.replace("fixSigXSecNominal_hypotest","upperlimit")) as data_file:
					print(">>> >>> Adding input file for upper limits")
					data = json.load(data_file)
					df = pd.DataFrame(data, columns=data[0].keys())
					df['fID'] = filename
					databaseListUpperLimit.append(df)
			except:
				pass

	database = pd.concat(databaseList, ignore_index=True)
	if len(databaseListTheoryUp):
		databaseTheoryUp   = pd.concat(databaseListTheoryUp, ignore_index=True)
		databaseTheoryDown = pd.concat(databaseListTheoryDown, ignore_index=True)
	if len(databaseListUpperLimit):
		databaseUpperLimit = pd.concat(databaseListUpperLimit, ignore_index=True)

	if args.debug:
		print(">>> Full database has length: %d"%len(database))

	# cleaning up the bad stuff!
	database = database[(database.CLsexp != 0) & database.failedstatus==0]

	try:
		listOfModels = database[args.modelDef.split(",")].drop_duplicates()
	except:
		print (">>> Problem! The model definition variables don't seem to exist! Quitting.")
		sys.exit(1)

	if args.debug:
		print (">>> ... List of Signal Models:")
		print (listOfModels)

	outputDB = doTheMuxing(args,database,listOfModels)

	# Handle the theory variation databses using the optimization from the nominal...

	if len(databaseListTheoryUp):
		outputDBTheoryUp   = doTheMuxing(args,databaseTheoryUp,listOfModels, outputDB)
		outputDBTheoryDown = doTheMuxing(args,databaseTheoryDown,listOfModels, outputDB)

	# Handle the theory variation databses using the optimization from the nominal...

	if len(databaseListUpperLimit):
		outputDBUpperLimit   = doTheMuxing(args,databaseUpperLimit,listOfModels, outputDB)

	return outputDB

def writeDBOutToFile(DB, filename):
	tmpDict = DB.to_dict(orient = "records")
	with open(filename, 'w') as f:
		f.write( json.dumps(tmpDict, indent=4) )

def doTheMuxing(args,database,listOfModels, nominalDatabase=0):

	outputDatabaseList = []

	for index, model in listOfModels.iterrows():

		tmpDatabase = database

		# Select for only those rows that correspond to my model point
		for item in list(model.index):
			tmpDatabase = tmpDatabase.loc[tmpDatabase[item] == model[item]]

		# Now I've filtered out so I'm only looking at statements for this specific model!

		# Use figure of merit to do optimization
		if type(nominalDatabase)==int:
			bestRow = tmpDatabase.loc[[tmpDatabase[args.figureOfMerit].idxmin()]]

		# Use nominalDatabase that has already been optimized for choosing best row
		else:

			tmpNominalDatabase = nominalDatabase
			# Select for only those rows that correspond to my model point
			for item in list(model.index):
				tmpNominalDatabase = tmpNominalDatabase.loc[tmpNominalDatabase[item] == model[item]]

			bestRow = tmpDatabase[(tmpDatabase.fID == tmpNominalDatabase.fID.iloc[0])]

		outputDatabaseList.append(bestRow)

	return pd.concat(outputDatabaseList)



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--inputFiles','-i', type=str, nargs='+', help='input json files', required=True)
	parser.add_argument("--outputFile","-o",  type=str, help="output json" , default = "outputJSON.json")
	parser.add_argument("--figureOfMerit","-f",  type=str, help="figure of merit", default = "CLsexp")
	parser.add_argument("--modelDef","-m",  type=str, help="comma separated list of variables that define a model", default = "mg,mlsp")
	parser.add_argument("--ignoreTheory","-t",      help = "ignore theory variation files", action="store_true", default=False)
	parser.add_argument("--ignoreUL",    "-u",      help = "ignore upper limit files", action="store_true", default=False)
	parser.add_argument("--debug","-d",      help = "debug", action="store_true", default=False)
	args = parser.parse_args()
	main(args)
