#!/usr/bin/env python

# harvestToContours.py #################
#
# An attempt to unify contour production so we don't have a million attempts to reinvent a broken wheel...
#
# Reads in a harvest list from the HF output in either text file or json formats (but c'mon.. use json.. what are you doing?)
# Then the CLs values are interpreted as p-values and converted to significances and then interpolated
# Interpolation of significance surface can be configured to e.g. nn, linear, cubic (options from mpl)
#
# By: Larry Lee - Dec 2017


# You'll have to have matplotlib and root setup side-by-side
# In an ATLAS environment, you can do that with:
#> localSetupSFT --cmtConfig=x86_64-slc6-gcc48-opt releases/LCG_79/pytools/1.9_python2.7,releases/LCG_79/pyanalysis/1.5_python2.7
#> lsetup root


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import ROOT, json, argparse, math, sys, os, pickle, copy

ROOT.gROOT.SetBatch()

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

def make_from_args(args):
	## If I need to use scipy, please let me have scipy. I'll even help you out!
	if args.useUpperLimit:

		args.noSig = True
		args.level = 1.0

		if args.ignoreUncertainty:
			listOfContours = ["upperLimit","expectedUpperLimit"]
		else:
			listOfContours = ["upperLimit","expectedUpperLimit","expectedUpperLimitPlus1Sig","expectedUpperLimitPlus2Sig","expectedUpperLimitMinus1Sig","expectedUpperLimitMinus2Sig"]
		listOfContours_OneSigma = ["expectedUpperLimitPlus1Sig","expectedUpperLimitMinus1Sig"]
		listOfContours_TwoSigma = ["expectedUpperLimitPlus2Sig","expectedUpperLimitMinus2Sig"]
		expectedContour = "expectedUpperLimit"
		observedContour = "upperLimit"
	else:
		if args.ignoreUncertainty:
			listOfContours = ["CLs","CLsexp","upperLimit","expectedUpperLimit"]
		else:
			listOfContours = ["CLs","CLsexp","clsu1s","clsu2s","clsd1s","clsd2s","upperLimit","expectedUpperLimit"]
		listOfContours_OneSigma = ["clsu1s","clsd1s"]
		listOfContours_TwoSigma = ["clsu2s","clsd2s"]
		expectedContour = "CLsexp"
		observedContour = "CLs"
	return listOfContours,expectedContour,observedContour,listOfContours_OneSigma,listOfContours_TwoSigma



def main(args,inputData):
	"""Main function for driving the whole thing..."""


	# Print out the settings
	for setting in dir(args):
		if not setting[0]=="_":
			pass


	if args.xMin == None:
		pass

	if args.forbiddenFunction == "x":
		pass

	outputs = processInputFile(args = args, inputData = inputData, label = "")

	return outputs

def processInputFile(args,inputData, label = ""):
	"""Do actual processing of a given input file"""

	listOfContours,expectedContour,observedContour,listOfContours_OneSigma,listOfContours_TwoSigma = make_from_args(args)

	############################################################
	# Step 1 - Read in harvest list in either text or json format and dump it into a dictionary

	resultsDict = harvestToDict( inputJSON =  inputData, args = args )

	if len(resultsDict)<3:
		print (">>> WARNING: You have fewer than three valid model points in your input. I can't interpolate that in 2D! You've given me %d valid points!"%( len(resultsDict) ) )
		return -1


	if label!="_UL":
		truncateSignificances( args = args, modelDict = resultsDict , sigmax = args.sigmax )

	############################################################
	# Step 1.5 - If there's a function for a kinematically forbidden region, add zeros to dictionary

	if "none" not in args.forbiddenFunction.lower():
		resultsDict = addValuesToDict(args = args, inputDict = resultsDict, function = args.forbiddenFunction, numberOfPoints=100 ,value = "mirror" )

	############################################################
	# Step 2 - Interpolate the fit results


	outputGraphs = interpolateSurface( args = args, modelDict = resultsDict , interpolationFunction =  args.interpolation , useROOT = args.useROOT , outputSurface=True if label=="" else False)

	############################################################
	# Step 4 - Make pretty curves (and bands) or try to...

	outputs = {}

	if not args.ignoreUncertainty and label=="":
		if len(outputGraphs[listOfContours_OneSigma[0] ])==0 and len(outputGraphs[listOfContours_OneSigma[1] ])>0:
			print (">>>")
			print (">>> WARNING: You don't have +1 sigma sensitivity,")
			print (">>> ... but you do have -1 sigma reach. Making a ")
			print (">>> ... +/-1 sigma band from only the -1 side.")
			print (">>> ")
			for icurve,curve1 in enumerate(outputGraphs[listOfContours_OneSigma[1] ]):
				band_1s = createBandFromContours(args, contour1 = curve1 )
				outputs["Band_1s_%d"%icurve] = band_1s
		for icurve,(curve1,curve2) in enumerate(zip(outputGraphs[listOfContours_OneSigma[0] ],outputGraphs[listOfContours_OneSigma[1] ]) ):
			band_1s = createBandFromContours(args, contour1 =  curve1, contour2 = curve2 )
			outputs["Band_1s_%d"%icurve] = band_1s
		for icurve,(curve1,curve2) in enumerate(zip(outputGraphs[listOfContours_TwoSigma[0] ],outputGraphs[listOfContours_TwoSigma[1] ]) ):
			band_2s = createBandFromContours(args, contour1 =  curve1, contour2 = curve2 )
			outputs["Band_2s_%d"%icurve] = band_2s


	for icurve,obsCurve in enumerate(outputGraphs[observedContour]):
		outputs["Obs_%s%s"%(icurve,label)] = obsCurve
	for icurve,expCurve in enumerate(outputGraphs[expectedContour]):
		outputs["Exp_%s%s"%(icurve,label)] = expCurve

	return outputs


def harvestToDict(args, inputJSON, tmpListOfContours = None ):
	listOfContours,expectedContour,observedContour,_,_ = make_from_args(args)
	tmpListOfContours = tmpListOfContours or listOfContours 

	"""This parses the input file into a dictionary object for simpler handling"""


	modelDict = {}

	for sample in inputJSON:

		## Allowing filtering of entries via a constraints file
		if args.fixedParamsFile:
			failConstraintCutList = [not sample[str(constrainedVar)]==constrainedVal for (constrainedVar, constrainedVal) in constraintsDict.iteritems()]
			if any(failConstraintCutList):
				continue

		try:
			sampleParams = (float(sample[args.xVariable]),float(sample[args.yVariable]) )
		except:
			print (">>> ... Error: %s or %s doesn't exist as an entry in the input file"%(args.xVariable,args.yVariable))
			print (">>> ... Use cmd line options -x and -y to point to variables that exist in the input")
			print (">>> Available variables are listed below:")
			print (">>> ")
			print (">>> "+"\n>>> ".join(sample.keys()))
			sys.exit(1)

		sampleParamsList = list(sampleParams)
		if args.logX:
			sampleParamsList[0] = math.log10(sampleParamsList[0])
		if args.logY:
			sampleParamsList[1] = math.log10(sampleParamsList[1])
		sampleParams = tuple(sampleParamsList)

		if not math.isinf(float(sample[expectedContour])) :
			tmpList = [float(sample["%s"%x]) if (args.noSig or x in ["upperLimit","expectedUpperLimit"]) else ROOT.RooStats.PValueToSignificance( float(sample["%s"%x]) ) for x in tmpListOfContours]
			modelDict[sampleParams] = dict(zip(tmpListOfContours,  tmpList ) )
			if "fID" in sample:
				modelDict[sampleParams]["fID"] = sample["fID"]
			else:
				modelDict[sampleParams]["fID"] = ""


		else:
			if not sampleParams in modelDict:
				modelDict[sampleParams] = dict(zip(tmpListOfContours,  [args.sigmax for x in tmpListOfContours] ) )
				modelDict[sampleParams]["fID"] = ""
		if(args.debug):
			print (sampleParams, float(sample[observedContour]), float(sample[expectedContour]) if args.noSig else ROOT.RooStats.PValueToSignificance( float(sample[observedContour])     ))



	return modelDict


def addValuesToDict(args, inputDict, function, numberOfPoints = 100, value = 0):
	"""This takes in a TF1 and dots zero points along that function, and adds to the dict"""
	listOfContours,_,_,_,_ = make_from_args(args)

	tmpListOfXValues = [entry[0] for entry in inputDict.keys()]
	lowerLimit = min(tmpListOfXValues)
	upperLimit = max(tmpListOfXValues)

	forbiddenFunction = ROOT.TF1("forbiddenFunction${}".format(value),args.forbiddenFunction,lowerLimit,upperLimit)

	if value == "mirror":
		from scipy.spatial.distance import cdist
		def closest_point(pt, others):
		    distances = cdist(pt, others)
		    return others[distances.argmin()]

		def rotate(origin, point, angle=math.pi):
		    ox, oy = origin
		    px, py = point

		    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
		    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
		    return qx, qy

		inputDictCopy = copy.deepcopy(inputDict)

		forbiddenLineArray = []
		for xValue in [lowerLimit + x*(upperLimit-lowerLimit)/float(numberOfPoints*100) for x in range(numberOfPoints*100)]:
			forbiddenLineArray.append( ( xValue,forbiddenFunction.Eval(xValue) ) )

		# now to loop over entries in the inputDict. rotate them about this closest point on the forbidden line
		for signalPoint in inputDict:
			closestPointOnLine = list(closest_point(np.array([signalPoint]),np.array(forbiddenLineArray) ) )
			inputDictCopy[tuple(closestPointOnLine)] = dict(zip(listOfContours,  [1 for x in listOfContours] ) )
			fakeMirroredSignalPoint = rotate(closestPointOnLine, signalPoint)
			tmpDict = copy.deepcopy(inputDictCopy[signalPoint])
			for key in tmpDict:
				if isinstance(tmpDict[key], (int, float)):
					tmpDict[key] *= -1*np.sign(tmpDict[key])
			inputDictCopy[fakeMirroredSignalPoint] = tmpDict

		inputDict = copy.deepcopy(inputDictCopy)

	else:
		for xValue in [lowerLimit + x*(upperLimit-lowerLimit)/float(numberOfPoints) for x in range(numberOfPoints)]:
			inputDict[(xValue,forbiddenFunction.Eval(xValue))] = dict(zip(listOfContours,  [value for x in listOfContours] ) )

	return inputDict

def interpolateSurface(args, modelDict = {}, interpolationFunction = "linear", useROOT=False, outputSurface=False, outputSurfaceTGraph=False, tmpListOfContours=None):
	"""The actual interpolation"""

	listOfContours,expectedContour,observedContour,_,_ = make_from_args(args)
	tmpListOfContours = tmpListOfContours or listOfContours 


	modelPoints = modelDict.keys()
	modelPointsValues = modelDict.values()
	x0 =   list( list(zip( *modelPoints ))[0] )
	y0 =   list( list(zip( *modelPoints ))[1] )

	zValues = {} # entry x points
	x={} # entry x points
	y={} # entry x points

	graphs = {}

	for whichContour in tmpListOfContours:
		zValues[whichContour] = [ tmpEntry[whichContour] for tmpEntry in modelPointsValues]
		x[whichContour]       = [ a for a in x0 ];
		y[whichContour]       = [ a for a in y0 ];

	# remove inf point in each entry
	for whichContour in tmpListOfContours:

		while any([math.isinf(tmp) or math.isnan(tmp) for tmp in zValues[whichContour]  ]):#  np.isinf( zValues[whichContour]  ).any():
			myindex = [math.isinf(tmp) or math.isnan(tmp) for tmp in zValues[whichContour] ].index(True)
			if (args.debug):
				print (">>> ... Remove Inf or NaN at i=%d x=%d y=%d" % (myindex,x[whichContour][myindex],y[whichContour][myindex]))
			x[whichContour].pop(myindex)
			y[whichContour].pop(myindex)
			zValues[whichContour].pop(myindex)
			pass;

		if any([math.isinf(tmp) or math.isnan(tmp) for tmp in zValues[whichContour]  ]):
			print (">>> ... Still infs or nans in %s!! This is a problem... Exiting." % whichContour)
			sys.exit(0)


	if True:

		xi = {}
		yi = {}
		zi = {}
		for whichContour in tmpListOfContours:

			# Convert everything to numpy arrays
			xArray = np.array(x[whichContour])
			yArray = np.array(y[whichContour])
			zArray = np.array( zValues[whichContour] )

			# this scaling here equalizes the axes such that using a radial basis function makes sense!
			yScaling = np.max(xArray)/np.max(yArray) if np.max(yArray) else 1
			yArray = yArray*yScaling

			# Creating some linspaces for interpolation
			xlinspace = np.linspace(xArray.min() if args.xMin == None else args.xMin,
									xArray.max() if args.xMax == None else args.xMax,
									args.xResolution)
			ylinspace = np.linspace(yArray.min() if args.yMin == None else args.yMin,
									yArray.max() if args.yMax == None else args.yMax,
									args.yResolution)

			# Creating meshgrid for interpolation
			xymeshgrid = np.meshgrid(xlinspace,ylinspace)


			# Optional smoothing given by -s option
			smoothingFactor = 0
			if args.smoothing:
				smoothingFactor = float(args.smoothing)

			try:
				# Actual interpolation done by RBF
				if args.interpolationEpsilon:
					rbf = scipy.interpolate.Rbf(xArray, yArray, zArray, function=interpolationFunction, smooth=smoothingFactor , epsilon=args.interpolationEpsilon)
				else:
					rbf = scipy.interpolate.Rbf(xArray, yArray, zArray, function=interpolationFunction, smooth=smoothingFactor )
			except:
				print (">>> Interpolation failing!!! Check to make sure there are no NANs or double defined points in your input JSON!")
				print (">>> Printing points we're trying to interpolate (x,y,z) triplets:")

				print (sorted( zip(xArray,yArray,zArray), key = lambda x: x[0]*x[1] ))
				sys.exit(1)


			ZI = rbf(xymeshgrid[0], xymeshgrid[1])

			# Undo the scaling from above to get back to original units
			xymeshgrid[1] = xymeshgrid[1] / yScaling

			# Turn this surface into contours!
			contourList = getContourPoints(xymeshgrid[0],xymeshgrid[1],ZI, args.level)

			graphs[whichContour] = []
			for contour in contourList:
				graph = ROOT.TGraph(len(contour[0]), contour[0].flatten('C'), contour[1].flatten('C') )
				if graph.Integral() > args.areaThreshold:
					graphs[whichContour].append(graph)

			# Let's sort output graphs by area so that the band construction later is more likely to get the right pairs
			graphs[whichContour] = sorted(graphs[whichContour], key=lambda g: g.Integral() , reverse=True)

		return graphs

def truncateSignificances(args,modelDict,sigmax=5):
	"""Truncates significance to sigmax option"""
	listOfContours,_,_,_,_ = make_from_args(args)
	for model in modelDict:
		for thing in listOfContours:
			if modelDict[model][thing] > sigmax:
				modelDict[model][thing] = sigmax

	return

def getContourPoints(xi,yi,zi,level ):

	plt.ioff()
	f,ax = plt.subplots(1,1)
	c = ax.contour(xi,yi,zi, [level])
	contour = c.collections[0]
	plt.ion()

	contourList = []

	for i in range( len(contour.get_paths() ) ):
		v = contour.get_paths()[i].vertices

		x = v[:,0]
		y = v[:,1]

		contourList.append( (x,y) )

	plt.close()
	del f
	return contourList
import ctypes
def createBandFromContours(args,contour1,contour2=None):

	if not contour2:
		outputGraph = contour1
	else:
		outputSize = contour1.GetN()+contour2.GetN()+1

		pointOffset = 0
		if args.closedBands:
			outputSize += 2
			pointOffset = 1

		outputGraph = ROOT.TGraph(outputSize)
		tmpx, tmpy = ctypes.c_double(), ctypes.c_double()
		for iPoint in range(contour2.GetN()):
			contour2.GetPoint(iPoint,tmpx,tmpy)
			outputGraph.SetPoint(iPoint,tmpx,tmpy)

		if args.closedBands:
			contour2.GetPoint(0,tmpx,tmpy)
			outputGraph.SetPoint(contour2.GetN()+1, tmpx,tmpy)

		for iPoint in range(contour1.GetN()):
			contour1.GetPoint(contour1.GetN()-1-iPoint,tmpx,tmpy)
			outputGraph.SetPoint(contour2.GetN()+pointOffset+iPoint,tmpx,tmpy)

		if args.closedBands:
			contour1.GetPoint(contour1.GetN()-1,tmpx,tmpy)
			outputGraph.SetPoint(contour1.GetN()+contour2.GetN(), tmpx,tmpy)

		contour2.GetPoint(0,tmpx,tmpy)
		outputGraph.SetPoint(contour1.GetN()+contour2.GetN()+pointOffset,tmpx,tmpy)

	return outputGraph

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--inputFile","-i",  type = str, help="input harvest file", default = "test.json")
	parser.add_argument("--outputFile","-o", type = str, help="output ROOT file", default = "outputGraphs.root")
	parser.add_argument("--interpolation",   type = str, help="type of interpolation for scipy (RBF). e.g. linear, cubic, gaussian, multiquadric.", default = "multiquadric")
	parser.add_argument("--interpolationEpsilon", type=float, help="scipy (RBF) epsilon parameter", default = 0)
	parser.add_argument("--level",           type = float, help="contour level output. Default to 95%% CL", default = 1.64485362695)
	parser.add_argument("--useROOT","-r",    help = "use the root interpolation engine instead of mpl", action="store_true", default=False)
	parser.add_argument("--debug","-d",      help = "print extra debugging info", action="store_true", default=False)
	parser.add_argument("--sigmax",          type = float, help="maximum significance in sigmas", default = 5.0)
	parser.add_argument("--xVariable","-x",  type = str, default = "mg" )
	parser.add_argument("--yVariable","-y",  type = str, default = "mlsp" )
	parser.add_argument("--xResolution",     type = int, default = 100 )
	parser.add_argument("--yResolution",     type = int, default = 100 )

	parser.add_argument("--xMin", type=float, default = None )
	parser.add_argument("--yMin", type=float, default = None )
	parser.add_argument("--xMax", type=float, default = None )
	parser.add_argument("--yMax", type=float, default = None )

	parser.add_argument("--logX", help="use log10 of x variable", action="store_true", default=False)
	parser.add_argument("--logY", help="use log10 of y variable", action="store_true", default=False)

	parser.add_argument("--fixedParamsFile","-f",   type=str, help="give a json file with key=variable and value=value. e.g. use for pinning down third parameter in harvest list", default="")
	parser.add_argument("--forbiddenFunction","-l", type=str, help="""a ROOT TF1 definition for a forbidden line e.g. kinematically forbidden regions. (defaults to diagonal, i.e. -l 'x'). Set to 'None' to turn off.""", default="x")
	parser.add_argument("--ignoreUncertainty","-u", help="""Don't care about uncertainty bands!""", action="store_true", default=False)

	parser.add_argument("--areaThreshold","-a",     type = float, help="Throw away contours with areas less than threshold", default=0)
	parser.add_argument("--smoothing",    "-s",     type = str, help="smoothing option. For ROOT, use {k5a, k5b, k3a}. For scipy, uses smoothing from RBF.", default="0.1")
	parser.add_argument("--noSig","-n",      help = "don't convert CLs to significance -- don't use this option unless you know what you're doing!", action="store_true", default=False)

	parser.add_argument("--nominalLabel",      help = "keyword in filename to look for nominal sig XS", type=str, default="Nominal")

	parser.add_argument("--useUpperLimit", help="use upper limit information instead of CLs. Automatically turns off significance transform.", action="store_true", default=False)

	parser.add_argument("--closedBands","-b",      help = "if contours are closed shapes in this space, this can help with stitching issues if you're seeing weird effects", action="store_true", default=False)

	args = parser.parse_args()
	main(args)
