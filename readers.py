#Author: Victor Mireles

#This module provides a set of readers for diverse gene expression formats. It also provides a wrapper method for ease of use.

#All of these readers return a m x n binary matrix, where m is the number of genes and n is the number of conditions for which the expression is measured.

import numpy as np
import csv


def reader(fileName,ftype,binarized=True,useReplicates=True):
	if ftype == 'expressionAtlas':
		return expressionAtlas(fileName,binarized=binarized,useReplicates=useReplicates)
	if ftype == 'Jaitin':
		p1 = fileName.split("/")[:-1];
		filePath = "/".join(p1)+"/"
		return jaitin2014(filePath,binarized=binarized)



def jaitin2014(filePath,binarized=True):

	#First we read the metadata
	metadataName = filePath+"GSE54006_experimental_design.txt";
	numSamples = 0;
	samplesPerType = dict();
	filePointer = open(metadataName,'r');
	cr = csv.reader(filePointer,delimiter='\t')
	toSkip = 1;
	for row in cr:
		toSkip -= 1;
		if toSkip >= 0:
			continue
		sampleId = row[14];
		if sampleId == "ommited":
			continue

		numSamples+=1
		sampleType = row[11];
		if sampleType in samplesPerType.keys():
			samplesPerType[sampleType].append(sampleId)
		else:
			samplesPerType[sampleType] = [sampleId];

	#Now we read the data and sort it into the different C matrices

	dataFileName = filePath + "GSE54006_umitab.txt";
	numGenes   = getNumLines(dataFileName) - 1;

	#Initialize the matrices
	allMatrices = dict()
	colsPerSample = dict();
	for sampleType,sampleIndices in samplesPerType.items():
		allMatrices[sampleType] = np.zeros((numGenes,len(sampleIndices)));
		colsPerSample[sampleType] = []

	filePointer = open(dataFileName,'r')
	rowNum = 0;
	cr = csv.reader(filePointer,delimiter='\t')
	goodGenes = np.zeros(numGenes);
	for row in cr:
		rowNum += 1;
		if rowNum == 1:
			numCols = len(row);
			for col in range(1,numCols):
				for sampleType,sampleIndices in samplesPerType.items():
					if row[col] in sampleIndices:
						colsPerSample[sampleType].append(int(col))


			continue
		if rowNum > numGenes:
			print(rowNum)
			continue
		for sampleType,colIndices in colsPerSample.items():
			vect = np.array([int(row[x]) for x in colIndices]);
			allMatrices[sampleType][rowNum-2,:] = vect
			if vect.sum() > 0:
				goodGenes[rowNum-2]  = 1;

	nz = np.nonzero(goodGenes)[0]
	for sampleType,C in allMatrices.items():
		if binarized:
			allMatrices[sampleType] = binarize(C[nz,:])
		else:
			allMatrices[sampleType] = C[nz,:]


	return allMatrices



def expressionAtlas(fileName,binarized=True,colsToIgnore=2,useReplicates=True):
	numGenes   = getNumLines(fileName)
	filePointer = open(fileName,'r')
	cr  = csv.reader(filePointer,delimiter='\t')
	toSkip = 1;
	numSamples = None
	skipped = 0;
	for row in cr:
		if row[0][0] == '#':
			skipped += 1;
			continue
		toSkip = toSkip - 1 if toSkip >= 0 else -1;
		if toSkip >= 0:
			skipped += 1;
			continue
		if numSamples == None:
			if useReplicates:
				numSamples = len(row) - colsToIgnore;
				numCols = numSamples;
			else:
				numSamples = 0;
				numCols = len(row) - colsToIgnore;
				for i in range(colsToIgnore,len(row)):
					numSamples+=len(row[i].split(','));
			C = np.zeros((numGenes,numSamples));
			numGene = 0;
		try:
			thisVect = [row[i] if row[i]!="" else 0 for i in range(colsToIgnore,numCols+colsToIgnore)]
			if useReplicates:
				thisVectF = [float(np.mean([float(xx) for xx in ff.split(',')])) if ',' in ff else float(ff) for ff in thisVect]
			else:
				thisVectF1 = [[float(xx) if xx != "" else 0 for xx in ff.split(',')] if ',' in ff else [float(ff)] for ff in thisVect]
				thisVectF = [item for sublist in thisVectF1 for item in sublist]

			C[numGene,:] = np.array(thisVectF)
		except ValueError:
			skipped+=1;
			continue
		numGene += 1;

	#If we skipped some lines, the matrix is probably too large.
	C = C[:numGenes-skipped,:];
	if binarized:
		return binarize(C);
	return C



def getNumLines(fileName):
	return sum(1 for line in open(fileName))

def binarize(D,threshold=0):
	C = np.zeros_like(D);
	C[np.nonzero(D>threshold)] = 1;
	return C
