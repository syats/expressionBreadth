#Author: Victor Mireles

#This module provides a set of readers for diverse gene expression formats. It also provides a wrapper method for ease of use.

#All of these readers return a m x n binary matrix, where m is the number of genes and n is the number of conditions for which the expression is measured.

import numpy as np
import csv


def reader(fileName,ftype,binarized=True):
	if ftype == 'expressionAtlas':
		return expressionAtlas(fileName,binarized)

def expressionAtlas(fileName,binarized=True):
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
			numSamples = len(row) - 2;
			C = np.zeros((numGenes,numSamples));
			numGene = 0;

		thisVect = [float(row[i+2]) if row[i+2]!="" else 0 for i in range(numSamples)]
		C[numGene,:] = np.array(thisVect)
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
