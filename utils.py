import os
import pickle
import numpy as np


def loadExperimentC(experimentName,outputDir,removeZeroesFromC=True):
	if 'Jaitin' in experimentName:
		dictFileName = outputDir+"JaitinMATRICES"+".pkl"
	else:
		dictFileName = outputDir+experimentName.split('_')[0]+"_MATRICES"+".pkl"
	d = pickle.load(open(dictFileName,'rb'));
	C = d[experimentName];

	if removeZeroesFromC:
		su = C.sum(axis=1);
		nz = np.nonzero(su)[0];
		C  = C[nz,:]

	return C

def getFilesFromDir(path,extension):
	fileList = [];
	for file in os.listdir(path):
	    if file.endswith(extension):
	        fileList.append(file)

	return fileList
