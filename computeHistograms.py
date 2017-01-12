import numpy as np
import matplotlib.pyplot as plt
from importlib import reload
import math
import os.path
import pickle

import readers as rs
import utils as us




#each directory is a tuple (path,extension,fileType,prefix). All files ending in said extension in said path will be passed to readers.reader with ftype=fileType. Therefore fileType must be one of the kown filetypes, see readers.py for a list of them. Prefix will be used internally to name the data sets.
dirs = [ ('../otherArticleWork/expressionAtlas_baseline/','tsv','expressionAtlas','ExAt_') ]
dirs = [ ('../otherArticleWork/expressionAtlas_baseline/small/','tsv','expressionAtlas','ExAtS_') ] #small for testing
dirs = [('../otherArticleWork/Jaitin2014/','umitab.txt','Jaitin','Jaitin')]
outputDir = './outData/'
forceRead = False;
useReplicates = True;
minNumTissues = 6;


sufix = "" if useReplicates else "_NR"
allCs = dict();
#Maybe we've already read some of the Files, so we take a look at outputDir
if not forceRead:
	for di in dirs:
		dictFileName = outputDir+di[3]+"MATRICES"+sufix+".pkl"
		if os.path.isfile(dictFileName):
				d = pickle.load(open(dictFileName,'rb'));
				for ke,C in d.items():
					allCs[ke] = C;


for di in dirs:
	count=-1;
	readCount = 0;
	fileList = us.getFilesFromDir(di[0],di[1]);
	thisCs = dict()
	for fi in fileList:
		expName = fi[:-len(di[1])-1].split('-')
		if len(expName) >= 3:
			orName = expName[1]+expName[2];
		else:
			orName = expName[0]
		count+=1;
		print("file "+str(count)+"/"+str(len(fileList))+" :"+fi),
		if di[3]+orName in allCs.keys():
			print('...loaded')
			thisCs[di[3]+orName] = allCs[di[3]+orName];
			continue
		C = rs.reader(di[0]+fi,di[2],useReplicates=useReplicates);
		readCount += 1;
		#print('....read '+str(C.shape))
		if isinstance(C,np.ndarray):
			if min(C.shape) < minNumTissues:
				continue
			thisCs[di[3]+orName]=C
			allCs[di[3]+orName]=C
		else:
			for ke,ma in C.items():
				if min(ma.shape) < minNumTissues:
					continue
				thisCs[di[3]+orName+ke] = ma
				allCs[di[3]+orName+ke]  = ma



	if readCount > 0:
		pickle.dump(thisCs,open(outputDir+di[3]+"MATRICES"+sufix+".pkl",'wb'))




numHists = len(allCs)
numSubplotCols = 1;
numSubplotRows = 1;
if numHists > 2:
	numSubplotCols = 2;
	numSubplotRows = math.ceil(numHists/2);

subPlotNum = 1;
for ke,C in allCs.items():
	su = C.sum(axis=1);
	numTissues = min(C.shape)
	hbins = [float(i)+0.5 for i in range(numTissues+1)]
	plt.subplot(numSubplotRows,numSubplotCols,subPlotNum)
	nA,nB,nC = plt.hist(su,hbins)

	plt.xlim([-0.5,numTissues+.5])

	plt.xticks(range(numTissues+1))
	plt.text(numTissues/3,max(nA)*0.5,ke)
	subPlotNum += 1;

plt.show();
