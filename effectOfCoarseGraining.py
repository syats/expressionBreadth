import numpy as np
import matplotlib.pyplot as plt
from importlib import reload
import math
import random as ran
import os.path
import pickle

import readers as rs
import utils as us
import rectangleFinder as rf
import coarseGraining as cg



outputDir = './outData/'
numTissues         = 12;
numReplicates 	   = 10;
numGenesReuse      = 180;
grainSizes = [1,2,4,6,8,10,20,30]
numReplicatesReuse = 20;
#For quick testing
grainSizes = [1,5,10,30,50]
numReplicatesReuse = 4;

analyzeModulePosibilities = False
removeZeroesFromC = True;
plotGeneReuseHisto = True;
sampleWithReplacement = True;
mixWhenSamplingFromSeveral = False


experiments = ['ExAt_MTAB2512'];
experiments += [['JaitinGSE54006NK_cell', 'JaitinGSE54006CD11c+', 'JaitinGSE54006GC B cell', 'JaitinGSE54006splenocyte'],
	['JaitinGSE54006CD8+pDC', 'JaitinGSE54006monocyte_or_neutrophil', 'JaitinGSE54006CD8-CD4+ESAM+'],
	[ 'JaitinGSE54006CD8+CD86+', 'JaitinGSE54006B cell', 'JaitinGSE54006CD8-pDC'],
	['JaitinGSE54006CD8+CD86-', 'JaitinGSE54006CD11c+(2hr_LPS)'], 'JaitinGSE54006hemato lineages']

experiments = ['JaitinGSE54006NK_cell', 'JaitinGSE54006CD11c+', 'JaitinGSE54006GC B cell', 'JaitinGSE54006splenocyte',
	'JaitinGSE54006CD8+pDC', 'JaitinGSE54006monocyte_or_neutrophil', 'JaitinGSE54006CD8-CD4+ESAM+',
	 'JaitinGSE54006CD8+CD86+', 'JaitinGSE54006B cell', 'JaitinGSE54006CD8-pDC',
	'JaitinGSE54006CD8+CD86-', 'JaitinGSE54006CD11c+(2hr_LPS)', 'JaitinGSE54006hemato lineages']


for experimentName in experiments:
	print(str(experimentName)+" ________ \n")
	fig = plt.figure()
	winTit = str(experimentName)+" "+("W_Rep" if sampleWithReplacement else "WO_Rep" )
	fig.canvas.set_window_title(winTit)


	if not isinstance(experimentName,list):
		C = us.loadExperimentC(experimentName,
						outputDir,
						removeZeroesFromC=removeZeroesFromC)
	else:
		print("Set --->",end="",flush=True)
		CsToConcat = [];
		for expName in experimentName:
			C = us.loadExperimentC(expName,
							outputDir,
							removeZeroesFromC=False)
			CsToConcat.append(C.copy())
			print(str(C.shape)+" ",end="",flush=True);


		totM = C.shape[0];
		totN = sum([CC.shape[1] for CC in CsToConcat])
		if any([CC.shape[0] != totM for CC in CsToConcat]):
			print("Matrices in "+str(experimentName)+" not of the same number of rows, skipping")
			continue
		C = np.zeros((totM,totN));
		print(" <--- "+str(C.shape))
		colBegin = 0;
		for CC in CsToConcat:
			numCols = CC.shape[1]
			C[:,colBegin:colBegin+numCols] = CC;
			colBegin = colBegin+numCols;


	if (not mixWhenSamplingFromSeveral) and isinstance(experimentName,list):
		colNumbers = [CC.shape[1] for CC in CsToConcat];
		rangesFromSamples=[]
		initCol = 0;
		for cn in colNumbers:
			rangesFromSamples.append([initCol,initCol+cn])
			initCol += cn
	else:
		rangesFromSamples = None

	allData = cg.coarseGrain(C,
							grainSizes,
							numReplicates,
							numTissues,
							sampleWithReplacement=sampleWithReplacement,
							rangesFromSamples=rangesFromSamples)

	if analyzeModulePosibilities:

		meanReusableModuleSize = dict()
		numberOfReusableModules = dict()
		print("computing module-wise reusabilities")
		for idx,grainSize in enumerate(grainSizes):
			print(str(idx)+" ("+str(grainSize)+") :")
			thisData=allData[grainSize]
			#For every possible number of tissues i \in [0,numTissues], dis[i] contains the average number of genes that are present in exactly i tissues in all the matrices generated with grain size grainSize;

			dis = [int(np.mean(xx)) for xx in thisData ]
			sumDis = sum(dis);
			dis = [dd/sumDis for dd in dis]  #This is an actual distribution (sums up to 1)

			#Now we normalized it to have the desired number of genes
			dis = [int(np.round(numGenesReuse*dd)) for dd in dis]

			#putativeCs is a list of matrices generated using dis as row sum distribution. Each is of size numTissues x C.shape[0].
			putatitiveCs = rf.generateMatrices(dis,numReps=numReplicatesReuse)


			#allH contains numReplicatesReuse matrices of size numTissues+1 x C.shape[0], one for each matrix in putatitiveCs. For each matrix in putatitiveCs a matrix M is generated that, such for every possible reuse r \in [0,numTissues], M[s,r] contains the number of different (maximal) sets of s rows of C2 that are all ones in r columns of C2.
			allH= rf.SRHistoForDistribution(putatitiveCs,maxSizeO=sum(dis))

			#For a given r \in [0,numTissues], rV[r] is the list given by [M.sum(axis=0)[r] for M in allH]
			rV = rf.fromHistosToAvgSizePerReuse(allH,dis)
			meanReusableModuleSize[grainSize] = rV;

			rV2 = rf.fromHistosToReuseV(allH)
			numberOfReusableModules[grainSize] = rV2;


	numHists = len(allData)
	extraRow = 1 if plotGeneReuseHisto else 0;
	numSubplotCols = 3 if analyzeModulePosibilities else 1;
	numSubplotRows = len(grainSizes)+extraRow;


	if plotGeneReuseHisto:
		su = C.sum(axis=1);
		nT = min(C.shape)
		nBinsT = nT if nT < 100 else 100
		binSize = nT/nBinsT;
		hbins = [float(i*binSize)+0.5 for i in range(nBinsT+1)]
		plt.subplot(numSubplotRows,numSubplotCols,1)
		nA,nB,nC = plt.hist(su,hbins)
		plt.xlim([-0.5,nT+.5])
		plt.xticks([int(i*binSize) for i in range(nBinsT+1)])

	for idx,grainSize in enumerate(grainSizes):

		#A violin plot of the element-wise reusability of all matrice of grain size grainSize
		thisData = allData[grainSize]
		varis = [np.std(a) for a in thisData]
		goodIndices = np.nonzero(varis)[0].tolist();
		plt.subplot(numSubplotRows,numSubplotCols,(extraRow+idx)*numSubplotCols+1)
		if len(goodIndices) > 0:
			plt.violinplot([thisData[i] for i in goodIndices ],[range(numTissues+1)[i] for i in goodIndices])
		plt.plot([np.mean(x) for x in thisData],'ro')
		plt.text(numTissues/3,16000,"Grain Size: "+str(grainSize) )
		plt.xlim([-0.5,numTissues+.5])
		plt.ylim([0,20000])
		if idx==0:
			plt.title('Expression Breadth')
		if idx < len(grainSizes)-1:
			plt.xticks([])
		else:
			plt.xlabel("Number of conditions a single \n gene is expressed in")

		if not analyzeModulePosibilities:
			continue

		#From the above violin plot, the mean (per replicate) element-wise reusability is computed, resulting in dis, a distribution of element-wise reusabilities.  The violin plot below shows, for each r \in [0,numTissues], how is the number of (non-singleton) modules of reusability r distributed for numReplicatesReuse matrices with element-wise reusability distributed according to dis.
		plt.subplot(numSubplotRows,numSubplotCols,(extraRow+idx)*numSubplotCols+2)
		thisSus = numberOfReusableModules[grainSize]
		varis = [np.std(a) for a in thisSus]
		goodIndices = np.nonzero(varis)[0].tolist()
		if len(goodIndices) > 0:
			plt.violinplot([thisSus[i]  for i in goodIndices],goodIndices);
		plt.plot([np.mean(x) for x in thisSus],'o--')
		plt.ylim([0,1.09])
		plt.xlim([-0.05,numTissues+0.05])
		if idx==0:
			plt.title('Number of reusable modules')
		if idx < len(grainSizes)-1:
			plt.xticks([])
		else:
			plt.xlabel("Number of conditions a module \n can be reused in")



		#Now, how is the size of (non-singleton) modules of reusability r distributed
		plt.subplot(numSubplotRows,numSubplotCols,(extraRow+idx)*numSubplotCols+3)
		thisSus = meanReusableModuleSize[grainSize]
		varis = [np.std(a) for a in thisSus]
		goodIndices = np.nonzero(varis)[0].tolist()
		if len(goodIndices) > 0:
			plt.violinplot([thisSus[i]  for i in goodIndices],goodIndices);
		plt.plot([np.mean(x) for x in thisSus],'ro--')
		#plt.ylim([0,1])
		plt.yticks([0,0.2,0.4,0.6,0.8])
		plt.xlim([-0.05,numTissues+0.05])
		if idx==0:
			plt.title('Mean size of reusable modules')
		if idx < len(grainSizes)-1:
			plt.xticks([])
		else:
			plt.xlabel("Number of conditions a module \n can be reused in")



plt.show()
