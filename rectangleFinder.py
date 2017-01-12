import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import random as ran
import math


def findRectangles(C,base):
	if C.shape[1] > 15:
		print("Too large matrix, aborting")
		return set()

	if (C.shape[1] < base) or (base <= 0):
		print("base size is out of range")
		return set()


	foundBases  = set();
	origIndices = [i for i in range(C.shape[1])]
	combis = it.combinations(origIndices,base)

	for comb in combis:
		subC = C[:,comb]
		su = subC.sum(axis=1);
		heig = np.nonzero(su==base)[0].tolist();
		if len(heig) > 1:
			heig.sort()
			foundBases.add(tuple(heig))

	return foundBases

def getMatrixFromDist(dist):
	n = len(dist)-1;
	m = sum(dist)
	C = np.zeros((m,n));
	nextRow = 0;
	for s,numS in enumerate([int(ss) for ss in dist]):
		for x in range(numS):
			if s == 0:
				nextRow+=1;
				continue
			colsToUse = ran.sample(range(n),s)
			C[nextRow,colsToUse] = 1;
			nextRow+=1;

	su = C.sum(axis=1);
	#hbins = [float(i)-0.5 for i in range(n+2)]
	#h1,h2 = np.histogram(su,hbins)
	#print("\t"+str(h1))
	return C


def generateMatrices(dis,numReps):
	matrices = []
	for i in range(numReps):
		C = getMatrixFromDist(dis)
		matrices.append(C)

	return matrices

def nCk(n,k):
	return math.factorial(n)/(math.factorial(k)*math.factorial(n-k))

def SRHistoForDistribution(putativeCs,maxSizeO = None,plot2D=False):

	numReps = len(putativeCs)
	allHistos = []
	totMax = -np.inf;
	totMin = np.inf;

	COld=None;
	for rn in range(numReps):

		# example----
		C = putativeCs[rn]
		m = C.shape[0];
		if maxSizeO == None:
			maxSize = m;
		else:
			maxSize = maxSizeO
		#print(str(C.sum()/(C.shape[0]*C.shape[1]))+" ",end="")
		#if COld != None:
		#	print("d:"+str(np.abs(C-COld).sum()  ),end=""  )

		SRHisto = np.zeros((maxSize,C.shape[1]+1))
		for r in range(1,C.shape[1]+1):
			fr = findRectangles(C,r);
			for s in [len(x) for x in fr if len(x) < maxSize]:
				SRHisto[s,r]+=1

		if SRHisto.max() >  totMax:
			totMax = SRHisto.max();

		if SRHisto.min() < totMin:
			totMin = SRHisto.min();

		allHistos.append(SRHisto)
		COld = C;

	'''
	hnz = [np.nonzero(hh>1) for hh in allHistos]

	mms = [(hnzz[0].min(),hnzz[0].max()) for hnzz in hnz]
	minX = min([mm[0] for mm in mms])
	maxX = max([mm[0] for mm in mms])


	for idx,SRHisto in enumerate(allHistos):
		plt.subplot(4,int(np.ceil(numReps/4)),idx+1)
		if plot2D:
			plt.imshow(SRHisto[minX:maxX,:],vmin=totMin,vmax=totMax,aspect='auto',interpolation='none')
		else:
			plt.plot(np.array(range(SRHisto.shape[1])),
			SRHisto.sum(axis=0)/SRHisto.sum() )
		if idx % int(np.ceil(numReps/4)) != 0:
			plt.yticks([])
		#plt.xticks([])

	plt.subplots_adjust(left=0.02,bottom=0.02,right=0.98,top=0.98,wspace=0.01,hspace=0.01)
	'''

	return allHistos

def fromHistosToAvgSizePerReuse(allHistos,dis):
	rs = allHistos[0].shape[1]
	reuseV = [[] for i in range(rs)]

	dS = sum([i*dis[i] for i in range(len(dis))])/(rs-1)
	print("\tds:"+str(dS))

	prevS = None
	for SRHisto in allHistos:
		if prevS != None:
			print("*"+str( np.abs(prevS-SRHisto).sum()  ))
		ixs = range(SRHisto.shape[0]);
		for r in range(rs):
			meanSizeR =(SRHisto[:,r]*ixs).sum() /  SRHisto[:,r].sum()
			reuseV[r].append(meanSizeR/dS)
		prevS = SRHisto

	return reuseV

#Of all the nCk sets of size k, we compute what is the probability of finding one which has reuse larger than 1.
def fromHistosToReuseV(allHistos,normalized=True):
	rs = allHistos[0].shape[1]
	reuseV = [[] for i in range(rs)]
	prevS = None
	for SRHisto in allHistos:

		su = SRHisto.sum(axis=0)
		#print(su)
		#print(np.array([nCk(rs,k) for k in range(rs)]))
		if normalized:
			su = su / np.array([nCk(rs-1,k) for k in range(rs)])

		for r in range(rs):
			reuseV[r].append(su[r])
		prevS = SRHisto

	return reuseV;
