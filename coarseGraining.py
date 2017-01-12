import random as ran
import numpy as np
import matplotlib.pyplot as plt

def coarseGrain(C,
				grainSizes,
				numReplicates,
				numTissues,
				sampleWithReplacement=True,
				rangesFromSamples=None):

	hbins = [float(i)-0.5 for i in range(numTissues+2)]
	allData = dict()
	for idx,grainSize in enumerate(grainSizes):  # !! --

		columnIndices = [i for i in range(C.shape[1])]
		subColums = [columnIndices]
		numTissuesPerSubM = numTissues
		if rangesFromSamples != None:
			subColums = [list(range(a[0],a[1])) for a in rangesFromSamples]
			numTissuesPerSubM = int(numTissues/ len(rangesFromSamples))
		if not sampleWithReplacement or (grainSize > len(subColums[0])):
			cids = columnIndices;
			while len(cids)<=len(grainSizes)*numTissues:
				cids += columnIndices;
			columnIndices = cids;
			for si in range(len(subColums)):
				sc = subColums[si];
				while len(sc) <= numTissuesPerSubM*grainSize:
					sc+=list(sc);
				subColums[si] = sc





		print("\n"+str(idx)+" ("+str(grainSize)+") :",end="", flush=True),
		dataThisGS = [[] for i in range(numTissues+1)]
		for replicateNum in range(numReplicates):  #!! --
			print("\t"+str(replicateNum)+" ",end="", flush=(replicateNum%4 == 0));
			ran.shuffle(columnIndices);
			#We generate a matrix
			newMat = np.zeros((C.shape[0],numTissues))
			for t in range(numTissues):
				if rangesFromSamples!=None:
					numSubM = int(t/numTissuesPerSubM)
					columnIndices = subColums[numSubM];
				if sampleWithReplacement:
					sa = ran.sample(columnIndices,grainSize);
				else:
					sa = [columnIndices[i] for i in range(t*grainSize,(t+1)*grainSize)]


				newMat[:,t] = C[:,sa].mean(axis=1);

			#We binarize it
			C2 = np.zeros_like(newMat);
			C2[np.nonzero(newMat)]=1;
			print(" "+str(C2.sum()/(C2.shape[0]*C2.shape[1])))
			#We compute its expression breadth
			su = C2.sum(axis=1);
			binCounts,nB,nC = plt.hist(su,hbins)
			for i in range(numTissues+1):
				dataThisGS[i].append(binCounts[i])
		print(".")
		allData[grainSize] = dataThisGS

	return allData
