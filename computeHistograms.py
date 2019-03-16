import numpy as np
import matplotlib.pyplot as plt

from importlib import reload
import math
import os.path
import pickle

import readers as rs
import utils as us

# each directory is a tuple (path,extension,fileType,prefix). All files ending in said extension in said path will be
#  passed to readers.reader with ftype=fileType. Therefore fileType must be one of the kown filetypes, see readers.py
# for a list of them. Prefix will be used internally to name the data sets.

dirs = [('./data/expressionAtlas/', 'tsv', 'expressionAtlas', 'ExAt_')]
toskip = ["MTAB2512"]
#dirs = [
#    ('../otherArticleWork/expressionAtlas_baseline/small/', 'tsv', 'expressionAtlas', 'ExAtS_')]  # small for testing
#dirs = [('../otherArticleWork/Jaitin2014/', 'umitab.txt', 'Jaitin', 'Jaitin')]
outputDir = './outData/'
forceRead = False;
useReplicates = True;
minNumTissues = 6;

sufix = "" if useReplicates else "_NR"
allCs = dict();
# Maybe we've already read some of the Files, so we take a look at outputDir
if not forceRead:
    for di in dirs:
        dictFileName = outputDir + di[3] + "MATRICES" + sufix + ".pkl"
        if os.path.isfile(dictFileName):
            d = pickle.load(open(dictFileName, 'rb'));
            for ke, C in d.items():
                allCs[ke] = C;

for di in dirs:
    count = -1;
    readCount = 0;
    fileList = us.getFilesFromDir(di[0], di[1]);
    thisCs = dict()
    for fi in fileList:
        expName = fi[:-len(di[1]) - 1].split('-')
        if len(expName) >= 3:
            orName = expName[1] + expName[2];
        else:
            orName = expName[0]
        count += 1;
        print("file " + str(count) + "/" + str(len(fileList)) + " : " + fi),
        if di[3] + orName in allCs.keys():
            print('...loaded')
            thisCs[di[3] + orName] = allCs[di[3] + orName];
            continue
        print("Skipped:  ",orName)
        continue
        C = rs.reader(di[0] + fi, di[2], useReplicates=useReplicates);
        readCount += 1;
        # print('....read '+str(C.shape))
        if isinstance(C, np.ndarray):
            if min(C.shape) < minNumTissues:
                continue
            thisCs[di[3] + orName] = C
            allCs[di[3] + orName] = C
        else:
            for ke, ma in C.items():
                if min(ma.shape) < minNumTissues:
                    continue
                thisCs[di[3] + orName + ke] = ma
                allCs[di[3] + orName + ke] = ma

    if readCount > 0:
        pickle.dump(thisCs, open(outputDir + di[3] + "MATRICES" + sufix + ".pkl", 'wb'))

numHists = len(allCs)
numSubplotCols = 1;
numSubplotRows = 1;
if numHists > 2:
    numSubplotCols = 2;
    numSubplotRows = math.ceil(numHists / 2);

subPlotNum = 1;
allkeys = list(allCs.keys())
allkeys.sort()
for ke in  allkeys:
    C = allCs[ke]
    if any([ke.endswith(x) for x in toskip]):
        print("Skipping ", ke)
        continue
    su = C.sum(axis=1);
    numTissues = min(C.shape)
    img = np.array([(su==i).sum()/float(max(C.shape)) for i in range(1, numTissues+1)])
    #plt.subplot(numSubplotRows, numSubplotCols, subPlotNum)

    fig = plt.figure(1)
    plt.subplot(numSubplotRows, numSubplotCols, subPlotNum)
    im = plt.imshow(img[np.newaxis,:],interpolation='none',aspect='auto')
    if  subPlotNum in [len(allkeys)-2, len(allkeys)-1]:
        if subPlotNum == len(allkeys)-1:
            off = float(numTissues)/20
        else:
            off = 0
        plt.xticks([off+i*(float(numTissues)/10) for i in range(10)],
                   ["{:.1f}".format(0.1+i*0.1)  for i in range(10)])
        plt.xlabel("Proportion of tissues",fontsize=18)
    else:
        plt.xticks([])
    plt.yticks([])
    h = plt.ylabel(ke, labelpad=45, rotation=0)
    fig.gca().yaxis.set_label_coords(-0.1,0)


    #
    # plt.figure(2)
    # plt.subplot(numSubplotRows, numSubplotCols, subPlotNum)
    # hbins = [float(i) + 0.5 for i in range(numTissues + 1)]
    # nA, nB, nC = plt.hist(su, hbins)
    # plt.xlim([-0.5, numTissues + .5])
    # plt.xticks(range(numTissues + 1))
    # plt.text(numTissues / 3, max(nA) * 0.5, ke)

    subPlotNum += 1;

cbar_ax = fig.add_axes([0.952, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax)
cbar_ax.set_title("Proportion\nof genes",fontsize=16)
plt.subplots_adjust(left=0.08,right=0.92,top=0.97,bottom=0.05,wspace=0.25,hspace=0.20)

plt.show();
