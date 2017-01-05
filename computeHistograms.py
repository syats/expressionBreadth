import numpy as np
import matplotlib.pyplot as plt
from importlib import reload

import readers as rs
import utils as us



#each directory is a tuple (path,extension,fileType,prefix). All files ending in said extension in said path will be passed to readers.reader with ftype=fileType. Therefore fileType must be one of the kown filetypes, see readers.py for a list of them. Prefix will be used internally to name the data sets.
dirs = [
	('../otherArticleWork/expressionAtlas_baseline/','tsv','expressionAtlas','ExAt')
]

allCs = dict();
for di in dirs:
	fileList = us.getFilesFromDir(di[0],di[1]);
	for fi in fileList:
		C = rs.reader(di[0]+fi,di[2]);
		orName = fi[:-len(di[1])-1].split('-')[-1];
		allCs[di[3]+orName]=C
