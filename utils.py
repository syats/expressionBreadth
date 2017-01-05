import os

def getFilesFromDir(path,extension):
	fileList = [];
	for file in os.listdir(path):
	    if file.endswith("."+extension):
	        fileList.append(file)

	return fileList
