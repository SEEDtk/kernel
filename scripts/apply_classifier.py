#!/usr/bin/env python

from __future__ import print_function  # <-- this lets python2 use python3's print function

import sys, time, os
from distutils import dir_util  # handle dirs

import numpy as np

from sklearn.externals import joblib

stime = time.time()

if (len(sys.argv) != 4):
    sys.stderr.write("usage: %s testDir trainDir classifierType\n" % (sys.argv[0]))
    sys.exit(-1)

testDir  = sys.argv[1]
trainDir = sys.argv[2]
clfType  = sys.argv[3]

if not os.path.isdir(trainDir):
    sys.stderr.write("not a valid trainDir: %s\n" % (trainDir))
    sys.exit(-1)

X = np.genfromtxt(testDir + '/X', delimiter='\t')
y = np.genfromtxt(testDir + '/y', delimiter='\t')

#...genfromtxt "helpfully" reformats a 1-row dataset as a vector, not an array,
#   so the number of dimensions must be tested, and reformatted if not 2D array.
if len(X.shape) != 2:
    X = X.reshape(1,-1)

if len(y.shape) != 1:
    # sys.stderr.write("y dimensions: %s\n", y.shape)
    y = np.reshape(y, (1,))
    # sys.stderr.write("New y dimensions: %s\n", y.shape)

clfFilename = trainDir + '/Classifiers/' + clfType + '/' + clfType
clf = joblib.load(clfFilename)

outDir = testDir + '/Classifiers/' + clfType
if not os.path.isdir(outDir):
    sys.stderr.write("making directory: %s\n" % (outDir))
    dir_util.mkpath(outDir)
outFilename = outDir + '/y.out'
outFile = open(outFilename,'w')

ncorrect = 0
nchecked = 0
for idx in range(X.shape[0]):
    nchecked += 1
    instance  = X[idx]
    instances = [instance]
    predictions = clf.predict(instances)
    prediction  = predictions[0]
    if prediction == y[idx]:
        ncorrect  += 1

    print(idx, y[idx], prediction, sep='\t')
    print(prediction, file=outFile)

outFile.close()

# <--- python2 wants pgmr to convert to float
sys.stderr.write("Accuracy: %0.2f\n"   % (float(ncorrect)/nchecked) )
sys.stderr.write("Total time: %0.2f\n" % (time.time()-stime) )
