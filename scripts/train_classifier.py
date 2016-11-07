#!/usr/bin/python

from __future__ import print_function

import sys, os, time
from distutils import dir_util  # handle dirs

import numpy as np

from sklearn.externals        import joblib


stime = time.time()


if (len(sys.argv) != 3):
    sys.stderr.write("usage: %s trainDir classifierType\n" % (sys.argv[0]))
    sys.exit(-1)

trainDir = sys.argv[1]
if not os.path.isdir(trainDir):
    sys.stderr.write("not a valid directory: %s\n" % (trainDir))
    sys.exit(-1)

clfType = sys.argv[2]
if clfType == "SVC":
    from sklearn.svm import SVC
    clf = SVC()
elif clfType == "DecisionTreeClassifier":
    from sklearn.tree import DecisionTreeClassifier
    clf = DecisionTreeClassifier()
elif clfType == "RandomForestClassifier":
    from sklearn.ensemble import RandomForestClassifier
    clf = RandomForestClassifier()
elif clfType == "ExtraTreesClassifier":
    from sklearn.ensemble import ExtraTreesClassifier
    clf = ExtraTreesClassifier()
elif clfType == "AdaBoostClassifier":
    from sklearn.ensemble import AdaBoostClassifier
    clf = AdaBoostClassifier()
elif clfType == "LogisticRegression":
    from sklearn.linear_model import LogisticRegression
    clf = LogisticRegression()
elif clfType == "MultinomialNB":
    from sklearn.naive_bayes import MultinomialNB
    clf = MultinomialNB()
elif clfType == "GaussianNB":
    from sklearn.naive_bayes import GaussianNB
    clf = GaussianNB()
elif clfType == "LinearDiscriminantAnalysis":
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    clf = LinearDiscriminantAnalysis()
elif clfType == "QuadraticDiscriminantAnalysis":
    from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
    clf = QuadraticDiscriminantAnalysis()
else:
    sys.stderr.write("unsupported classifierType: %s\n" % clfType)
    sys.exit(-1)

clfDir = trainDir + '/Classifiers/' + clfType
if not os.path.isdir(clfDir):
    sys.stderr.write("making dir: %s\n" % (clfDir))
    dir_util.mkpath(clfDir)


X = np.genfromtxt(trainDir + '/X', delimiter='\t')
y = np.genfromtxt(trainDir + '/y', delimiter='\t')

clf.fit(X,y)

sys.stderr.write("joblib.dump-ing clf\n")
clfFilename = clfDir + '/%s' % (clfType) 
joblib.dump(clf,clfFilename)

sys.stderr.write("total time: %0.2f\n" % (time.time()-stime) )
