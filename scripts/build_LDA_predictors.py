import numpy as np
import sys, os, time
import shutil
import argparse
from sklearn.model_selection import KFold
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib

def make_predictor(n_col):
        role_ID = col_names[n_col][1]

        pred_dir = args.probDir + "/Predictors/" + role_ID
        if not os.path.isdir(pred_dir):
                os.mkdir(pred_dir)
        clf_dir = pred_dir + "/Classifiers"
        if not os.path.isdir(clf_dir):
                os.mkdir(clf_dir)
        pkl_dir = clf_dir + '/' + args.classifier + 'Classifier'
        if not os.path.isdir(pkl_dir):
                os.mkdir(pkl_dir)

        y  = X_all[:, n_col]
        X  = np.delete(X_all, n_col, axis=1)
        kf = KFold(n_splits=n_fold, shuffle=True, random_state=30258509)
        accuracies = np.zeros(n_fold)
        n_cross = 0

        for training, testing in kf.split(X):
                clf.fit(X[training], y[training])
                prediction = clf.predict(X[testing])
                accuracies[n_cross] = np.mean(prediction == y[testing])*100.00
                n_cross += 1

        summary = [[args.classifier + 'Classifier']]
        summary[0].append(np.round(np.mean(accuracies),2))
        Q = np.percentile(accuracies, np.arange(0,125,25))
        Q = np.round(Q, 2)
        summary[0].extend(Q.tolist())
        summary[0].append(np.round(0.25*Q[1] + 0.5*Q[2] + 0.25*Q[3],2))
        summary[0].append(np.round(Q[3]-Q[1],2))
        summary = np.array(summary)

        clf.fit(X, y)
        clfFilename = pkl_dir + '/' + args.classifier + 'Classifier'
        joblib.dump(clf, clfFilename)
        acc_file = pkl_dir + "/accuracy"
        np.savetxt(acc_file, summary, fmt="%s", delimiter='\t')
        print("Completed %d: %s predictor." % (n_col, role_ID))



stime = time.time()

parser = argparse.ArgumentParser(description='Build role predictors')
parser.add_argument("probDir", help="Directory which contains a role matrix",
                    type=str)
parser.add_argument("-f", "--fraction", dest="fraction", default=0.2, type=float, help="Fraction of data to use in testing")
parser.add_argument("-c", "--classifier", dest="classifier", default="RandomForest", help="Type of sklearn classifier to use")
parser.add_argument("-n", "--n_jobs", dest="n_jobs", default=16, help="Number of parallel jobs to run")
parser.add_argument("--clear", action="store_true",
                  help="Clear input probDir Predictors directory")
args = parser.parse_args()

try:
        X_all = np.genfromtxt(args.probDir + '/X', delimiter='\t')
except:
        sys.stderr.write("Matrix X not found in %s!\n" % args.probDir)
        sys.exit(-1)

try:
        col_names = np.genfromtxt(args.probDir + '/col.h', delimiter='\t', dtype=str)
except:
        sys.stderr.write("Column names file col.h not found in %s!\n" % args.probDir)
        sys.exit(-1)

pred_completed = 0
pred_to_run = []

if os.path.isdir(args.probDir + "/Predictors"):
    if __name__ == '__main__':
        if args.clear:
                print("Clearing " + args.probDir + "...")
                shutil.rmtree(args.probDir + "/Predictors")
                os.mkdir(args.probDir + "/Predictors")
                pred_to_run = range(col_names.shape[0])

        else:
                for i in range(col_names.shape[0]):
                        acc_file = args.probDir + "/Predictors/" + col_names[i,1] + "/Classifiers/" + args.classifier + "Classifier/accuracy"
                        if os.path.isfile(acc_file):
                                pred_completed += 1
                        else:
                                pred_to_run.append(i)

else:
        os.mkdir(args.probDir + "/Predictors")
        pred_to_run = range(col_names.shape[0])

if __name__ == '__main__':
        print("%d out of %d functions already processed." % (pred_completed, X_all.shape[1]))

n_fold = int(1./args.fraction)

clfType = args.classifier
if clfType == "SVC":
        from sklearn.svm import SVC
        clf = SVC()
elif clfType == "DecisionTree":
        from sklearn.tree import DecisionTreeClassifier
        clf = DecisionTreeClassifier()
elif clfType == "RandomForest":
        from sklearn.ensemble import RandomForestClassifier
        clf = RandomForestClassifier(n_estimators=100, max_features=0.2, criterion='entropy', random_state=4721359)
elif clfType == "ExtraTrees":
        from sklearn.ensemble import ExtraTreesClassifier
        clf = ExtraTreesClassifier()
elif clfType == "AdaBoost":
        from sklearn.ensemble import AdaBoostClassifier
        clf = AdaBoostClassifier()
elif clfType == "LogisticRegression":
        from sklearn.linear_model import LogisticRegression
#	clf = LogisticRegression(solver="saga", penalty="l1", C=0.1, multi_class="multinomial", max_iter=50)
        clf = LogisticRegression(solver="newton-cg", penalty="l2", C=1.0, multi_class="ovr", max_iter=50)
elif clfType == "MultinomialNB":
        from sklearn.naive_bayes import MultinomialNB
        clf = MultinomialNB()
elif clfType == "GaussianNB":
        from sklearn.naive_bayes import GaussianNB
        clf = GaussianNB()
elif clfType == "LDA":
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    clf = LinearDiscriminantAnalysis()#solver="eigen"
else:
        sys.stderr.write("unsupported classifierType: %s\n" % clfType)
        sys.exit(-1)

if __name__ == '__main__':
        err_file = open(args.probDir + '/train.err', 'w')
        sys.stderr = err_file
        results = joblib.Parallel(n_jobs=args.n_jobs)(joblib.delayed(make_predictor)(n_col) for n_col in pred_to_run)
        err_file.close()
        print("Finished in %0.3f seconds." % (time.time() - stime))
