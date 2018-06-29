from __future__ import print_function  # <-- this lets python2 use python3's print function

import sys, time, os, subprocess
from distutils import dir_util  # handle dirs
import argparse
import shutil
import numpy as np
from sklearn.externals import joblib

def run_predictor(n_col):
        y = X_all[:, n_col]
        X = np.delete(X_all, n_col, axis=1)
        role_ID = col_names[n_col, 1]

        clfDir = args.trainDir + "/Predictors/" + role_ID + "/Classifiers/" + args.clfType
        clfFile = clfDir + "/"  + args.clfType
        ldaFile = clfDir + "/LDA_vars"

        if args.LDA:
                vars_to_use = np.genfromtxt(ldaFile, delimiter='\t', dtype=int)
                X_feat = X[:,vars_to_use]
                if vars_to_use.shape[0] == 1:
                        X_feat = X_feat.reshape(-1, 1)

        else:
                X_feat = X

        clf = joblib.load(clfFile)
        results = clf.predict(X_feat)
        results = results.tolist()
        results.append(n_col)
        print("Completed %d: %s predictor." % (n_col, role_ID))
        return results


stime = time.time()

parser = argparse.ArgumentParser(description='Evaluate a list of gtos for consistency.')
parser.add_argument("trainDir", help="Directory which contains built predictors",
                    type=str)
parser.add_argument("testDir", help="Directory for output files",
                    type=str)
parser.add_argument("gtoList", help="Plain text list of gto files to evaluate",
                    type=str)
parser.add_argument("roles_in_subsystems", help="List of roles in subsystems",
                    type=str)
parser.add_argument("roles_to_use", help="List of roles to use in evaluations",
                    type=str)
parser.add_argument("-c", "--classifier", dest="clfType", default="RandomForestClassifier", help="Type of sklearn classifier to use")
parser.add_argument("--clear", action="store_true",
                  help="Clear files from testDir so that new output can be stored")
parser.add_argument("--LDA", action="store_true",
                  help="Use Linear Discriminant Analysis")
args = parser.parse_args()
if __name__ == '__main__':
    print("Using predictors from " + args.trainDir + ".")
    print("Output will be in " + args.testDir + ".")
    if not os.path.isdir(args.testDir):
        sys.stderr.write("not a valid testDir: %s\n" % (args.testDir))
        sys.exit(-1)

    if args.clear:
            print("Clearing " + args.testDir + "/summaries...")
            if os.path.isdir(args.testDir + "/summaries"):
                    shutil.rmtree(args.testDir + "/summaries")
            if os.path.isdir(args.testDir + "/predictions"):
                    shutil.rmtree(args.testDir + "/predictions")
            retcode = subprocess.call(["gtos_to_matrix", args.gtoList, args.testDir, args.roles_in_subsystems, args.roles_to_use, "--clear"])

            if retcode:
                    sys.exit(-1)
            print("Generation of matrix finished in %0.2f seconds." % (time.time() - stime))

X_all = np.genfromtxt(args.trainDir + '/X', delimiter='\t')
col_names = np.genfromtxt(args.trainDir + '/col.h', delimiter='\t', dtype=str)
genomes = np.genfromtxt(args.trainDir + '/row.h', delimiter='\t', dtype=str)


#...genfromtxt reformats a 1-row dataset as a vector, not an array,
#   so the number of dimensions must be reformatted if not 2D array.
if len(X_all.shape) != 2:
        X_all = X_all.reshape(1,-1)
if len(genomes.shape) != 2:
        genomes = genomes.reshape(1, -1)

#X_all[X_all > 5.] = 6.

if __name__ == '__main__':

    predictions = joblib.Parallel(n_jobs=32)(joblib.delayed(run_predictor)(n_col) for n_col in range(X_all.shape[1]))

    predictions = np.asarray(predictions)
    predictions = np.transpose(predictions)
    col_ind = np.asarray(predictions[-1], dtype=int)
    predictions = np.delete(predictions, -1, axis=0)
    predictions = predictions[:,col_ind]

    real_present = X_all > 0
    pred_present = predictions > 0

    coarse_const = 100.0*np.mean(real_present == pred_present, axis=1)
    fine_const = 100.0*np.mean(X_all == predictions, axis=1)
    score_table = genomes[:,1].reshape(-1, 1)
    score_table = np.hstack((score_table, coarse_const.round(2).reshape(-1,1)))
    score_table = np.hstack((score_table, fine_const.round(2).reshape(-1,1)))

    if not os.path.isdir(args.testDir + "/summaries"):
            os.mkdir(args.testDir + "/summaries")
    #if not os.path.isdir(args.testDir + "/predictions"):
    #	os.mkdir(args.testDir + "/predictions")

    for n_row in range(X_all.shape[0]):
            gtoID = genomes[n_row, 1]
            gto_sum_file = args.testDir + "/summaries/" + gtoID + ".out"
            summary = ["Coarse Consistency: " + str(np.round(coarse_const[n_row], decimals = 1))]
            summary.append("Fine Consistency: " + str(np.round(fine_const[n_row], decimals = 1)))
            for n_col in range(X_all.shape[1]):
                    n_pred = predictions[n_row, n_col]
                    n_real = X_all[n_row, n_col]

                    if n_pred != n_real:
                            summary.append(col_names[n_col, 1] + "\t" + str(X_all[n_row, n_col]) + "\t" + str(predictions[n_row, n_col]))

            np.savetxt(gto_sum_file, summary, fmt="%s", delimiter='\t')

    pred_file = args.testDir + "/predictions"
    score_file = args.testDir + "/scores"
    np.savetxt(score_file, score_table, fmt="%s", delimiter='\t')
    np.savetxt(pred_file, predictions, fmt="%d", delimiter='\t')
    print("Finished %d evaluations with %d roles in %0.2f seconds." % (predictions.shape[0], predictions.shape[1], time.time()-stime))
