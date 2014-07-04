from numpy import *
from scipy.sparse import csr_matrix

from svmlight_loader import dump_svmlight_file

positives = loadtxt('positives_imputed.csv',delimiter=',')
negatives = loadtxt('negatives_imputed.csv', delimiter=',')

positives = positives[:,1:]
negatives = negatives[:,1:]

positives_samples = positives.shape[0]
negatives_samples = negatives.shape[0]

y = ones(positives_samples+negatives_samples, dtype='int64')
y[0:negatives_samples] = -1

X = vstack((negatives,positives))
X = csr_matrix(X)

dump_svmlight_file(X, y, "Jul2_training.svmlight", zero_based=False)
