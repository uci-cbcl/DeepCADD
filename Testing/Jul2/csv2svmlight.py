from numpy import *
from scipy.sparse import csr_matrix

from svmlight_loader import dump_svmlight_file

clinvar = loadtxt('clinvar_imputed.csv',delimiter=',')
ESP6500 = loadtxt('ESP6500_imputed.csv', delimiter=',')

clinvar = clinvar[:,1:]
ESP6500 = ESP6500[:,1:]

clinvar_samples = clinvar.shape[0]
ESP6500_samples = ESP6500.shape[0]

y = ones(clinvar_samples+ESP6500_samples, dtype='int64')
y[0:ESP6500_samples] = -1

X = vstack((ESP6500,clinvar))
X = csr_matrix(X)

dump_svmlight_file(X, y, "Jul2_testing.svmlight", zero_based=False)
