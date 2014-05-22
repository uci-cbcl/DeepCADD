INSTALL
=======

Required
--------

* [svmlight-loader](https://github.com/mblondel/svmlight-loader). For converting svmlight files to sparse CSR matrices very efficiently.

* [Neural Network Toolbox](https://github.com/IssamLaradji/NeuralNetworks). Contains deep learning Python source code. Meant to mimic [scikit-learn](http://scikit-learn.org/stable/).

All other required libraries should be installed with Enthought Canopy already.

USAGE
=====

Download the four .gz files from [here](http://krishna.gs.washington.edu/martin/download/cadd_training/). Unzip them with the 'gzip -d' command.

Generating svm-light files
--------------------------

You need to first extract features. You need to use the impute2svmlight.py script. This script was supplied by the original CADD authors. From the uncompressed tsv files you downloaded, randomly sample 13,141,299 SNVs, 627,071 insertions and 926,968 deletions from both the simulated variant (Simulated) and observed (HCanc_diff) variant data sets. Hopefully the filenames make it obvious where these variants are located. You can determine the type of variant from the sixth column. You can randomly sample using a combination of the awk and shuf UNIX commands. Put all of the variants into a file for storage. Make sure the file has observed variants first and then simulated variants.

```
$ awk '$6=="SNV"' HCanc_diff_SNVs.tsv > Observed_SNVS.tsv #Store all Observed SNVS
$ awk '$6=="INS"' HCanc_diff_InDels.tsv > Observed_INS.tsv #Store all Observed INS
$ awk '$6=="DEL"' HCanc_diff_InDels.tsv > Observed_DEL.tsv #Store all Observed DEL
$ awk '$6=="SNV"' Simulated_SNVs.tsv > Simulated_SNVS.tsv #Store all Simulated SNVS
$ awk '$6=="INS"' Simulated_InDels.tsv > Simulated_INS.tsv #Store all Simulated INS
$ awk '$6=="DEL"' Simulated_InDels.tsv > Simulated_DEL.tsv #Store all Simulated DEL

$ shuf -n 13141299 Observed_SNVs.tsv > dataset.tsv #sample samples into dataset.tsv
$ shuf -n 627071 Observed_INSs.tsv >> dataset.tsv
$ shuf -n 926968 Observed_DELs.tsv >> dataset.tsv
$ shuf -n 13141299 Simul_SNVs.tsv >> dataset.tsv
$ shuf -n 627071 Simul_INSs.tsv >> dataset.tsv
$ shuf -n 926968 Simul_DELs.tsv >> dataset.tsv

$ cat dataset.tsv | python impute2svmlight.py --noMap --noSegDup -p 0.99 #Convert dataset.tsv into two svmlight files. 99% of the data are stored as training in training.svmlight. The rest are stored in testing.svmlight.
```

To store the data in Python sparse csr_matrix format, you need to write a script or do it from iPython. I will do it from iPython for demonstration purposes:

```
from numpy import *
from scipy.sparse import csr_matrix
from svmlight_loader import load_svmlight_file

#will load the training file. The same can be done for the testing file.
X, y = load_svmlight_file('training.svmlight')
#X is a csr_matrix of the features and y is an array of the labels (-1 for observed, 1 for simulated)
#To save the csr_matrix so that you do not have to load it again, you need to use the numpy save function
#on the three arrays that comprise a csr_matrix

save('training_data.npy',X.data)
save('training_indices.npy',X.indices)
save('training_indptr.npy',X.indptr)

#We can now load back the csr_matrix much faster than we loaded it using svmlight_loader
data = load('training_data.npy')
indices = load('training_indices.npy')
indptr = load('training_indptr.npy')

X = csr_matrix((data,indices,indptr))
```

Running neural network code
---------------------------

Looking into the help to see how to run the Python neural network codes. I wrote a logistic regression and a multilayer perceptron trainer for both scikit-learn and theano. Make sure you have the csr array files for testing and training stored in a dataset folder in '../data/'. For example, suppose I am in the folder containing the Python code, and I have a folder '../data/Impute1/' containing 'training_data.npy', 'training_indices.npy', 'training_indptr.npy', 'testing_data.npy', 'testing_indices.npy', 'testing_indptr.npy'. Then, to run the sklearn logistic regression classifier:

```
$ python sklearn_CADD_sgd.py -scale Impute1
```

Note that the codes use stochastic gradient descent for training. Sklearn has L2-regularization for logistic regression while theano does not, hence the difference in their results. The scale flag tells the code to use the scikit-learn StandardScaler to scale the data before training. This is because training is sensitive to the variance of the data, so we need to make sure the features have unit variance first.
