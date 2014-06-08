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
X, y = load_svmlight_file('training.svmlight', dtype='float32')
#X is a csr_matrix of the features and y is an array of the labels (-1 for observed, 1 for simulated)
#To save the csr_matrix so that you do not have to load it again, you need to use the numpy save function
#on the three arrays that comprise a csr_matrix
#It is recommended to save as float32 because GPUs only work with float32. It will save you memory too.

save('training.data.npy',X.data)
save('training.indices.npy',X.indices)
save('training.indptr.npy',X.indptr)
save('training.labels.npy',y)

#We can now load back the csr_matrix much faster than we loaded it using svmlight_loader
data = load('training.data.npy')
indices = load('training.indices.npy')
indptr = load('training.indptr.npy')
y = load('training.labels.npy')

X = csr_matrix((data,indices,indptr))
```

For the deepnet package, you cannot save the sparse matrix as three separate numpy files. You have to save all three of these numpy arrays in a zipped archive. Because the normal savez command in numpy only works for files up to 4 GB, you have to do a little hacking and use a hidden package in numpy as follows:

```
from numpy.lib import npyio
npyio.savez('TrainingFeatures.npz', data=X.data, indices=X.indices, indptr=X.indptr,shape=array(list(X.shape)))
``` 

Running neural network code
---------------------------

Looking into the help to see how to run the Python neural network codes. I wrote a logistic regression and a multilayer perceptron trainer for both scikit-learn and theano. Make sure you have the csr array files for testing and training stored in a dataset folder in '../data/'. For example, suppose I am in the folder containing the Python code, and I have a folder '../data/Impute1/' containing 'training.data.npy', 'training.indices.npy', 'training.indptr.npy', 'training.labels.npy', 'testing.data.npy', 'testing.indices.npy', 'testing.indptr.npy', 'testing.labels.npy'. Then, to run the sklearn logistic regression classifier:

```
$ python sklearn_CADD_sgd.py -scale Impute1
```

Note that the codes use stochastic gradient descent for training. Sklearn has L2-regularization for logistic regression while theano does not, hence the difference in their results. The scale flag tells the code to use the scikit-learn StandardScaler to scale the data before training. This is because training is sensitive to the variance of the data, so we need to make sure the features have unit variance first.

Here is example on how to scale the data, using the same X and y from above:

```
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
scaler.fit(X)  # Don't cheat - fit only on training data
X = scaler.transform(X)
X_test = scaler.transform(X_test)  # apply same transformation to test data
```

deepnet
=======

For the deepnet package, you cannot simply use the same files as before. The features must be saved into a single .npz file, not as three separate matrices for the CSR format. Unfortunately, there is a nasty bug in Python that prevents you from saving zip files bigger than 4GB. If you try to do this directly, you will get an error that says: 'L' format requires 0 <= number <= 4294967295. To fix the bug, you must first find where zipfile.py is located. Here is a sample of how I did it from IPython:

```
In [1]: import zipfile

In [2]: print zipfile
<module 'zipfile' from '/data/users/dxquang/Canopy/appdata/canopy-1.3.0.1715.rh5-x86_64/lib/python2.7/zipfile.py'>
```

Once you located zipfile.py, you must patch. Download the patch file [here] (http://bugs.python.org/file18685/zipfile_zip64_header.patch). Put the patch file in the same folder as zipfile.py, cd into that folder, and type the following:

```
patch zipfile.py < zipfile_zip64_header.patch
```

Now you are ready to save CSR matrices in .npz format. Using the same 'X' from above:

```
In [1]: from numpy.lib import npyio

In [2]: npyio.savez('train.npz', data=X.data, indices=X.indices, indptr=X.indptr, shape=array(list(X.shape)))

```

You will need to do this for validation and test sets as well. Try to do some array slicing to get a validation set. I personally use training/validation/testing sizes of 28800000/296000/294000. Save the labels in appropriately named .npy files as well. I also recommend doing the shuffling ahead of time (sklearn.utils.shuffle can help). 
The labels also have to be preprocessed. Unfortunately, deepnet is not smart enough to differentiate unique labels. It assumes each label is an array index. Due to the way the CADD authors save their data, one of the labels end up being '-1', which of course is an invalid array index. Using the same 'y' as before, replace '-1' labels with the following code:

```
In [1]: from numpy import save

In [2]: y[y == -1] = 0#replace -1 values with 0

In [3]: y = y.astype('int64')#for consistency with the deepnet author data, but I think float32 and int32 will work fine too

In [4]: save('train.newlabels.npy',y)#save the data

```
Put the data into a folder. You will need to modify 4 .pbtxt files. I included all of them in the deepnet folder. cadd.pbtxt specifies the locations of the data and labels and specifies how much cpu/gpu memory to use. eval.pbtxt specifies the batch size used for evaluating the test and validation set. eval.pbtxt can also set the options for the stop condition based on evaluating the testing/validation sets, but I have not explored this. model.pbtxt specifies the modify, such as number of hidden units and layers, the objective function, and the activation units to use. train.pbtxt specifies the training algorithm to use, the stop condition to use, the minibatch size, when to evaluate the test/validation sets, and how often to save the results. Each step evaluates one minibatch, so you will have to do some math to determine how many steps corresponds to an epoch. train.pbtxt also references the proto file cadd.pbtxt. When everything looks good, run with './runall.sh'. 


LIBOCAS
=======

This is the SVM algorithm the original CADD authors used. Download the source code from [here] (http://cmp.felk.cvut.cz/~xfrancv/ocas/html/). Hopefully it is straightforward for you to compile and run their test examples. You can use the training.svmlight and testing.svmlight files you generated earlier to train this. The original authors used C parameter of 0.0025. You can also run it with multiple threads. I also included my LIBOCAS output for reference.

```
./svmocas -c 0.0025 training.svmlight cadd.model
./linclassif -e -o myoutput.pred testing.svmlight cadd.model
```

The first line trains a model and saves it to cadd.model. The second line uses the model saved in cadd.model to predict the testing data and outputs the results to myoutput.pred.

ROC Analysis
============

I modified two files in the deepnet library: util.py and neuralnet.py. I uploaded the modified python scripts onto the github. Replace the files in your deepnet/deepnet folder with these two files. To make deepnet output the posterior probabilities, you need to set the compute_MAP boolean to true in the output layer's performance states in the model.pbtxt file. Your model.pbtxt file should have these lines now:

```
performance_stats {
  compute_cross_entropy: true
  compute_correct_preds: true
  compute_MAP: true  
}
```

I modified util.py and neuralnet.py so that deepnet will not actually calculate the MAP, but instead output the posterior probabilities and the correct labels in the checkpoint directory. It will do this for both the validation set and testing set. The files will have the same name as their corresponding checkpoint file, plus "_valid_preds.npy", "_valid_targets.npy", "_test_preds.npy", "_test_targets.npy" appended at the end of the filename. the "_preds.npy" files contain the posterior probabilities of each sample while the "_targets.npy" files contain the correct labels for each sample. I've also included a script that will accept two .npy files and plot a ROC curve for you. Go [here] (http://scikit-learn.org/stable/auto_examples/plot_roc.html) if you would like to learn more about generating ROC curves.

Here is an example of me using the ROC plotting script on the last checkpoint test file arrays. In theory, this should work for the best output as well, but I have not checked this. 

```
$ python ROC_plotter.py cadd_3layer_relu_train_op_LAST_test_preds.npy cadd_3layer_relu_train_op_LAST_test_targets.npy
Area under the ROC curve : 0.673855 
```

Hopefully you can get an area better than what I got.
