#!/usr/bin/env python
from argparse import ArgumentParser
from sklearn.utils import shuffle
from sklearn.linear_model import SGDClassifier
from sklearn.preprocessing import StandardScaler

import numpy
import scipy.sparse

def apply_minibatch_sgd(datasets, minibatch, epoch=5, cores=1, seed=1):
    ''' Applies the logistic regression sgd method

    :type datasets: list
    :param datasets: List containing training/testing data
    
    :type minibatch: int
    :param minibatch: minibatch size
        
    :type cores: int
    :param cores: Number of cores
    
    :type seed: int
    :param seed: Random seed
    '''
    print 'Applying mini-batch SGD with mini-batch size of ', minibatch
    training_X, training_y = datasets[0]
    testing_X, testing_y = datasets[1]
    print 'Shuffling training data'
    training_X, training_y = shuffle(training_X, training_y, random_state = seed)
    clf = SGDClassifier(loss="log", random_state=seed, n_iter=epoch, verbose=0, n_jobs=cores)
    classes = numpy.unique([-1, 1])
    minibatches = training_X.shape[0]/minibatch + 1
    samples = training_X.shape[0]
    for i in xrange(epoch):
        print "Epoch ", i+1
        for j in xrange(minibatches):
            clf.partial_fit(training_X[j*minibatch:min(samples,(j+1)*minibatch)], training_y[j*minibatch:min(samples,(j+1)*minibatch)], classes=classes)
        print "Accuracy on testing data:", clf.score(testing_X, testing_y)
        #training_X, training_y = shuffle(training_X, training_y)
    
def apply_sgd(datasets, epoch=5, cores=1, seed=1):
    ''' Applies the logistic regression sgd method

    :type datasets: list
    :param datasets: List containing training/testing data
    
    :type cores: int
    :param cores: Number of cores
    
    :type seed: int
    :param seed: Random seed
    '''
    print 'Applying SGD'
    training_X, training_y = datasets[0]
    testing_X, testing_y = datasets[1]
    clf = SGDClassifier(loss="log", random_state=seed, shuffle=True, n_iter=epoch, verbose=1, n_jobs=cores)
    print 'Fitting...'
    clf.fit(training_X, training_y)
    print "Accuracy on testing data:", clf.score(testing_X, testing_y)


def load_data(dataset, scale=False):
    ''' Loads the dataset

    :type dataset: string
    :param dataset: The folder in ../data/ containing the training/testing numpy arrays
    '''

    print '... loading data'
    path = "../data/" + dataset + "/"
    
    #training set
    trainingData = numpy.load(path + "training.data.npy") 
    trainingIndices = numpy.load(path + "training.indices.npy")
    trainingIndptr = numpy.load(path + "training.indptr.npy")
    training_y = numpy.load(path + "training.labels.npy")
    training_X = scipy.sparse.csr_matrix((trainingData, trainingIndices, trainingIndptr))

    #testing set
    testingData = numpy.load(path + "testing.data.npy") 
    testingIndices = numpy.load(path + "testing.indices.npy")
    testingIndptr = numpy.load(path + "testing.indptr.npy")
    testing_y = numpy.load(path + "testing.labels.npy")
    testing_X = scipy.sparse.csr_matrix((testingData, testingIndices, testingIndptr))

    #scale the data 
    if scale:
        print "..training scaler"
        scaler = StandardScaler(with_mean=False)
        scaler.fit(training_X)
        print "..scaling features"
        training_X = scaler.transform(training_X)
        testing_X = scaler.transform(testing_X)
    
    return [(training_X, training_y),(testing_X, testing_y)]

def main():
    description = "Apply stochastic gradient descent to CADD training data"
    parser = ArgumentParser(description=description)
    parser.add_argument('dataset', metavar='dataset', help='Folder in data/ containing training/testing set')
    parser.add_argument("-p", "--parallel", dest="parallel", help="Number of cores to use", type=int, default=1)
    parser.add_argument("-s", "--seed", dest="seed", help="Random seed", type=int, default=1)
    parser.add_argument("-m", "--minibatch", dest="minibatch", help="Mini-batch size", type=int, default=1)
    parser.add_argument("-e", "--epoch", dest="epoch", help="Number of epochs to run", type=int, default=5)    
    parser.add_argument("-scale", "--scale", dest="scale", help="If specified, scale the data", action='store_true')
    args = parser.parse_args()
    dataset = args.dataset
    datasets = load_data(dataset, args.scale)
    m = int(args.minibatch)
    p = int(args.parallel)
    s = int(args.seed)
    e = int(args.epoch)
    if m < 2:
        apply_sgd(datasets, e, p, s)
    else:
        apply_minibatch_sgd(datasets, m, e, p, s)
    
if __name__=='__main__':
    main()