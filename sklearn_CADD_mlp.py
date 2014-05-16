#!/usr/bin/env python
from argparse import ArgumentParser
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler

from multilayer_perceptron  import MultilayerPerceptronClassifier

import numpy
import scipy.sparse

    
def apply_mlp(datasets, algorithm="l-bfgs", h=40, minibatchsize=200, epoch=5, seed=1):
    ''' Applies the MLP

    :type datasets: list
    :param datasets: List containing training/testing data
    
    :type cores: int
    :param cores: Number of cores
    
    :type seed: int
    :param seed: Random seed
    '''
    print 'Applying', algorithm
    training_X, training_y = datasets[0]
    testing_X, testing_y = datasets[1]
    #clf = SGDClassifier(loss="log", random_state=seed, shuffle=True, n_iter=epoch, verbose=1, n_jobs=cores)
    mlp = MultilayerPerceptronClassifier(n_hidden = h, activation="tanh", batch_size=minibatchsize, algorithm=algorithm, max_iter = epoch,shuffle=True,verbose=True,random_state=seed)
    print 'Fitting...'
    #clf.fit(training_X, training_y)
    mlp.fit(training_X,training_y)
    #predict_y = clf.predict(testing_X)
    print "Accuracy on testing data:", mlp.score(testing_X, testing_y)



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
    parser.add_argument("-s", "--seed", dest="seed", help="Random seed", type=int, default=1)
    parser.add_argument("-m", "--minibatch", dest="minibatch", help="Mini-batch size. Only applied with SGD. Default: 200", type=int, default=200)
    parser.add_argument("-e", "--epoch", dest="epoch", help="Number of epochs to run", type=int, default=5)    
    parser.add_argument("-hidden", "--hidden", dest="hidden", help="Number of hidden units. Default 40", type=int, default=40)    
    parser.add_argument("-scale", "--scale", dest="scale", help="If specified, scale the data", action='store_true')
    parser.add_argument("-a", "--algorithm", dest="algorithm", help="Algorithm to use. sgd or l-bgfs (default)", default="l-bfgs")
    args = parser.parse_args()
    dataset = args.dataset
    datasets = load_data(dataset, args.scale)
    m = int(args.minibatch)
    h = int(args.hidden)
    s = int(args.seed)
    e = int(args.epoch)
    a = args.algorithm
    apply_mlp(datasets, a, h, m, e, s)

    
if __name__=='__main__':
    main()