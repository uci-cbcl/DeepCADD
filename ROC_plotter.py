import numpy as np
import pylab as pl
import sys
from sklearn.utils import shuffle
from sklearn.metrics import roc_curve, auc
probas = np.load(sys.argv[1])
labels = np.load(sys.argv[2])
fpr, tpr, thresholds = roc_curve(labels[:,1], probas[:, 1])
roc_auc = auc(fpr, tpr)
print("Area under the ROC curve : %f" % roc_auc)

# Plot ROC curve
pl.clf()
pl.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
pl.plot([0, 1], [0, 1], 'k--')
pl.xlim([0.0, 1.0])
pl.ylim([0.0, 1.0])
pl.xlabel('False Positive Rate')
pl.ylabel('True Positive Rate')
pl.title('Receiver operating characteristic example')
pl.legend(loc="lower right")
pl.show()
