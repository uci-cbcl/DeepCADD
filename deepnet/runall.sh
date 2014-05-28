#!/bin/bash
# Trains a feed forward net on MNIST.
train_deepnet='python /data/users/dxquang/deepnet/deepnet/trainer.py'
${train_deepnet} model.pbtxt train.pbtxt eval.pbtxt
