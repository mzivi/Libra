# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 15:49:38 2015

@author: Michele
"""

import numpy as np
import random as rnd

class LogisticRegression(object):
    """Logistic regression based on sigmoid function and trained via gradient
    descent method.
    """

    def __init__(self, theta=None):
        """Provide the initial values for the parameters. If theta is not given,
        it will be inferred from the design matrix.
        """
        self._theta = theta

    def _sigmoid(self, x):
        return 1 / (1 + np.exp(-x))
    
    @property
    def theta(self):
        return self._theta
    
    def train(self, x, y, alpha, stop_thr, max_iter = 100):
        
        if not self._theta:
            # TODO: this should actually init to random epsilons
            self._theta = np.array([rnd.random() * 0.1 for i in xrange(x.shape[1])])
        # TODO: elsif check theta is consistent with x

        m = x.shape[0]
        h = self._sigmoid(np.dot(x, self._theta))
        cost_fn = np.empty(max_iter)
        cost_fn[0] = np.sum(- y * np.log(h) - (1 - y) * np.log(1 - h)) / m
        for i in xrange(max_iter):
            cost_grad = np.sum(x * (h-y)[:, np.newaxis], axis=0) / m
            self._theta = self._theta - alpha * cost_grad
            h = self._sigmoid(np.dot(x, self._theta))
            cost_fn[i+1] = np.sum(- y * np.log(h) - (1 - y) * np.log(1 - h)) / m
            if np.abs(cost_fn[i] - cost_fn[i+1]) < stop_thr:
                cost_fn = np.resize(cost_fn, i + 1)
                break

        return cost_fn
    
    def predict(self, x, thr=0.5):
        
        if self._theta is None:
            return
            
        return self._sigmoid(np.dot(x, self._theta)) > thr

    
if __name__ == '__main__':
    
    from matplotlib.pyplot import plot
    l = LogisticRegression()
    x = np.array([[1., 1.],[1., 2.], [1., 1.1],[1.,1.9]], dtype=float)
    y = np.array([1, 0, 1, 0], dtype=float).T
    rnd.seed(17)
    cost_fn = l.train(x, y, alpha=0.5, stop_thr=0.0001, max_iter=10000)
    plot(cost_fn)
    print l.theta
    print l.predict(x)