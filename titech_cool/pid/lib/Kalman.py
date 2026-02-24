#/usr/bin/python3

# code based on:
# https://qiita.com/hanon/items/7f03621414c59f06d7ca

import pandas as pd
import numpy as np

import math
import time

np.random.seed( int( time.time() ) )
pd.set_option('display.max_columns', None)

class KalmanSmoother:

    def __init__( self, y, w, v ):
    
        self.G = np.array([[1]])
        self.F = np.array([[1]])
        self.W = np.array([[w]])
        self.V = np.array([[v]])
        self.size = len(y)
        self.K = 0
        
        # w is the process noise, assumed to be known
        # w = np.random.multivariate_normal(np.zeros(1), self.W, size+K)
        
        # v is the observation noise, assumed to be known
        # v = np.random.multivariate_normal(np.zeros(1), self.V, size+K)
        
        # zero array with the size of (size+K) = 105
        self.y = y
        
        # generate true series x and observation series y
        # using the process and observable noises w, v
        # x is not observable
        
        self.m0 = np.array([[0]])
        self.C0 = np.array([[1e7]])
        
        
        self.m = np.zeros((self.size, 1))
        self.C = np.zeros((self.size, 1, 1))
        self.s = np.zeros((self.size, 1))
        self.S = np.zeros((self.size, 1, 1))


    def kalman_filter(self, m, C, y):
        a = self.G @ m
        R = self.G @ C @ self.G.T + self.W
        f = self.F @ a
        Q = self.F @ R @ self.F.T + self.V
        K = (np.linalg.solve(Q.T, self.F @ R.T)).T
        m = a + K @ (y - f)
        C = R - K @ self.F @ R
        return m, C

    def kalman_prediction(self, a, R):
        a = self.G @ a
        R = self.G @ R @ self.G.T + self.W
        return a, R

    def kalman_smoothing(self, s, S, m, C):
        a = self.G @ m
        R = self.G @ C @ self.G.T + self.W
        A = np.linalg.solve(R, C @ self.G.T)
        s = m + A @ (s - a)
        S = C + A @ (S - R) @ A.T
        return s, S


    def process(self):
        
        for t in range((self.size)):
            if t == 0:
                self.m[t], self.C[t] = self.kalman_filter(self.m0, self.C0, self.y[t:t+1])
            else:
                self.m[t], self.C[t] = self.kalman_filter(self.m[t-1:t], self.C[t-1:t], self.y[t:t+1])

        for t in range((self.size)):
            t = (self.size) - t - 1
            if t == (self.size) - 1:
                self.s[t] = self.m[t]
                self.S[t] = self.C[t]
            else:
                self.s[t], self.S[t] = self.kalman_smoothing(self.s[t+1], self.S[t+1], self.m[t], self.C[t])
                
        
        return [ k[0] for k in self.s ]


if __name__ == '__main__':
    w = np.random.multivariate_normal(np.zeros(1), [[0.25]], 100)
    v = np.random.multivariate_normal(np.zeros(1), [[1]], 100)
    x = np.zeros(100)
    y = np.zeros(100)
    x0 = 0
    
    x[0] = x0 + w[0]
    y[0] = x[0] + v[0]
    for t in range(100):
        x[t] = x[t-1] + w[t]
        y[t] = x[t] + v[t]

    out = KalmanSmoother( list(y), 0.25, 1 ).process()
    print( list(y) )
    print( out )
