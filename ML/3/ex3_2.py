################################################################################
# Machine Learning
# Exercise 3.2 - Logistic Regression
################################################################################

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import sys

# load data from single file
def loadData1(inputFile):
    x = np.loadtxt(inputFile, usecols = (0,1))
    tmp, y = np.loadtxt(inputFile, usecols =(1,2), unpack=True)
    return x, np.transpose(np.matrix(y))

# load data from two files
def loadData2(inputFile1, inputFile2):
    x = np.loadtxt(inputFile1)
    y = np.loadtxt(inputFile2, unpack=True)
    return x, np.transpose(np.matrix(y))    
    
# compute linear features
def linfeatures(x):
    X = np.insert(x, 0, 1, axis=1)
    return X
 
# compute quadratic features
def quadfeatures(x):
    n = 0
    for row in x:
        tmp = []
        for w in row:
            for z in row:
                tmp = np.append(tmp,w*z)
                
        row = np.insert(row, 0, 1)
        row = np.insert(row, 3, tmp)
        if n == 0:
            X = row
            n = 1
        else:
            X = np.vstack([X,row])
    return X
    
# calculate features    
def calculateFeatures(x,func):
    return func(x)

# generate Test points
def generateTestPoints(n, dimensions, minimum, maximum):
    retVal = []
    for i in range(n):
        row = []
        for j in range(dimensions):
            column =  random.uniform(minimum, maximum)
            row = np.insert(row, j, column)
        retVal = np.insert(x, i, row, axis=0)
        
    return retVal
    
# calculate optimal Beta
def calcOptBeta(X, y, lam, beta0reg=False):
    # get dimension of x to create identity matrix
    L = np.identity(X.shape[1]) * lam
    if not beta0reg:
        L[0,0] = 0
   
    beta = np.linalg.inv(X.T.dot(X)+L).dot(X.T).dot(y)
    return beta

# calculate Neg-log-Likelihood
def calcMeanNegLogLikelihood(X, y, lam, beta):
    n = len(y)
    for i in range(0, n-1):
        condProb = (1/(np.exp(-X[i].dot(beta))+1))
        negLogLike = y[i]*np.log10(condProb)+(1-y[i])*np.log10(1-condProb)
    
    negLogLike = negLogLike - lam*np.power(np.linalg.norm(beta),2)
    return (1.0/float(n))*negLogLike.item(0)

# plot data
def plot(xs, ys, zs, beta, feature, title):
    fig = plt.figure()
    fig.clf()
    ax = Axes3D(fig)
    
    ax.set_xlabel('X1')
    ax.set_ylabel('X2')
    ax.set_zlabel('Y')
    
    # Plot Data
    for i in range(0,len(xs)-1):
        col = 'r'
        if zs[i] == 1.0:
            col = 'r'
        if zs[i] == 0.0:
            col = 'b'
        ax.scatter(xs[i], ys[i], zs[i], c=col)
    
    # Plot probability function
    x,y = np.meshgrid(np.linspace(xs.min()-.5, xs.max()+.5, 100), np.linspace(ys.min()-.5, ys.max()+.5, 100))
    
    if( feature == "lin" ):
        tmp = beta.item(0) - beta.item(1)*x + beta.item(2)*y
    else:
        tmp = beta.item(0) - beta.item(1)*x + beta.item(2)*y + beta.item(3)*np.power(x,2) + beta.item(4)*x*y + beta.item(5)*x*y + beta.item(6)*np.power(y,2)
        
    z = (1/(np.exp(-tmp)+1))
    ax.plot_wireframe(x, y, z, rstride=10, cstride=10)
    
    plt.title(title)
    plt.show()
    
# Main
if __name__ == "__main__":
    offset = 0
    if(len(sys.argv) == 5):
        offset = 1
        x,y = loadData2(sys.argv[1], sys.argv[2])
    else:
        x,y = loadData1(sys.argv[1])
    
    if sys.argv[(2+offset)] == "lin":
        X = calculateFeatures(x,linfeatures)
    if sys.argv[(2+offset)] == "quad":
        X = calculateFeatures(x,quadfeatures)
        
    optBeta = calcOptBeta(X,y,float(sys.argv[3+offset]))
    meanLogLike = calcMeanNegLogLikelihood(X,y,float(sys.argv[3+offset]),optBeta)
    
    print "Optimal Beta:"
    print optBeta
    print "Mean Neg-Log-Likelihood:"
    print meanLogLike
        
    # samples to plot values for calculated beta
    # only for ex 3.2
    if offset == 0:
        plot(x[:,0],x[:,1],y, optBeta, sys.argv[2+offset], "Class Data")
    
        testx = generateTestPoints(100,2,-3,3)
        if sys.argv[2+offset] == "lin":
            testX = calculateFeatures(testx, linfeatures)
        if sys.argv[2+offset] == "quad":
            testX = calculateFeatures(testx, quadfeatures)
        testX = np.matrix(testX)
        testy = testX.dot(optBeta)
        
        plot(testx[:,0],testx[:,1],testy, optBeta, sys.argv[2+offset], "Test Data")
    