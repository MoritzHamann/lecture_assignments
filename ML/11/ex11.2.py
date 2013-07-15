################################################################################
# Machine Learning
# Exercise 11.2 - Gaussians and Singular Value Decomposition
################################################################################

import numpy as np
import matplotlib.pyplot as plt

def readData(filename):
    X = []
    numOfPoints = 0
    f = open(filename)
    
    for line in f.readlines():
        line = line.strip().split()
        #convert string into float
        line = [float(value) for value in line]
        X.append(line)
        numOfPoints += 1
    
    f.close()
     
    return np.array(X), numOfPoints 

def centerData(X, n, mu):
    
    ones = np.ones((n,1))
    return X - ones*mu

if __name__ == "__main__":
    # read data into lists
    X,n  = readData("./gauss.txt")
    
    # compute mean
    mu = np.array(np.mean(X,axis=0))
    
    # center Data
    X_c = centerData(X,n,mu)
    
    # Covariance Matrix
    C = (1.0/n)*(X_c.T.dot(X_c))
    C_uncentered = (1.0/n)*(X.T.dot(X))-mu.dot(mu.T)
    
    # Singular Value decomposition
    # see: http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.svd.html
    U,s,V = np.linalg.svd(C)
    D = np.diag(s)
    
    # Eigenvalues:
    eigenvalues = []
    for i in range(0, len(D)):
        eigenvalues.append(D[i][i])
        
    # Eigenvektors:
    eigenvectors = []
    for i in range(0, len(V)):
        eigenvectors.append(V[i][:])
    eigenvectors = np.array(eigenvectors)
        
    # Output
    print "Covariance matrix with centered data:"
    print C
    print 
    print "Covariance matrix with uncentered data:"
    print C_uncentered
    print "--"
    print "U:"
    print U
    print "D:"
    print D
    print "V:"
    print V
    print "--"
    print "Eigenvalues:"
    print eigenvalues    
    print "Eigenvectors:"
    for i in range(0, len(eigenvectors)):
        print eigenvectors[i]
        
    # Plot
    # data
    plt.scatter(X[:,0],X[:,1])
    # first segment
    tmp1 = [mu[0], mu[0] + np.sqrt(eigenvalues[0])*eigenvectors[0][0].T]
    tmp2 = [mu[1], mu[1] + np.sqrt(eigenvalues[0])*eigenvectors[0][1].T]
    plt.plot(tmp1,tmp2,color='r')
    
    # second segment
    tmp1 = [mu[0], mu[0] + np.sqrt(eigenvalues[1])*eigenvectors[1][0].T]
    tmp2 = [mu[1], mu[1] + np.sqrt(eigenvalues[1])*eigenvectors[1][1].T]
    plt.plot(tmp1,tmp2,color='g')
    plt.show()
    