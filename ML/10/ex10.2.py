################################################################################
# Machine Learning
# Exercise 10.2 - Sampling from a Gaussian
################################################################################

import numpy as np
import random

#draw gaussian random variable x
def draw_x():
    return np.random.normal(0.0,1.0);
    
#draw binary variable y
def draw_y(x):
    list = [];
    
    if( x > 0):
        list = [1]*90 + [0]*10;
    else:
        list = [1]*10 + [0]*90;
        
    return random.choice(list);
    
 # rejection sampling for "obs" with "K" number of samples
def rej_sampling(obs, K):
    samples = [];
    
    while(len(samples) < K):
        X = draw_x();
        Y = draw_y(X);
        
        if(Y == obs):
            samples.append([X,Y]);
            #print [X,Y];
        
    return samples;
    
# Main
if __name__ == "__main__":

    num_samples = 1000;
    obs = 1;
    
    samples = rej_sampling(obs,num_samples);
    
    p_x = 0.0;
    for i in range(0, num_samples-1):
        sampleMember = samples[i];
        p_x +=  sampleMember[0];

    mean = p_x / num_samples;
    
    print "Mean: %f" % mean;