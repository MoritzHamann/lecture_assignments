################################################################################
# Machine Learning
# Exercise 10.1 - Sampling
################################################################################

#import numpy as np
import random

# Propability of Battery being bad (=0)
def P_Battery():
    list = [0] * 2 + [1] * 98;
    return random.choice(list);

# Propability of Fuel being empty (=0)
def P_Fuel():
    list = [0] * 5 + [1] * 95;
    return random.choice(list);
    
# Propability of Gauge showing empty (=0) given P_Battery & P_Fuel    
def P_Gauge(B, F):
    list = [];
    if(B == 1 and F == 1):
        list = [0] * 4 + [1] * 96;
    if(B == 1 and F == 0):
        list = [0] * 97 + [1] * 3;
    if(B == 0 and F == 1):
        list = [0] * 10 + [1] * 90;
    if(B == 0 and F == 0):
        list = [0] * 99 + [1] * 1;
    
    return random.choice(list);
    
# Propability of TurnOver (no=0) given P_Battery
def P_TurnOver(B):
    list = [];
    if(B == 1):
        list = [0] * 3 + [1] * 97;
    if(B == 0):
        list = [0] * 98 + [1] * 2;
    
    return random.choice(list);
    
# Propability of Start (no=0) given P_TurnOver and P_Fuel
def P_Start(T, F):
    list = [];
    if(T == 1 and F == 1):
        list = [0] * 1 + [1] * 99;
    if(T == 1 and F == 0):
        list = [0] * 92 + [1] * 8;
    if(T == 0 and F == 1):
        list = [0] * 100;
    if(T == 0 and F == 0):
        list = [0] * 100;
    
    return random.choice(list);
    
# rejection sampling for "obs" with "K" number of samples
def rej_sampling(obs, K):
    samples = [];
    
    #print "S = %d" % obs;
    #print "[B,F,G,T]"
    while(len(samples) < K):
        B = P_Battery();
        F = P_Fuel();
        G = P_Gauge(B,F);
        T = P_TurnOver(B);
        S = P_Start(T,F);
        
        if(S == obs):
            samples.append([B,F,G,T,S]);
            #print [B,F,G,T,S];
        
    return samples;
    
# importance sampling with "K" number of samples
def imp_sampling( K):
    samples = [];
    
    #print "S = %d" % obs;
    #print "[B,F,G,T]"
    while(len(samples) < K):
        B = P_Battery();
        F = P_Fuel();
        G = P_Gauge(B,F);
        T = P_TurnOver(B);
        S = P_Start(T,F);
        
        if(S == 0):
            if(T == 1 and F == 1):
                W = 0.01;
            if(T == 1 and F == 0):
                W = 0.92;
            if(T == 0 and F == 1):
                W = 1.0;
            if(T == 0 and F == 0):
                W = 1.0;
        else:
            if(T == 1 and F == 1):
                W = 0.99;
            if(T == 1 and F == 0):
                W = 0.08;
            if(T == 0 and F == 1):
                W = 0.0;
            if(T == 0 and F == 0):
                W = 0.0;
            
        samples.append([W,B,F,G,T,S])
        
    return samples;
    
# Main
if __name__ == "__main__":

    num_samples = 1000;
    obs = 0;
    
    # compute P(F|S=0) with rejection Sampling 
    samples = rej_sampling(obs, num_samples);
    p_empty = 0.0;
    p_not_empty = 0.0;
    for i in range(0,(len(samples)-1)):
        sampleMember = samples[i];
        if (sampleMember[1] == 0):
            p_empty += 1.0;
        else:
            p_not_empty += 1.0;
            
    p_empty = p_empty / len(samples);
    p_not_empty = p_not_empty / len(samples);
    
    print "Rejection Sampling:"
    print "P(F=empty|S=no) = %f" % p_empty;
    print "P(F=not_empty|S=no) = %f" % p_not_empty;
    
    # compute P(F|S=0) with importance Sampling 
    samples = imp_sampling( num_samples);
    p_empty = 0.0;
    p_not_empty = 0.0;
    for i in range(0,(len(samples)-1)):
        sampleMember = samples[i];
        if(sampleMember[len(sampleMember)-1] == obs):
            sum_weight = sampleMember[0]
            if (sampleMember[2] == 0):
                p_empty += sampleMember[0];
            else:
                p_not_empty += sampleMember[0];
            
    p_empty = p_empty / sum_weight;
    p_not_empty = p_not_empty / sum_weight;
    
    print "Importance Sampling:"
    print "P(F=empty|S=no) = %f" % p_empty;
    print "P(F=not_empty|S=no) = %f" % p_not_empty;