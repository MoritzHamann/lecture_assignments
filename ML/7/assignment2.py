def possible_rolls(D1,n):
    """All possible rolls where the first dice is d1 and the sum is n"""
    possibilities = []
    for D2 in range(1,7):
        for D3 in range(1,7):
            if D1+D2+D3 == n:
                possibilities.append((D1,D2,D3))
    return possibilities


if __name__ == "__main__":

    # P(S|D1) should be the same as P(S,D1)/P(D1)
    # where P(S,D1) is the probability of all rolls
    # where the sum is S and the first dice is D1


    # possibility tables
    d1 = [1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0]
    d2 = [1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0]
    d3 = [1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0]
    
    # get the both values
    D1 = raw_input("Enter first dice: ")
    S = raw_input("Enter Sum of dice: ")
    D1=int(D1)
    S=int(S)
    
    # get all rolls, for which the sum is S
    # and the first dice is D1
    rolls = possible_rolls(D1, S)
    
    # calculate possibility for every roll and sum them up
    p=0
    for roll in rolls:
        p+= d1[roll[0]-1]*d2[roll[1]-1]*d3[roll[2]-1]
    
    print("")
    print("Possible rolls are:")
    print(rolls)
    print("====================")
    pd1=0
    # calculate P(D1)
    for D2 in range(5):
        for D3 in range(5):
            pd1+= d1[D1-1]*d2[D2]*d3[D3]
    
    print("Probability is: "+str(p/pd1))
