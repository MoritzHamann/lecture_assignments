def independent(x,y,p,dims):
    # checks for a specific entry (x,y) if row * column = p(x,y)
    sum_over_x = 0
    sum_over_y = 0
    
    for i in range(dims[0]):
        sum_over_x += p[i][y]
        
    for i in range(dims[1]):
        sum_over_y += p[x][i]
    
    # use round as a dirty hack, to avoid issues with
    # floating point precision
    return round(sum_over_x * sum_over_y,4) == p[x][y]


if __name__ == "__main__":
    dims = (2,3)
    p = [[0.08, 0.12, 0.2],[0.12, 0.18, 0.3]]
    i = True
    
    # if for all entries in the table the values are independent,
    # both probability variables are independent
    for x in range(dims[0]):
        for y in range(dims[1]):
            if not independent(x,y,p,dims):
                print("X: "+str(x)+" Y: "+str(y)+" are not independent")
                i = False
                
    if i:
        print("Variables are independent")
    else:
        print("Variables are not independent")
