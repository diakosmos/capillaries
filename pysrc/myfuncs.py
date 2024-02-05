def myfrange(a, b, n):
    L = [0.0] * n        # <--- allocates memory for list; faster than appending/realloc
    nm1 = n - 1
    nm1inv = 1.0 / nm1   # <---do float division only once, since divide is slower than multiply
    for i in range(n):
        L[i] = nm1inv * (a*(nm1 - i) + b*i) #<--- symmetric under swap (a <--> b), (i + 1 <--> n - i)
    return L

def rtbis(func,x1,x2,xacc):
    # root finding by bisection, straight out of numerical recipes
    jmax = 40
    fmid = func(x2)
    f = func(x1)
#    if(f*fmid > 0) error
    if (f < 0):
        rtbis = x1
        dx = x2-x1
    else:
        rtbis = x2
        dx = x1-x2
    for j in range(1,jmax):
        dx = dx*0.5
        xmid = rtbis + dx
        fmid = func(xmid)
        if (fmid < 0):
            rtbis = xmid
        if ((abs(dx) < xacc) or (fmid ==0)):
            return rtbis
            break
         
