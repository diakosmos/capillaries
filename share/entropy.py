import math

def SJent(f):
    ii = range(len(f))
    e = 0
    for i in ii:
        ff = f[i]
        e += -ff * math.log( ff )
    return e
