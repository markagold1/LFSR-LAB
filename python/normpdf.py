import numpy as np

def normpdf(x,mu=0,sigma=1):
    z = (x - mu) / sigma
    sq2pi = 2.5066282746310005024157652848110
    p = np.exp(-z**2 / 2) / (sigma*sq2pi)

    #p[(x == mu) and (sigma == 0)] = np.inf
    #p[np.isinf(z)] = 0
    #p[np.isnan(x) or np.isnan(mu) or (sigma < 0)] = np.nan

    return p

