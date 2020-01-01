import numpy as np

# Binary vector to decimal scalar converter, good to 64-bits.
def bi2de(bin,flag = 'left-msb'):
    if bin.ndim == 2:
        d = np.empty(bin.shape[1])
        for i in xrange(0, bin.shape[1]):
            d[i] = bi2de1(bin[:,i],flag)
    else:
        d = bi2de1(bin,flag)
    return d


def bi2de1(bin,flag = 'left-msb'):
    if flag.lower() == 'right-msb':
        bin = bin[ : : -1]
    return bin.dot(1 << np.arange(bin.shape[-1] - 1, -1, -1,dtype=np.uint64))

