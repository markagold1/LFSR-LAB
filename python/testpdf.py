# Test generating AWGN via lfsr + masks + box-muller

from matplotlib.pyplot import *
import scipy.io as mat
import numpy as np
import lfsr as l
import boxmuller as bm
import normpdf

# Degree 64 generator polynomial: poly = [64,63,61,60,0]
# Q = Number of bits to approximate distribution
# jump by 1 year (at fs = 50 MHz): jump = np.uint64(365.25 * 86400 * 50e6)
def process(poly = [64,63,61,60,0], Q = 10, jump = 1e10, verbose=1):
    ion()

    # lfsr seq length to generate
    N = 2**18

    # jump, masks and seqs
    jump = np.uint64(jump)
    mask = np.zeros(Q,dtype=np.uint64)
    seq = np.zeros([Q,N],dtype=np.uint64)
    for u in range(0,Q):
        mask[u] = l.jump2mask(jump,poly)
        (seq_, fill_) = l.ssrg_mask(N,poly,1,mask[u])
        seq[u] = np.array(seq_)
        jump += jump

    de = bi2de(seq,'left-msb')
    de += 0.5
    (u,n) = bm.boxmuller(de)
    if verbose:
        rr = np.arange(-10.0,10.0,0.05)
        hist,bins = np.histogram(n,rr)
        delta = bins[1] - bins[0]
        area = np.sum(hist) * delta
        bar(bins[:-1]-delta/2,hist/area,align='edge',width=delta)
        hold
        y = normpdf.normpdf(rr,0,1)
        plot(rr,y,'r')
        ax = axis([-6, 6, 0, 0.42])
        show()
        tstr = 'True and estimated Gaussian PDF '
        tstr += '(' + str(Q) + 'bits)'
        title(tstr)
        legend(('True','Estimated'))

    return(u,n,de,seq)    

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

#mat.savemat('seq.mat',mdict={'seq':seq})
