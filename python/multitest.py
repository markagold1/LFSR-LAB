from __future__ import print_function
from matplotlib.pyplot import *
import scipy.io as mat
import numpy as np
import lfsr as l
import bi2de as bd
import boxmuller as bm
import normpdf
import testpdf as test

ion()
sb = 221
rr = np.arange(-10.0,10.0,0.05)
N = list()
poly = [64,63,61,60,0]
jump = 1e10
for Q0 in range(4,20,4):
    (u,n,de,seq) = test.process(poly, Q0, jump, verbose=0)
    N.append(n) 
    hist,bins = np.histogram(n,rr)
    delta = bins[1] - bins[0]
    area = np.sum(hist) * delta
    subplot(sb)
    bar(bins[:-1]-delta/2,hist/area,align='edge',width=delta)
    y = normpdf.normpdf(rr,0,1)
    plot(rr,y,'r')
    ax = axis([-6, 6, 0, 0.42])
    title('Q = ' + str(Q0)) 
    show()
    sb += 1
    print("Q = {0}".format(Q0))
