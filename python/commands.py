from matplotlib.pyplot import *
import scipy.io as mat
import numpy as np
import lfsr as l
import bi2de as bd
import boxmuller as bm
import normpdf
ion()

N = 2**18
poly = [64,63,61,60,0]
jump = np.uint64(365.25 * 86400 * 50e6)
mask = l.jump2mask(jump,poly)
sfill = 1
sfill = l.ssrgmask2ssrg(poly,sfill,mask)
(seq, fill) = l.ssrg(N,poly,sfill)
Q = 12
ss = seq[3*(Q-1)::]
L = ss.size
smat = ss
for i in range(Q-2,-1,-1):
    smat = np.vstack((smat,seq[3*i:3*i+L]))

#smat2 = np.flipud(smat)
smat2 = smat
u0 = bd.bi2de(smat2)
(u,n) = bm.boxmuller(u0)

figure(1)
rr = np.arange(-10.0,10.0,0.05)
hist,bins = np.histogram(n,rr)
delta = bins[1] - bins[0]
area = np.sum(hist) * delta
bar(bins[:-1]-delta/2,hist/area,align='edge',width=delta)
y = normpdf.normpdf(rr,0,1)
plot(rr,y,'r')
ax = axis([-6, 6, 0, 0.42])
show()

# uniform deviates
figure(2)
rru = np.arange(0.0,1.0,0.005)
histu,binsu = np.histogram(u,rru)
delta_u = binsu[1] -binsu[0]
area_u = np.sum(histu) * delta_u
bar(binsu[:-1]-delta_u/2,histu/area_u,align='edge',width=delta_u)
axu = axis([0, 1, 0, 1.2])


# test ssrgmask2ssrg
N = 2**18
poly = [64,63,61,60,0]
jump = np.uint64(365.25 * 86400 * 50e6)
mask = l.jump2mask(jump,poly)
mfill = 1
(seq, fill) = l.ssrg_mask(N,poly,1,mask)
sfill = l.ssrgmask2ssrg(poly,1,mask)
(sseq, fill) = l.ssrg(N,poly,sfill)
any(seq != sseq)


