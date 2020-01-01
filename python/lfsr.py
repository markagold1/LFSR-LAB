'''
LFSR generator utilities.
Simple shift register (SSRG) aka Fibonacci generator model:

   +--------(+)<-----(+)-----(+)<-------+
   | r       ^ r-1    ^ r-2   ^         |
   |z        |z       |z      |z        |1
   +-->|r-1|-+->|r-2|-+- ... -+->| 0 |--+------> seq

Modular shift register (MSRG) aka Galois generator model:

   +<----------+-----------+----------+----------+
   |  r        |  r-1      |  r-2     |          ^
   | z         v z         v z        v z        | 1
   +-->|r-1|->(+)->|r-2|->(+)- ...-->(+)->| 0 |--+---> seq

Terminology:

  poly - Generator polynomial, a degree r polynomial 
         in z whose coefficients represent the non-
         zero tap connections to the xor feedback
         network.  In this package poly is a vector 
         of powers in z with non-zero coefficients.
         For example poly = [5,3,0] represents
         generator polynomial: g(z) = z^5 + z^3 + 1
  taps - A binary vector whose non-zero elements 
         represent the coefficients of the generator
         polynomial.  For example the taps vector
         corresponding to poly = [5,3,0] is given
         by taps = [1 0 1 0 0 1]
  fill - Shift register contents, in this package
         fill is represented as a shift-right integer
         whose left-msb binary equivalent represents
         the bit values in the shift register.  For
         example fill=1 corresponds to a "1" in
         the rightmost shift register stage while
         fill=8 corresponds to a "1" two stages to 
         the left of the rightmost stage.
  mask - Specifies the shift register stages to
         combine in order to effect a code phase change.
         In this package mask is an integer whose
         left-msb binary equivalent represents the 
         shift register stages combined.  For example
         if mask=7 the contents of the three rightmost
         stages are modulo-2 summed to form the LFSR 
         output sequence.
  num -  An input parameter to the generator functions
         specifying the length of sequence to generate.

Masked SSRG model:

   +--------(+)<-----(+)-----(+)<-------+
   | r       ^ r-1    ^ r-2   ^         |
   |z        |z       |z      |z        |1
   +-->|r-1|-+->|r-2|-+- ... -+->| 0 |--+
             |        |       |         |
             v        v       v         v
 A=AND gate +-+      +-+     +-+       +-+
            |A|      |A|     |A|       |A|
       r-1/ +-+ r-2/ +-+  1/ +-+    0/ +-+
 mask ---+---|----+---|---+---|-----+   |
             |        |       |         | 
             v        v       v         v
          +--------------------------------+
          |             XOR                |
          +--------------+-----------------+
                         |
                         +------> seq

Masked MSRG model:

    +<----------+-----------+----------+----------+
    |  r        |  r-1      |  r-2     |          ^
    | z         v z         v z        v z        |1
    +-->|r-1|->(+)->|r-2|->(+)- ...-->(+)->| 0 |--+
             |           |          |             |
             v           v          v             v
 A=AND gate +-+         +-+        +-+           +-+
            |A|         |A|        |A|           |A|
       r-1/ +-+    r-2/ +-+     1/ +-+        0/ +-+
 mask ---+---|-------+---|------+---|---------+   |
             |           |          |             | 
             v           v          v             v
          +------------------------------------------+
          |                   XOR                    |
          +--------------------+---------------------+
                               |
                               +------> seq
'''
import numpy as np

# LFSR generator using simple shift register (SSRG) 
# aka Fibonacci structure
def ssrg(num=16,poly=[4,3,0],ifill=8):
    poly = np.asarray(poly)
    degree = poly[0]
    seq = np.empty(num,dtype=np.uint64)
    taps = np.zeros(degree + 1,dtype=np.uint64)
    taps[degree - poly] = int(1)
    sr = de2bi(ifill,degree,'left-msb')

    for nn in range(num):
        seq[nn] = sr[-1]
        parity = (taps[1:] & sr).sum() % np.uint64(2)
        sr = np.hstack((np.uint64(parity),sr[0:-1]))

    # final fill
    fill = bi2de(sr,'left-msb')
    return (seq,fill)


# Masked LFSR generator using simple shift register (SSRG)
# aka Fibonacci structure
def ssrg_mask(num=16,poly=[4,3,0],ifill=8,mask=15):
    poly = np.asarray(poly)
    degree = poly[0]
    seq = np.empty(num,dtype=np.uint64)
    taps = np.zeros(degree + 1,dtype=np.uint64)
    taps[degree - poly] = int(1)

    sr = de2bi(ifill,degree,'left-msb')
    ma = de2bi( mask,degree,'left-msb')

    for nn in xrange(num):
        seq[nn] = (sr & ma).sum() % np.uint64(2)
        parity = (taps[1:] & sr).sum() % np.uint64(2)
        sr = np.hstack((np.uint64(parity),sr[0:-1]))

    # final fill
    fill = bi2de(sr,'left-msb')
    return (seq,fill)

# LFSR generator using modular shift register (MSRG) 
# aka Galois structure
def msrg(num=16,poly=[4,3,0],ifill=8):
    poly = np.asarray(poly)
    degree = poly[0]
    seq = np.empty(num,dtype=np.uint64)
    taps = np.zeros(degree + 1,dtype=np.uint64)
    taps[degree - poly] = int(1)
    sr = de2bi(ifill,degree,'left-msb')

    for nn in range(num):
        seq[nn] = sr[-1]
        if sr[-1]:
            sr = np.hstack((np.uint64(1),xor(sr[0:-1],taps[1:-1]))) 
        else:
            sr = np.hstack((np.uint64(0),sr[0:-1]))

    # final fill
    fill = bi2de(sr,'left-msb')
    return (seq,fill)


# Masked LFSR generator using modular shift register (MSRG)
# aka Galois structure
def msrg_mask(num=16,poly=[4,3,0],ifill=8,mask=15):
    poly = np.asarray(poly)
    degree = poly[0]
    seq = np.empty(num,dtype=np.uint64)
    taps = np.zeros(degree + 1,dtype=np.uint64)
    taps[degree - poly] = int(1)
    sr = de2bi(ifill,degree,'left-msb')
    ma = de2bi( mask,degree,'left-msb')

    for nn in range(num):
        seq[nn] = (sr & ma).sum() % np.uint64(2)
        if sr[-1]:
            sr = np.hstack((np.uint64(1),xor(sr[0:-1],taps[1:-1]))) 
        else:
            sr = np.hstack((np.uint64(0),sr[0:-1]))

    # final fill
    fill = bi2de(sr,'left-msb')
    return (seq,fill)

# State space generator formulation of MSRG
def ssgm(num,poly,ifill):
    poly = np.asarray(poly)
    degree = poly[0]
    seq = np.empty(num,dtype=np.uint64)
    taps = np.zeros(degree + 1,dtype=np.uint64)
    taps[degree - poly] = int(1)
    sr = de2bi(ifill,degree,'left-msb')
    # form msrg characteristic matrix Tm
    pm = taps[-2::-1]
    pm = pm.reshape(pm.size,1)
    Tm = np.vstack((np.eye(degree-1,dtype=np.uint64),np.zeros((1,degree-1),dtype=np.uint64)))
    Tm = np.flipud(np.fliplr(np.hstack((pm,Tm))))

    for nn in xrange(num):
        seq[nn] = sr[-1]
        sr = mvMultGF2(Tm,sr)

    # final fill
    fill = bi2de(sr,'left-msb')
    return (seq,fill,Tm)

# State space generator formulation of SSRG
def ssgs(num,poly,ifill):
    poly = np.asarray(poly)
    degree = poly[0]
    seq = np.empty(num,dtype=np.uint64)
    taps = np.zeros(degree + 1,dtype=np.uint64)
    taps[degree - poly] = int(1)
    sr = de2bi(ifill,degree,'left-msb')
    # form ssrg characteristic matrix Ts
    ps = taps[-1:0:-1]
    Ts = np.hstack((np.zeros((degree-1,1),dtype=np.uint64),np.eye(degree-1,dtype=np.uint64)))
    Ts = np.flipud(np.fliplr(np.vstack((Ts,ps))))

    for nn in xrange(num):
        seq[nn] = sr[-1]
        sr = mvMultGF2(Ts,sr)

    # final fill
    fill = bi2de(sr,'left-msb')
    return (seq,fill,Ts)

# Convert MSRG polynomial and initial fill value to SSRG.
def msrg2ssrg(mpoly,mfill):
    mpoly = np.asarray(mpoly)
    degree = mpoly[0]
    spoly= (degree - mpoly)[::-1]
    (seq,fill) = msrg(degree,mpoly,mfill)
    seq = seq[::-1]
    sfill = bi2de(seq,'left-msb')
    return (spoly,sfill)

# Convert SSRG polynomial and initial fill value to MSRG.
def ssrg2msrg(spoly,sfill):
    spoly = np.asarray(spoly)
    degree = spoly[0]
    staps = np.zeros(degree + 1,dtype=np.uint64)
    staps[degree - spoly] = int(1)

    mpoly= (degree - spoly)[::-1]
    mtaps = np.zeros(degree + 1,dtype=np.uint64)
    mtaps[degree - mpoly] = int(1)

    ssr = de2bi(sfill,degree,'left-msb')
    msr = np.zeros(degree,dtype=np.uint64)
    for nn in range(degree-1,-1,-1):
        if ssr[nn]:
            msr = np.hstack((np.uint64(1),xor(msr[0:-1],mtaps[1:-1]))) 
        else:
            msr = np.hstack((np.uint64(0),msr[0:-1]))
    
    pm = mtaps[-2::-1]
    Tm = np.vstack((np.eye(degree-1,dtype=np.uint64),np.zeros((1,degree-1),dtype=np.uint64)))
    Tm = np.flipud(np.fliplr(np.hstack((pm.reshape(degree,1),Tm))))
    invTm = matInvGF2(Tm)
    for nn in range(0,degree,1):
        msr = mvMultGF2(invTm,msr)
    mfill = bi2de(msr,'left-msb')
    return (mpoly,mfill)

# Given a jump value and tap polynomial, compute
# the code phase jump mask for an MSRG
def msrg_jump2mask(jump=6,poly=[4,3,0]):
    degree = poly[0]
    ifill = 2**(degree - 1)
    (fill,Tm) = ssgm_jump(-jump, poly, ifill)
    maskb, mfill = msrg(degree, poly, fill)
    mask = bi2de(maskb,'left-msb')
    return mask

# Compute mask corresponding to a code phase change 
# of jump states. jump > 0 delays the code phase, 
# while jump < 0 advances it. The resulting mask
# is for an SSRG (aka Fibonacci) generator.
# This implementation uses a fast O(log(N)) 
# algorithm where N is the jump value.
def jump2mask(jump=6,poly=[4,3,0]):
    (mask,Tm) = ssgm_jump(jump,poly,1)
    return mask

# Compute mask corresponding to a code phase change
# of shift states. Only shift > 0 are accepted.
# The resulting mask is for an SSRG generator.
# This is a slower and less flexible version of 
# jump2mask and was used to test correctness of the 
# former.  This implementation uses a brute-force
# approach requiring O(N) where N is the shift. 
def shift2mask(shift=6,poly=[4,3,0]):
    (seq,mask,Tm) = ssgm(shift,poly,1)
    return mask

# Compute the MSRG fill value and transition 
# matrix corresponding to a code phase change
# of jump states. Jump can be positive or
# negative corresponding respectively to a
# code phase delay or advance.
def ssgm_jump(jump=6,poly=[4,3,0],ifill=8):
    poly = np.asarray(poly)
    degree = poly[0]
    taps = np.zeros(degree + 1,dtype=np.uint64)
    taps[degree - poly] = int(1)
    sr = de2bi(ifill,degree,'left-msb')
    # form msrg characteristic matrix T
    p = taps[-2::-1]
    p = p.reshape(p.size,1)
    T = np.vstack((np.eye(degree-1,dtype=np.uint64),np.zeros((1,degree-1),dtype=np.uint64)))
    T = np.flipud(np.fliplr(np.hstack((p,T))))
    if jump < 0:
        T = matInvGF2(T)
        jump = -jump;
    sr = jumpSR(jump, T, sr)

    # final fill
    fill = bi2de(sr,'left-msb')
    return (fill,T)

# Compute the SSRG fill value and transition 
# matrix corresponding to a code phase change
# of jump states. Jump can be positive or
# negative corresponding respectively to a
# code phase delay or advance.
def ssgs_jump(jump=6,poly=[4,3,0],ifill=8):
    poly = np.asarray(poly)
    degree = poly[0]
    taps = np.zeros(degree + 1,dtype=np.uint64)
    taps[degree - poly] = int(1)
    sr = de2bi(ifill,degree,'left-msb')
    # form ssrg characteristic matrix Ts
    ps = taps[-1:0:-1]
    Ts = np.hstack((np.zeros((degree-1,1),dtype=np.uint64),np.eye(degree-1,dtype=np.uint64)))
    Ts = np.flipud(np.fliplr(np.vstack((Ts,ps))))
    if jump < 0:
        Ts = matInvGF2(Ts)
        jump = -jump;
    sr = jumpSR(jump, Ts, sr)

    # final fill
    fill = bi2de(sr,'left-msb')
    return (fill,Ts)

# Given a masked SSRG and its initial fill compute
# the equivalent fill for an maskless SSRG
def ssrgmask2ssrg(poly=[4,3,0],ifill=1,mask=15):
  (seq, mfill) = ssrg_mask(poly[0],poly,ifill,mask)
  sfill = bi2de(seq,'right-msb')
  return sfill


##################
# Helper functions
##################

# Binary vector to decimal scalar converter, good to 64-bits.
def bi2de(bin,flag = 'left-msb'):
    if flag.lower() == 'right-msb':
        bin = bin[ : : -1]
    return bin.dot(1 << np.arange(bin.shape[-1] - 1, -1, -1,dtype=np.uint64))

# Decimal scalar to binary vector converter, good to 64-bits.
def de2bi(din,nb,flag = 'left-msb'):
    bstr = np.binary_repr(din,nb)
    if flag.lower() == 'right-msb':
        bstr = bstr[ : : -1]
    return np.array([np.uint64(i) for i in bstr],dtype=np.uint64)

# Given an initial SSRG fill value (sr), one-step 
# transition matrix (TC), and a jump value (jump),
# compute the final SSRG fill value corresponding
# to a code phase change of jump states.  Jump 
# can be positive or negative, corresponding 
# respectively to a code phase delay or advance.
def jumpSR(jump, TC, sr):
    if jump == 0:
        return sr
    N = np.uint64(np.ceil(np.log2(jump + 0.5)))
    jumpb = de2bi(jump,N,'right-msb')
    TM = np.uint64(np.zeros((N,TC.shape[0],TC.shape[1])))
    TM[0] = TC
    for pp in range(1,N,1):
        T = TM[pp-1]
        TM[pp] = matSqrGF2(T)
    #
    for pp in range(N-np.uint64(1),-1,-1):
        if jumpb[pp] == 1:
            T = TM[pp]
            sr = mvMultGF2(T,sr)
        #
    #
    return sr

# Square of a matrix in GF2
def matSqrGF2(M):
    return np.uint64(np.mod(np.int64(np.dot(M,M)),2))

# Matrix-vector product in GF2
def mvMultGF2(M,v):
    return np.uint64(np.mod(np.int64(np.dot(M,v)),2))

# Invert matrix in GF2
def matInvGF2(M):
    return np.uint64(np.mod(np.int64(np.linalg.inv(M)),2))

# Elementwise XOR
def xor(a,b):
    return np.uint64(a != b)
