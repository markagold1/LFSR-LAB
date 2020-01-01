# Test lfsr.jump2mask()

import matplotlib.pyplot as plt
import numpy as np
import lfsr as l

# Expected results:
expected_fill = 15495097962036131017
expected_mask = 12263066945517737910
expected_fillj = 13932952120993008214

# Degree 64 generator polynomial:
poly = [64,63,61,60,0]

# Generate the first 32k states from initial fill=1:
(seq,fill)=l.ssrg(2**15,poly,1)
if fill == expected_fill:
  print "FILL = ", fill, "OK"
else:
  print fill, "FAIL, expected", expected_fill

# jump by almost 1 million code states
jump = 990005

# Compute a mask which, when applied to an SSRG,
# will delay the code sequence by our jump value.
mask = l.jump2mask(jump,poly)
if mask == expected_mask:
  print "MASK = ", mask, "OK"
else:
  print mask, "FAIL, expected", expected_mask

# Now generate a length 1M sequence using
# a masked generator with the same initial
# fill=1 and our mask. This sequence is the
# original sequence delayed by jump code
# states. Therefore the original sequence can
# be found within this new sequence starting
# at index=jump=990005.
(seqj,fillj)=l.ssrg_mask(2**20,poly,1,mask)
if fillj == expected_fillj:
  print "FILL after jump = ", fillj, "OK"
else:
  print fillj, "FAIL, expected", expected_fillj

# Plot the difference between the original
# sequence and the delayed sequence evaluated
# starting at index=jump. The two should be
# equal and therefore the difference of 0 is
# reflected in the plot.
plt.plot(seq[1:2**15]-seqj[jump+1:jump+2**15])
plt.show()



# fill is 15495097962036131017
# mask is 12263066945517737910
# fillj is 13932952120993008214
