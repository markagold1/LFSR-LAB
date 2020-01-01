function [spoly,sfill] = lfsr_msrg2ssrg(mpoly,mfill)

degree = mpoly(1);
spoly = fliplr(degree - mpoly);
seq = fliplr(lfsr_msrg(degree,mpoly,mfill));
sfill = bi2de(seq,'left-msb');
