function mask=lfsr_shift2mask(shift,poly)
% Usage: mask=lfsr_shift2mask(shift,poly)
%
% Given a shift value and polynomial compute a shift mask for an SSRG.
%

[seq,mask,Tm] = lfsr_ssgm(shift,poly,1);

