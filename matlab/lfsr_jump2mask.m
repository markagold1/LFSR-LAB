function mask=lfsr_jump2mask(jump,poly)
% Usage: mask=lfsr_jump2mask(jump,poly)
%
% Given a jump value and polynomial compute a phase jump mask for an SSRG.
%
% JUMP...number of code states to jump. JUMP > 0
%        delays the code sequence, jump < 0 advances
% POLY...SSRG generator polynomial; numeric vector 
%        containing the exponents of z for the 
%        nonzero terms of the polynomial in 
%        descending order of powers of z
% MASK...Phase shift mask corresponding to a code
%        sequence shift of JUMP states
%

[mask,Tm] = lfsr_ssgm_jump(jump,poly,1);


