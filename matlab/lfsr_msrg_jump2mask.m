function mask=lfsr_msrg_jump2mask(jump, poly)
% Usage: mask=lfsr_marg_jump2mask(jump, poly)
%
% Given a jump value and polynomial compute the code
% phase jump mask for an MSRG.
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
% Reference: Lee & Miller, "CDMA Systems Engineering Handbook"
%            Artech House, 1998, section 6.3.3
%

degree = poly(1);
ifill = 2^(degree - 1);
fill = lfsr_ssgm_jump(-jump, poly, ifill);
maskb = lfsr_msrg(poly(1), poly, fill);
mask = bi2de(maskb, 'left-msb');
