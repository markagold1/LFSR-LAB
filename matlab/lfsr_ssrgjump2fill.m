function fill=lfsr_ssrgjump2fill(jump,poly,ifill)
% Usage: fill=lfsr_ssrgjump2fill(jump,poly,ifill)
%
% +-----------------------------------------------------+
% | THIS FILE IS DEPRECATED. USE LFSR_SSGS_JUMP INSTEAD |
% +-----------------------------------------------------+
%
% Given a jump value, polynomial, and initial fill 
% compute the fill value corresponding to a code
% phase shift of jump states
% 
% JUMP...number of code states to jump. JUMP > 0
%        delays the code sequence, jump < 0 advances
% POLY...SSRG generator polynomial; numeric vector 
%        containing the exponents of z for the 
%        nonzero terms of the polynomial in 
%        descending order of powers of z
% IFILL..initial fill before JUMP
% FILL...fill value after JUMP
%
% All binary vectors use 'left-msb' orientation
%

warning('lfsr_ssrgjump2fill is deprecated. Use lfsr_ssgs_jump instead.');
jmask = lfsr_jump2mask(jump, poly);
fill = lfsr_ssrgmask2ssrg(poly, ifill, jmask);
