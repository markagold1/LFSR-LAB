function fill = lfsr_ssrgmask2ssrg(poly,ifill,mask)
% Usage: fill = lfsr_ssrgmask2ssrg(poly,ifill,mask)
%
% Given a masked SSRG and its initial fill compute
% the equivalent fill for a maskless SSRG.
%
% poly...numeric vector containing the exponents of z 
%        for the nonzero terms of the polynomial in 
%        descending order of powers
% ifill..scalar, initial shift register state
% mask...scalar, mask applied to shift register to
%        effect a code phase shift
% fill...equivalent fill for maskless SSRG
%
% Masked SSRG:
%
%   +--------(+)<-----(+)-----(+)<-------+
%   | r       ^ r-1    ^ r-2   ^         |
%   |z        |z       |z      |z        |1
%   +-->|r-1|-+->|r-2|-+- ... -+->| 0 |--+------> seq
%             |        |       |         |
%             v        v       v         v
%          +--------------------------------+
% mask --->|             AND                |
%          +--------------++----------------+
%                         ||
%                         \/
%          +--------------------------------+
%          |             XOR                |
%          +--------------+-----------------+
%                         |
%                         +------> seq
%
% Maskless SSRG:
%
%   +--------(+)<-----(+)-----(+)<-------+
%   | r       ^ r-1    ^ r-2   ^         |
%   |z        |z       |z      |z        |1
%   +-->|r-1|-+->|r-2|-+- ... -+->| 0 |--+------> seq
%
% Example:
%
%  fill=lfsr_ssrgmask2ssrg([5,3,0],1,14);
%
% All binary vectors use 'left-msb' orientation
%

num = poly(1); % number of code bits to generate = degree
[seq,fill] = lfsr_ssrg_mask(num,poly,ifill,mask);
fill = bi2de(seq, 'right-msb');
