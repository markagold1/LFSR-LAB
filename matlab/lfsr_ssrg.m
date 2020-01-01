function [seq,fill] = lfsr_ssrg(num,poly,ifill)
% Usage: [seq,fill] = lfsr_ssrg(num,poly,ifill)
%
% poly...numeric vector containing the exponents of z 
%        for the nonzero terms of the polynomial in 
%        descending order of powers
% ifill..scalar, initial shift register state
%
%   +--------(+)<-----(+)-----(+)<-------+
%   | r       ^ r-1    ^ r-2   ^         |
%   |z        |z       |z      |z        |1
%   +-->|r-1|-+->|r-2|-+- ... -+->| 0 |--+------> seq
%
% Example:
%
%  [seq,fill]=lfsr_ssrg(32,[5,3,0],1);
%
% All binary vectors use 'left-msb' orientation
%

seq = NaN(1,num);
degree = poly(1);
taps(1+degree-poly) = 1;
sr = de2bi(ifill,degree,'left-msb');
for nn = 1:num
  seq(nn) = sr(end);
  parity = mod(sum(and(taps(2:end),sr)),2);
  sr = [parity sr(1:end-1)];
end;

% final fill
fill = bi2de(sr,'left-msb');
