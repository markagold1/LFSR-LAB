function [seq,fill] = lfsr_msrg_mask(num,poly,ifill,mask)
% Usage: [seq,fill] = lfsr_msrg_mask(num,poly,ifill,mask)
%
% poly...numeric vector containing the exponents of z 
%        for the nonzero terms of the polynomial in 
%        descending order of powers
% ifill..scalar, initial shift register state
% mask...scalar, mask applied to shift register
%
%   +<----------+-----------+----------+----------+
%   |  r        |  r-1      |  r-2     |          ^
%   | z         v z         v z        v z        | 1
%   +-->|r-1|->(+)->|r-2|->(+)- ...-->(+)->| 0 |--+---> seq
%            |           |               |       |
%            v           v               v       v
%          +----------------------------------------+
% mask --->|                  AND                   |
%          +--------------++------------------------+
%                             ||
%                             \/
%              +--------------------------------+
%              |             XOR                |
%              +--------------+-----------------+
%                             |
%                             +------> seq
%
% Example:
%
%  [seq,fill]=lfsr_msrg_mask(32,[5,3,0],1,3);
%
% All binary vectors use 'left-msb' orientation
%

seq = NaN(1,num);
degree = poly(1);
taps(1+degree-poly) = 1;
sr = de2bi(ifill,degree,'left-msb');
ma = de2bi(mask,degree,'left-msb');
for nn = 1:num
  seq(nn) = mod(sum(and(sr,ma)),2);
  if sr(end)
    sr = [1 xor(sr(1:end-1),taps(2:end-1))];
  else
    sr = [0 sr(1:end-1)];
  end;
end;

% final fill
fill = bi2de(sr,'left-msb');
