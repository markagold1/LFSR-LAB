function [seq,fill] = lfsr_ssrg_mask(num,poly,ifill,mask)
% Usage: [seq,fill] = lfsr_ssrg_mask(num,poly,ifill,mask)
%
% poly...numeric vector containing the exponents of z 
%        for the nonzero terms of the polynomial in 
%        descending order of powers
% ifill..scalar, initial shift register state
% mask...scalar, mask applied to shift register
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
% Example:
%
%  [seq,fill]=lfsr_ssrg_mask(32,[5,3,0],1,3);
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
  parity = mod(sum(and(taps(2:end),sr)),2);
  sr = [parity sr(1:end-1)];
end;

% final fill
fill = bi2de(sr,'left-msb');

% TODO: If de2bi overflows use the following:
function b = d2b(d,N,flag)
% flag = 'left-msb'
d = uint64(d);
ds = sprintf('%.16x\n',d);
dsh = ds(1:8);
dsl = ds(9:16);
bh = de2bi(hex2dec(dsh),32,2,'left-msb');
bl = de2bi(hex2dec(dsl),32,2,'left-msb');
b = [bh bl];
