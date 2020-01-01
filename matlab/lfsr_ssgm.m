function [seq,fill,Tm] = lfsr_ssgm(num,poly,ifill)
% Usage: [seq,fill,Tm] = lfsr_ssgm(num,poly,ifill)
%
% State space generator formulation of MSRG.
%
% poly...numeric vector containing the exponents of z 
%        for the nonzero terms of the polynomial in 
%        descending order of powers
% ifill..scalar, initial shift register state
%
%
%   +<----------+-----------+----------+----------+
%   |  r        |  r-1      |  r-2     |          ^
%   | z         v z         v z        v z        | 1
%   +-->|r-1|->(+)->|r-2|->(+)- ...-->(+)->| 0 |--+---> seq
%
% Example:
%
%  [seq,Tm]=lfsr_ssgm(32,[5,3,0],1);
%
% All binary vectors use 'left-msb' orientation
%

seq = NaN(1,num);
degree = poly(1);
taps(1+degree-poly) = 1;
sr = de2bi(ifill,degree,'left-msb').';
% form msrg characteristic matrix
p = fliplr(taps(1:end-1));
T = [eye(degree-1);zeros(1,degree-1)];
T = flipud(fliplr([p(:) T]));

for nn = 1:num
  seq(nn) = sr(end);
  sr = mod(T*sr,2);
end;

% final fill
fill = bi2de(sr.','left-msb');
Tm = T;
