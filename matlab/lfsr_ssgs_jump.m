function [fill,Ts] = lfsr_ssgs_jump(jump,poly,ifill)
% Usage: [fill,Ts] = lfsr_ssgs_jump(jump,poly,ifill)
%
% Compute the SSRG state space generator transition matrix and 
% final SSRG fill state given a jump value and initial
% SSRG fill state.
%
% JUMP...number of code states to jump. JUMP > 0
%        delays the code sequence, jump < 0 advances
% POLY...SSRG generator polynomial; numeric vector 
%        containing the exponents of z for the 
%        nonzero terms of the polynomial in 
%        descending order of powers of z
% IFILL..scalar, initial shift register state
%
% FILL...Final SSRG shift register fill state
%        representing JUMP states from the initial
%        fill state
% TS.....State space generator transition matrix
%        which applies JUMP state transition 
%        each application. TM = TC^JUMP where
%        TC is the characteristic transition matrix
%
%
%   +--------(+)<-----(+)-----(+)<-------+
%   | r       ^ r-1    ^ r-2   ^         |
%   |z        |z       |z      |z        |1
%   +-->|r-1|-+->|r-2|-+- ... -+->| 0 |--+------> seq
%
% Example:
%
%  [fill,Ts]=lfsr_ssgs_jump(15,[5,3,0],1);
%
% All binary vectors use 'left-msb' orientation
%


degree = poly(1);
taps(1+degree-poly) = 1;
sr = de2bi(ifill,degree,'left-msb').';
% form ssrg characteristic matrix
p = fliplr(taps(2:end));
T = [zeros(degree-1,1) eye(degree-1)];
T = flipud(fliplr([T;p]));

if jump < 0
  T = mod(T^-1,2);
  jump = -jump;
end

% Compute jumped fill state
sr = jumpSR(jump, T, sr);

% final fill
fill = b2d(sr.','left-msb');
Ts = T;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sr = jumpSR(jump, TC, sr)

if jump == 0
  return
end
N = ceil(log2(jump + 0.5));
jump2 = de2bi(jump,N,2,'right-msb');
TM = zeros(N,size(TC,1),size(TC,2));
TM(1,:,:) = TC;
for kk = 2:N
  T = squeeze(TM(kk-1,:,:));
  TM(kk,:,:) = mod(T^2,2);
end

for kk = N:-1:1
  if jump2(kk) == 1
    T = squeeze(TM(kk,:,:));
    sr = mod(T*sr,2);
  end
end
  
function d = b2d(b,flag)

N = length(b);
d = uint64(0);
for kk=N:-1:1
  d = d + bitshift(uint64(b(kk)),N-kk);
end
d = double(d);

