function P = oct2poly(x, flip, format)
% Usage: P = oct2poly(x, flip, format)
%
% Convert octal value to polynomial vector suitable
% for processing by LFSR
%
%   x.........scalar octal integer
%   flip......optional, non-zero returns the reciprical polynomial (default=0)
%   format....optional, non-zero returns a binary vector of coefficents
%   P.........vector of polynomial exponents of terms
%             with non-zero coefficients
%
% Example: 
%  Polynomial:  X^10 + X^9 + X^8 + X^6 + X^4 + X^2 + 1
%  Octal input: 3525
%  Output (default format): [10 9 8 6 4 2 0]
%  Output  (binary format): [ 1 1 1 0 1 0 1 0 1 0 1]
%

if nargin == 1 || nargin == 2 || format == 0
  % vector of exponents by default
  format = 0;
else
  % binary vector of gf(2) coefficients
  format = 1;
end

if nargin == 1 || flip == 0
  flip = 0;
else
  % return reciprical polynomial
  flip = 1;
end

if x == 0
  P = x;
  return
end

% ensure input is integer
x = floor(x);
xs = num2str(x);

% nn=number of octal digits
nn = length(xs);
Pv = nan(1,3*nn);
for kk = 0:nn-1
  xi = xs(kk+1);
  Pv(3*kk + (1:3)) = lut(xi);
end
Pv = Pv(end:-1:1);
P = fliplr(find(Pv) - 1);
if flip
  P = fliplr(max(P) - P);
end
if format
  P = fliplr(Pv);
  idx = find(P);
  P = P(idx(1):end);
end
%P = fliplr(find(Pv(end:-1:1) ~= 0) -1);


function y = lut(xc)
  tab = [ 0 0 0 
          0 0 1
          0 1 0
          0 1 1
          1 0 0
          1 0 1
          1 1 0
          1 1 1];

y = tab(str2num(xc)+1, :);
% return
