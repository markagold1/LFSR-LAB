function l = a2l(a, invert)
% Usage: l = a2l(a, invert)
%
% Analytic to logical conversion.
%
% A2L converts from antipodal values to logical bits.
% Conversion rule is:
%
%    (default)
%    invert=0        invert=1
%   ----------      ---------
%    >0 ->  1        >0 ->  0
%   <=0 ->  0       <=0 ->  1
%

if nargin == 1 || invert == 0
  invert = 0;
elseif invert
  invert = 1;
end

l = logical(0.5 * (sign(a) + 1));
if invert
  l = ~l;
end

