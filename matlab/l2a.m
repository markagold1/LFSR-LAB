function a = l2a(l, invert)
% Usage: a = l2a(l, invert)
%
% Logical to analytic conversion.
%
% L2A converts from logical bits to antipodal values.
% Conversion rule is:
%
%    (default)
%    invert=0      invert=1
%    --------      --------
%    0 -> -1       0 -> +1
%    1 -> +1       1 -> -1
%

if nargin == 1 || invert == 0
  invert = 0;
elseif invert
  invert = 1;
end

a = round(2 * (double(l) - 0.5));
if invert
  a = -1*a;
end
