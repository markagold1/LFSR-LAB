function b = de2bi(d, n, p, f)
% Usage: b = de2bi(d, n, p, f)
%
% Convert decimal numbers to binary vectors.
%
% d..........decimal number to convert
% n..........number of binary digits to create
% p..........optional base to convert to (default=2)
% f..........flag defining the target orientation
%            one of : 'left-msb', 'right-msb'
%            default is 'right'msb'
%

  if (nargin == 1)
    p = 2;
    n = floor ( log (max (max (d), 1)) ./ log (p) ) + 1;
    f = "right-msb";
  elseif (nargin == 2)
    p = 2;
    f = "right-msb";
  elseif (nargin == 3)
    if (ischar (p))
      f = p;
      p = 2;
    else
      f = "right-msb";
    end
  elseif (nargin == 4)
    if (ischar (p))
      tmp = f;
      f = p;
      p = tmp;
    end
  else
    help de2bi
    return
  end

  d = d(:);
  if (~ (all (d == fix (d)) && all (d >= 0)))
    error ("de2bi: all elements of D must be non-negative integers");
  end

  if (isempty (n))
    n = floor ( log (max (max (d), 1)) ./ log (p) ) + 1;
  end

  power = ones (length (d), 1) * (p .^ [0 : n-1] );
  d = d * ones (1, n);
  b = floor (rem (d, p*power) ./ power);

  if (strcmp (f, "left-msb"))
    b = b(:,columns(b):-1:1);
  elseif (~strcmp (f, "right-msb"))
    error ("de2bi: invalid option '%s'", f);
  end

function y = columns(x)
  y = size(x,2);
