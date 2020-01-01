function d = bi2de(b, p, f)
% Usage: b = bi2de(b, p, f)
%
% Convert binary vectors to decimal numbers.
%
% b..........binary vector to convert
% p..........optional base of vector to convert from (default=2)
% f..........flag defining the target orientation
%            one of : 'left-msb', 'right-msb'
%            default is 'right'msb'
%

  switch (nargin)
    case 1
      p = 2;
      f = "right-msb";
    case 2
      if (ischar (p))
        f = p;
        p = 2;
      else
        f = "right-msb";
      end
    case 3
      if (ischar (p))
        tmp = f;
        f = p;
        p = tmp;
      end
    otherwise
      help bi2de
      return
  end

  if (~ (all (b(:) == fix (b(:))) && all (b(:) >= 0) && all (b(:) < p)))
    error ("bi2de: all elements of B must be integers in the range [0,P-1]");
  end

  if (strcmp (f, "left-msb"))
    b = b(:,size (b, 2):-1:1);
  elseif (~strcmp (f, "right-msb"))
    error ("bi2de: invalid option '%s'", f);
  end

  if (length (b) == 0)
    d = [];
  else
    d = b * (p .^ [0:(columns (b) - 1)]');
  end

function y = columns(x)
  y = size(x,2);
