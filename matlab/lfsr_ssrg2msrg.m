function [mpoly,mfill]=lfsr_ssrg2msrg(spoly,sfill)

% 1. Compute mpoly directly from spoly
% 2. Open the MSRG feedback output and clock
%    sfill bit by bit into the MSRG feedback input
% 3. Use inv(Tm) to propagate in reverse by
%    length(sfill) steps

% 0. Initialize degree and staps
degree = spoly(1);
staps(1+degree-spoly) = 1;

% 1. Compute mpoly directly from spoly
mpoly = fliplr(degree - spoly);
mtaps(1+degree-mpoly) = 1;

% 2. Open the MSRG feedback output and clock
%    sfill bit by bit into the MSRG feedback input
ssr = de2bi(sfill,degree,'left-msb');
msr = zeros(1,degree);
for nn = degree:-1:1
  if ssr(nn)
    msr = [1 xor(msr(1:end-1),mtaps(2:end-1))];
  else
    msr = [0 msr(1:end-1)];
  end;
end;

% 3. Use inv(Tm) to propagate in reverse by
%    length(sfill) steps
% form msrg characteristic matrix
pm = fliplr(mtaps(1:end-1));
Tm = [eye(degree-1);zeros(1,degree-1)];
Tm = flipud(fliplr([pm(:) Tm]));
invTm = mod(round(det(Tm)*inv(Tm)),2);
for nn = 1:degree
  msr = mod(invTm*msr(:),2);
end;

mfill = bi2de(msr.','left-msb');

