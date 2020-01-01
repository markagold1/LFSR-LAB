function Sdln = psc(code_group, primary_code)
% Usage: Sdln = psc(code_group, primary_code)
%
% WCDMA primary scrambling code generator
%
% Reference: TS 25.213 section 5.2.2:

% code_group range: 0..63, 
% primary_code range: 0..7
code_index = 16 * 8 * code_group + 16 * primary_code;
assert(code_group < 64 && code_group >= 0, 'Scrambling code group must be in the range [0,64)');
assert(primary_code < 8 && primary_code >= 0, 'Primary code in code_group must be in the range [0,8)');

% Code generators
xpoly = [18, 7, 0];        % gx = x^18 + x^7 + 1
ypoly = [18, 10, 7, 5, 0]; % gy = y^18 + y^10 + y^7 + y^5 + 1

% Initial states
xfill = 1;
yfill = 2^18 - 1;

% Construct nth Gold code sequence for I channel
xIjump =  -code_index;
xIFillShift = lfsr_ssgs_jump(-xIjump, xpoly, xfill);
xIseq = lfsr_ssrg(38400, xpoly, xIFillShift);
yIseq = lfsr_ssrg(38400, ypoly, yfill);
znI = mod(xIseq + yIseq, 2);

% Construct nth Gold code sequence for Q channel
xQjump = xIjump - 131072;
xQFillShift = lfsr_ssgs_jump(-xQjump, xpoly, xfill);
xQseq = lfsr_ssrg(38400, xpoly, xQFillShift);
yQjump = -131072;
yQFillShift = lfsr_ssgs_jump(-yQjump, ypoly, yfill);
yQseq = lfsr_ssrg(38400, ypoly, yQFillShift);
znQ = mod(xQseq + yQseq, 2);

% Logical to analytic (0->+1, 1->-1)
% nth complex scrambling code sequence Sdln
Sdln = l2a(not(znI)) + j*l2a(not(znQ));
