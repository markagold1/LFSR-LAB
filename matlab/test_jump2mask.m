% Some polynomials from
% http://courses.cse.tamu.edu/walker/csce680/lfsr_table.pdf
% https://www.xilinx.com/support/documentation/application_notes/xapp210.pdf
if 0
poly = [64,63,61,60,0];  % order 64
poly = [32,30,26,25,0];  % order 32
poly = [65,64,62,61,0];  % order 65
end

if 1
%desired_jump = 1294281470000000; % 299d + 14h + 27m + 9.4s @fs
%mask = lfsr_jump2mask(desired_jump,poly);
%mask = 3197494096723161336
%mask = 2c5f c7f0 9e97 b0f8
poly = [64,63,61,60,0];  % order 64 (matlab breaks when order > 63)
poly = [42,40,37,35,0];  % order 42
jump = 990005;
[seq,fill]=lfsr_ssrg(2^15,poly,1);
mask = lfsr_jump2mask(jump,poly);
[seqj,fillj]=lfsr_ssrg_mask(2^20,poly,1,mask);
plot(seq(1:1000)-seqj(jump+(1:1000)));gg
end

if 0
% From C.S0002
polyI=[15,13,9,8,7,5,0];
polyQ=[15,12,11,10,6,5,4,3,0];

% Generate 2 periods of each sequence
[seqI,fill] = lfsr_ssrg(65536,polyI,1);bpseqI=1-2*seqI;
[seqQ,fill] = lfsr_ssrg(65536,polyQ,1);bpseqQ=1-2*seqQ;

% Autocorrelation sequence for each sequence
[xcI,lagI]=xcorr(bpseqI);
[xcQ,lagQ]=xcorr(bpseqQ);

% Get jump mask
jump = 1069;
maskI = lfsr_jump2mask(jump,polyI);
[seqIjump,fill] = lfsr_ssrg_mask(65536,polyI,1,maskI);

bpseqIjump = 1-2*seqIjump;

% Cross correlation
[xcIJ,lagIJ]=xcorr(bpseqI,bpseqIjump);
abs_xc = abs(xcIJ);

% Plots
idx = find(abs_xc == max(abs_xc));
rng = 5;
plot(lagIJ(idx-rng:idx+rng), abs_xc(idx-rng:idx+rng));gg
%plot(lagIJ,abs(xcIJ));gg
%plot(lagI,abs(xcI));gg
%plot(lagQ,abs(xcQ));gg

end


