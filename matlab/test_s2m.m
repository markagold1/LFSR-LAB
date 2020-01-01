[mpoly,mfill]=lfsr_ssrg2msrg([6,5,0],1);
[seq,fill]=lfsr_ssrg(64,[6,5,0],1);
[mseq,mfill]=lfsr_msrg(64,mpoly,mfill);
any(mseq~=seq)

[mpoly,mfill]=lfsr_ssrg2msrg([15,14,0],1);
[seq,fill]=lfsr_ssrg(2^15 - 1,[15,14,0],1);
[mseq,mfill]=lfsr_msrg(2^15 - 1,mpoly,mfill);
any(mseq~=seq)


