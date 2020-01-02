/* Harness for testing functions in lfsr.h */

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <cstdlib>
#include "lfsr.h"

void print_gf2_vec(std::vector<int> &v, int wrap=16);

int main()
{
        int num = 34;
        uint64_t taps;
        uint64_t ifill;
        uint64_t staps;
        uint64_t sfill;
        uint64_t mtaps;
        uint64_t mfill;
        uint64_t mask;
        lfsr_info info;
        int shift;
        std::vector<int> seq;

        // lfsr_msrg()
        {
          ifill = 1;
          mfill = ifill;
          mtaps = 25;
          lfsr msrg(mtaps, mfill, "msrg");
          seq = msrg(num);
          std::cout << "lfsr_msrg():" << std::endl;
          std::cout << "MSRG taps: " << msrg.get_taps() << std::endl;
          std::cout << "MSRG initial fill: " << mfill << std::endl;
          std::cout << "MSRG final fill: " << msrg.get_fill() << std::endl;
          std::cout << "MSRG seq:" << std::endl;
          print_gf2_vec(seq, 15);
          std::cout << std::endl;
        }

        // lfsr_msrg2ssrg()
        {
          ifill = 1;
          mfill = ifill;
          mtaps = 25;
          lfsr msrg(mtaps, mfill, "msrg");
          info = msrg.convert();
          std::cout << "lfsr_msrg2ssrg():" << std::endl;
          std::cout << "MSRG taps: " << msrg.get_taps() << std::endl;
          std::cout << "MSRG initial fill: " << msrg.get_fill() << std::endl;
          std::cout << "SSRG taps: " << info.taps << std::endl;
          std::cout << "SSRG initial fill: " << info.fill << std::endl;
          lfsr ssrg(info.taps, info.fill, "ssrg");
          seq = ssrg(num);
          std::cout << "SSRG seq:" << std::endl;
          print_gf2_vec(seq, 15);
          std::cout << std::endl;
        }

        // lfsr_ssrg()
        {
          ifill = 9;
          sfill = ifill;
          staps = 19;
          lfsr ssrg(staps, sfill, "ssrg");
          std::cout << "lfsr_ssrg():" << std::endl;
          std::cout << "SSRG taps: " << ssrg.get_taps() << std::endl;
          std::cout << "SSRG initial fill: " << ssrg.get_fill() << std::endl;
          seq = ssrg(num);
          std::cout << "SSRG final fill: " << ssrg.get_fill() << std::endl;
          std::cout << "SSRG seq:" << std::endl;
          print_gf2_vec(seq, 15);
          std::cout << std::endl;
        }

        // lfsr_ssrg2msrg()
        {
          ifill = 9;
          sfill = ifill;
          staps = 19;
          lfsr ssrg(staps, sfill, "ssrg");
          info = ssrg.convert();
          std::cout << "lfsr_ssrg2msrg():" << std::endl;
          std::cout << "SSRG taps: " << ssrg.get_taps() << std::endl;
          std::cout << "SSRG initial fill: " << ssrg.get_fill() << std::endl;
          std::cout << "MSRG taps: " << info.taps << std::endl;
          std::cout << "MSRG initial fill: " << info.fill << std::endl;
          lfsr msrg(info.taps, info.fill, "msrg");
          seq = msrg(num);
          std::cout << "MSRG seq:" << std::endl;
          print_gf2_vec(seq, 15);
          std::cout << std::endl;
        }

        // lfsr_ssgm()
        ifill = 1;
        mfill = ifill;
        mtaps = 25;
        seq = lfsr_ssgm(num, mtaps, ifill);
        std::cout << "lfsr_ssgm():" << std::endl;
        std::cout << "SSGM taps: " << mtaps << std::endl;
        std::cout << "SSGM initial fill: " << mfill << std::endl;
        std::cout << "SSGM final fill: " << ifill << std::endl;
        std::cout << "SSGM seq:" << std::endl;
        print_gf2_vec(seq, 15);
        std::cout << std::endl;

        // lfsr_ssgs()
        ifill = 9;
        sfill = ifill;
        staps = 19;
        seq = lfsr_ssgs(num, staps, ifill);
        std::cout << "lfsr_ssgs():" << std::endl;
        std::cout << "SSGS taps: " << staps << std::endl;
        std::cout << "SSGS initial fill: " << sfill << std::endl;
        std::cout << "SSGS final fill: " << ifill << std::endl;
        std::cout << "SSGS seq:" << std::endl;
        print_gf2_vec(seq, 15);
        std::cout << std::endl;

        // lfsr_shift2mask() (ssrg shift to mask)
        {
          shift = 8;
          staps = 19;
          mask = lfsr_shift2mask(shift, staps);
          std::cout << "lfsr_shift2mask():" << std::endl;
          std::cout << "SSRG taps: " << staps << std::endl;
          std::cout << "SSRG shift: " << shift << std::endl;
          std::cout << "SSRG initial fill: 1 (assumed)" << std::endl;
          std::cout << "SSRG mask: " << mask << std::endl;
          std::cout << std::endl;
        }

        // lfsr_ssrg_mask()
        {
          ifill = 1;
          sfill = ifill;
          staps = 19;
          mask = 11;
          lfsr ssrg_mask(staps, sfill, "ssrg", mask);
          std::cout << "lfsr_ssrg_mask():" << std::endl;
          std::cout << "Masked SSRG taps: " << ssrg_mask.get_taps() << std::endl;
          std::cout << "Masked SSRG initial fill: " << ssrg_mask.get_fill() << std::endl;
          seq = ssrg_mask(num);
          std::cout << "Masked SSRG final fill: " << ssrg_mask.get_fill() << std::endl;
          std::cout << "Masked SSRG mask: " << ssrg_mask.get_mask() << std::endl;
          std::cout << "Masked SSRG seq (delayed by 8 relative to unmasked):" << std::endl;
          print_gf2_vec(seq, 15);
          std::cout << "lfsr_ssrg(): (unmasked)" << std::endl;
          lfsr ssrg(staps, sfill, "ssrg");
          std::cout << "Unasked SSRG taps: " << ssrg.get_taps() << std::endl;
          std::cout << "Unasked SSRG initial fill: " << ssrg.get_fill() << std::endl;
          seq = ssrg(num);
          std::cout << "Unasked SSRG final fill: " << ssrg.get_fill() << std::endl;
          std::cout << "Unnasked SSRG seq:" << std::endl;
          print_gf2_vec(seq, 15);
          std::cout << std::endl;
        }

        // lfsr_jump2mask() (ssrg jump to mask)
        {
          shift = 8;
          sfill = 11111;
          staps = 19;
          lfsr ssrg(staps, sfill, "ssrg");
          mask = ssrg.calculate_mask(shift);
          std::cout << "lfsr_jump2mask():" << std::endl;
          std::cout << "SSRG taps: " << ssrg.get_taps() << std::endl;
          std::cout << "SSRG shift: " << shift << std::endl;
          std::cout << "SSRG initial fill: 1 (assumed)" << std::endl;
          std::cout << "SSRG mask: " << mask << std::endl;
          std::cout << std::endl;
        }

        // lfsr_msrg_mask()
        {
          mfill = 1;
          mtaps = 25;
          mask = 11;
          lfsr msrg_mask(mtaps, mfill, "msrg", mask);
          std::cout << "lfsr_msrg_mask():" << std::endl;
          std::cout << "Masked MSRG taps: " << msrg_mask.get_taps() << std::endl;
          std::cout << "Masked MSRG initial fill: " << msrg_mask.get_fill() << std::endl;
          seq = msrg_mask(num);
          std::cout << "Masked MSRG final fill: " << msrg_mask.get_fill() << std::endl;
          std::cout << "Masked MSRG mask: " << msrg_mask.get_mask() << std::endl;
          std::cout << "Masked MSRG seq (delayed by 6 relative to unmasked):" << std::endl;
          print_gf2_vec(seq, 15);
          std::cout << "lfsr_msrg(): (unmasked)" << std::endl;
          lfsr msrg(mtaps, mfill, "msrg");
          std::cout << "Unasked MSRG taps: " << msrg.get_taps() << std::endl;
          std::cout << "Unasked MSRG initial fill: " << msrg.get_fill() << std::endl;
          seq = msrg(num);
          std::cout << "Unasked MSRG final fill: " << msrg.get_fill() << std::endl;
          std::cout << "Unnasked MSRG seq:" << std::endl;
          print_gf2_vec(seq, 15);
          std::cout << std::endl;
        }

        // lfsr_ssgs_jump()
        {
          int jump = 8064;
          taps = 0x40081; // 262273;
          ifill = 1;
          lfsr ssrg(taps, ifill, "ssrg");
          uint64_t fill = ssrg.jump(jump);
          //uint64_t fill = lfsr_ssgs_jump(jump, taps, ifill);
          std::cout << "lfsr_ssgs_jump():" << std::endl;
          std::cout << "SSRG jump: " << jump << std::endl;
          std::cout << "SSRG taps: " << ssrg.get_taps() << std::endl;
          std::cout << "Fill for jump: " << fill << std::endl;
          std::cout << std::endl;
        }

        // get_t()
        {
          std::cout << "get_t():" << std::endl;
          taps = 0x29;
          ifill = 1;
          lfsr ssrg(taps, ifill, "ssrg");
          info = ssrg.convert();
          lfsr msrg(info.taps, info.fill, "msrg");

          std::cout << "Ts: " << std::endl;
          ssrg.get_T().PrintGF2();
          std::cout << std::endl;
          std::cout << "Ts: " << std::endl;
          msrg.get_T().PrintGF2();
          std::cout << std::endl;
        }

        return 0;
}

void print_gf2_vec(std::vector<int> &v, int wrap)
{
        for (size_t ii = 0; ii < v.size(); ++ii) {
        	std::cout << " " << v[ii];
        	if (ii % wrap == (wrap-1)) {
        		std::cout << std::endl;
        	}
        }
        std::cout << std::endl;
}
