/*
 * LFSR demonstration: Generate W-CDMA primary scrambling code
 */
#include <iostream>
#include <cstdlib>
#include <vector>
#include "lfsr.h"

int main(int argc, char *argv[])
{
  int code_group = 63;
  int primary_code = 0;
  int seq_length = 38400;

  switch (argc) {
    case 3:
      code_group = atoi(argv[1]);
      primary_code = atoi(argv[2]);
      seq_length = 38400;
      break;
    case 4:
      code_group = atoi(argv[1]);
      primary_code = atoi(argv[2]);
      seq_length = atoi(argv[3]);
      break;
    default:
      std::cerr << "Usage: psc <code_group 0..63>"
                << " <primary_code 0..7>"
                << " [<seq_length 1..2^18-1>]"
                << std::endl;
      std::cerr << std::endl
                << "\t* If unsepecified seq_length defaults to 38400."
                << std::endl;
      std::cerr << "\t* psc outputs its results to stdout"
                << std::endl;
      std::cerr << "\t* psc outputs two columns: I(:) Q(:)"
                << std::endl;
      std::cerr << "\t* row format: '0' -> +1,  '1' -> -1"
                << std::endl;
      exit(EXIT_FAILURE);
  }

#if 0
  std::cout << "argc: " << argc << std::endl;
  std::cout << "code_group: " << code_group << std::endl;
  std::cout << "primary_code: " << primary_code << std::endl;
  std::cout << "seq_length: " << seq_length << std::endl;
#endif

  int code_index = 16*8*code_group + 16*primary_code;
  uint64_t xtaps = (1<<18) + (1<<7) + 1;
  uint64_t ytaps = (1<<18) + (1<<10) + (1<<7) + (1<<5) + 1;
  uint64_t xfill = 1;
  uint64_t yfill = (1<<18) - 1;

  // X and Y m-sequence generators
  lfsr ssrgX(xtaps, xfill, "ssrg");
  lfsr ssrgY(ytaps, yfill, "ssrg");

  // Construct nth Gold code sequence for I channel
  int xIjump = -code_index;
  uint64_t  xIFillShift = ssrgX.jump(-xIjump);
  std::vector<int> xIseq = ssrgX(seq_length);
  std::vector<int> yIseq = ssrgY(seq_length);
  std::vector<int> znI;
  // map "0" to +1 and "1" to -1
  for (size_t ii = 0; ii < xIseq.size(); ++ii) {
    znI.push_back((xIseq[ii] + yIseq[ii]) % 2 ? -1 : 1);
  }

  // Propagate m-sequence generators back to their initial state
  ssrgX.jump(-seq_length + xIjump);
  ssrgY.jump(-seq_length);

  // Construct nth Gold code sequence for Q channel
  int xQjump = xIjump - 131072;
  uint64_t xQFillShift = ssrgX.jump(-xQjump);
  std::vector<int> xQseq = ssrgX(seq_length);
  int yQjump = -131072;
  uint64_t yQFillShift = ssrgY.jump(-yQjump);
  std::vector<int> yQseq = ssrgY(seq_length);
  std::vector<int> znQ;
  // map "0" to +1 and "1" to -1
  for (size_t ii = 0; ii < xQseq.size(); ++ii) {
    znQ.push_back((xQseq[ii] + yQseq[ii]) % 2 ? -1 : 1);
  }

  // Write psc to stdout
  for (int ii = 0; ii < znI.size(); ++ii) {
    printf("%2d %2d\n", znI[ii], znQ[ii]);
  }

  return 0;
}
