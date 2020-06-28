/* LFSR class and library
 *
 * LFSR generator utilities.
 * Simple shift register (SSRG) aka Fibonacci generator model:
 * 
 *    +--------(+)<-----(+)-----(+)<-------+
 *    | r       ^ r-1    ^ r-2   ^         |
 *    |z        |z       |z      |z        |1
 *    +-->|r-1|-+->|r-2|-+- ... -+->| 0 |--+------> seq
 * 
 * Modular shift register (MSRG) aka Galois generator model:
 * 
 *    +<----------+-----------+----------+----------+
 *    |  r        |  r-1      |  r-2     |          ^
 *    | z         v z         v z        v z        | 1
 *    +-->|r-1|->(+)->|r-2|->(+)- ...-->(+)->| 0 |--+---> seq
 * 
 * Terminology:
 * 
 *   poly - Generator polynomial, a degree r polynomial 
 *          in z whose coefficients represent the non-
 *          zero tap connections to the xor feedback
 *          network.  In this package poly is a vector 
 *          of powers in z with non-zero coefficients.
 *          For example poly = [5,3,0] represents
 *          generator polynomial: g(z) = z^5 + z^3 + 1
 *          Poly is a carry-over from the matlab and python
 *          models. In the c/c++ implementation tap
 *          polynomials are represented exclusively by
 *          left-msb integers (see taps definition).
 *   taps - An unsigned 64-bit integer whose bits
 *          represent the coefficients of the generator
 *          polynomial.  For example a tap value of 41
 *          decimal is equivalent to binary 101001b
 *          which represents poly = [5,3,0] because
 *          bits 5, 3, and 0 are set.
 *   fill - Shift register contents, in this package
 *          fill is represented as a shift-right integer
 *          whose left-msb binary equivalent represents
 *          the bit values in the shift register.  For
 *          example fill=1 corresponds to a "1" in
 *          the rightmost shift register stage while
 *          fill=8 corresponds to a "1" two stages to 
 *          the left of the rightmost stage.
 *   mask - Specifies the shift register stages to
 *          combine in order to effect a code phase change.
 *          In this package mask is an integer whose
 *          left-msb binary equivalent represents the 
 *          shift register stages combined.  For example
 *          if mask=7 the contents of the three rightmost
 *          stages are modulo-2 summed to form the LFSR 
 *          output sequence.
 *   num -  An input parameter to the generator functions
 *          specifying the length of sequence to generate.
 * 
 * Masked SSRG model:
 * 
 *    +--------(+)<-----(+)-----(+)<-------+
 *    | r       ^ r-1    ^ r-2   ^         |
 *    |z        |z       |z      |z        |1
 *    +-->|r-1|-+->|r-2|-+- ... -+->| 0 |--+
 *              |        |       |         |
 *              v        v       v         v
 *  A=AND gate +-+      +-+     +-+       +-+
 *             |A|      |A|     |A|       |A|
 *        r-1/ +-+ r-2/ +-+  1/ +-+    0/ +-+
 *  mask ---+---|----+---|---+---|-----+   |
 *              |        |       |         | 
 *              v        v       v         v
 *           +--------------------------------+
 *           |       Modulo-2 sum             |
 *           +--------------+-----------------+
 *                          |
 *                          +------> seq
 * 
 * Masked MSRG model:
 * 
 *     +<----------+-----------+----------+----------+
 *     |  r        |  r-1      |  r-2     |          ^
 *     | z         v z         v z        v z        |1
 *     +-->|r-1|->(+)->|r-2|->(+)- ...-->(+)->| 0 |--+
 *              |           |          |             |
 *              v           v          v             v
 *  A=AND gate +-+         +-+        +-+           +-+
 *             |A|         |A|        |A|           |A|
 *        r-1/ +-+    r-2/ +-+     1/ +-+        0/ +-+
 *  mask ---+---|-------+---|------+---|---------+   |
 *              |           |          |             | 
 *              v           v          v             v
 *           +------------------------------------------+
 *           |             Modulo-2 sum                 |
 *           +--------------------+---------------------+
 *                                |
 *                                +------> seq
 *
 * State-space generator (SSG) model:
 *  Ref [1] Section 6.2.3 
 *
 * SSG models are useful for fast propagation of LFSR state.
 *
 *  T   : Characteristic matrix, N-by-N, where N is the degree
 *        of the generator (tap) polynomial
 *
 *  v   : Shift register state (fill) at time step n, an N-by-1 vector
 *   n
 *
 *
 *    v    = Tv
 *     n+1     n
 *
 *
 *            m
 *   v     = T v
 *    n+m       n
 *
 *
 *  where m>0 propagates forward and m<0 propagates backward.
 *
 *  For SSRG, T is denoted Ts and is defined in matlab notation as:
 *
 *            Ts = [ [     c(N:-1:1)       ]
 *                   [eye(N-1) zeros(N-1,1)] ]
 *
 *  where c is a binary row vector of feeback tap coefficients,
 *  eye is an identity matrix, and zeros is an all zeros vector.
 *
 *
 *  For MSRG, T is denoted Tm and is defined in matlab notation as:
 *
 *            Tm = [ [zeros(1,N-1); eye(N-1)] c(N:-1:1).']
 *
 *  The T matrix formulations above are slightly different from, 
 *  but equivalent to [1] in order to accomodate efficient use of
 *  unsigned integer bit-wise arithmetic.
 *
 * % Reference: 
 *    [1] Lee & Miller, "CDMA Systems Engineering Handbook"
 *        Artech House, 1998, Chapter 6
 */

#include <cmath>
#include <algorithm>
#include <cctype>
#include <string>
#include "matrix.h"

/**************************
 * Helper Functions
 **************************/

void error(const std::string msg)
{
  std::cerr << msg << std::endl;
  exit(EXIT_FAILURE);
}

std::string lower(std::string s)
{
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  return s;
}

// Map logical input {0, 1} to one of the following formats:
//       type    output
//       input   format
//       -----   ---------------
//         0     logical (nop)
//         1     0 -> +1, 1 -> -1
//         2     0 -> -1, 1 -> +1
//        else   logical invert
int output(bool in, int type)
{
  int out;
  if (type == 0) {
    // nop
    out = (int) in;
  }
  else if (type == 1) {
    // 0 -> +1,  1 -> -1
    out = in ? -1 : 1;
  }
  else if (type == 2) {
    // 0 -> -1, 1 -> +1
    out = in ? 1 : -1;
  }
  else {
    // invert: 0 -> 1,  1 -> 0
    out = in ? 0 : 1;
  }
  return out;
}

// Compute sum modulo-2 of set bits in an unsigned integer
int parity(uint64_t in, int degree = 64) {
  int sum = 0;
  for (int d = degree - 1; d >= 0; --d) {
    sum += (in >> d) & 1;
  }
  return sum % 2;
}

// Bit reverse an unsigned integer of up to 64 bits long
uint64_t bitreverse(uint64_t in, int wordlength = 0) {
  uint64_t out = 0;
  if (wordlength == 0) {
    wordlength = int(std::log2(in) + 1);
  }
  while (wordlength--) {
    int bit = in & 1;
    out = (out << 1) | bit;
    in >>= 1;
  }
  return out;
}

// Convert a vector of binary digits to its decimal integer equivalent
uint64_t bi2de(std::vector<int> & seq, std::string flag = "left-msb")
{
  size_t sz = seq.size();
  uint64_t out = 0;
  if (sz > 64) {
    error("bi2de: Input vector length exceeds max of 64.");
  }
  if (lower(flag) == "left-msb") {
    for (size_t ii = 0; ii < sz; ++ii) {
      out = (out<<1) | seq[ii];
    }
  }
  else if (lower(flag) == "right-msb") {
    for (int ii = sz; ii > 0; --ii) {
      out = (out << 1) | seq[ii - 1];
    }
  }
  else {
    error("bi2de: unknown flag value");
  }
  return out;
}

// Convert the bits in a uint64_t to elements of a column vector
// MSB is placed at the first element (top of) the column vector
Matrix uint2vec(uint64_t in, int wordlength = 0, std::string flag = "left-msb")
{
  if (wordlength == 0) {
    wordlength = int(std::log2(in) + 1);
  }
  if (lower(flag) == "left-msb") {
    in = bitreverse(in, wordlength);
  }
  Matrix V(wordlength, 1);
  for (int ii = 1; ii <= wordlength; ++ii) {
    V(ii, 1) = double(in & 1);
    in = in >> 1;
  }
  return V;
}

// Put the elements of a vector into bit positions of an integer
uint64_t vec2uint(const Matrix & V, std::string vflag = "msb-first", std::string iflag = "left-msb")
{
  Matrix Vr = V; // operate on a row vector
  int rows = V.Size(1);
  int cols = V.Size(2);
  int wordlength = rows > cols ? rows : cols;
  if (rows > cols) {
    Vr = Transpose(V);
  }
  if (lower(vflag) == "msb-last") {
    Vr = FlipLR(Vr);
  }
  uint64_t out = 0;
  for (int ii = 1; ii <= wordlength; ++ii) {
    out = (out << 1) | int(Vr.get(1, ii));
  }
  if (lower(iflag) == "right-msb") {
    out = bitreverse(out, wordlength);
  }
  return out;
}

// Form MSRG characteristic matrix
// taps is a left-msb integer
Matrix msrg_char_mtx(uint64_t taps, int degree)
{
  uint64_t p = taps >> 1;
  Matrix Vp = Matrix(degree, 1);
  p = bitreverse(p);
  Vp = uint2vec(p, degree);
  Matrix IO = Vcat(Eye(degree - 1), Zeros(1, degree - 1));
  Matrix T = ToGF2(FlipUD(FlipLR(Hcat(Vp, IO))));
  return T;
}

// Form SSRG characteristic matrix
// taps is a left-msb integer
Matrix ssrg_char_mtx(uint64_t taps, int degree)
{
  uint64_t p = taps & ((1ull << degree) - 1);
  Matrix Vp = Transpose(uint2vec(p, degree));
  Matrix IO = Hcat(Eye(degree - 1), Zeros(degree - 1, 1));
  Matrix T = ToGF2(Vcat(Vp, IO));
  return T;
}

// Invert matrix in GF2
Matrix matInvGF2(const Matrix& A)
{
  return ToGF2(Inv(A));
}

// Multiply a vector by a Matrix in GF2
Matrix mvMultGF2(const Matrix& A, const Matrix & v)
{
  // Input dimensions enforced by matrix multiply
  int degree = A.Size(1);
  Matrix C = A * v;
  return ToGF2(C);
}

// Matrix square (A^2) in GF2
Matrix matSqrGF2(const Matrix& A)
{
  // Input dimensions enforced by matrix multiply
  return ToGF2(A * A);
}

struct lfsr_info {
  std::string type; // lfsr type: "ssrg" or "msrg"
  uint64_t taps;
  uint64_t fill;
  uint64_t mask;    // mask=0 means maskless lfsr
  int outtype;      // 0=logical{0,1} , 1=analytic{1,-1}, 2=analytic{-1,1}
};


/**************************
* LFSR Functions
**************************/

// LFSR generator using simple shift register (SSRG) aka Fibonacci structure
std::vector<int> lfsr_ssrg(int num, uint64_t taps, uint64_t & fill, int outtype = 0)
{
  int degree = int(std::log2(taps));
  uint64_t sr = fill;
  taps &= ((1ull << degree) - 1);
  std::vector<int> seq;
  for (int nn = 0; nn < num; ++nn) {
    seq.push_back(output(sr & 1, outtype));
    int parbit = parity(sr & taps , degree);
    sr = ((uint64_t) parbit << (degree - 1)) | (sr >> 1);
  }
  fill = sr;
  return seq;
}

// LFSR generator using modular shift register (SSRG) aka Galois structure
std::vector<int> lfsr_msrg(int num, uint64_t taps, uint64_t & fill, int outtype = 0)
{
  int degree = int(std::log2(taps));
  uint64_t sr = fill;
  taps &= ((1ull << degree) - 1);
  std::vector<int> seq;
  for (int nn = 0; nn < num; ++nn) {
    seq.push_back(output(sr & 1, outtype));
    if (sr & 1) {
      sr = long(1ull << (degree - 1)) | ((sr >> 1) ^ (taps >> 1));
    }
    else {
      sr = (0ull << (degree - 1)) | (sr >> 1);
    }
  }
  fill = sr;
  return seq;
}

// Masked LFSR generator using modular shift register (MSRG) aka Galois structure
std::vector<int> lfsr_msrg_mask(int num, uint64_t taps, uint64_t & fill, uint64_t mask, int outtype = 0)
{
  int degree = int(std::log2(taps));
  uint64_t sr = fill;
  taps &= ((1ull << degree) - 1);
  std::vector<int> seq;
  for (int nn = 0; nn < num; ++nn) {
    seq.push_back(output(parity(sr & mask, degree) != 0, outtype));
    if (sr & 1) {
      sr = long(1ull << (degree - 1)) | ((sr >> 1) ^ (taps >> 1));
    }
    else {
      sr = (0ull << (degree - 1)) | (sr >> 1);
    }
  }
  fill = sr;
  return seq;
}

// Masked LFSR generator using simple shift register (SSRG) aka Fibonacci structure
std::vector<int> lfsr_ssrg_mask(int num, uint64_t taps, uint64_t & fill, uint64_t mask, int outtype = 0)
{
  int degree = int(std::log2(taps));
  uint64_t sr = fill;
  taps &= ((1ull << degree) - 1);
  std::vector<int> seq;
  for (int nn = 0; nn < num; ++nn) {
    seq.push_back(output(parity(sr & mask, degree) != 0, outtype));
    int parbit = parity(sr & taps, degree);
    sr = ((uint64_t) parbit << (degree - 1)) | (sr >> 1);
  }
  fill = sr;
  return seq;
}

// Convert MSRG polynomial and fill to the
// equivalent SSRG polynomial and fill
uint64_t lfsr_msrg2ssrg(uint64_t mtaps, uint64_t mfill, uint64_t & sfill)
{
  int degree = int(std::log2(mtaps));
  uint64_t staps = bitreverse(mtaps);
  std::vector<int> seq = lfsr_msrg(degree, mtaps, mfill);
  sfill = bi2de(seq, "right-msb");
  return staps;
}

// Given a masked SSRG and its fill compute
// the equivalent fill for a maskless SSRG
uint64_t lfsr_ssrgmask2ssrg(uint64_t taps, uint64_t ifill, uint64_t mask)
{
  int num = int(std::log2(taps));
  std::vector<int> seq = lfsr_ssrg_mask(num, taps, ifill, mask);
  uint64_t sfill = bi2de(seq, "right-msb");
  return sfill;
}

// State Space generator formulation for MSRG
std::vector<int> lfsr_ssgm(int num, uint64_t taps, uint64_t & fill, int outtype = 0)
{
  int degree = int(std::log2(taps));
  uint64_t sr = fill;
  // Vsr is a column vector
  Matrix Vsr = uint2vec(sr, degree);
  // form msrg characteristic matrix Tm
  Matrix Tm = msrg_char_mtx(taps, degree);
  std::vector<int> seq;
  for (int ii = 0; ii < num; ++ii) {
    int outbit = int(Vsr.get(degree, 1));
    seq.push_back(output(outbit != 0, outtype));
    Vsr = mvMultGF2(Tm, Vsr);
  }
  fill = vec2uint(Vsr, "msb-first", "left-msb");
  return seq;
}

// State Space generator formulation for SSRG
std::vector<int> lfsr_ssgs(int num, uint64_t taps, uint64_t & fill, int outtype = 0)
{
  int degree = int(std::log2(taps));
  uint64_t sr = fill;
  // Vsr is a column vector
  Matrix Vsr = uint2vec(sr, degree); // MSB is in first element of Vsr
  // form ssrg characteristic matrix Ts
  Matrix Ts = ssrg_char_mtx(taps, degree);
  std::vector<int> seq;
  for (int ii = 0; ii < num; ++ii) {
    int outbit = int(Vsr.get(degree, 1));
    seq.push_back(output(outbit != 0, outtype));
    Vsr = mvMultGF2(Ts, Vsr);
  }
  fill = vec2uint(Vsr, "msb-first", "left-msb");
  return seq;
}

// Compute SSRG mask corresponding to a code phase change
// of "shift" states. Only shifts >= 0 are accepted.
// This function uses a brute-force approach that is O(n)
// in shifts. Use jump2mask for a more efficient and
// flexible algorithm. This function was used to validate
// jump2mask.
uint64_t lfsr_shift2mask(int shift, uint64_t taps)
{
  if (shift < 0) {
    error("lfsr_shift2mask: shift value must be >= 0.");
  }
  uint64_t fill = 1;
  lfsr_ssgm(shift, taps, fill);
  // final fill value is the mask
  return fill;
}

// Convert SSRG tap polynomial and fill value to
// the equivalent MSRG polynomial and fill value
uint64_t lfsr_ssrg2msrg(uint64_t staps, uint64_t sfill, uint64_t & mfill)
{
  uint64_t ssr = sfill;
  int degree = int(std::log2(staps));
  uint64_t mtaps = bitreverse(staps);
  // form msrg characteristic matrix Tm
  Matrix Tm = msrg_char_mtx(mtaps, degree);
  // mtaps_nolsb because nn loop uses all but the lsb of mtaps
  uint64_t mtaps_nolsb = (mtaps & ((1ull << degree) - 1)) >> 1;
  uint64_t msr = 0;
  // inject ssrg shift reg contents into msrg feedback path
  for (int nn = 0; nn < degree; ++nn) {
    if (ssr & 1) {
      msr = long(1ull << (degree - 1)) | ((msr >> 1) ^ mtaps_nolsb);
    }
    else {
      msr = (0ull << (degree - 1)) | (msr >> 1);
    }
    ssr >>= 1;
  }
  // invert Tm to propagate msrg state backwards
  Matrix invTm = matInvGF2(Tm);
  // Vmsr is a column vector of msrg shift reg contents
  Matrix Vmsr = uint2vec(msr, degree, "left-msb");
  // propagate msrg backwards
  for (int nn = 0; nn < degree; ++nn) {
    Vmsr = mvMultGF2(invTm, Vmsr);
  }
  msr = vec2uint(Vmsr, "msb-first", "left-msb");
  mfill = msr;
  return mtaps;
}

// Given a fill value "sr", one-step transition matrix "T"
// (aka characteristic matrix), and jump value "jump", compute the
// fill value corresponding to a code phase change of "jump" states.
// This function can be used for both SSRG and MSRG LFSRs.
// Jump must be a non-negative integer.
uint64_t jumpSR(int jump, Matrix T, uint64_t sr)
{
  if (jump < 0) {
    error("jumpSR: jump value must be >= 0.");
  }
  int N = int(std::log2(jump) + 1);
  int wordlength = T.Size(2);
  // Vsr is a column vector of shift reg contents
  Matrix Vsr = uint2vec(sr, wordlength);
  while (jump != 0) {
    if (jump & 1) {
      Vsr = mvMultGF2(T, Vsr); // propagtate state
    }
    T = matSqrGF2(T); // propagation "stride" *= 2
    jump >>= 1;
  }
  return vec2uint(Vsr, "msb-first", "left-msb");
}

// Compute the MSRG fill value corresponding to a code phase
// change of jump states. Jump can be positive or negative, 
// corresponding respectively to a code phase delay or advance.
uint64_t lfsr_ssgm_jump(int jump, uint64_t taps, uint64_t ifill)
{
  int degree = int(std::log2(taps));
  // form msrg characteristic matrix Tm
  Matrix Tm = msrg_char_mtx(taps, degree);
  if (jump < 0) {
    Tm = matInvGF2(Tm);
    jump = -jump;
  }
  uint64_t fill = jumpSR(jump, Tm, ifill); // propagate
  return fill;
}

// Compute the SSRG fill value corresponding to a code phase
// change of jump states. Jump can be positive or negative, 
// corresponding respectively to a code phase delay or advance.
uint64_t lfsr_ssgs_jump(int jump, uint64_t taps, uint64_t ifill)
{
  int degree = int(std::log2(taps));
  // form ssrg characteristic matrix Tm
  Matrix Tm = ssrg_char_mtx(taps, degree);
  if (jump < 0) {
    Tm = matInvGF2(Tm);
    jump = -jump;
  }
  uint64_t fill = jumpSR(jump, Tm, ifill); // propagate
  return fill;
}

// Given a jump value and SSRG tap polynomial compute
// mask value to effect a phase shift of jump states
// in a masked SSRG
uint64_t lfsr_jump2mask(int jump, uint64_t taps)
{
  return lfsr_ssgm_jump(jump, taps, 1);
}

// Given a jump value and MSRG tap polynomial compute
// mask value to effect a phase shift of jump states
// in a masked MSRG
uint64_t lfsr_msrg_jump2mask(int jump, uint64_t taps)
{
  int degree = int(std::log2(taps));
  uint64_t ifill = 1ull << (degree - 1);
  uint64_t fill = lfsr_ssgm_jump(-jump, taps, ifill);
  std::vector<int> seq = lfsr_msrg(degree, taps, fill);
  uint64_t mask = bi2de(seq, "left-msb");
  return mask;
}

/**************************
 *  LFSR class
 **************************/
// TODO: add "reset" member function
class lfsr {
  uint64_t taps_;
  uint64_t fill_;
  uint64_t mask_;
  Matrix T_;
  int outtype_;
  std::string type_;
public:
  explicit lfsr(uint64_t taps    = 0x19
              , uint64_t fill    = 1
              , std::string type = "ssrg"
              , uint64_t mask    = 0
              , int outtype      = 0)
  : taps_(taps), fill_(fill), type_(type), mask_(mask), outtype_(outtype)
  { }
  // Generator methods
  // Generate "num" code bits
  // Note: This method changes the state of the current LFSR object
  std::vector<int> operator()(const int num)
  {
    uint64_t tempfill = fill_;
    std::vector<int> seq;
    if (type_ == "ssrg") {
      if (mask_ == 0) {
        seq = lfsr_ssrg(num, taps_, tempfill, outtype_);
      } else {
        seq = lfsr_ssrg_mask(num, taps_, tempfill, mask_, outtype_);
      }
    } else  if (type_ == "msrg"){
      if (mask_ == 0) {
        seq = lfsr_msrg(num, taps_, tempfill, outtype_);
      } else {
        seq = lfsr_msrg_mask(num, taps_, tempfill, mask_, outtype_);
      }
    }
    fill_ = tempfill;
    return seq;
  }
  // Convert between LFSR types:
  // If the current LFSR object is SSRG then return
  // the equivalent MSRG polynomial and fill
  // If the current LFSR object is MSRG then return
  // the equivalent SSRG polynomial and fill
  // Note: This method does not change the state of
  // the current LFSR object but merely computes
  // the polynomial and fill for its dual
  lfsr_info convert() const
  {
    uint64_t newfill;
    uint64_t taps;
    std::string infotype;
    lfsr_info info;
    if (type_ == "ssrg") {
      taps = lfsr_ssrg2msrg(taps_, fill_, newfill);
      infotype = "msrg";
    } else if (type_ == "msrg") {
      taps = lfsr_msrg2ssrg(taps_, fill_, newfill);
      infotype = "ssrg";
    } else {
      error("convert: unknown generator type.");
    }
    return get_info(infotype, taps, newfill);
  }
  // Propagate the LFSR state by an amount
  // corresponding to "num" states
  // num>0 forward propagates the state relative to
  // num = 0 thereby advancing the output sequence
  // num<0 propagates backwards thereby delaying
  // the output sequence
  // Note: This method changes the state of the
  // LFSR object by updating its fill value
  uint64_t jump(int num)
  {
    uint64_t newfill;
    if (type_ == "ssrg") {
      newfill = lfsr_ssgs_jump(num, taps_, fill_);
    } else if (type_ == "msrg") {
      newfill = lfsr_ssgm_jump(num, taps_, fill_);
    } else {
      error("jump: unknown generator type.");
    }
    // update state per requested jump
    fill_ = newfill;
    return newfill;
  }
  // Calculate the LFSR mask corresponding to an
  // advance or delay of "num" states
  // num>0 produces a mask that delays the output 
  // sequence num bits relative to num=0
  // num<0 produces a mask that advances it
  // Note: This method does not change the
  // state of the LFSR object but merely computes
  // the mask value to use with a masked LFSR
  uint64_t calculate_mask(int num) const
  {
    uint64_t mask;
    if (type_ == "ssrg") {
      mask = lfsr_jump2mask(num, taps_);
    } else if (type_ == "msrg") {
      mask = lfsr_msrg_jump2mask(num, taps_);
    } else {
      error("calculate_mask: unknown generator type.");
    }
    return mask;
  }
  // Return shift register fill for the current LFSR object
  uint64_t get_fill() const
  {
    return fill_;
  }
  // Return generator polynomial for the current LFSR object
  uint64_t get_taps() const
  {
    return taps_;
  }
  // Return mask for the current LFSR object
  uint64_t get_mask() const
  {
    return mask_;
  }
  // Return LFSR type for the current LFSR object
  std::string get_type() const
  {
    return type_;
  }
  // Return characteristic matrix T
  Matrix get_T() const
  {
    Matrix T;
    int degree = int(std::log2(taps_));
    if (type_ == "ssrg") {
      T = msrg_char_mtx(taps_, degree);
    } else if (type_ == "msrg") {
      T = ssrg_char_mtx(taps_, degree);
    }
    return T;
  }
  // Return info structure for the current LFSR object
  lfsr_info get_info() const
  {
    lfsr_info info;
    info.type = type_;
    info.taps = taps_;
    info.fill = fill_;
    info.mask = mask_;
    info.outtype = outtype_;
    return info;
  }
  lfsr_info get_info(std::string type
                    , uint64_t taps
                    , uint64_t fill
                    , uint64_t mask = 0
                    , uint64_t outtype = 0) const
  {
    lfsr_info info;
    info.type = type;
    info.taps = taps;
    info.fill = fill;
    info.mask = mask;
    info.outtype = outtype;
    return info;
  }
  // Print info for the current LFSR object
  void PrintInfo() const
  {
    lfsr_info info = get_info();
    printf("LFSR type:\t\t%s\n", info.type.c_str());
    printf("LFSR taps:\t\t0x%lx\n", info.taps);
    printf("LFSR fill:\t\t0x%lx\n", info.fill);
    printf("LFSR mask:\t\t0x%lx\n", info.mask);
    printf("LFSR output type:\t%d\n", info.outtype);
  }
  void PrintInfo(lfsr_info info) const
  {
    printf("LFSR type:\t\t%s\n", info.type.c_str());
    printf("LFSR taps:\t\t0x%lx\n", info.taps);
    printf("LFSR fill:\t\t0x%lx\n", info.fill);
    printf("LFSR mask:\t\t0x%lx\n", info.mask);
    printf("LFSR output type:\t%d\n", info.outtype);
  }

};

