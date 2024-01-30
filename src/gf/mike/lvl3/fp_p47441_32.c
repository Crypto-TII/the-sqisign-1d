#ifdef RADIX_32

#include <stdint.h>
#include <stdio.h>

#include <x86intrin.h>

#include <stdbool.h>
#include <fp.h>

#define uspint uint32_t
#define sspint int32_t
#define spint uint32_t
#define dpint uint64_t

// propagate carries - return sign
static __attribute__((always_inline)) sspint prop(spint *n) {
  spint d, mask = ((spint)1 << 28) - 1;
  sspint carry = (sspint)n[0] >> 28;
  n[0] &= mask;
  for (int i = 1; i < 13; i++) {
    d = n[i] + carry;
    n[i] = d & mask;
    carry = (sspint)d >> 28;
  }
  n[13] += carry;
  return ((sspint)n[13] >> 31);
}

// propagate carries and add p if negative, propagate carries again
static __attribute__((noinline)) void flatten(spint *n) {
  spint q = ((spint)1 << 28);
  sspint carry = prop(n);
  n[0] -= 1 & carry;
  n[3] += (1 * (spint)0x74c2000) & carry;
  n[4] += (1 * (spint)0x4684c61) & carry;
  n[5] += (1 * (spint)0x69356ea) & carry;
  n[6] += (1 * (spint)0x1c722f6) & carry;
  n[7] += (1 * (spint)0x90aeb75) & carry;
  n[8] += (1 * (spint)0x5bc2e0a) & carry;
  n[9] += (1 * (spint)0xd10ad66) & carry;
  n[10] += (1 * (spint)0xe604a45) & carry;
  n[11] += (1 * (spint)0x71a2c6a) & carry;
  n[12] += (1 * (spint)0xeeeab08) & carry;
  n[13] += (1 * (spint)0x3df6) & carry;
  prop(n);
}

// Montgomery final subtract
void __attribute__((noinline)) modfsb(spint *n) {
  spint q = ((spint)1 << 28);
  n[0] += 1;
  n[3] -= 1 * (spint)0x74c2000;
  n[4] -= 1 * (spint)0x4684c61;
  n[5] -= 1 * (spint)0x69356ea;
  n[6] -= 1 * (spint)0x1c722f6;
  n[7] -= 1 * (spint)0x90aeb75;
  n[8] -= 1 * (spint)0x5bc2e0a;
  n[9] -= 1 * (spint)0xd10ad66;
  n[10] -= 1 * (spint)0xe604a45;
  n[11] -= 1 * (spint)0x71a2c6a;
  n[12] -= 1 * (spint)0xeeeab08;
  n[13] -= 1 * (spint)0x3df6;
  flatten(n);
}

// Modular addition - reduce less than 2p
void __attribute__((noinline)) modadd(spint *a, spint *b, spint *n) {
  spint q = ((spint)1 << 28);
  sspint carry;
  n[0] = a[0] + b[0];
  n[1] = a[1] + b[1];
  n[2] = a[2] + b[2];
  n[3] = a[3] + b[3];
  n[4] = a[4] + b[4];
  n[5] = a[5] + b[5];
  n[6] = a[6] + b[6];
  n[7] = a[7] + b[7];
  n[8] = a[8] + b[8];
  n[9] = a[9] + b[9];
  n[10] = a[10] + b[10];
  n[11] = a[11] + b[11];
  n[12] = a[12] + b[12];
  n[13] = a[13] + b[13];
  n[0] += 2;
  n[3] -= 2 * (spint)0x74c2000;
  n[4] -= 2 * (spint)0x4684c61;
  n[5] -= 2 * (spint)0x69356ea;
  n[6] -= 2 * (spint)0x1c722f6;
  n[7] -= 2 * (spint)0x90aeb75;
  n[8] -= 2 * (spint)0x5bc2e0a;
  n[9] -= 2 * (spint)0xd10ad66;
  n[10] -= 2 * (spint)0xe604a45;
  n[11] -= 2 * (spint)0x71a2c6a;
  n[12] -= 2 * (spint)0xeeeab08;
  n[13] -= 2 * (spint)0x3df6;
  carry = prop(n);
  n[0] -= 2 & carry;
  n[3] += (2 * (spint)0x74c2000) & carry;
  n[4] += (2 * (spint)0x4684c61) & carry;
  n[5] += (2 * (spint)0x69356ea) & carry;
  n[6] += (2 * (spint)0x1c722f6) & carry;
  n[7] += (2 * (spint)0x90aeb75) & carry;
  n[8] += (2 * (spint)0x5bc2e0a) & carry;
  n[9] += (2 * (spint)0xd10ad66) & carry;
  n[10] += (2 * (spint)0xe604a45) & carry;
  n[11] += (2 * (spint)0x71a2c6a) & carry;
  n[12] += (2 * (spint)0xeeeab08) & carry;
  n[13] += (2 * (spint)0x3df6) & carry;
  prop(n);
}

// Modular subtraction - reduce less than 2p
void __attribute__((noinline)) modsub(spint *a, spint *b, spint *n) {
  spint q = ((spint)1 << 28);
  sspint carry;
  n[0] = a[0] - b[0];
  n[1] = a[1] - b[1];
  n[2] = a[2] - b[2];
  n[3] = a[3] - b[3];
  n[4] = a[4] - b[4];
  n[5] = a[5] - b[5];
  n[6] = a[6] - b[6];
  n[7] = a[7] - b[7];
  n[8] = a[8] - b[8];
  n[9] = a[9] - b[9];
  n[10] = a[10] - b[10];
  n[11] = a[11] - b[11];
  n[12] = a[12] - b[12];
  n[13] = a[13] - b[13];
  carry = prop(n);
  n[0] -= 2 & carry;
  n[3] += (2 * (spint)0x74c2000) & carry;
  n[4] += (2 * (spint)0x4684c61) & carry;
  n[5] += (2 * (spint)0x69356ea) & carry;
  n[6] += (2 * (spint)0x1c722f6) & carry;
  n[7] += (2 * (spint)0x90aeb75) & carry;
  n[8] += (2 * (spint)0x5bc2e0a) & carry;
  n[9] += (2 * (spint)0xd10ad66) & carry;
  n[10] += (2 * (spint)0xe604a45) & carry;
  n[11] += (2 * (spint)0x71a2c6a) & carry;
  n[12] += (2 * (spint)0xeeeab08) & carry;
  n[13] += (2 * (spint)0x3df6) & carry;
  prop(n);
}

// Overflow limit   = 18446744073709551616
// maximum possible = 1388517001176019771
// Modular multiplication, c=a*b mod 2p
void __attribute__((noinline)) modmul(spint *a, spint *b, spint *c) {
  spint v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13;
  dpint t = 0;
  spint p3 = 0x74c2000;
  spint p4 = 0x4684c61;
  spint p5 = 0x69356ea;
  spint p6 = 0x1c722f6;
  spint p7 = 0x90aeb75;
  spint p8 = 0x5bc2e0a;
  spint p9 = 0xd10ad66;
  spint p10 = 0xe604a45;
  spint p11 = 0x71a2c6a;
  spint p12 = 0xeeeab08;
  spint p13 = 0x3df6;
  spint s, q = ((spint)1 << 28); // q is unsaturated radix
  spint mask = q - 1;
  t += (dpint)a[0] * b[0];
  v0 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[1];
  t += (dpint)a[1] * b[0];
  v1 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[2];
  t += (dpint)a[1] * b[1];
  t += (dpint)a[2] * b[0];
  v2 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[3];
  t += (dpint)a[1] * b[2];
  t += (dpint)a[2] * b[1];
  t += (dpint)a[3] * b[0];
  t += (dpint)v0 * p3;
  v3 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[4];
  t += (dpint)a[1] * b[3];
  t += (dpint)a[2] * b[2];
  t += (dpint)a[3] * b[1];
  t += (dpint)a[4] * b[0];
  t += (dpint)v0 * p4;
  t += (dpint)v1 * p3;
  v4 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[5];
  t += (dpint)a[1] * b[4];
  t += (dpint)a[2] * b[3];
  t += (dpint)a[3] * b[2];
  t += (dpint)a[4] * b[1];
  t += (dpint)a[5] * b[0];
  t += (dpint)v0 * p5;
  t += (dpint)v1 * p4;
  t += (dpint)v2 * p3;
  v5 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[6];
  t += (dpint)a[1] * b[5];
  t += (dpint)a[2] * b[4];
  t += (dpint)a[3] * b[3];
  t += (dpint)a[4] * b[2];
  t += (dpint)a[5] * b[1];
  t += (dpint)a[6] * b[0];
  t += (dpint)v0 * p6;
  t += (dpint)v1 * p5;
  t += (dpint)v2 * p4;
  t += (dpint)v3 * p3;
  v6 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[7];
  t += (dpint)a[1] * b[6];
  t += (dpint)a[2] * b[5];
  t += (dpint)a[3] * b[4];
  t += (dpint)a[4] * b[3];
  t += (dpint)a[5] * b[2];
  t += (dpint)a[6] * b[1];
  t += (dpint)a[7] * b[0];
  t += (dpint)v0 * p7;
  t += (dpint)v1 * p6;
  t += (dpint)v2 * p5;
  t += (dpint)v3 * p4;
  t += (dpint)v4 * p3;
  v7 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[8];
  t += (dpint)a[1] * b[7];
  t += (dpint)a[2] * b[6];
  t += (dpint)a[3] * b[5];
  t += (dpint)a[4] * b[4];
  t += (dpint)a[5] * b[3];
  t += (dpint)a[6] * b[2];
  t += (dpint)a[7] * b[1];
  t += (dpint)a[8] * b[0];
  t += (dpint)v0 * p8;
  t += (dpint)v1 * p7;
  t += (dpint)v2 * p6;
  t += (dpint)v3 * p5;
  t += (dpint)v4 * p4;
  t += (dpint)v5 * p3;
  v8 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[9];
  t += (dpint)a[1] * b[8];
  t += (dpint)a[2] * b[7];
  t += (dpint)a[3] * b[6];
  t += (dpint)a[4] * b[5];
  t += (dpint)a[5] * b[4];
  t += (dpint)a[6] * b[3];
  t += (dpint)a[7] * b[2];
  t += (dpint)a[8] * b[1];
  t += (dpint)a[9] * b[0];
  t += (dpint)v0 * p9;
  t += (dpint)v1 * p8;
  t += (dpint)v2 * p7;
  t += (dpint)v3 * p6;
  t += (dpint)v4 * p5;
  t += (dpint)v5 * p4;
  t += (dpint)v6 * p3;
  v9 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[10];
  t += (dpint)a[1] * b[9];
  t += (dpint)a[2] * b[8];
  t += (dpint)a[3] * b[7];
  t += (dpint)a[4] * b[6];
  t += (dpint)a[5] * b[5];
  t += (dpint)a[6] * b[4];
  t += (dpint)a[7] * b[3];
  t += (dpint)a[8] * b[2];
  t += (dpint)a[9] * b[1];
  t += (dpint)a[10] * b[0];
  t += (dpint)v0 * p10;
  t += (dpint)v1 * p9;
  t += (dpint)v2 * p8;
  t += (dpint)v3 * p7;
  t += (dpint)v4 * p6;
  t += (dpint)v5 * p5;
  t += (dpint)v6 * p4;
  t += (dpint)v7 * p3;
  v10 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[11];
  t += (dpint)a[1] * b[10];
  t += (dpint)a[2] * b[9];
  t += (dpint)a[3] * b[8];
  t += (dpint)a[4] * b[7];
  t += (dpint)a[5] * b[6];
  t += (dpint)a[6] * b[5];
  t += (dpint)a[7] * b[4];
  t += (dpint)a[8] * b[3];
  t += (dpint)a[9] * b[2];
  t += (dpint)a[10] * b[1];
  t += (dpint)a[11] * b[0];
  t += (dpint)v0 * p11;
  t += (dpint)v1 * p10;
  t += (dpint)v2 * p9;
  t += (dpint)v3 * p8;
  t += (dpint)v4 * p7;
  t += (dpint)v5 * p6;
  t += (dpint)v6 * p5;
  t += (dpint)v7 * p4;
  t += (dpint)v8 * p3;
  v11 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[12];
  t += (dpint)a[1] * b[11];
  t += (dpint)a[2] * b[10];
  t += (dpint)a[3] * b[9];
  t += (dpint)a[4] * b[8];
  t += (dpint)a[5] * b[7];
  t += (dpint)a[6] * b[6];
  t += (dpint)a[7] * b[5];
  t += (dpint)a[8] * b[4];
  t += (dpint)a[9] * b[3];
  t += (dpint)a[10] * b[2];
  t += (dpint)a[11] * b[1];
  t += (dpint)a[12] * b[0];
  t += (dpint)v0 * p12;
  t += (dpint)v1 * p11;
  t += (dpint)v2 * p10;
  t += (dpint)v3 * p9;
  t += (dpint)v4 * p8;
  t += (dpint)v5 * p7;
  t += (dpint)v6 * p6;
  t += (dpint)v7 * p5;
  t += (dpint)v8 * p4;
  t += (dpint)v9 * p3;
  v12 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[13];
  t += (dpint)a[1] * b[12];
  t += (dpint)a[2] * b[11];
  t += (dpint)a[3] * b[10];
  t += (dpint)a[4] * b[9];
  t += (dpint)a[5] * b[8];
  t += (dpint)a[6] * b[7];
  t += (dpint)a[7] * b[6];
  t += (dpint)a[8] * b[5];
  t += (dpint)a[9] * b[4];
  t += (dpint)a[10] * b[3];
  t += (dpint)a[11] * b[2];
  t += (dpint)a[12] * b[1];
  t += (dpint)a[13] * b[0];
  t += (dpint)v0 * p13;
  t += (dpint)v1 * p12;
  t += (dpint)v2 * p11;
  t += (dpint)v3 * p10;
  t += (dpint)v4 * p9;
  t += (dpint)v5 * p8;
  t += (dpint)v6 * p7;
  t += (dpint)v7 * p6;
  t += (dpint)v8 * p5;
  t += (dpint)v9 * p4;
  t += (dpint)v10 * p3;
  v13 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[1] * b[13];
  t += (dpint)a[2] * b[12];
  t += (dpint)a[3] * b[11];
  t += (dpint)a[4] * b[10];
  t += (dpint)a[5] * b[9];
  t += (dpint)a[6] * b[8];
  t += (dpint)a[7] * b[7];
  t += (dpint)a[8] * b[6];
  t += (dpint)a[9] * b[5];
  t += (dpint)a[10] * b[4];
  t += (dpint)a[11] * b[3];
  t += (dpint)a[12] * b[2];
  t += (dpint)a[13] * b[1];
  t += (dpint)v1 * p13;
  t += (dpint)v2 * p12;
  t += (dpint)v3 * p11;
  t += (dpint)v4 * p10;
  t += (dpint)v5 * p9;
  t += (dpint)v6 * p8;
  t += (dpint)v7 * p7;
  t += (dpint)v8 * p6;
  t += (dpint)v9 * p5;
  t += (dpint)v10 * p4;
  t += (dpint)v11 * p3;
  c[0] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[2] * b[13];
  t += (dpint)a[3] * b[12];
  t += (dpint)a[4] * b[11];
  t += (dpint)a[5] * b[10];
  t += (dpint)a[6] * b[9];
  t += (dpint)a[7] * b[8];
  t += (dpint)a[8] * b[7];
  t += (dpint)a[9] * b[6];
  t += (dpint)a[10] * b[5];
  t += (dpint)a[11] * b[4];
  t += (dpint)a[12] * b[3];
  t += (dpint)a[13] * b[2];
  t += (dpint)v2 * p13;
  t += (dpint)v3 * p12;
  t += (dpint)v4 * p11;
  t += (dpint)v5 * p10;
  t += (dpint)v6 * p9;
  t += (dpint)v7 * p8;
  t += (dpint)v8 * p7;
  t += (dpint)v9 * p6;
  t += (dpint)v10 * p5;
  t += (dpint)v11 * p4;
  t += (dpint)v12 * p3;
  c[1] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[3] * b[13];
  t += (dpint)a[4] * b[12];
  t += (dpint)a[5] * b[11];
  t += (dpint)a[6] * b[10];
  t += (dpint)a[7] * b[9];
  t += (dpint)a[8] * b[8];
  t += (dpint)a[9] * b[7];
  t += (dpint)a[10] * b[6];
  t += (dpint)a[11] * b[5];
  t += (dpint)a[12] * b[4];
  t += (dpint)a[13] * b[3];
  t += (dpint)v3 * p13;
  t += (dpint)v4 * p12;
  t += (dpint)v5 * p11;
  t += (dpint)v6 * p10;
  t += (dpint)v7 * p9;
  t += (dpint)v8 * p8;
  t += (dpint)v9 * p7;
  t += (dpint)v10 * p6;
  t += (dpint)v11 * p5;
  t += (dpint)v12 * p4;
  t += (dpint)v13 * p3;
  c[2] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[4] * b[13];
  t += (dpint)a[5] * b[12];
  t += (dpint)a[6] * b[11];
  t += (dpint)a[7] * b[10];
  t += (dpint)a[8] * b[9];
  t += (dpint)a[9] * b[8];
  t += (dpint)a[10] * b[7];
  t += (dpint)a[11] * b[6];
  t += (dpint)a[12] * b[5];
  t += (dpint)a[13] * b[4];
  t += (dpint)v4 * p13;
  t += (dpint)v5 * p12;
  t += (dpint)v6 * p11;
  t += (dpint)v7 * p10;
  t += (dpint)v8 * p9;
  t += (dpint)v9 * p8;
  t += (dpint)v10 * p7;
  t += (dpint)v11 * p6;
  t += (dpint)v12 * p5;
  t += (dpint)v13 * p4;
  c[3] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[5] * b[13];
  t += (dpint)a[6] * b[12];
  t += (dpint)a[7] * b[11];
  t += (dpint)a[8] * b[10];
  t += (dpint)a[9] * b[9];
  t += (dpint)a[10] * b[8];
  t += (dpint)a[11] * b[7];
  t += (dpint)a[12] * b[6];
  t += (dpint)a[13] * b[5];
  t += (dpint)v5 * p13;
  t += (dpint)v6 * p12;
  t += (dpint)v7 * p11;
  t += (dpint)v8 * p10;
  t += (dpint)v9 * p9;
  t += (dpint)v10 * p8;
  t += (dpint)v11 * p7;
  t += (dpint)v12 * p6;
  t += (dpint)v13 * p5;
  c[4] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[6] * b[13];
  t += (dpint)a[7] * b[12];
  t += (dpint)a[8] * b[11];
  t += (dpint)a[9] * b[10];
  t += (dpint)a[10] * b[9];
  t += (dpint)a[11] * b[8];
  t += (dpint)a[12] * b[7];
  t += (dpint)a[13] * b[6];
  t += (dpint)v6 * p13;
  t += (dpint)v7 * p12;
  t += (dpint)v8 * p11;
  t += (dpint)v9 * p10;
  t += (dpint)v10 * p9;
  t += (dpint)v11 * p8;
  t += (dpint)v12 * p7;
  t += (dpint)v13 * p6;
  c[5] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[7] * b[13];
  t += (dpint)a[8] * b[12];
  t += (dpint)a[9] * b[11];
  t += (dpint)a[10] * b[10];
  t += (dpint)a[11] * b[9];
  t += (dpint)a[12] * b[8];
  t += (dpint)a[13] * b[7];
  t += (dpint)v7 * p13;
  t += (dpint)v8 * p12;
  t += (dpint)v9 * p11;
  t += (dpint)v10 * p10;
  t += (dpint)v11 * p9;
  t += (dpint)v12 * p8;
  t += (dpint)v13 * p7;
  c[6] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[8] * b[13];
  t += (dpint)a[9] * b[12];
  t += (dpint)a[10] * b[11];
  t += (dpint)a[11] * b[10];
  t += (dpint)a[12] * b[9];
  t += (dpint)a[13] * b[8];
  t += (dpint)v8 * p13;
  t += (dpint)v9 * p12;
  t += (dpint)v10 * p11;
  t += (dpint)v11 * p10;
  t += (dpint)v12 * p9;
  t += (dpint)v13 * p8;
  c[7] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[9] * b[13];
  t += (dpint)a[10] * b[12];
  t += (dpint)a[11] * b[11];
  t += (dpint)a[12] * b[10];
  t += (dpint)a[13] * b[9];
  t += (dpint)v9 * p13;
  t += (dpint)v10 * p12;
  t += (dpint)v11 * p11;
  t += (dpint)v12 * p10;
  t += (dpint)v13 * p9;
  c[8] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[10] * b[13];
  t += (dpint)a[11] * b[12];
  t += (dpint)a[12] * b[11];
  t += (dpint)a[13] * b[10];
  t += (dpint)v10 * p13;
  t += (dpint)v11 * p12;
  t += (dpint)v12 * p11;
  t += (dpint)v13 * p10;
  c[9] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[11] * b[13];
  t += (dpint)a[12] * b[12];
  t += (dpint)a[13] * b[11];
  t += (dpint)v11 * p13;
  t += (dpint)v12 * p12;
  t += (dpint)v13 * p11;
  c[10] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[12] * b[13];
  t += (dpint)a[13] * b[12];
  t += (dpint)v12 * p13;
  t += (dpint)v13 * p12;
  c[11] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[13] * b[13];
  t += (dpint)v13 * p13;
  c[12] = (spint)((uspint)t & mask);
  t >>= 28;
  c[13] = (spint)t;
}

// Modular squaring, c=a*a  mod 2p
void __attribute__((noinline)) modsqr(spint *a, spint *c) {
  spint v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13;
  dpint tot, t = 0;
  spint p3 = 0x74c2000;
  spint p4 = 0x4684c61;
  spint p5 = 0x69356ea;
  spint p6 = 0x1c722f6;
  spint p7 = 0x90aeb75;
  spint p8 = 0x5bc2e0a;
  spint p9 = 0xd10ad66;
  spint p10 = 0xe604a45;
  spint p11 = 0x71a2c6a;
  spint p12 = 0xeeeab08;
  spint p13 = 0x3df6;
  spint s, q = ((spint)1 << 28); // q is unsaturated radix
  spint mask = q - 1;
  tot = (dpint)a[0] * a[0];
  t = tot;
  v0 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[1];
  tot *= 2;
  t += tot;
  v1 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[2];
  tot *= 2;
  tot += (dpint)a[1] * a[1];
  t += tot;
  v2 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[3];
  tot += (dpint)a[1] * a[2];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p3;
  v3 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[4];
  tot += (dpint)a[1] * a[3];
  tot *= 2;
  tot += (dpint)a[2] * a[2];
  t += tot;
  t += (dpint)v0 * p4;
  t += (dpint)v1 * p3;
  v4 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[5];
  tot += (dpint)a[1] * a[4];
  tot += (dpint)a[2] * a[3];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p5;
  t += (dpint)v1 * p4;
  t += (dpint)v2 * p3;
  v5 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[6];
  tot += (dpint)a[1] * a[5];
  tot += (dpint)a[2] * a[4];
  tot *= 2;
  tot += (dpint)a[3] * a[3];
  t += tot;
  t += (dpint)v0 * p6;
  t += (dpint)v1 * p5;
  t += (dpint)v2 * p4;
  t += (dpint)v3 * p3;
  v6 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[7];
  tot += (dpint)a[1] * a[6];
  tot += (dpint)a[2] * a[5];
  tot += (dpint)a[3] * a[4];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p7;
  t += (dpint)v1 * p6;
  t += (dpint)v2 * p5;
  t += (dpint)v3 * p4;
  t += (dpint)v4 * p3;
  v7 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[8];
  tot += (dpint)a[1] * a[7];
  tot += (dpint)a[2] * a[6];
  tot += (dpint)a[3] * a[5];
  tot *= 2;
  tot += (dpint)a[4] * a[4];
  t += tot;
  t += (dpint)v0 * p8;
  t += (dpint)v1 * p7;
  t += (dpint)v2 * p6;
  t += (dpint)v3 * p5;
  t += (dpint)v4 * p4;
  t += (dpint)v5 * p3;
  v8 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[9];
  tot += (dpint)a[1] * a[8];
  tot += (dpint)a[2] * a[7];
  tot += (dpint)a[3] * a[6];
  tot += (dpint)a[4] * a[5];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p9;
  t += (dpint)v1 * p8;
  t += (dpint)v2 * p7;
  t += (dpint)v3 * p6;
  t += (dpint)v4 * p5;
  t += (dpint)v5 * p4;
  t += (dpint)v6 * p3;
  v9 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[10];
  tot += (dpint)a[1] * a[9];
  tot += (dpint)a[2] * a[8];
  tot += (dpint)a[3] * a[7];
  tot += (dpint)a[4] * a[6];
  tot *= 2;
  tot += (dpint)a[5] * a[5];
  t += tot;
  t += (dpint)v0 * p10;
  t += (dpint)v1 * p9;
  t += (dpint)v2 * p8;
  t += (dpint)v3 * p7;
  t += (dpint)v4 * p6;
  t += (dpint)v5 * p5;
  t += (dpint)v6 * p4;
  t += (dpint)v7 * p3;
  v10 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[11];
  tot += (dpint)a[1] * a[10];
  tot += (dpint)a[2] * a[9];
  tot += (dpint)a[3] * a[8];
  tot += (dpint)a[4] * a[7];
  tot += (dpint)a[5] * a[6];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p11;
  t += (dpint)v1 * p10;
  t += (dpint)v2 * p9;
  t += (dpint)v3 * p8;
  t += (dpint)v4 * p7;
  t += (dpint)v5 * p6;
  t += (dpint)v6 * p5;
  t += (dpint)v7 * p4;
  t += (dpint)v8 * p3;
  v11 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[12];
  tot += (dpint)a[1] * a[11];
  tot += (dpint)a[2] * a[10];
  tot += (dpint)a[3] * a[9];
  tot += (dpint)a[4] * a[8];
  tot += (dpint)a[5] * a[7];
  tot *= 2;
  tot += (dpint)a[6] * a[6];
  t += tot;
  t += (dpint)v0 * p12;
  t += (dpint)v1 * p11;
  t += (dpint)v2 * p10;
  t += (dpint)v3 * p9;
  t += (dpint)v4 * p8;
  t += (dpint)v5 * p7;
  t += (dpint)v6 * p6;
  t += (dpint)v7 * p5;
  t += (dpint)v8 * p4;
  t += (dpint)v9 * p3;
  v12 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[13];
  tot += (dpint)a[1] * a[12];
  tot += (dpint)a[2] * a[11];
  tot += (dpint)a[3] * a[10];
  tot += (dpint)a[4] * a[9];
  tot += (dpint)a[5] * a[8];
  tot += (dpint)a[6] * a[7];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p13;
  t += (dpint)v1 * p12;
  t += (dpint)v2 * p11;
  t += (dpint)v3 * p10;
  t += (dpint)v4 * p9;
  t += (dpint)v5 * p8;
  t += (dpint)v6 * p7;
  t += (dpint)v7 * p6;
  t += (dpint)v8 * p5;
  t += (dpint)v9 * p4;
  t += (dpint)v10 * p3;
  v13 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[1] * a[13];
  tot += (dpint)a[2] * a[12];
  tot += (dpint)a[3] * a[11];
  tot += (dpint)a[4] * a[10];
  tot += (dpint)a[5] * a[9];
  tot += (dpint)a[6] * a[8];
  tot *= 2;
  tot += (dpint)a[7] * a[7];
  t += tot;
  t += (dpint)v1 * p13;
  t += (dpint)v2 * p12;
  t += (dpint)v3 * p11;
  t += (dpint)v4 * p10;
  t += (dpint)v5 * p9;
  t += (dpint)v6 * p8;
  t += (dpint)v7 * p7;
  t += (dpint)v8 * p6;
  t += (dpint)v9 * p5;
  t += (dpint)v10 * p4;
  t += (dpint)v11 * p3;
  c[0] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[2] * a[13];
  tot += (dpint)a[3] * a[12];
  tot += (dpint)a[4] * a[11];
  tot += (dpint)a[5] * a[10];
  tot += (dpint)a[6] * a[9];
  tot += (dpint)a[7] * a[8];
  tot *= 2;
  t += tot;
  t += (dpint)v2 * p13;
  t += (dpint)v3 * p12;
  t += (dpint)v4 * p11;
  t += (dpint)v5 * p10;
  t += (dpint)v6 * p9;
  t += (dpint)v7 * p8;
  t += (dpint)v8 * p7;
  t += (dpint)v9 * p6;
  t += (dpint)v10 * p5;
  t += (dpint)v11 * p4;
  t += (dpint)v12 * p3;
  c[1] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[3] * a[13];
  tot += (dpint)a[4] * a[12];
  tot += (dpint)a[5] * a[11];
  tot += (dpint)a[6] * a[10];
  tot += (dpint)a[7] * a[9];
  tot *= 2;
  tot += (dpint)a[8] * a[8];
  t += tot;
  t += (dpint)v3 * p13;
  t += (dpint)v4 * p12;
  t += (dpint)v5 * p11;
  t += (dpint)v6 * p10;
  t += (dpint)v7 * p9;
  t += (dpint)v8 * p8;
  t += (dpint)v9 * p7;
  t += (dpint)v10 * p6;
  t += (dpint)v11 * p5;
  t += (dpint)v12 * p4;
  t += (dpint)v13 * p3;
  c[2] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[4] * a[13];
  tot += (dpint)a[5] * a[12];
  tot += (dpint)a[6] * a[11];
  tot += (dpint)a[7] * a[10];
  tot += (dpint)a[8] * a[9];
  tot *= 2;
  t += tot;
  t += (dpint)v4 * p13;
  t += (dpint)v5 * p12;
  t += (dpint)v6 * p11;
  t += (dpint)v7 * p10;
  t += (dpint)v8 * p9;
  t += (dpint)v9 * p8;
  t += (dpint)v10 * p7;
  t += (dpint)v11 * p6;
  t += (dpint)v12 * p5;
  t += (dpint)v13 * p4;
  c[3] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[5] * a[13];
  tot += (dpint)a[6] * a[12];
  tot += (dpint)a[7] * a[11];
  tot += (dpint)a[8] * a[10];
  tot *= 2;
  tot += (dpint)a[9] * a[9];
  t += tot;
  t += (dpint)v5 * p13;
  t += (dpint)v6 * p12;
  t += (dpint)v7 * p11;
  t += (dpint)v8 * p10;
  t += (dpint)v9 * p9;
  t += (dpint)v10 * p8;
  t += (dpint)v11 * p7;
  t += (dpint)v12 * p6;
  t += (dpint)v13 * p5;
  c[4] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[6] * a[13];
  tot += (dpint)a[7] * a[12];
  tot += (dpint)a[8] * a[11];
  tot += (dpint)a[9] * a[10];
  tot *= 2;
  t += tot;
  t += (dpint)v6 * p13;
  t += (dpint)v7 * p12;
  t += (dpint)v8 * p11;
  t += (dpint)v9 * p10;
  t += (dpint)v10 * p9;
  t += (dpint)v11 * p8;
  t += (dpint)v12 * p7;
  t += (dpint)v13 * p6;
  c[5] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[7] * a[13];
  tot += (dpint)a[8] * a[12];
  tot += (dpint)a[9] * a[11];
  tot *= 2;
  tot += (dpint)a[10] * a[10];
  t += tot;
  t += (dpint)v7 * p13;
  t += (dpint)v8 * p12;
  t += (dpint)v9 * p11;
  t += (dpint)v10 * p10;
  t += (dpint)v11 * p9;
  t += (dpint)v12 * p8;
  t += (dpint)v13 * p7;
  c[6] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[8] * a[13];
  tot += (dpint)a[9] * a[12];
  tot += (dpint)a[10] * a[11];
  tot *= 2;
  t += tot;
  t += (dpint)v8 * p13;
  t += (dpint)v9 * p12;
  t += (dpint)v10 * p11;
  t += (dpint)v11 * p10;
  t += (dpint)v12 * p9;
  t += (dpint)v13 * p8;
  c[7] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[9] * a[13];
  tot += (dpint)a[10] * a[12];
  tot *= 2;
  tot += (dpint)a[11] * a[11];
  t += tot;
  t += (dpint)v9 * p13;
  t += (dpint)v10 * p12;
  t += (dpint)v11 * p11;
  t += (dpint)v12 * p10;
  t += (dpint)v13 * p9;
  c[8] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[10] * a[13];
  tot += (dpint)a[11] * a[12];
  tot *= 2;
  t += tot;
  t += (dpint)v10 * p13;
  t += (dpint)v11 * p12;
  t += (dpint)v12 * p11;
  t += (dpint)v13 * p10;
  c[9] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[11] * a[13];
  tot *= 2;
  tot += (dpint)a[12] * a[12];
  t += tot;
  t += (dpint)v11 * p13;
  t += (dpint)v12 * p12;
  t += (dpint)v13 * p11;
  c[10] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[12] * a[13];
  tot *= 2;
  t += tot;
  t += (dpint)v12 * p13;
  t += (dpint)v13 * p12;
  c[11] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[13] * a[13];
  t += tot;
  t += (dpint)v13 * p13;
  c[12] = (spint)((uspint)t & mask);
  t >>= 28;
  c[13] = (spint)t;
}

// copy
void __attribute__((always_inline)) modcpy(spint *a, spint *c) {
  int i;
  for (i = 0; i < 14; i++) {
    c[i] = a[i];
  }
}

// Calculate progenitor
void modpro(spint *w, spint *z) {
  int i;
  spint x[14];
  spint t0[14], t1[14], t2[14], t3[14], t4[14], t5[14], t6[14], t7[14], t8[14],
      t9[14], t10[14], t11[14], t12[14], t13[14], t14[14], t15[14], t16[14],
      t17[14], t18[14], t19[14], t20[14], t21[14], t22[14], t23[14], t24[14];
  modcpy(w, x);
  modsqr(x, t0);
  modmul(x, t0, t1);
  modmul(x, t1, t8);
  modmul(t0, t8, t6);
  modmul(t8, t6, t20);
  modmul(t8, t20, t2);
  modmul(t0, t2, t4);
  modmul(t8, t4, z);
  modmul(t2, z, t5);
  modmul(t1, t5, t19);
  modmul(t8, t19, t7);
  modmul(t20, t7, t16);
  modmul(z, t16, t10);
  modmul(t20, t10, t18);
  modmul(t8, t18, t22);
  modmul(t0, t22, t13);
  modmul(t6, t13, t2);
  modmul(t8, t2, t21);
  modmul(t0, t21, t3);
  modmul(t5, t3, t15);
  modmul(t6, t15, t9);
  modmul(t6, t9, t5);
  modmul(t4, t5, t4);
  modmul(t0, t4, t11);
  modmul(t20, t11, t17);
  modmul(t20, t17, t12);
  modmul(t8, t12, t23);
  modmul(t6, t23, t14);
  modmul(z, t14, t6);
  modmul(t8, t6, t8);
  modmul(t0, t8, t0);
  modmul(t20, t8, t20);
  modmul(z, t20, z);
  for (i = 0; i < 2; i++) {
    modsqr(z, z);
  }
  modcpy(z, t24);
  for (i = 0; i < 5; i++) {
    modsqr(t24, t24);
  }
  modmul(z, t24, z);
  modmul(t11, z, z);
  modsqr(t24, t24);
  modmul(t0, t24, t24);
  modsqr(z, t0);
  modmul(x, t0, t0);
  for (i = 0; i < 8; i++) {
    modsqr(t24, t24);
  }
  modmul(t23, t24, t23);
  for (i = 0; i < 7; i++) {
    modsqr(t23, t23);
  }
  modmul(t22, t23, t22);
  for (i = 0; i < 8; i++) {
    modsqr(t22, t22);
  }
  modmul(t21, t22, t21);
  for (i = 0; i < 12; i++) {
    modsqr(t21, t21);
  }
  modmul(t20, t21, t20);
  for (i = 0; i < 9; i++) {
    modsqr(t20, t20);
  }
  modmul(t9, t20, t20);
  for (i = 0; i < 11; i++) {
    modsqr(t20, t20);
  }
  modmul(t6, t20, t20);
  for (i = 0; i < 6; i++) {
    modsqr(t20, t20);
  }
  modmul(t16, t20, t20);
  for (i = 0; i < 12; i++) {
    modsqr(t20, t20);
  }
  modmul(t19, t20, t19);
  for (i = 0; i < 10; i++) {
    modsqr(t19, t19);
  }
  modmul(t9, t19, t19);
  for (i = 0; i < 7; i++) {
    modsqr(t19, t19);
  }
  modmul(t18, t19, t18);
  for (i = 0; i < 12; i++) {
    modsqr(t18, t18);
  }
  modmul(t17, t18, t17);
  for (i = 0; i < 7; i++) {
    modsqr(t17, t17);
  }
  modmul(t16, t17, t16);
  for (i = 0; i < 10; i++) {
    modsqr(t16, t16);
  }
  modmul(t12, t16, t16);
  for (i = 0; i < 8; i++) {
    modsqr(t16, t16);
  }
  modmul(t15, t16, t15);
  for (i = 0; i < 8; i++) {
    modsqr(t15, t15);
  }
  modmul(t14, t15, t14);
  for (i = 0; i < 7; i++) {
    modsqr(t14, t14);
  }
  modmul(t7, t14, t14);
  for (i = 0; i < 11; i++) {
    modsqr(t14, t14);
  }
  modmul(t13, t14, t13);
  for (i = 0; i < 9; i++) {
    modsqr(t13, t13);
  }
  modmul(t12, t13, t12);
  for (i = 0; i < 9; i++) {
    modsqr(t12, t12);
  }
  modmul(t11, t12, t11);
  for (i = 0; i < 7; i++) {
    modsqr(t11, t11);
  }
  modmul(t10, t11, t10);
  for (i = 0; i < 10; i++) {
    modsqr(t10, t10);
  }
  modmul(t9, t10, t9);
  for (i = 0; i < 8; i++) {
    modsqr(t9, t9);
  }
  modmul(t8, t9, t8);
  for (i = 0; i < 6; i++) {
    modsqr(t8, t8);
  }
  modmul(t7, t8, t7);
  for (i = 0; i < 10; i++) {
    modsqr(t7, t7);
  }
  modmul(t6, t7, t6);
  for (i = 0; i < 7; i++) {
    modsqr(t6, t6);
  }
  modmul(t2, t6, t6);
  for (i = 0; i < 9; i++) {
    modsqr(t6, t6);
  }
  modmul(t5, t6, t5);
  for (i = 0; i < 8; i++) {
    modsqr(t5, t5);
  }
  modmul(t4, t5, t4);
  for (i = 0; i < 9; i++) {
    modsqr(t4, t4);
  }
  modmul(t3, t4, t3);
  for (i = 0; i < 11; i++) {
    modsqr(t3, t3);
  }
  modmul(t2, t3, t2);
  for (i = 0; i < 4; i++) {
    modsqr(t2, t2);
  }
  modmul(t1, t2, t1);
  for (i = 0; i < 21; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 16; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 16; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 16; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 16; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t0);
  for (i = 0; i < 15; i++) {
    modsqr(t0, t0);
  }
  modmul(z, t0, z);
}

// calculate inverse, provide progenitor h if available
void modinv(spint *x, spint *h, spint *z) {
  int i;
  spint s[14], t[14];
  if (h == NULL) {
    modpro(x, t);
  } else {
    modcpy(h, t);
  }
  modcpy(x, s);
  for (i = 0; i <= 1; i++) {
    modsqr(t, t);
  }
  modmul(s, t, z);
}

// Convert m to n-residue form, n=nres(m)
void nres(spint *m, spint *n) {
  spint c[14] = {0x17ee9b4, 0x2686646, 0xc47b3e8, 0x157ca1a, 0x5365e7b,
                 0xd3dc7f5, 0xdb4bdee, 0xbe00e06, 0x7ac68b0, 0x3771313,
                 0xc212c80, 0xfdb63da, 0x3db3011, 0x23cf};
  modmul(m, c, n);
}

// Convert n back to normal form, m=redc(n)
void redc(spint *n, spint *m) {
  spint c[14];
  c[0] = 1;
  for (int i = 1; i < 14; i++)
    c[i] = 0;
  modmul(n, c, m);
  modfsb(m);
}

// is unity?
int modis1(spint *a) {
  spint c[14];
  sspint c0, d = 0;
  redc(a, c);
  for (int i = 1; i < 14; i++) {
    d |= c[i];
  }
  c0 = (sspint)c[0];
  return (1 & ((d - 1) >> 28) & (((c0 ^ 1) - 1) >> 28));
}

// is zero?
int modis0(spint *a) {
  sspint d = 0;
  for (int i = 0; i < 14; i++) {
    d |= a[i];
  }
  return (1 & ((d - 1) >> 28));
}

// set to zero
void modzer(spint *a) {
  for (int i = 0; i < 14; i++)
    a[i] = 0;
}

// set to one
void modone(spint *a) {
  a[0] = 1;
  for (int i = 1; i < 14; i++)
    a[i] = 0;
  nres(a, a);
}

// Test for quadratic residue
int modqr(spint *h, spint *x) {
  int i;
  spint r[14];
  if (h == NULL) {
    modpro(x, r);
  } else {
    modcpy(h, r);
  }
  modsqr(r, r);
  modmul(r, x, r);
  return modis1(r);
}

// Modular square root, provide progenitor h if available, NULL if not
void modsqrt(spint *x, spint *h, spint *r) {
  spint s[14], y[14];
  if (h == NULL) {
    modpro(x, y);
  } else {
    modcpy(h, y);
  }
  modmul(y, x, s);
  modcpy(s, r);
}

// shift left by less than a word
void modshl(int n, spint *a) {
  a[14 - 1] = ((a[14 - 1] << n)) | (a[14 - 2] >> (28 - n));
  for (int i = 14 - 2; i > 0; i--) {
    a[i] = ((a[i] << n) & 0xfffffff) | (a[i - 1] >> (28 - n));
  }
  a[0] = (a[0] << n) & 0xfffffff;
}

// shift right by less than a word. Return shifted out part
int modshr(int n, spint *a) {
  spint r = a[0] & (((spint)1 << n) - 1);
  for (int i = 0; i < 14 - 1; i++) {
    a[i] = (a[i] >> n) | ((a[i + 1] << (28 - n)) & 0xfffffff);
  }
  a[14 - 1] = a[14 - 1] >> n;
  return r;
}

/* API functions calling generated code */

const digit_t p[NWORDS_ORDER] =  { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x4C6174C1, 0x356EA468, 0xC722F669, 0x90AEB751, 0x65BC2E0A, 0x45D10AD6, 0xC6AE604A, 0xAB0871A2, 0x03DF6EEE };

bool fp_is_zero(const digit_t* a) {
    return (bool) modis0(a);
}

void fp_copy(digit_t* out, const digit_t* a) {
    modcpy(a, out);
}

void fp_add(digit_t* out, const digit_t* a, const digit_t* b) {
    modadd(a, b, out);
    modfsb(out);
}

void fp_sub(digit_t* out, const digit_t* a, const digit_t* b) {
    modsub(a, b, out);
    modfsb(out);
}

void fp_neg(digit_t* out, const digit_t* a) {
    spint zero[NWORDS_FIELD];
    modzer(zero);
    modsub(zero, a, out);
    modfsb(out);
}

void fp_sqr(digit_t* out, const digit_t* a) {
    modsqr(a, out);
    modfsb(out);
}

void fp_mul(digit_t* out, const digit_t* a, const digit_t* b) {
    modmul(a, b, out);
    modfsb(out);
}

void fp_inv(digit_t* a) {
    modinv(a, NULL, a);
}

bool fp_is_square(const digit_t* a) {
    return (bool) modqr(NULL, a);
}

void fp_sqrt(digit_t* a) {
    modsqrt(a, NULL, a);
}

void fp_tomont(digit_t* out, const digit_t* a) {
    nres(a, out);
}

void fp_frommont(digit_t* out, const digit_t* a) {
    redc(a, out);
}

void fp_mont_setone(digit_t* out) {
    modone(out);
}

void fp_to_digit_array(digit_t* out, const digit_t* a) {
    digit_t x[NWORDS_FIELD];
    modcpy(a, x);
    for (int i = 0; i < NWORDS_ORDER; i++) {
        out[i] = 0;
    }
    for (int i = 0; i < 48; i++) {
        ((char *) out)[i] = x[0] & 0xff;
        modshr(8, x);
    }
}

void fp_from_digit_array(digit_t* out, const digit_t* a) {
    for (int i = 0; i < NWORDS_FIELD; i++) {
        out[i] = 0;
    }
    for (int i = 48 - 1; i >= 0; i--) {
        modshl(8, out);
        out[0] += (digit_t)((unsigned char *) a)[i];
    }
}

#endif
