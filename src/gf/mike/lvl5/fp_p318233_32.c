#if defined(RADIX_32) || defined(PQM4)

#include <stdint.h>
#include <stdio.h>

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
  for (int i = 1; i < 17; i++) {
    d = n[i] + carry;
    n[i] = d & mask;
    carry = (sspint)d >> 28;
  }
  n[17] += carry;
  return ((sspint)n[17] >> 31);
}

// propagate carries and add p if negative, propagate carries again
static __attribute__((noinline)) void flatten(spint *n) {
  spint q = ((spint)1 << 28);
  sspint carry = prop(n);
  n[0] -= 1 & carry;
  n[5] += (1 * (spint)0xada6e20) & carry;
  n[6] += (1 * (spint)0xe994c68) & carry;
  n[7] += (1 * (spint)0x781974c) & carry;
  n[8] += (1 * (spint)0xaf0a29a) & carry;
  n[9] += (1 * (spint)0xdea65f) & carry;
  n[10] += (1 * (spint)0xac5904a) & carry;
  n[11] += (1 * (spint)0x7d01fe3) & carry;
  n[12] += (1 * (spint)0xe632650) & carry;
  n[13] += (1 * (spint)0x9202bdb) & carry;
  n[14] += (1 * (spint)0x36936e7) & carry;
  n[15] += (1 * (spint)0x8c15b00) & carry;
  n[16] += (1 * (spint)0x8869bc6) & carry;
  n[17] += (1 * (spint)0x255946a) & carry;
  prop(n);
}

// Montgomery final subtract
void __attribute__((noinline)) modfsb(spint *n) {
  spint q = ((spint)1 << 28);
  n[0] += 1;
  n[5] -= 1 * (spint)0xada6e20;
  n[6] -= 1 * (spint)0xe994c68;
  n[7] -= 1 * (spint)0x781974c;
  n[8] -= 1 * (spint)0xaf0a29a;
  n[9] -= 1 * (spint)0xdea65f;
  n[10] -= 1 * (spint)0xac5904a;
  n[11] -= 1 * (spint)0x7d01fe3;
  n[12] -= 1 * (spint)0xe632650;
  n[13] -= 1 * (spint)0x9202bdb;
  n[14] -= 1 * (spint)0x36936e7;
  n[15] -= 1 * (spint)0x8c15b00;
  n[16] -= 1 * (spint)0x8869bc6;
  n[17] -= 1 * (spint)0x255946a;
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
  n[14] = a[14] + b[14];
  n[15] = a[15] + b[15];
  n[16] = a[16] + b[16];
  n[17] = a[17] + b[17];
  n[0] += 2;
  n[5] -= 2 * (spint)0xada6e20;
  n[6] -= 2 * (spint)0xe994c68;
  n[7] -= 2 * (spint)0x781974c;
  n[8] -= 2 * (spint)0xaf0a29a;
  n[9] -= 2 * (spint)0xdea65f;
  n[10] -= 2 * (spint)0xac5904a;
  n[11] -= 2 * (spint)0x7d01fe3;
  n[12] -= 2 * (spint)0xe632650;
  n[13] -= 2 * (spint)0x9202bdb;
  n[14] -= 2 * (spint)0x36936e7;
  n[15] -= 2 * (spint)0x8c15b00;
  n[16] -= 2 * (spint)0x8869bc6;
  n[17] -= 2 * (spint)0x255946a;
  carry = prop(n);
  n[0] -= 2 & carry;
  n[5] += (2 * (spint)0xada6e20) & carry;
  n[6] += (2 * (spint)0xe994c68) & carry;
  n[7] += (2 * (spint)0x781974c) & carry;
  n[8] += (2 * (spint)0xaf0a29a) & carry;
  n[9] += (2 * (spint)0xdea65f) & carry;
  n[10] += (2 * (spint)0xac5904a) & carry;
  n[11] += (2 * (spint)0x7d01fe3) & carry;
  n[12] += (2 * (spint)0xe632650) & carry;
  n[13] += (2 * (spint)0x9202bdb) & carry;
  n[14] += (2 * (spint)0x36936e7) & carry;
  n[15] += (2 * (spint)0x8c15b00) & carry;
  n[16] += (2 * (spint)0x8869bc6) & carry;
  n[17] += (2 * (spint)0x255946a) & carry;
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
  n[14] = a[14] - b[14];
  n[15] = a[15] - b[15];
  n[16] = a[16] - b[16];
  n[17] = a[17] - b[17];
  carry = prop(n);
  n[0] -= 2 & carry;
  n[5] += (2 * (spint)0xada6e20) & carry;
  n[6] += (2 * (spint)0xe994c68) & carry;
  n[7] += (2 * (spint)0x781974c) & carry;
  n[8] += (2 * (spint)0xaf0a29a) & carry;
  n[9] += (2 * (spint)0xdea65f) & carry;
  n[10] += (2 * (spint)0xac5904a) & carry;
  n[11] += (2 * (spint)0x7d01fe3) & carry;
  n[12] += (2 * (spint)0xe632650) & carry;
  n[13] += (2 * (spint)0x9202bdb) & carry;
  n[14] += (2 * (spint)0x36936e7) & carry;
  n[15] += (2 * (spint)0x8c15b00) & carry;
  n[16] += (2 * (spint)0x8869bc6) & carry;
  n[17] += (2 * (spint)0x255946a) & carry;
  prop(n);
}

// Overflow limit   = 18446744073709551616
// maximum possible = 1791946672152748246
// Modular multiplication, c=a*b mod 2p
void __attribute__((noinline)) modmul(spint *a, spint *b, spint *c) {
  spint v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15,
      v16, v17;
  dpint t = 0;
  spint p5 = 0xada6e20;
  spint p6 = 0xe994c68;
  spint p7 = 0x781974c;
  spint p8 = 0xaf0a29a;
  spint p9 = 0xdea65f;
  spint p10 = 0xac5904a;
  spint p11 = 0x7d01fe3;
  spint p12 = 0xe632650;
  spint p13 = 0x9202bdb;
  spint p14 = 0x36936e7;
  spint p15 = 0x8c15b00;
  spint p16 = 0x8869bc6;
  spint p17 = 0x255946a;
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
  v3 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[4];
  t += (dpint)a[1] * b[3];
  t += (dpint)a[2] * b[2];
  t += (dpint)a[3] * b[1];
  t += (dpint)a[4] * b[0];
  v4 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[5];
  t += (dpint)a[1] * b[4];
  t += (dpint)a[2] * b[3];
  t += (dpint)a[3] * b[2];
  t += (dpint)a[4] * b[1];
  t += (dpint)a[5] * b[0];
  t += (dpint)v0 * p5;
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
  v13 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[14];
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
  t += (dpint)a[14] * b[0];
  t += (dpint)v0 * p14;
  t += (dpint)v1 * p13;
  t += (dpint)v2 * p12;
  t += (dpint)v3 * p11;
  t += (dpint)v4 * p10;
  t += (dpint)v5 * p9;
  t += (dpint)v6 * p8;
  t += (dpint)v7 * p7;
  t += (dpint)v8 * p6;
  t += (dpint)v9 * p5;
  v14 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[15];
  t += (dpint)a[1] * b[14];
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
  t += (dpint)a[14] * b[1];
  t += (dpint)a[15] * b[0];
  t += (dpint)v0 * p15;
  t += (dpint)v1 * p14;
  t += (dpint)v2 * p13;
  t += (dpint)v3 * p12;
  t += (dpint)v4 * p11;
  t += (dpint)v5 * p10;
  t += (dpint)v6 * p9;
  t += (dpint)v7 * p8;
  t += (dpint)v8 * p7;
  t += (dpint)v9 * p6;
  t += (dpint)v10 * p5;
  v15 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[16];
  t += (dpint)a[1] * b[15];
  t += (dpint)a[2] * b[14];
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
  t += (dpint)a[14] * b[2];
  t += (dpint)a[15] * b[1];
  t += (dpint)a[16] * b[0];
  t += (dpint)v0 * p16;
  t += (dpint)v1 * p15;
  t += (dpint)v2 * p14;
  t += (dpint)v3 * p13;
  t += (dpint)v4 * p12;
  t += (dpint)v5 * p11;
  t += (dpint)v6 * p10;
  t += (dpint)v7 * p9;
  t += (dpint)v8 * p8;
  t += (dpint)v9 * p7;
  t += (dpint)v10 * p6;
  t += (dpint)v11 * p5;
  v16 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[0] * b[17];
  t += (dpint)a[1] * b[16];
  t += (dpint)a[2] * b[15];
  t += (dpint)a[3] * b[14];
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
  t += (dpint)a[14] * b[3];
  t += (dpint)a[15] * b[2];
  t += (dpint)a[16] * b[1];
  t += (dpint)a[17] * b[0];
  t += (dpint)v0 * p17;
  t += (dpint)v1 * p16;
  t += (dpint)v2 * p15;
  t += (dpint)v3 * p14;
  t += (dpint)v4 * p13;
  t += (dpint)v5 * p12;
  t += (dpint)v6 * p11;
  t += (dpint)v7 * p10;
  t += (dpint)v8 * p9;
  t += (dpint)v9 * p8;
  t += (dpint)v10 * p7;
  t += (dpint)v11 * p6;
  t += (dpint)v12 * p5;
  v17 = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[1] * b[17];
  t += (dpint)a[2] * b[16];
  t += (dpint)a[3] * b[15];
  t += (dpint)a[4] * b[14];
  t += (dpint)a[5] * b[13];
  t += (dpint)a[6] * b[12];
  t += (dpint)a[7] * b[11];
  t += (dpint)a[8] * b[10];
  t += (dpint)a[9] * b[9];
  t += (dpint)a[10] * b[8];
  t += (dpint)a[11] * b[7];
  t += (dpint)a[12] * b[6];
  t += (dpint)a[13] * b[5];
  t += (dpint)a[14] * b[4];
  t += (dpint)a[15] * b[3];
  t += (dpint)a[16] * b[2];
  t += (dpint)a[17] * b[1];
  t += (dpint)v1 * p17;
  t += (dpint)v2 * p16;
  t += (dpint)v3 * p15;
  t += (dpint)v4 * p14;
  t += (dpint)v5 * p13;
  t += (dpint)v6 * p12;
  t += (dpint)v7 * p11;
  t += (dpint)v8 * p10;
  t += (dpint)v9 * p9;
  t += (dpint)v10 * p8;
  t += (dpint)v11 * p7;
  t += (dpint)v12 * p6;
  t += (dpint)v13 * p5;
  c[0] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[2] * b[17];
  t += (dpint)a[3] * b[16];
  t += (dpint)a[4] * b[15];
  t += (dpint)a[5] * b[14];
  t += (dpint)a[6] * b[13];
  t += (dpint)a[7] * b[12];
  t += (dpint)a[8] * b[11];
  t += (dpint)a[9] * b[10];
  t += (dpint)a[10] * b[9];
  t += (dpint)a[11] * b[8];
  t += (dpint)a[12] * b[7];
  t += (dpint)a[13] * b[6];
  t += (dpint)a[14] * b[5];
  t += (dpint)a[15] * b[4];
  t += (dpint)a[16] * b[3];
  t += (dpint)a[17] * b[2];
  t += (dpint)v2 * p17;
  t += (dpint)v3 * p16;
  t += (dpint)v4 * p15;
  t += (dpint)v5 * p14;
  t += (dpint)v6 * p13;
  t += (dpint)v7 * p12;
  t += (dpint)v8 * p11;
  t += (dpint)v9 * p10;
  t += (dpint)v10 * p9;
  t += (dpint)v11 * p8;
  t += (dpint)v12 * p7;
  t += (dpint)v13 * p6;
  t += (dpint)v14 * p5;
  c[1] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[3] * b[17];
  t += (dpint)a[4] * b[16];
  t += (dpint)a[5] * b[15];
  t += (dpint)a[6] * b[14];
  t += (dpint)a[7] * b[13];
  t += (dpint)a[8] * b[12];
  t += (dpint)a[9] * b[11];
  t += (dpint)a[10] * b[10];
  t += (dpint)a[11] * b[9];
  t += (dpint)a[12] * b[8];
  t += (dpint)a[13] * b[7];
  t += (dpint)a[14] * b[6];
  t += (dpint)a[15] * b[5];
  t += (dpint)a[16] * b[4];
  t += (dpint)a[17] * b[3];
  t += (dpint)v3 * p17;
  t += (dpint)v4 * p16;
  t += (dpint)v5 * p15;
  t += (dpint)v6 * p14;
  t += (dpint)v7 * p13;
  t += (dpint)v8 * p12;
  t += (dpint)v9 * p11;
  t += (dpint)v10 * p10;
  t += (dpint)v11 * p9;
  t += (dpint)v12 * p8;
  t += (dpint)v13 * p7;
  t += (dpint)v14 * p6;
  t += (dpint)v15 * p5;
  c[2] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[4] * b[17];
  t += (dpint)a[5] * b[16];
  t += (dpint)a[6] * b[15];
  t += (dpint)a[7] * b[14];
  t += (dpint)a[8] * b[13];
  t += (dpint)a[9] * b[12];
  t += (dpint)a[10] * b[11];
  t += (dpint)a[11] * b[10];
  t += (dpint)a[12] * b[9];
  t += (dpint)a[13] * b[8];
  t += (dpint)a[14] * b[7];
  t += (dpint)a[15] * b[6];
  t += (dpint)a[16] * b[5];
  t += (dpint)a[17] * b[4];
  t += (dpint)v4 * p17;
  t += (dpint)v5 * p16;
  t += (dpint)v6 * p15;
  t += (dpint)v7 * p14;
  t += (dpint)v8 * p13;
  t += (dpint)v9 * p12;
  t += (dpint)v10 * p11;
  t += (dpint)v11 * p10;
  t += (dpint)v12 * p9;
  t += (dpint)v13 * p8;
  t += (dpint)v14 * p7;
  t += (dpint)v15 * p6;
  t += (dpint)v16 * p5;
  c[3] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[5] * b[17];
  t += (dpint)a[6] * b[16];
  t += (dpint)a[7] * b[15];
  t += (dpint)a[8] * b[14];
  t += (dpint)a[9] * b[13];
  t += (dpint)a[10] * b[12];
  t += (dpint)a[11] * b[11];
  t += (dpint)a[12] * b[10];
  t += (dpint)a[13] * b[9];
  t += (dpint)a[14] * b[8];
  t += (dpint)a[15] * b[7];
  t += (dpint)a[16] * b[6];
  t += (dpint)a[17] * b[5];
  t += (dpint)v5 * p17;
  t += (dpint)v6 * p16;
  t += (dpint)v7 * p15;
  t += (dpint)v8 * p14;
  t += (dpint)v9 * p13;
  t += (dpint)v10 * p12;
  t += (dpint)v11 * p11;
  t += (dpint)v12 * p10;
  t += (dpint)v13 * p9;
  t += (dpint)v14 * p8;
  t += (dpint)v15 * p7;
  t += (dpint)v16 * p6;
  t += (dpint)v17 * p5;
  c[4] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[6] * b[17];
  t += (dpint)a[7] * b[16];
  t += (dpint)a[8] * b[15];
  t += (dpint)a[9] * b[14];
  t += (dpint)a[10] * b[13];
  t += (dpint)a[11] * b[12];
  t += (dpint)a[12] * b[11];
  t += (dpint)a[13] * b[10];
  t += (dpint)a[14] * b[9];
  t += (dpint)a[15] * b[8];
  t += (dpint)a[16] * b[7];
  t += (dpint)a[17] * b[6];
  t += (dpint)v6 * p17;
  t += (dpint)v7 * p16;
  t += (dpint)v8 * p15;
  t += (dpint)v9 * p14;
  t += (dpint)v10 * p13;
  t += (dpint)v11 * p12;
  t += (dpint)v12 * p11;
  t += (dpint)v13 * p10;
  t += (dpint)v14 * p9;
  t += (dpint)v15 * p8;
  t += (dpint)v16 * p7;
  t += (dpint)v17 * p6;
  c[5] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[7] * b[17];
  t += (dpint)a[8] * b[16];
  t += (dpint)a[9] * b[15];
  t += (dpint)a[10] * b[14];
  t += (dpint)a[11] * b[13];
  t += (dpint)a[12] * b[12];
  t += (dpint)a[13] * b[11];
  t += (dpint)a[14] * b[10];
  t += (dpint)a[15] * b[9];
  t += (dpint)a[16] * b[8];
  t += (dpint)a[17] * b[7];
  t += (dpint)v7 * p17;
  t += (dpint)v8 * p16;
  t += (dpint)v9 * p15;
  t += (dpint)v10 * p14;
  t += (dpint)v11 * p13;
  t += (dpint)v12 * p12;
  t += (dpint)v13 * p11;
  t += (dpint)v14 * p10;
  t += (dpint)v15 * p9;
  t += (dpint)v16 * p8;
  t += (dpint)v17 * p7;
  c[6] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[8] * b[17];
  t += (dpint)a[9] * b[16];
  t += (dpint)a[10] * b[15];
  t += (dpint)a[11] * b[14];
  t += (dpint)a[12] * b[13];
  t += (dpint)a[13] * b[12];
  t += (dpint)a[14] * b[11];
  t += (dpint)a[15] * b[10];
  t += (dpint)a[16] * b[9];
  t += (dpint)a[17] * b[8];
  t += (dpint)v8 * p17;
  t += (dpint)v9 * p16;
  t += (dpint)v10 * p15;
  t += (dpint)v11 * p14;
  t += (dpint)v12 * p13;
  t += (dpint)v13 * p12;
  t += (dpint)v14 * p11;
  t += (dpint)v15 * p10;
  t += (dpint)v16 * p9;
  t += (dpint)v17 * p8;
  c[7] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[9] * b[17];
  t += (dpint)a[10] * b[16];
  t += (dpint)a[11] * b[15];
  t += (dpint)a[12] * b[14];
  t += (dpint)a[13] * b[13];
  t += (dpint)a[14] * b[12];
  t += (dpint)a[15] * b[11];
  t += (dpint)a[16] * b[10];
  t += (dpint)a[17] * b[9];
  t += (dpint)v9 * p17;
  t += (dpint)v10 * p16;
  t += (dpint)v11 * p15;
  t += (dpint)v12 * p14;
  t += (dpint)v13 * p13;
  t += (dpint)v14 * p12;
  t += (dpint)v15 * p11;
  t += (dpint)v16 * p10;
  t += (dpint)v17 * p9;
  c[8] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[10] * b[17];
  t += (dpint)a[11] * b[16];
  t += (dpint)a[12] * b[15];
  t += (dpint)a[13] * b[14];
  t += (dpint)a[14] * b[13];
  t += (dpint)a[15] * b[12];
  t += (dpint)a[16] * b[11];
  t += (dpint)a[17] * b[10];
  t += (dpint)v10 * p17;
  t += (dpint)v11 * p16;
  t += (dpint)v12 * p15;
  t += (dpint)v13 * p14;
  t += (dpint)v14 * p13;
  t += (dpint)v15 * p12;
  t += (dpint)v16 * p11;
  t += (dpint)v17 * p10;
  c[9] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[11] * b[17];
  t += (dpint)a[12] * b[16];
  t += (dpint)a[13] * b[15];
  t += (dpint)a[14] * b[14];
  t += (dpint)a[15] * b[13];
  t += (dpint)a[16] * b[12];
  t += (dpint)a[17] * b[11];
  t += (dpint)v11 * p17;
  t += (dpint)v12 * p16;
  t += (dpint)v13 * p15;
  t += (dpint)v14 * p14;
  t += (dpint)v15 * p13;
  t += (dpint)v16 * p12;
  t += (dpint)v17 * p11;
  c[10] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[12] * b[17];
  t += (dpint)a[13] * b[16];
  t += (dpint)a[14] * b[15];
  t += (dpint)a[15] * b[14];
  t += (dpint)a[16] * b[13];
  t += (dpint)a[17] * b[12];
  t += (dpint)v12 * p17;
  t += (dpint)v13 * p16;
  t += (dpint)v14 * p15;
  t += (dpint)v15 * p14;
  t += (dpint)v16 * p13;
  t += (dpint)v17 * p12;
  c[11] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[13] * b[17];
  t += (dpint)a[14] * b[16];
  t += (dpint)a[15] * b[15];
  t += (dpint)a[16] * b[14];
  t += (dpint)a[17] * b[13];
  t += (dpint)v13 * p17;
  t += (dpint)v14 * p16;
  t += (dpint)v15 * p15;
  t += (dpint)v16 * p14;
  t += (dpint)v17 * p13;
  c[12] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[14] * b[17];
  t += (dpint)a[15] * b[16];
  t += (dpint)a[16] * b[15];
  t += (dpint)a[17] * b[14];
  t += (dpint)v14 * p17;
  t += (dpint)v15 * p16;
  t += (dpint)v16 * p15;
  t += (dpint)v17 * p14;
  c[13] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[15] * b[17];
  t += (dpint)a[16] * b[16];
  t += (dpint)a[17] * b[15];
  t += (dpint)v15 * p17;
  t += (dpint)v16 * p16;
  t += (dpint)v17 * p15;
  c[14] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[16] * b[17];
  t += (dpint)a[17] * b[16];
  t += (dpint)v16 * p17;
  t += (dpint)v17 * p16;
  c[15] = (spint)((uspint)t & mask);
  t >>= 28;
  t += (dpint)a[17] * b[17];
  t += (dpint)v17 * p17;
  c[16] = (spint)((uspint)t & mask);
  t >>= 28;
  c[17] = (spint)t;
}

// Modular squaring, c=a*a  mod 2p
void __attribute__((noinline)) modsqr(spint *a, spint *c) {
  spint v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15,
      v16, v17;
  dpint tot, t = 0;
  spint p5 = 0xada6e20;
  spint p6 = 0xe994c68;
  spint p7 = 0x781974c;
  spint p8 = 0xaf0a29a;
  spint p9 = 0xdea65f;
  spint p10 = 0xac5904a;
  spint p11 = 0x7d01fe3;
  spint p12 = 0xe632650;
  spint p13 = 0x9202bdb;
  spint p14 = 0x36936e7;
  spint p15 = 0x8c15b00;
  spint p16 = 0x8869bc6;
  spint p17 = 0x255946a;
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
  v3 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[4];
  tot += (dpint)a[1] * a[3];
  tot *= 2;
  tot += (dpint)a[2] * a[2];
  t += tot;
  v4 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[5];
  tot += (dpint)a[1] * a[4];
  tot += (dpint)a[2] * a[3];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p5;
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
  v13 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[14];
  tot += (dpint)a[1] * a[13];
  tot += (dpint)a[2] * a[12];
  tot += (dpint)a[3] * a[11];
  tot += (dpint)a[4] * a[10];
  tot += (dpint)a[5] * a[9];
  tot += (dpint)a[6] * a[8];
  tot *= 2;
  tot += (dpint)a[7] * a[7];
  t += tot;
  t += (dpint)v0 * p14;
  t += (dpint)v1 * p13;
  t += (dpint)v2 * p12;
  t += (dpint)v3 * p11;
  t += (dpint)v4 * p10;
  t += (dpint)v5 * p9;
  t += (dpint)v6 * p8;
  t += (dpint)v7 * p7;
  t += (dpint)v8 * p6;
  t += (dpint)v9 * p5;
  v14 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[15];
  tot += (dpint)a[1] * a[14];
  tot += (dpint)a[2] * a[13];
  tot += (dpint)a[3] * a[12];
  tot += (dpint)a[4] * a[11];
  tot += (dpint)a[5] * a[10];
  tot += (dpint)a[6] * a[9];
  tot += (dpint)a[7] * a[8];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p15;
  t += (dpint)v1 * p14;
  t += (dpint)v2 * p13;
  t += (dpint)v3 * p12;
  t += (dpint)v4 * p11;
  t += (dpint)v5 * p10;
  t += (dpint)v6 * p9;
  t += (dpint)v7 * p8;
  t += (dpint)v8 * p7;
  t += (dpint)v9 * p6;
  t += (dpint)v10 * p5;
  v15 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[16];
  tot += (dpint)a[1] * a[15];
  tot += (dpint)a[2] * a[14];
  tot += (dpint)a[3] * a[13];
  tot += (dpint)a[4] * a[12];
  tot += (dpint)a[5] * a[11];
  tot += (dpint)a[6] * a[10];
  tot += (dpint)a[7] * a[9];
  tot *= 2;
  tot += (dpint)a[8] * a[8];
  t += tot;
  t += (dpint)v0 * p16;
  t += (dpint)v1 * p15;
  t += (dpint)v2 * p14;
  t += (dpint)v3 * p13;
  t += (dpint)v4 * p12;
  t += (dpint)v5 * p11;
  t += (dpint)v6 * p10;
  t += (dpint)v7 * p9;
  t += (dpint)v8 * p8;
  t += (dpint)v9 * p7;
  t += (dpint)v10 * p6;
  t += (dpint)v11 * p5;
  v16 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[0] * a[17];
  tot += (dpint)a[1] * a[16];
  tot += (dpint)a[2] * a[15];
  tot += (dpint)a[3] * a[14];
  tot += (dpint)a[4] * a[13];
  tot += (dpint)a[5] * a[12];
  tot += (dpint)a[6] * a[11];
  tot += (dpint)a[7] * a[10];
  tot += (dpint)a[8] * a[9];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p17;
  t += (dpint)v1 * p16;
  t += (dpint)v2 * p15;
  t += (dpint)v3 * p14;
  t += (dpint)v4 * p13;
  t += (dpint)v5 * p12;
  t += (dpint)v6 * p11;
  t += (dpint)v7 * p10;
  t += (dpint)v8 * p9;
  t += (dpint)v9 * p8;
  t += (dpint)v10 * p7;
  t += (dpint)v11 * p6;
  t += (dpint)v12 * p5;
  v17 = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[1] * a[17];
  tot += (dpint)a[2] * a[16];
  tot += (dpint)a[3] * a[15];
  tot += (dpint)a[4] * a[14];
  tot += (dpint)a[5] * a[13];
  tot += (dpint)a[6] * a[12];
  tot += (dpint)a[7] * a[11];
  tot += (dpint)a[8] * a[10];
  tot *= 2;
  tot += (dpint)a[9] * a[9];
  t += tot;
  t += (dpint)v1 * p17;
  t += (dpint)v2 * p16;
  t += (dpint)v3 * p15;
  t += (dpint)v4 * p14;
  t += (dpint)v5 * p13;
  t += (dpint)v6 * p12;
  t += (dpint)v7 * p11;
  t += (dpint)v8 * p10;
  t += (dpint)v9 * p9;
  t += (dpint)v10 * p8;
  t += (dpint)v11 * p7;
  t += (dpint)v12 * p6;
  t += (dpint)v13 * p5;
  c[0] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[2] * a[17];
  tot += (dpint)a[3] * a[16];
  tot += (dpint)a[4] * a[15];
  tot += (dpint)a[5] * a[14];
  tot += (dpint)a[6] * a[13];
  tot += (dpint)a[7] * a[12];
  tot += (dpint)a[8] * a[11];
  tot += (dpint)a[9] * a[10];
  tot *= 2;
  t += tot;
  t += (dpint)v2 * p17;
  t += (dpint)v3 * p16;
  t += (dpint)v4 * p15;
  t += (dpint)v5 * p14;
  t += (dpint)v6 * p13;
  t += (dpint)v7 * p12;
  t += (dpint)v8 * p11;
  t += (dpint)v9 * p10;
  t += (dpint)v10 * p9;
  t += (dpint)v11 * p8;
  t += (dpint)v12 * p7;
  t += (dpint)v13 * p6;
  t += (dpint)v14 * p5;
  c[1] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[3] * a[17];
  tot += (dpint)a[4] * a[16];
  tot += (dpint)a[5] * a[15];
  tot += (dpint)a[6] * a[14];
  tot += (dpint)a[7] * a[13];
  tot += (dpint)a[8] * a[12];
  tot += (dpint)a[9] * a[11];
  tot *= 2;
  tot += (dpint)a[10] * a[10];
  t += tot;
  t += (dpint)v3 * p17;
  t += (dpint)v4 * p16;
  t += (dpint)v5 * p15;
  t += (dpint)v6 * p14;
  t += (dpint)v7 * p13;
  t += (dpint)v8 * p12;
  t += (dpint)v9 * p11;
  t += (dpint)v10 * p10;
  t += (dpint)v11 * p9;
  t += (dpint)v12 * p8;
  t += (dpint)v13 * p7;
  t += (dpint)v14 * p6;
  t += (dpint)v15 * p5;
  c[2] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[4] * a[17];
  tot += (dpint)a[5] * a[16];
  tot += (dpint)a[6] * a[15];
  tot += (dpint)a[7] * a[14];
  tot += (dpint)a[8] * a[13];
  tot += (dpint)a[9] * a[12];
  tot += (dpint)a[10] * a[11];
  tot *= 2;
  t += tot;
  t += (dpint)v4 * p17;
  t += (dpint)v5 * p16;
  t += (dpint)v6 * p15;
  t += (dpint)v7 * p14;
  t += (dpint)v8 * p13;
  t += (dpint)v9 * p12;
  t += (dpint)v10 * p11;
  t += (dpint)v11 * p10;
  t += (dpint)v12 * p9;
  t += (dpint)v13 * p8;
  t += (dpint)v14 * p7;
  t += (dpint)v15 * p6;
  t += (dpint)v16 * p5;
  c[3] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[5] * a[17];
  tot += (dpint)a[6] * a[16];
  tot += (dpint)a[7] * a[15];
  tot += (dpint)a[8] * a[14];
  tot += (dpint)a[9] * a[13];
  tot += (dpint)a[10] * a[12];
  tot *= 2;
  tot += (dpint)a[11] * a[11];
  t += tot;
  t += (dpint)v5 * p17;
  t += (dpint)v6 * p16;
  t += (dpint)v7 * p15;
  t += (dpint)v8 * p14;
  t += (dpint)v9 * p13;
  t += (dpint)v10 * p12;
  t += (dpint)v11 * p11;
  t += (dpint)v12 * p10;
  t += (dpint)v13 * p9;
  t += (dpint)v14 * p8;
  t += (dpint)v15 * p7;
  t += (dpint)v16 * p6;
  t += (dpint)v17 * p5;
  c[4] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[6] * a[17];
  tot += (dpint)a[7] * a[16];
  tot += (dpint)a[8] * a[15];
  tot += (dpint)a[9] * a[14];
  tot += (dpint)a[10] * a[13];
  tot += (dpint)a[11] * a[12];
  tot *= 2;
  t += tot;
  t += (dpint)v6 * p17;
  t += (dpint)v7 * p16;
  t += (dpint)v8 * p15;
  t += (dpint)v9 * p14;
  t += (dpint)v10 * p13;
  t += (dpint)v11 * p12;
  t += (dpint)v12 * p11;
  t += (dpint)v13 * p10;
  t += (dpint)v14 * p9;
  t += (dpint)v15 * p8;
  t += (dpint)v16 * p7;
  t += (dpint)v17 * p6;
  c[5] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[7] * a[17];
  tot += (dpint)a[8] * a[16];
  tot += (dpint)a[9] * a[15];
  tot += (dpint)a[10] * a[14];
  tot += (dpint)a[11] * a[13];
  tot *= 2;
  tot += (dpint)a[12] * a[12];
  t += tot;
  t += (dpint)v7 * p17;
  t += (dpint)v8 * p16;
  t += (dpint)v9 * p15;
  t += (dpint)v10 * p14;
  t += (dpint)v11 * p13;
  t += (dpint)v12 * p12;
  t += (dpint)v13 * p11;
  t += (dpint)v14 * p10;
  t += (dpint)v15 * p9;
  t += (dpint)v16 * p8;
  t += (dpint)v17 * p7;
  c[6] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[8] * a[17];
  tot += (dpint)a[9] * a[16];
  tot += (dpint)a[10] * a[15];
  tot += (dpint)a[11] * a[14];
  tot += (dpint)a[12] * a[13];
  tot *= 2;
  t += tot;
  t += (dpint)v8 * p17;
  t += (dpint)v9 * p16;
  t += (dpint)v10 * p15;
  t += (dpint)v11 * p14;
  t += (dpint)v12 * p13;
  t += (dpint)v13 * p12;
  t += (dpint)v14 * p11;
  t += (dpint)v15 * p10;
  t += (dpint)v16 * p9;
  t += (dpint)v17 * p8;
  c[7] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[9] * a[17];
  tot += (dpint)a[10] * a[16];
  tot += (dpint)a[11] * a[15];
  tot += (dpint)a[12] * a[14];
  tot *= 2;
  tot += (dpint)a[13] * a[13];
  t += tot;
  t += (dpint)v9 * p17;
  t += (dpint)v10 * p16;
  t += (dpint)v11 * p15;
  t += (dpint)v12 * p14;
  t += (dpint)v13 * p13;
  t += (dpint)v14 * p12;
  t += (dpint)v15 * p11;
  t += (dpint)v16 * p10;
  t += (dpint)v17 * p9;
  c[8] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[10] * a[17];
  tot += (dpint)a[11] * a[16];
  tot += (dpint)a[12] * a[15];
  tot += (dpint)a[13] * a[14];
  tot *= 2;
  t += tot;
  t += (dpint)v10 * p17;
  t += (dpint)v11 * p16;
  t += (dpint)v12 * p15;
  t += (dpint)v13 * p14;
  t += (dpint)v14 * p13;
  t += (dpint)v15 * p12;
  t += (dpint)v16 * p11;
  t += (dpint)v17 * p10;
  c[9] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[11] * a[17];
  tot += (dpint)a[12] * a[16];
  tot += (dpint)a[13] * a[15];
  tot *= 2;
  tot += (dpint)a[14] * a[14];
  t += tot;
  t += (dpint)v11 * p17;
  t += (dpint)v12 * p16;
  t += (dpint)v13 * p15;
  t += (dpint)v14 * p14;
  t += (dpint)v15 * p13;
  t += (dpint)v16 * p12;
  t += (dpint)v17 * p11;
  c[10] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[12] * a[17];
  tot += (dpint)a[13] * a[16];
  tot += (dpint)a[14] * a[15];
  tot *= 2;
  t += tot;
  t += (dpint)v12 * p17;
  t += (dpint)v13 * p16;
  t += (dpint)v14 * p15;
  t += (dpint)v15 * p14;
  t += (dpint)v16 * p13;
  t += (dpint)v17 * p12;
  c[11] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[13] * a[17];
  tot += (dpint)a[14] * a[16];
  tot *= 2;
  tot += (dpint)a[15] * a[15];
  t += tot;
  t += (dpint)v13 * p17;
  t += (dpint)v14 * p16;
  t += (dpint)v15 * p15;
  t += (dpint)v16 * p14;
  t += (dpint)v17 * p13;
  c[12] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[14] * a[17];
  tot += (dpint)a[15] * a[16];
  tot *= 2;
  t += tot;
  t += (dpint)v14 * p17;
  t += (dpint)v15 * p16;
  t += (dpint)v16 * p15;
  t += (dpint)v17 * p14;
  c[13] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[15] * a[17];
  tot *= 2;
  tot += (dpint)a[16] * a[16];
  t += tot;
  t += (dpint)v15 * p17;
  t += (dpint)v16 * p16;
  t += (dpint)v17 * p15;
  c[14] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[16] * a[17];
  tot *= 2;
  t += tot;
  t += (dpint)v16 * p17;
  t += (dpint)v17 * p16;
  c[15] = (spint)((uspint)t & mask);
  t >>= 28;
  tot = (dpint)a[17] * a[17];
  t += tot;
  t += (dpint)v17 * p17;
  c[16] = (spint)((uspint)t & mask);
  t >>= 28;
  c[17] = (spint)t;
}

// copy
void __attribute__((always_inline)) modcpy(spint *a, spint *c) {
  int i;
  for (i = 0; i < 18; i++) {
    c[i] = a[i];
  }
}

// Calculate progenitor
void modpro(spint *w, spint *z) {
  int i;
  spint x[18];
  spint t0[18], t1[18], t2[18], t3[18], t4[18], t5[18], t6[18], t7[18], t8[18],
      t9[18], t10[18], t11[18], t12[18], t13[18], t14[18], t15[18], t16[18],
      t17[18], t18[18], t19[18], t20[18], t21[18], t22[18], t23[18], t24[18],
      t25[18], t26[18], t27[18], t28[18], t29[18], t30[18];
  modcpy(w, x);
  modsqr(x, t16);
  modsqr(t16, t19);
  modmul(t16, t19, z);
  modmul(t16, z, t6);
  modmul(x, t6, t1);
  modmul(x, t1, t20);
  modmul(t16, t20, t8);
  modmul(z, t1, t9);
  modmul(z, t8, t13);
  modmul(t16, t13, t0);
  modmul(t1, t0, t17);
  modmul(z, t17, t26);
  modmul(t16, t26, t5);
  modmul(t20, t5, t12);
  modmul(t6, t12, t1);
  modmul(t8, t1, t27);
  modmul(t16, t27, t3);
  modmul(t19, t3, t22);
  modmul(t19, t22, t10);
  modmul(t19, t10, t11);
  modmul(t16, t11, t7);
  modmul(z, t7, t15);
  modmul(t20, t15, t4);
  modmul(t16, t4, t29);
  modmul(t6, t29, t2);
  modmul(t16, t2, t14);
  modmul(t13, t2, z);
  modmul(t19, z, t18);
  modmul(t13, t18, t21);
  modmul(t19, t21, t19);
  modmul(t16, t19, t24);
  modmul(t13, t24, t25);
  modmul(t16, t25, t21);
  modmul(t13, t25, t13);
  modmul(t8, t13, t8);
  modmul(t20, t8, t28);
  modmul(t20, t28, t20);
  modmul(t6, t20, t23);
  modmul(t16, t23, t6);
  modmul(t16, t6, t16);
  modmul(t0, t16, t0);
  modmul(t27, t23, t30);
  for (i = 0; i < 7; i++) {
    modsqr(t30, t30);
  }
  modmul(t29, t30, t29);
  for (i = 0; i < 11; i++) {
    modsqr(t29, t29);
  }
  modmul(t28, t29, t28);
  for (i = 0; i < 10; i++) {
    modsqr(t28, t28);
  }
  modmul(t27, t28, t27);
  for (i = 0; i < 9; i++) {
    modsqr(t27, t27);
  }
  modmul(t24, t27, t27);
  for (i = 0; i < 7; i++) {
    modsqr(t27, t27);
  }
  modmul(t4, t27, t27);
  for (i = 0; i < 7; i++) {
    modsqr(t27, t27);
  }
  modmul(t26, t27, t26);
  for (i = 0; i < 13; i++) {
    modsqr(t26, t26);
  }
  modmul(t25, t26, t25);
  modsqr(t25, t25);
  modmul(x, t25, t25);
  for (i = 0; i < 17; i++) {
    modsqr(t25, t25);
  }
  modmul(t2, t25, t25);
  for (i = 0; i < 10; i++) {
    modsqr(t25, t25);
  }
  modmul(t24, t25, t24);
  for (i = 0; i < 9; i++) {
    modsqr(t24, t24);
  }
  modmul(t23, t24, t23);
  for (i = 0; i < 7; i++) {
    modsqr(t23, t23);
  }
  modmul(t22, t23, t22);
  for (i = 0; i < 15; i++) {
    modsqr(t22, t22);
  }
  modmul(t21, t22, t21);
  for (i = 0; i < 9; i++) {
    modsqr(t21, t21);
  }
  modmul(t20, t21, t20);
  for (i = 0; i < 9; i++) {
    modsqr(t20, t20);
  }
  modmul(t4, t20, t20);
  for (i = 0; i < 10; i++) {
    modsqr(t20, t20);
  }
  modmul(t19, t20, t19);
  for (i = 0; i < 9; i++) {
    modsqr(t19, t19);
  }
  modmul(t18, t19, t18);
  for (i = 0; i < 5; i++) {
    modsqr(t18, t18);
  }
  modmul(t17, t18, t17);
  for (i = 0; i < 15; i++) {
    modsqr(t17, t17);
  }
  modmul(t0, t17, t17);
  for (i = 0; i < 11; i++) {
    modsqr(t17, t17);
  }
  modmul(t16, t17, t16);
  for (i = 0; i < 10; i++) {
    modsqr(t16, t16);
  }
  modmul(t15, t16, t15);
  for (i = 0; i < 11; i++) {
    modsqr(t15, t15);
  }
  modmul(t5, t15, t15);
  for (i = 0; i < 12; i++) {
    modsqr(t15, t15);
  }
  modmul(t14, t15, t14);
  for (i = 0; i < 8; i++) {
    modsqr(t14, t14);
  }
  modmul(t7, t14, t14);
  for (i = 0; i < 10; i++) {
    modsqr(t14, t14);
  }
  modmul(t13, t14, t13);
  for (i = 0; i < 7; i++) {
    modsqr(t13, t13);
  }
  modmul(t12, t13, t12);
  for (i = 0; i < 11; i++) {
    modsqr(t12, t12);
  }
  modmul(t11, t12, t11);
  for (i = 0; i < 8; i++) {
    modsqr(t11, t11);
  }
  modmul(t10, t11, t10);
  for (i = 0; i < 6; i++) {
    modsqr(t10, t10);
  }
  modmul(t9, t10, t9);
  for (i = 0; i < 14; i++) {
    modsqr(t9, t9);
  }
  modmul(t8, t9, t8);
  for (i = 0; i < 7; i++) {
    modsqr(t8, t8);
  }
  modmul(t7, t8, t7);
  for (i = 0; i < 10; i++) {
    modsqr(t7, t7);
  }
  modmul(t6, t7, t6);
  for (i = 0; i < 6; i++) {
    modsqr(t6, t6);
  }
  modmul(t5, t6, t5);
  for (i = 0; i < 9; i++) {
    modsqr(t5, t5);
  }
  modmul(t4, t5, t4);
  for (i = 0; i < 8; i++) {
    modsqr(t4, t4);
  }
  modmul(t3, t4, t3);
  for (i = 0; i < 8; i++) {
    modsqr(t3, t3);
  }
  modmul(t2, t3, t2);
  for (i = 0; i < 8; i++) {
    modsqr(t2, t2);
  }
  modmul(t1, t2, t1);
  for (i = 0; i < 12; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 8; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t0);
  for (i = 0; i < 7; i++) {
    modsqr(t0, t0);
  }
  modmul(z, t0, z);
}

// calculate inverse, provide progenitor h if available
void modinv(spint *x, spint *h, spint *z) {
  int i;
  spint s[18], t[18];
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
  spint c[18] = {0x8a0c754, 0xa546e4e, 0x943e89e, 0xcb993b5, 0xbf555c8,
                 0x8dda851, 0x328de11, 0xc717daa, 0xb090199, 0x307cf4e,
                 0xd096114, 0x8537d92, 0x2dce4ee, 0x96f0b3d, 0x5bc94fd,
                 0x73aacac, 0x943ddd4, 0x16de108};
  modmul(m, c, n);
}

// Convert n back to normal form, m=redc(n)
void redc(spint *n, spint *m) {
  spint c[18];
  c[0] = 1;
  for (int i = 1; i < 18; i++)
    c[i] = 0;
  modmul(n, c, m);
  modfsb(m);
}

// is unity?
int modis1(spint *a) {
  spint c[18];
  sspint c0, d = 0;
  redc(a, c);
  for (int i = 1; i < 18; i++) {
    d |= c[i];
  }
  c0 = (sspint)c[0];
  return (1 & ((d - 1) >> 28) & (((c0 ^ 1) - 1) >> 28));
}

// is zero?
int modis0(spint *a) {
  sspint d = 0;
  for (int i = 0; i < 18; i++) {
    d |= a[i];
  }
  return (1 & ((d - 1) >> 28));
}

// set to zero
void modzer(spint *a) {
  for (int i = 0; i < 18; i++)
    a[i] = 0;
}

// set to one
void modone(spint *a) {
  a[0] = 1;
  for (int i = 1; i < 18; i++)
    a[i] = 0;
  nres(a, a);
}

// Test for quadratic residue
int modqr(spint *h, spint *x) {
  int i;
  spint r[18];
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
  spint s[18], y[18];
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
  a[18 - 1] = ((a[18 - 1] << n)) | (a[18 - 2] >> (28 - n));
  for (int i = 18 - 2; i > 0; i--) {
    a[i] = ((a[i] << n) & 0xfffffff) | (a[i - 1] >> (28 - n));
  }
  a[0] = (a[0] << n) & 0xfffffff;
}

// shift right by less than a word. Return shifted out part
int modshr(int n, spint *a) {
  spint r = a[0] & (((spint)1 << n) - 1);
  for (int i = 0; i < 18 - 1; i++) {
    a[i] = (a[i] >> n) | ((a[i + 1] << (28 - n)) & 0xfffffff);
  }
  a[18 - 1] = a[18 - 1] >> n;
  return r;
}

/* API functions calling generated code */

const digit_t p[NWORDS_ORDER] =  { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xA6E1FFFF, 0x994C68AD, 0x781974CE, 0xFAF0A29A, 0x4A0DEA65, 0xFE3AC590, 0x26507D01, 0x02BDBE63, 0x6936E792, 0x8C15B003, 0xA8869BC6, 0x00255946 };

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
    for (int i = 0; i < 63; i++) {
        ((char *) out)[i] = x[0] & 0xff;
        modshr(8, x);
    }
}

void fp_from_digit_array(digit_t* out, const digit_t* a) {
    for (int i = 0; i < NWORDS_FIELD; i++) {
        out[i] = 0;
    }
    for (int i = 63 - 1; i >= 0; i--) {
        modshl(8, out);
        out[0] += (digit_t)((unsigned char *) a)[i];
    }
}

#endif
