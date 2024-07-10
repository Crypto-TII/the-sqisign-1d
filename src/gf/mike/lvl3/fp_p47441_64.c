#ifdef RADIX_64

#include <stdint.h>
#include <stdio.h>

#include <stdbool.h>
#include <fp.h>

#define uspint uint64_t
#define sspint int64_t
#define spint uint64_t
#define dpint __uint128_t

// propagate carries - return sign
static __attribute__((always_inline)) sspint prop(spint *n) {
  spint d, mask = ((spint)1 << 55) - 1;
  sspint carry = (sspint)n[0] >> 55;
  n[0] &= mask;
  for (int i = 1; i < 6; i++) {
    d = n[i] + carry;
    n[i] = d & mask;
    carry = (sspint)d >> 55;
  }
  n[6] += carry;
  return ((sspint)n[6] >> 63);
}

// propagate carries and add p if negative, propagate carries again
static __attribute__((always_inline)) void flatten(spint *n) {
  spint q = ((spint)1 << 55);
  sspint carry = prop(n);
  n[0] -= 1 & carry;
  n[1] += (1 * (spint)0x69840000000000) & carry;
  n[2] += (1 * (spint)0x24d5ba91a13185) & carry;
  n[3] += (1 * (spint)0x575ba8e3917b3) & carry;
  n[4] += (1 * (spint)0x10ad665bc2e0a9) & carry;
  n[5] += (1 * (spint)0x3458d5cc0948ba) & carry;
  n[6] += (1 * (spint)0xf7dbbbaac21c) & carry;
  prop(n);
}

// Montgomery final subtract
void __attribute__((always_inline)) modfsb(spint *n) {
  spint q = ((spint)1 << 55);
  n[0] += 1;
  n[1] -= 1 * (spint)0x69840000000000;
  n[2] -= 1 * (spint)0x24d5ba91a13185;
  n[3] -= 1 * (spint)0x575ba8e3917b3;
  n[4] -= 1 * (spint)0x10ad665bc2e0a9;
  n[5] -= 1 * (spint)0x3458d5cc0948ba;
  n[6] -= 1 * (spint)0xf7dbbbaac21c;
  flatten(n);
}

// Modular addition - reduce less than 2p
void __attribute__((always_inline)) modadd(spint *a, spint *b, spint *n) {
  spint q = ((spint)1 << 55);
  sspint carry;
  n[0] = a[0] + b[0];
  n[1] = a[1] + b[1];
  n[2] = a[2] + b[2];
  n[3] = a[3] + b[3];
  n[4] = a[4] + b[4];
  n[5] = a[5] + b[5];
  n[6] = a[6] + b[6];
  n[0] += 2;
  n[1] -= 2 * (spint)0x69840000000000;
  n[2] -= 2 * (spint)0x24d5ba91a13185;
  n[3] -= 2 * (spint)0x575ba8e3917b3;
  n[4] -= 2 * (spint)0x10ad665bc2e0a9;
  n[5] -= 2 * (spint)0x3458d5cc0948ba;
  n[6] -= 2 * (spint)0xf7dbbbaac21c;
  carry = prop(n);
  n[0] -= 2 & carry;
  n[1] += (2 * (spint)0x69840000000000) & carry;
  n[2] += (2 * (spint)0x24d5ba91a13185) & carry;
  n[3] += (2 * (spint)0x575ba8e3917b3) & carry;
  n[4] += (2 * (spint)0x10ad665bc2e0a9) & carry;
  n[5] += (2 * (spint)0x3458d5cc0948ba) & carry;
  n[6] += (2 * (spint)0xf7dbbbaac21c) & carry;
  prop(n);
}

// Modular subtraction - reduce less than 2p
void __attribute__((always_inline)) modsub(spint *a, spint *b, spint *n) {
  spint q = ((spint)1 << 55);
  sspint carry;
  n[0] = a[0] - b[0];
  n[1] = a[1] - b[1];
  n[2] = a[2] - b[2];
  n[3] = a[3] - b[3];
  n[4] = a[4] - b[4];
  n[5] = a[5] - b[5];
  n[6] = a[6] - b[6];
  carry = prop(n);
  n[0] -= 2 & carry;
  n[1] += (2 * (spint)0x69840000000000) & carry;
  n[2] += (2 * (spint)0x24d5ba91a13185) & carry;
  n[3] += (2 * (spint)0x575ba8e3917b3) & carry;
  n[4] += (2 * (spint)0x10ad665bc2e0a9) & carry;
  n[5] += (2 * (spint)0x3458d5cc0948ba) & carry;
  n[6] += (2 * (spint)0xf7dbbbaac21c) & carry;
  prop(n);
}

// Overflow limit   = 340282366920938463463374607431768211456
// maximum possible = 11295303651797634038207769553587024
// Modular multiplication, c=a*b mod 2p
void __attribute__((always_inline)) modmul(spint *a, spint *b, spint *c) {
  spint v0, v1, v2, v3, v4, v5, v6;
  dpint t = 0;
  spint p1 = 0x69840000000000;
  spint p2 = 0x24d5ba91a13185;
  spint p3 = 0x575ba8e3917b3;
  spint p4 = 0x10ad665bc2e0a9;
  spint p5 = 0x3458d5cc0948ba;
  spint p6 = 0xf7dbbbaac21c;
  spint s, q = ((spint)1 << 55); // q is unsaturated radix
  spint mask = q - 1;
  t += (dpint)a[0] * b[0];
  v0 = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[1];
  t += (dpint)a[1] * b[0];
  t += (dpint)v0 * p1;
  v1 = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[2];
  t += (dpint)a[1] * b[1];
  t += (dpint)a[2] * b[0];
  t += (dpint)v0 * p2;
  t += (dpint)v1 * p1;
  v2 = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[3];
  t += (dpint)a[1] * b[2];
  t += (dpint)a[2] * b[1];
  t += (dpint)a[3] * b[0];
  t += (dpint)v0 * p3;
  t += (dpint)v1 * p2;
  t += (dpint)v2 * p1;
  v3 = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[4];
  t += (dpint)a[1] * b[3];
  t += (dpint)a[2] * b[2];
  t += (dpint)a[3] * b[1];
  t += (dpint)a[4] * b[0];
  t += (dpint)v0 * p4;
  t += (dpint)v1 * p3;
  t += (dpint)v2 * p2;
  t += (dpint)v3 * p1;
  v4 = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[5];
  t += (dpint)a[1] * b[4];
  t += (dpint)a[2] * b[3];
  t += (dpint)a[3] * b[2];
  t += (dpint)a[4] * b[1];
  t += (dpint)a[5] * b[0];
  t += (dpint)v0 * p5;
  t += (dpint)v1 * p4;
  t += (dpint)v2 * p3;
  t += (dpint)v3 * p2;
  t += (dpint)v4 * p1;
  v5 = (spint)((uspint)t & mask);
  t >>= 55;
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
  t += (dpint)v4 * p2;
  t += (dpint)v5 * p1;
  v6 = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[1] * b[6];
  t += (dpint)a[2] * b[5];
  t += (dpint)a[3] * b[4];
  t += (dpint)a[4] * b[3];
  t += (dpint)a[5] * b[2];
  t += (dpint)a[6] * b[1];
  t += (dpint)v1 * p6;
  t += (dpint)v2 * p5;
  t += (dpint)v3 * p4;
  t += (dpint)v4 * p3;
  t += (dpint)v5 * p2;
  t += (dpint)v6 * p1;
  c[0] = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[2] * b[6];
  t += (dpint)a[3] * b[5];
  t += (dpint)a[4] * b[4];
  t += (dpint)a[5] * b[3];
  t += (dpint)a[6] * b[2];
  t += (dpint)v2 * p6;
  t += (dpint)v3 * p5;
  t += (dpint)v4 * p4;
  t += (dpint)v5 * p3;
  t += (dpint)v6 * p2;
  c[1] = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[3] * b[6];
  t += (dpint)a[4] * b[5];
  t += (dpint)a[5] * b[4];
  t += (dpint)a[6] * b[3];
  t += (dpint)v3 * p6;
  t += (dpint)v4 * p5;
  t += (dpint)v5 * p4;
  t += (dpint)v6 * p3;
  c[2] = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[4] * b[6];
  t += (dpint)a[5] * b[5];
  t += (dpint)a[6] * b[4];
  t += (dpint)v4 * p6;
  t += (dpint)v5 * p5;
  t += (dpint)v6 * p4;
  c[3] = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[5] * b[6];
  t += (dpint)a[6] * b[5];
  t += (dpint)v5 * p6;
  t += (dpint)v6 * p5;
  c[4] = (spint)((uspint)t & mask);
  t >>= 55;
  t += (dpint)a[6] * b[6];
  t += (dpint)v6 * p6;
  c[5] = (spint)((uspint)t & mask);
  t >>= 55;
  c[6] = (spint)t;
}

// Modular squaring, c=a*a  mod 2p
void __attribute__((always_inline)) modsqr(spint *a, spint *c) {
  spint v0, v1, v2, v3, v4, v5, v6;
  dpint tot, t = 0;
  spint p1 = 0x69840000000000;
  spint p2 = 0x24d5ba91a13185;
  spint p3 = 0x575ba8e3917b3;
  spint p4 = 0x10ad665bc2e0a9;
  spint p5 = 0x3458d5cc0948ba;
  spint p6 = 0xf7dbbbaac21c;
  spint s, q = ((spint)1 << 55); // q is unsaturated radix
  spint mask = q - 1;
  tot = (dpint)a[0] * a[0];
  t = tot;
  v0 = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[0] * a[1];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p1;
  v1 = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[0] * a[2];
  tot *= 2;
  tot += (dpint)a[1] * a[1];
  t += tot;
  t += (dpint)v0 * p2;
  t += (dpint)v1 * p1;
  v2 = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[0] * a[3];
  tot += (dpint)a[1] * a[2];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p3;
  t += (dpint)v1 * p2;
  t += (dpint)v2 * p1;
  v3 = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[0] * a[4];
  tot += (dpint)a[1] * a[3];
  tot *= 2;
  tot += (dpint)a[2] * a[2];
  t += tot;
  t += (dpint)v0 * p4;
  t += (dpint)v1 * p3;
  t += (dpint)v2 * p2;
  t += (dpint)v3 * p1;
  v4 = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[0] * a[5];
  tot += (dpint)a[1] * a[4];
  tot += (dpint)a[2] * a[3];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p5;
  t += (dpint)v1 * p4;
  t += (dpint)v2 * p3;
  t += (dpint)v3 * p2;
  t += (dpint)v4 * p1;
  v5 = (spint)((uspint)t & mask);
  t >>= 55;
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
  t += (dpint)v4 * p2;
  t += (dpint)v5 * p1;
  v6 = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[1] * a[6];
  tot += (dpint)a[2] * a[5];
  tot += (dpint)a[3] * a[4];
  tot *= 2;
  t += tot;
  t += (dpint)v1 * p6;
  t += (dpint)v2 * p5;
  t += (dpint)v3 * p4;
  t += (dpint)v4 * p3;
  t += (dpint)v5 * p2;
  t += (dpint)v6 * p1;
  c[0] = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[2] * a[6];
  tot += (dpint)a[3] * a[5];
  tot *= 2;
  tot += (dpint)a[4] * a[4];
  t += tot;
  t += (dpint)v2 * p6;
  t += (dpint)v3 * p5;
  t += (dpint)v4 * p4;
  t += (dpint)v5 * p3;
  t += (dpint)v6 * p2;
  c[1] = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[3] * a[6];
  tot += (dpint)a[4] * a[5];
  tot *= 2;
  t += tot;
  t += (dpint)v3 * p6;
  t += (dpint)v4 * p5;
  t += (dpint)v5 * p4;
  t += (dpint)v6 * p3;
  c[2] = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[4] * a[6];
  tot *= 2;
  tot += (dpint)a[5] * a[5];
  t += tot;
  t += (dpint)v4 * p6;
  t += (dpint)v5 * p5;
  t += (dpint)v6 * p4;
  c[3] = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[5] * a[6];
  tot *= 2;
  t += tot;
  t += (dpint)v5 * p6;
  t += (dpint)v6 * p5;
  c[4] = (spint)((uspint)t & mask);
  t >>= 55;
  tot = (dpint)a[6] * a[6];
  t += tot;
  t += (dpint)v6 * p6;
  c[5] = (spint)((uspint)t & mask);
  t >>= 55;
  c[6] = (spint)t;
}

// copy
void __attribute__((always_inline)) modcpy(spint *a, spint *c) {
  int i;
  for (i = 0; i < 7; i++) {
    c[i] = a[i];
  }
}

// Calculate progenitor
void modpro(spint *w, spint *z) {
  int i;
  spint x[7];
  spint t0[7], t1[7], t2[7], t3[7], t4[7], t5[7], t6[7], t7[7], t8[7], t9[7],
      t10[7], t11[7], t12[7], t13[7], t14[7], t15[7], t16[7], t17[7], t18[7],
      t19[7], t20[7], t21[7], t22[7], t23[7], t24[7];
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
  spint s[7], t[7];
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
  spint c[7] = {0x4fa09a199185fb, 0x75a1ff250d623d, 0x20a5c269fa2c7a,
                0x294e86bbda4a59, 0x47bd98ef8d71e,  0x1a7afdeed9dd4a,
                0xa1841b2d0fcc};
  modmul(m, c, n);
}

// Convert n back to normal form, m=redc(n)
void redc(spint *n, spint *m) {
  spint c[7];
  c[0] = 1;
  for (int i = 1; i < 7; i++)
    c[i] = 0;
  modmul(n, c, m);
  modfsb(m);
}

// is unity?
int modis1(spint *a) {
  spint c[7];
  sspint c0, d = 0;
  redc(a, c);
  for (int i = 1; i < 7; i++) {
    d |= c[i];
  }
  c0 = (sspint)c[0];
  return (1 & ((d - 1) >> 55) & (((c0 ^ 1) - 1) >> 55));
}

// is zero?
int modis0(spint *a) {
  sspint d = 0;
  for (int i = 0; i < 7; i++) {
    d |= a[i];
  }
  return (1 & ((d - 1) >> 55));
}

// set to zero
void modzer(spint *a) {
  for (int i = 0; i < 7; i++)
    a[i] = 0;
}

// set to one
void modone(spint *a) {
  a[0] = 1;
  for (int i = 1; i < 7; i++)
    a[i] = 0;
  nres(a, a);
}

// Test for quadratic residue
int modqr(spint *h, spint *x) {
  int i;
  spint r[7];
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
  spint s[7], y[7];
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
  a[7 - 1] = ((a[7 - 1] << n)) | (a[7 - 2] >> (55 - n));
  for (int i = 7 - 2; i > 0; i--) {
    a[i] = ((a[i] << n) & 0x7fffffffffffff) | (a[i - 1] >> (55 - n));
  }
  a[0] = (a[0] << n) & 0x7fffffffffffff;
}

// shift right by less than a word. Return shifted out part
int modshr(int n, spint *a) {
  spint r = a[0] & (((spint)1 << n) - 1);
  for (int i = 0; i < 7 - 1; i++) {
    a[i] = (a[i] >> n) | ((a[i + 1] << (55 - n)) & 0x7fffffffffffff);
  }
  a[7 - 1] = a[7 - 1] >> n;
  return r;
}

/* API functions calling generated code */

const digit_t p[NWORDS_ORDER] =  { 0xFFFFFFFFFFFFFFFF, 0x4C6174C1FFFFFFFF, 0xC722F669356EA468, 0x65BC2E0A90AEB751, 0xC6AE604A45D10AD6, 0x03DF6EEEAB0871A2 };

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

void fp_exp3div4(digit_t* out, const digit_t* a) {
    modpro(a, out);
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
