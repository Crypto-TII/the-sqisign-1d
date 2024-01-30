#ifdef RADIX_64

#include <stdint.h>
#include <stdio.h>

#include <x86intrin.h>

#include <stdbool.h>
#include <fp.h>

#define uspint uint64_t
#define sspint int64_t
#define spint uint64_t
#define dpint __uint128_t

// propagate carries - return sign
static __attribute__((always_inline)) sspint prop(spint *n) {
  spint d, mask = ((spint)1 << 56) - 1;
  sspint carry = (sspint)n[0] >> 56;
  n[0] &= mask;
  for (int i = 1; i < 8; i++) {
    d = n[i] + carry;
    n[i] = d & mask;
    carry = (sspint)d >> 56;
  }
  n[8] += carry;
  return ((sspint)n[8] >> 63);
}

// propagate carries and add p if negative, propagate carries again
static __attribute__((always_inline)) void flatten(spint *n) {
  spint q = ((spint)1 << 56);
  sspint carry = prop(n);
  n[0] -= 1 & carry;
  n[2] += (1 * (spint)0xada6e200000000) & carry;
  n[3] += (1 * (spint)0x781974ce994c68) & carry;
  n[4] += (1 * (spint)0xdea65faf0a29a) & carry;
  n[5] += (1 * (spint)0x7d01fe3ac5904a) & carry;
  n[6] += (1 * (spint)0x9202bdbe632650) & carry;
  n[7] += (1 * (spint)0x8c15b0036936e7) & carry;
  n[8] += (1 * (spint)0x255946a8869bc6) & carry;
  prop(n);
}

// Montgomery final subtract
void __attribute__((always_inline)) modfsb(spint *n) {
  spint q = ((spint)1 << 56);
  n[0] += 1;
  n[2] -= 1 * (spint)0xada6e200000000;
  n[3] -= 1 * (spint)0x781974ce994c68;
  n[4] -= 1 * (spint)0xdea65faf0a29a;
  n[5] -= 1 * (spint)0x7d01fe3ac5904a;
  n[6] -= 1 * (spint)0x9202bdbe632650;
  n[7] -= 1 * (spint)0x8c15b0036936e7;
  n[8] -= 1 * (spint)0x255946a8869bc6;
  flatten(n);
}

// Modular addition - reduce less than 2p
void __attribute__((always_inline)) modadd(spint *a, spint *b, spint *n) {
  spint q = ((spint)1 << 56);
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
  n[0] += 2;
  n[2] -= 2 * (spint)0xada6e200000000;
  n[3] -= 2 * (spint)0x781974ce994c68;
  n[4] -= 2 * (spint)0xdea65faf0a29a;
  n[5] -= 2 * (spint)0x7d01fe3ac5904a;
  n[6] -= 2 * (spint)0x9202bdbe632650;
  n[7] -= 2 * (spint)0x8c15b0036936e7;
  n[8] -= 2 * (spint)0x255946a8869bc6;
  carry = prop(n);
  n[0] -= 2 & carry;
  n[2] += (2 * (spint)0xada6e200000000) & carry;
  n[3] += (2 * (spint)0x781974ce994c68) & carry;
  n[4] += (2 * (spint)0xdea65faf0a29a) & carry;
  n[5] += (2 * (spint)0x7d01fe3ac5904a) & carry;
  n[6] += (2 * (spint)0x9202bdbe632650) & carry;
  n[7] += (2 * (spint)0x8c15b0036936e7) & carry;
  n[8] += (2 * (spint)0x255946a8869bc6) & carry;
  prop(n);
}

// Modular subtraction - reduce less than 2p
void __attribute__((always_inline)) modsub(spint *a, spint *b, spint *n) {
  spint q = ((spint)1 << 56);
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
  carry = prop(n);
  n[0] -= 2 & carry;
  n[2] += (2 * (spint)0xada6e200000000) & carry;
  n[3] += (2 * (spint)0x781974ce994c68) & carry;
  n[4] += (2 * (spint)0xdea65faf0a29a) & carry;
  n[5] += (2 * (spint)0x7d01fe3ac5904a) & carry;
  n[6] += (2 * (spint)0x9202bdbe632650) & carry;
  n[7] += (2 * (spint)0x8c15b0036936e7) & carry;
  n[8] += (2 * (spint)0x255946a8869bc6) & carry;
  prop(n);
}

// Overflow limit   = 340282366920938463463374607431768211456
// maximum possible = 62066584718639432850653955984558016
// Modular multiplication, c=a*b mod 2p
void __attribute__((always_inline)) modmul(spint *a, spint *b, spint *c) {
  spint v0, v1, v2, v3, v4, v5, v6, v7, v8;
  dpint t = 0;
  spint p2 = 0xada6e200000000;
  spint p3 = 0x781974ce994c68;
  spint p4 = 0xdea65faf0a29a;
  spint p5 = 0x7d01fe3ac5904a;
  spint p6 = 0x9202bdbe632650;
  spint p7 = 0x8c15b0036936e7;
  spint p8 = 0x255946a8869bc6;
  spint s, q = ((spint)1 << 56); // q is unsaturated radix
  spint mask = q - 1;
  t += (dpint)a[0] * b[0];
  v0 = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[0] * b[1];
  t += (dpint)a[1] * b[0];
  v1 = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[0] * b[2];
  t += (dpint)a[1] * b[1];
  t += (dpint)a[2] * b[0];
  t += (dpint)v0 * p2;
  v2 = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[0] * b[3];
  t += (dpint)a[1] * b[2];
  t += (dpint)a[2] * b[1];
  t += (dpint)a[3] * b[0];
  t += (dpint)v0 * p3;
  t += (dpint)v1 * p2;
  v3 = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[0] * b[4];
  t += (dpint)a[1] * b[3];
  t += (dpint)a[2] * b[2];
  t += (dpint)a[3] * b[1];
  t += (dpint)a[4] * b[0];
  t += (dpint)v0 * p4;
  t += (dpint)v1 * p3;
  t += (dpint)v2 * p2;
  v4 = (spint)((uspint)t & mask);
  t >>= 56;
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
  v5 = (spint)((uspint)t & mask);
  t >>= 56;
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
  v6 = (spint)((uspint)t & mask);
  t >>= 56;
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
  t += (dpint)v5 * p2;
  v7 = (spint)((uspint)t & mask);
  t >>= 56;
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
  t += (dpint)v6 * p2;
  v8 = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[1] * b[8];
  t += (dpint)a[2] * b[7];
  t += (dpint)a[3] * b[6];
  t += (dpint)a[4] * b[5];
  t += (dpint)a[5] * b[4];
  t += (dpint)a[6] * b[3];
  t += (dpint)a[7] * b[2];
  t += (dpint)a[8] * b[1];
  t += (dpint)v1 * p8;
  t += (dpint)v2 * p7;
  t += (dpint)v3 * p6;
  t += (dpint)v4 * p5;
  t += (dpint)v5 * p4;
  t += (dpint)v6 * p3;
  t += (dpint)v7 * p2;
  c[0] = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[2] * b[8];
  t += (dpint)a[3] * b[7];
  t += (dpint)a[4] * b[6];
  t += (dpint)a[5] * b[5];
  t += (dpint)a[6] * b[4];
  t += (dpint)a[7] * b[3];
  t += (dpint)a[8] * b[2];
  t += (dpint)v2 * p8;
  t += (dpint)v3 * p7;
  t += (dpint)v4 * p6;
  t += (dpint)v5 * p5;
  t += (dpint)v6 * p4;
  t += (dpint)v7 * p3;
  t += (dpint)v8 * p2;
  c[1] = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[3] * b[8];
  t += (dpint)a[4] * b[7];
  t += (dpint)a[5] * b[6];
  t += (dpint)a[6] * b[5];
  t += (dpint)a[7] * b[4];
  t += (dpint)a[8] * b[3];
  t += (dpint)v3 * p8;
  t += (dpint)v4 * p7;
  t += (dpint)v5 * p6;
  t += (dpint)v6 * p5;
  t += (dpint)v7 * p4;
  t += (dpint)v8 * p3;
  c[2] = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[4] * b[8];
  t += (dpint)a[5] * b[7];
  t += (dpint)a[6] * b[6];
  t += (dpint)a[7] * b[5];
  t += (dpint)a[8] * b[4];
  t += (dpint)v4 * p8;
  t += (dpint)v5 * p7;
  t += (dpint)v6 * p6;
  t += (dpint)v7 * p5;
  t += (dpint)v8 * p4;
  c[3] = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[5] * b[8];
  t += (dpint)a[6] * b[7];
  t += (dpint)a[7] * b[6];
  t += (dpint)a[8] * b[5];
  t += (dpint)v5 * p8;
  t += (dpint)v6 * p7;
  t += (dpint)v7 * p6;
  t += (dpint)v8 * p5;
  c[4] = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[6] * b[8];
  t += (dpint)a[7] * b[7];
  t += (dpint)a[8] * b[6];
  t += (dpint)v6 * p8;
  t += (dpint)v7 * p7;
  t += (dpint)v8 * p6;
  c[5] = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[7] * b[8];
  t += (dpint)a[8] * b[7];
  t += (dpint)v7 * p8;
  t += (dpint)v8 * p7;
  c[6] = (spint)((uspint)t & mask);
  t >>= 56;
  t += (dpint)a[8] * b[8];
  t += (dpint)v8 * p8;
  c[7] = (spint)((uspint)t & mask);
  t >>= 56;
  c[8] = (spint)t;
}

// Modular squaring, c=a*a  mod 2p
void __attribute__((always_inline)) modsqr(spint *a, spint *c) {
  spint v0, v1, v2, v3, v4, v5, v6, v7, v8;
  dpint tot, t = 0;
  spint p2 = 0xada6e200000000;
  spint p3 = 0x781974ce994c68;
  spint p4 = 0xdea65faf0a29a;
  spint p5 = 0x7d01fe3ac5904a;
  spint p6 = 0x9202bdbe632650;
  spint p7 = 0x8c15b0036936e7;
  spint p8 = 0x255946a8869bc6;
  spint s, q = ((spint)1 << 56); // q is unsaturated radix
  spint mask = q - 1;
  tot = (dpint)a[0] * a[0];
  t = tot;
  v0 = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[0] * a[1];
  tot *= 2;
  t += tot;
  v1 = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[0] * a[2];
  tot *= 2;
  tot += (dpint)a[1] * a[1];
  t += tot;
  t += (dpint)v0 * p2;
  v2 = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[0] * a[3];
  tot += (dpint)a[1] * a[2];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p3;
  t += (dpint)v1 * p2;
  v3 = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[0] * a[4];
  tot += (dpint)a[1] * a[3];
  tot *= 2;
  tot += (dpint)a[2] * a[2];
  t += tot;
  t += (dpint)v0 * p4;
  t += (dpint)v1 * p3;
  t += (dpint)v2 * p2;
  v4 = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[0] * a[5];
  tot += (dpint)a[1] * a[4];
  tot += (dpint)a[2] * a[3];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p5;
  t += (dpint)v1 * p4;
  t += (dpint)v2 * p3;
  t += (dpint)v3 * p2;
  v5 = (spint)((uspint)t & mask);
  t >>= 56;
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
  v6 = (spint)((uspint)t & mask);
  t >>= 56;
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
  t += (dpint)v5 * p2;
  v7 = (spint)((uspint)t & mask);
  t >>= 56;
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
  t += (dpint)v6 * p2;
  v8 = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[1] * a[8];
  tot += (dpint)a[2] * a[7];
  tot += (dpint)a[3] * a[6];
  tot += (dpint)a[4] * a[5];
  tot *= 2;
  t += tot;
  t += (dpint)v1 * p8;
  t += (dpint)v2 * p7;
  t += (dpint)v3 * p6;
  t += (dpint)v4 * p5;
  t += (dpint)v5 * p4;
  t += (dpint)v6 * p3;
  t += (dpint)v7 * p2;
  c[0] = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[2] * a[8];
  tot += (dpint)a[3] * a[7];
  tot += (dpint)a[4] * a[6];
  tot *= 2;
  tot += (dpint)a[5] * a[5];
  t += tot;
  t += (dpint)v2 * p8;
  t += (dpint)v3 * p7;
  t += (dpint)v4 * p6;
  t += (dpint)v5 * p5;
  t += (dpint)v6 * p4;
  t += (dpint)v7 * p3;
  t += (dpint)v8 * p2;
  c[1] = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[3] * a[8];
  tot += (dpint)a[4] * a[7];
  tot += (dpint)a[5] * a[6];
  tot *= 2;
  t += tot;
  t += (dpint)v3 * p8;
  t += (dpint)v4 * p7;
  t += (dpint)v5 * p6;
  t += (dpint)v6 * p5;
  t += (dpint)v7 * p4;
  t += (dpint)v8 * p3;
  c[2] = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[4] * a[8];
  tot += (dpint)a[5] * a[7];
  tot *= 2;
  tot += (dpint)a[6] * a[6];
  t += tot;
  t += (dpint)v4 * p8;
  t += (dpint)v5 * p7;
  t += (dpint)v6 * p6;
  t += (dpint)v7 * p5;
  t += (dpint)v8 * p4;
  c[3] = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[5] * a[8];
  tot += (dpint)a[6] * a[7];
  tot *= 2;
  t += tot;
  t += (dpint)v5 * p8;
  t += (dpint)v6 * p7;
  t += (dpint)v7 * p6;
  t += (dpint)v8 * p5;
  c[4] = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[6] * a[8];
  tot *= 2;
  tot += (dpint)a[7] * a[7];
  t += tot;
  t += (dpint)v6 * p8;
  t += (dpint)v7 * p7;
  t += (dpint)v8 * p6;
  c[5] = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[7] * a[8];
  tot *= 2;
  t += tot;
  t += (dpint)v7 * p8;
  t += (dpint)v8 * p7;
  c[6] = (spint)((uspint)t & mask);
  t >>= 56;
  tot = (dpint)a[8] * a[8];
  t += tot;
  t += (dpint)v8 * p8;
  c[7] = (spint)((uspint)t & mask);
  t >>= 56;
  c[8] = (spint)t;
}

// copy
void __attribute__((always_inline)) modcpy(spint *a, spint *c) {
  int i;
  for (i = 0; i < 9; i++) {
    c[i] = a[i];
  }
}

// Calculate progenitor
void modpro(spint *w, spint *z) {
  int i;
  spint x[9];
  spint t0[9], t1[9], t2[9], t3[9], t4[9], t5[9], t6[9], t7[9], t8[9], t9[9],
      t10[9], t11[9], t12[9], t13[9], t14[9], t15[9], t16[9], t17[9], t18[9],
      t19[9], t20[9], t21[9], t22[9], t23[9], t24[9], t25[9], t26[9], t27[9],
      t28[9], t29[9], t30[9];
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
  spint s[9], t[9];
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
  spint c[9] = {0xa546e4e8a0c754, 0xcb993b5943e89e, 0x8dda851bf555c8,
                0xc717daa328de11, 0x307cf4eb090199, 0x8537d92d096114,
                0x96f0b3d2dce4ee, 0x73aacac5bc94fd, 0x16de108943ddd4};
  modmul(m, c, n);
}

// Convert n back to normal form, m=redc(n)
void redc(spint *n, spint *m) {
  spint c[9];
  c[0] = 1;
  for (int i = 1; i < 9; i++)
    c[i] = 0;
  modmul(n, c, m);
  modfsb(m);
}

// is unity?
int modis1(spint *a) {
  spint c[9];
  sspint c0, d = 0;
  redc(a, c);
  for (int i = 1; i < 9; i++) {
    d |= c[i];
  }
  c0 = (sspint)c[0];
  return (1 & ((d - 1) >> 56) & (((c0 ^ 1) - 1) >> 56));
}

// is zero?
int modis0(spint *a) {
  sspint d = 0;
  for (int i = 0; i < 9; i++) {
    d |= a[i];
  }
  return (1 & ((d - 1) >> 56));
}

// set to zero
void modzer(spint *a) {
  for (int i = 0; i < 9; i++)
    a[i] = 0;
}

// set to one
void modone(spint *a) {
  a[0] = 1;
  for (int i = 1; i < 9; i++)
    a[i] = 0;
  nres(a, a);
}

// Test for quadratic residue
int modqr(spint *h, spint *x) {
  int i;
  spint r[9];
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
  spint s[9], y[9];
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
  a[9 - 1] = ((a[9 - 1] << n)) | (a[9 - 2] >> (56 - n));
  for (int i = 9 - 2; i > 0; i--) {
    a[i] = ((a[i] << n) & 0xffffffffffffff) | (a[i - 1] >> (56 - n));
  }
  a[0] = (a[0] << n) & 0xffffffffffffff;
}

// shift right by less than a word. Return shifted out part
int modshr(int n, spint *a) {
  spint r = a[0] & (((spint)1 << n) - 1);
  for (int i = 0; i < 9 - 1; i++) {
    a[i] = (a[i] >> n) | ((a[i + 1] << (56 - n)) & 0xffffffffffffff);
  }
  a[9 - 1] = a[9 - 1] >> n;
  return r;
}

/* API functions calling generated code */

const digit_t p[NWORDS_ORDER] =  { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x994C68ADA6E1FFFF, 0xFAF0A29A781974CE, 0xFE3AC5904A0DEA65, 0x02BDBE6326507D01, 0x8C15B0036936E792, 0x00255946A8869BC6 };

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
