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
  spint d, mask = ((spint)1 << 52) - 1;
  sspint carry = (sspint)n[0] >> 52;
  n[0] &= mask;
  for (int i = 1; i < 4; i++) {
    d = n[i] + carry;
    n[i] = d & mask;
    carry = (sspint)d >> 52;
  }
  n[4] += carry;
  return ((sspint)n[4] >> 63);
}

// propagate carries and add p if negative, propagate carries again
static __attribute__((always_inline)) void flatten(spint *n) {
  spint q = ((spint)1 << 52);
  sspint carry = prop(n);
  n[0] -= 1 & carry;
  n[1] += (1 * (spint)0x4935514800000) & carry;
  n[2] += (1 * (spint)0x7407437252c9e) & carry;
  n[3] += (1 * (spint)0xd98c33a6a8658) & carry;
  n[4] += (1 * (spint)0x34e29e286b95) & carry;
  prop(n);
}

// Montgomery final subtract
void __attribute__((always_inline)) modfsb(spint *n) {
  spint q = ((spint)1 << 52);
  n[0] += 1;
  n[1] -= 1 * (spint)0x4935514800000;
  n[2] -= 1 * (spint)0x7407437252c9e;
  n[3] -= 1 * (spint)0xd98c33a6a8658;
  n[4] -= 1 * (spint)0x34e29e286b95;
  flatten(n);
}

// Modular addition - reduce less than 2p
void __attribute__((always_inline)) modadd(spint *a, spint *b, spint *n) {
  spint q = ((spint)1 << 52);
  sspint carry;
  n[0] = a[0] + b[0];
  n[1] = a[1] + b[1];
  n[2] = a[2] + b[2];
  n[3] = a[3] + b[3];
  n[4] = a[4] + b[4];
  n[0] += 2;
  n[1] -= 2 * (spint)0x4935514800000;
  n[2] -= 2 * (spint)0x7407437252c9e;
  n[3] -= 2 * (spint)0xd98c33a6a8658;
  n[4] -= 2 * (spint)0x34e29e286b95;
  carry = prop(n);
  n[0] -= 2 & carry;
  n[1] += (2 * (spint)0x4935514800000) & carry;
  n[2] += (2 * (spint)0x7407437252c9e) & carry;
  n[3] += (2 * (spint)0xd98c33a6a8658) & carry;
  n[4] += (2 * (spint)0x34e29e286b95) & carry;
  prop(n);
}

// Modular subtraction - reduce less than 2p
void __attribute__((always_inline)) modsub(spint *a, spint *b, spint *n) {
  spint q = ((spint)1 << 52);
  sspint carry;
  n[0] = a[0] - b[0];
  n[1] = a[1] - b[1];
  n[2] = a[2] - b[2];
  n[3] = a[3] - b[3];
  n[4] = a[4] - b[4];
  carry = prop(n);
  n[0] -= 2 & carry;
  n[1] += (2 * (spint)0x4935514800000) & carry;
  n[2] += (2 * (spint)0x7407437252c9e) & carry;
  n[3] += (2 * (spint)0xd98c33a6a8658) & carry;
  n[4] += (2 * (spint)0x34e29e286b95) & carry;
  prop(n);
}

// Overflow limit   = 340282366920938463463374607431768211456
// maximum possible = 133902696273723909579593872826746
// Modular multiplication, c=a*b mod 2p
void __attribute__((always_inline)) modmul(spint *a, spint *b, spint *c) {
  spint v0, v1, v2, v3, v4;
  dpint t = 0;
  spint p1 = 0x4935514800000;
  spint p2 = 0x7407437252c9e;
  spint p3 = 0xd98c33a6a8658;
  spint p4 = 0x34e29e286b95;
  spint s, q = ((spint)1 << 52); // q is unsaturated radix
  spint mask = q - 1;
  t += (dpint)a[0] * b[0];
  v0 = (spint)((uspint)t & mask);
  t >>= 52;
  t += (dpint)a[0] * b[1];
  t += (dpint)a[1] * b[0];
  t += (dpint)v0 * p1;
  v1 = (spint)((uspint)t & mask);
  t >>= 52;
  t += (dpint)a[0] * b[2];
  t += (dpint)a[1] * b[1];
  t += (dpint)a[2] * b[0];
  t += (dpint)v0 * p2;
  t += (dpint)v1 * p1;
  v2 = (spint)((uspint)t & mask);
  t >>= 52;
  t += (dpint)a[0] * b[3];
  t += (dpint)a[1] * b[2];
  t += (dpint)a[2] * b[1];
  t += (dpint)a[3] * b[0];
  t += (dpint)v0 * p3;
  t += (dpint)v1 * p2;
  t += (dpint)v2 * p1;
  v3 = (spint)((uspint)t & mask);
  t >>= 52;
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
  t >>= 52;
  t += (dpint)a[1] * b[4];
  t += (dpint)a[2] * b[3];
  t += (dpint)a[3] * b[2];
  t += (dpint)a[4] * b[1];
  t += (dpint)v1 * p4;
  t += (dpint)v2 * p3;
  t += (dpint)v3 * p2;
  t += (dpint)v4 * p1;
  c[0] = (spint)((uspint)t & mask);
  t >>= 52;
  t += (dpint)a[2] * b[4];
  t += (dpint)a[3] * b[3];
  t += (dpint)a[4] * b[2];
  t += (dpint)v2 * p4;
  t += (dpint)v3 * p3;
  t += (dpint)v4 * p2;
  c[1] = (spint)((uspint)t & mask);
  t >>= 52;
  t += (dpint)a[3] * b[4];
  t += (dpint)a[4] * b[3];
  t += (dpint)v3 * p4;
  t += (dpint)v4 * p3;
  c[2] = (spint)((uspint)t & mask);
  t >>= 52;
  t += (dpint)a[4] * b[4];
  t += (dpint)v4 * p4;
  c[3] = (spint)((uspint)t & mask);
  t >>= 52;
  c[4] = (spint)t;
}

// Modular squaring, c=a*a  mod 2p
void __attribute__((always_inline)) modsqr(spint *a, spint *c) {
  spint v0, v1, v2, v3, v4;
  dpint tot, t = 0;
  spint p1 = 0x4935514800000;
  spint p2 = 0x7407437252c9e;
  spint p3 = 0xd98c33a6a8658;
  spint p4 = 0x34e29e286b95;
  spint s, q = ((spint)1 << 52); // q is unsaturated radix
  spint mask = q - 1;
  tot = (dpint)a[0] * a[0];
  t = tot;
  v0 = (spint)((uspint)t & mask);
  t >>= 52;
  tot = (dpint)a[0] * a[1];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p1;
  v1 = (spint)((uspint)t & mask);
  t >>= 52;
  tot = (dpint)a[0] * a[2];
  tot *= 2;
  tot += (dpint)a[1] * a[1];
  t += tot;
  t += (dpint)v0 * p2;
  t += (dpint)v1 * p1;
  v2 = (spint)((uspint)t & mask);
  t >>= 52;
  tot = (dpint)a[0] * a[3];
  tot += (dpint)a[1] * a[2];
  tot *= 2;
  t += tot;
  t += (dpint)v0 * p3;
  t += (dpint)v1 * p2;
  t += (dpint)v2 * p1;
  v3 = (spint)((uspint)t & mask);
  t >>= 52;
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
  t >>= 52;
  tot = (dpint)a[1] * a[4];
  tot += (dpint)a[2] * a[3];
  tot *= 2;
  t += tot;
  t += (dpint)v1 * p4;
  t += (dpint)v2 * p3;
  t += (dpint)v3 * p2;
  t += (dpint)v4 * p1;
  c[0] = (spint)((uspint)t & mask);
  t >>= 52;
  tot = (dpint)a[2] * a[4];
  tot *= 2;
  tot += (dpint)a[3] * a[3];
  t += tot;
  t += (dpint)v2 * p4;
  t += (dpint)v3 * p3;
  t += (dpint)v4 * p2;
  c[1] = (spint)((uspint)t & mask);
  t >>= 52;
  tot = (dpint)a[3] * a[4];
  tot *= 2;
  t += tot;
  t += (dpint)v3 * p4;
  t += (dpint)v4 * p3;
  c[2] = (spint)((uspint)t & mask);
  t >>= 52;
  tot = (dpint)a[4] * a[4];
  t += tot;
  t += (dpint)v4 * p4;
  c[3] = (spint)((uspint)t & mask);
  t >>= 52;
  c[4] = (spint)t;
}

// copy
void __attribute__((always_inline)) modcpy(spint *a, spint *c) {
  int i;
  for (i = 0; i < 5; i++) {
    c[i] = a[i];
  }
}

// Calculate progenitor
void modpro(spint *w, spint *z) {
  int i;
  spint x[5];
  spint t0[5], t1[5], t2[5], t3[5], t4[5], t5[5], t6[5], t7[5], t8[5], t9[5],
      t10[5], t11[5], t12[5], t13[5], t14[5], t15[5], t16[5], t17[5], t18[5];
  modcpy(w, x);
  modsqr(x, t13);
  modsqr(t13, t11);
  modmul(x, t11, t0);
  modsqr(t11, t1);
  modmul(t11, t1, t4);
  modmul(t13, t4, t9);
  modmul(t11, t9, t3);
  modsqr(t4, z);
  modmul(t13, z, t2);
  modmul(t0, z, t7);
  modmul(x, t7, t16);
  modmul(t1, t7, t5);
  modmul(t3, t5, t6);
  modmul(t4, t6, t15);
  modmul(t13, t15, t17);
  modmul(t1, t15, t8);
  modmul(t13, t8, t10);
  modmul(t1, t10, t1);
  modmul(t11, t1, t4);
  modmul(t9, t1, t12);
  modmul(t11, t12, t11);
  modmul(t3, t11, t3);
  modmul(t2, t3, t2);
  modmul(t9, t2, t9);
  modmul(z, t9, t14);
  modmul(t13, t14, t13);
  modmul(z, t13, z);
  modmul(t16, z, t16);
  modsqr(z, t18);
  modmul(t4, t18, z);
  for (i = 0; i < 6; i++) {
    modsqr(t18, t18);
  }
  modmul(t17, t18, t17);
  for (i = 0; i < 10; i++) {
    modsqr(t17, t17);
  }
  modmul(t16, t17, t16);
  for (i = 0; i < 8; i++) {
    modsqr(t16, t16);
  }
  modmul(t15, t16, t15);
  for (i = 0; i < 9; i++) {
    modsqr(t15, t15);
  }
  modmul(t14, t15, t14);
  for (i = 0; i < 9; i++) {
    modsqr(t14, t14);
  }
  modmul(t13, t14, t13);
  for (i = 0; i < 9; i++) {
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
  for (i = 0; i < 9; i++) {
    modsqr(t10, t10);
  }
  modmul(t9, t10, t9);
  for (i = 0; i < 7; i++) {
    modsqr(t9, t9);
  }
  modmul(t8, t9, t8);
  for (i = 0; i < 9; i++) {
    modsqr(t8, t8);
  }
  modmul(t7, t8, t8);
  for (i = 0; i < 12; i++) {
    modsqr(t8, t8);
  }
  modmul(t7, t8, t7);
  for (i = 0; i < 10; i++) {
    modsqr(t7, t7);
  }
  modmul(t6, t7, t6);
  for (i = 0; i < 8; i++) {
    modsqr(t6, t6);
  }
  modmul(t5, t6, t5);
  for (i = 0; i < 9; i++) {
    modsqr(t5, t5);
  }
  modmul(t4, t5, t4);
  for (i = 0; i < 9; i++) {
    modsqr(t4, t4);
  }
  modmul(t3, t4, t3);
  for (i = 0; i < 10; i++) {
    modsqr(t3, t3);
  }
  modmul(t2, t3, t2);
  for (i = 0; i < 8; i++) {
    modsqr(t2, t2);
  }
  modmul(t1, t2, t1);
  for (i = 0; i < 6; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t0);
  for (i = 0; i < 12; i++) {
    modsqr(t0, t0);
  }
  modmul(z, t0, t0);
  for (i = 0; i < 9; i++) {
    modsqr(t0, t0);
  }
  modmul(z, t0, t0);
  for (i = 0; i < 9; i++) {
    modsqr(t0, t0);
  }
  modmul(z, t0, t0);
  for (i = 0; i < 9; i++) {
    modsqr(t0, t0);
  }
  modmul(z, t0, t0);
  for (i = 0; i < 9; i++) {
    modsqr(t0, t0);
  }
  modmul(z, t0, t0);
  for (i = 0; i < 9; i++) {
    modsqr(t0, t0);
  }
  modmul(z, t0, t0);
  for (i = 0; i < 9; i++) {
    modsqr(t0, t0);
  }
  modmul(z, t0, t0);
  modsqr(t0, t0);
  modmul(x, t0, t0);
  for (i = 0; i < 9; i++) {
    modsqr(t0, t0);
  }
  modmul(z, t0, z);
}

// calculate inverse, provide progenitor h if available
void modinv(spint *x, spint *h, spint *z) {
  int i;
  spint s[5], t[5];
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
  spint c[5] = {0x5ae400674d441, 0x6bd078e623362, 0x907db203f8290,
                0x978af0e941c, 0x55dae38dfc1};
  modmul(m, c, n);
}

// Convert n back to normal form, m=redc(n)
void redc(spint *n, spint *m) {
  spint c[5];
  c[0] = 1;
  for (int i = 1; i < 5; i++)
    c[i] = 0;
  modmul(n, c, m);
  modfsb(m);
}

// is unity?
int modis1(spint *a) {
  spint c[5];
  sspint c0, d = 0;
  redc(a, c);
  for (int i = 1; i < 5; i++) {
    d |= c[i];
  }
  c0 = (sspint)c[0];
  return (1 & ((d - 1) >> 52) & (((c0 ^ 1) - 1) >> 52));
}

// is zero?
int modis0(spint *a) {
  sspint d = 0;
  for (int i = 0; i < 5; i++) {
    d |= a[i];
  }
  return (1 & ((d - 1) >> 52));
}

// set to zero
void modzer(spint *a) {
  for (int i = 0; i < 5; i++)
    a[i] = 0;
}

// set to one
void modone(spint *a) {
  a[0] = 1;
  for (int i = 1; i < 5; i++)
    a[i] = 0;
  nres(a, a);
}

// Test for quadratic residue
int modqr(spint *h, spint *x) {
  int i;
  spint r[5];
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
  spint s[5], y[5];
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
  a[5 - 1] = ((a[5 - 1] << n)) | (a[5 - 2] >> (52 - n));
  for (int i = 5 - 2; i > 0; i--) {
    a[i] = ((a[i] << n) & 0xfffffffffffff) | (a[i - 1] >> (52 - n));
  }
  a[0] = (a[0] << n) & 0xfffffffffffff;
}

// shift right by less than a word. Return shifted out part
int modshr(int n, digit_t *a) {
  digit_t r = a[0] & (((digit_t)1 << n) - 1);
  for (int i = 0; i < 5 - 1; i++) {
    a[i] = (a[i] >> n) | ((a[i + 1] << (52 - n)) & 0xfffffffffffff);
  }
  a[5 - 1] = a[5 - 1] >> n;
  return r;
}

/* API functions calling generated code */

const digit_t p[NWORDS_ORDER] =  { 0xFFFFFFFFFFFFFFFF, 0x252C9E49355147FF, 0x33A6A86587407437, 0x34E29E286B95D98C };

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
    for (int i = 0; i < 32; i++) {
        ((char *) out)[i] = x[0] & 0xff;
        modshr(8, x);
    }
}

void fp_from_digit_array(digit_t* out, const digit_t* a) {
    for (int i = 0; i < NWORDS_FIELD; i++) {
        out[i] = 0;
    }
    for (int i = 32 - 1; i >= 0; i--) {
        modshl(8, out);
        out[0] += (digit_t)((unsigned char *) a)[i];
    }
}

#endif
