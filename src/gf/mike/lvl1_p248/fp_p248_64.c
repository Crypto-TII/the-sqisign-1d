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
  spint d, mask = ((spint)1 << 51) - 1;
  sspint carry = (sspint)n[0] >> 51;
  n[0] &= mask;
  for (int i = 1; i < 4; i++) {
    d = n[i] + carry;
    n[i] = d & mask;
    carry = (sspint)d >> 51;
  }
  n[4] += carry;
  return ((sspint)n[4] >> 63);
}

// propagate carries and add p if negative, propagate carries again
static __attribute__((always_inline)) void flatten(spint *n) {
  spint q = ((spint)1 << 51);
  sspint carry = prop(n);
  n[0] -= 1 & carry;
  n[4] += (1 * (spint)0x500000000000) & carry;
  prop(n);
}

// Montgomery final subtract
void __attribute__((always_inline)) modfsb(spint *n) {
  spint q = ((spint)1 << 51);
  n[0] += 1;
  n[4] -= 1 * (spint)0x500000000000;
  flatten(n);
}

// Modular addition - reduce less than 2p
void __attribute__((always_inline)) modadd(spint *a, spint *b, spint *n) {
  spint q = ((spint)1 << 51);
  sspint carry;
  n[0] = a[0] + b[0];
  n[1] = a[1] + b[1];
  n[2] = a[2] + b[2];
  n[3] = a[3] + b[3];
  n[4] = a[4] + b[4];
  n[0] += 2;
  n[4] -= 2 * (spint)0x500000000000;
  carry = prop(n);
  n[0] -= 2 & carry;
  n[4] += (2 * (spint)0x500000000000) & carry;
  prop(n);
}

// Modular subtraction - reduce less than 2p
void __attribute__((always_inline)) modsub(spint *a, spint *b, spint *n) {
  spint q = ((spint)1 << 51);
  sspint carry;
  n[0] = a[0] - b[0];
  n[1] = a[1] - b[1];
  n[2] = a[2] - b[2];
  n[3] = a[3] - b[3];
  n[4] = a[4] - b[4];
  carry = prop(n);
  n[0] -= 2 & carry;
  n[4] += (2 * (spint)0x500000000000) & carry;
  prop(n);
}

// Overflow limit   = 340282366920938463463374607431768211456
// maximum possible = 25551082561965953719787503747077
// Modular multiplication, c=a*b mod 2p
void __attribute__((always_inline)) modmul(spint *a, spint *b, spint *c) {
  spint v0, v1, v2, v3, v4;
  dpint t = 0;
  spint p4 = 0x500000000000;
  spint s, q = ((spint)1 << 51); // q is unsaturated radix
  spint mask = q - 1;
  t += (dpint)a[0] * b[0];
  v0 = (spint)((uspint)t & mask);
  t >>= 51;
  t += (dpint)a[0] * b[1];
  t += (dpint)a[1] * b[0];
  v1 = (spint)((uspint)t & mask);
  t >>= 51;
  t += (dpint)a[0] * b[2];
  t += (dpint)a[1] * b[1];
  t += (dpint)a[2] * b[0];
  v2 = (spint)((uspint)t & mask);
  t >>= 51;
  t += (dpint)a[0] * b[3];
  t += (dpint)a[1] * b[2];
  t += (dpint)a[2] * b[1];
  t += (dpint)a[3] * b[0];
  v3 = (spint)((uspint)t & mask);
  t >>= 51;
  t += (dpint)a[0] * b[4];
  t += (dpint)a[1] * b[3];
  t += (dpint)a[2] * b[2];
  t += (dpint)a[3] * b[1];
  t += (dpint)a[4] * b[0];
  t += (dpint)v0 * p4;
  v4 = (spint)((uspint)t & mask);
  t >>= 51;
  t += (dpint)a[1] * b[4];
  t += (dpint)a[2] * b[3];
  t += (dpint)a[3] * b[2];
  t += (dpint)a[4] * b[1];
  t += (dpint)v1 * p4;
  c[0] = (spint)((uspint)t & mask);
  t >>= 51;
  t += (dpint)a[2] * b[4];
  t += (dpint)a[3] * b[3];
  t += (dpint)a[4] * b[2];
  t += (dpint)v2 * p4;
  c[1] = (spint)((uspint)t & mask);
  t >>= 51;
  t += (dpint)a[3] * b[4];
  t += (dpint)a[4] * b[3];
  t += (dpint)v3 * p4;
  c[2] = (spint)((uspint)t & mask);
  t >>= 51;
  t += (dpint)a[4] * b[4];
  t += (dpint)v4 * p4;
  c[3] = (spint)((uspint)t & mask);
  t >>= 51;
  c[4] = (spint)t;
}

// Modular squaring, c=a*a  mod 2p
void __attribute__((always_inline)) modsqr(spint *a, spint *c) {
  spint v0, v1, v2, v3, v4;
  dpint tot, t = 0;
  spint p4 = 0x500000000000;
  spint s, q = ((spint)1 << 51); // q is unsaturated radix
  spint mask = q - 1;
  tot = (dpint)a[0] * a[0];
  t = tot;
  v0 = (spint)((uspint)t & mask);
  t >>= 51;
  tot = (dpint)a[0] * a[1];
  tot *= 2;
  t += tot;
  v1 = (spint)((uspint)t & mask);
  t >>= 51;
  tot = (dpint)a[0] * a[2];
  tot *= 2;
  tot += (dpint)a[1] * a[1];
  t += tot;
  v2 = (spint)((uspint)t & mask);
  t >>= 51;
  tot = (dpint)a[0] * a[3];
  tot += (dpint)a[1] * a[2];
  tot *= 2;
  t += tot;
  v3 = (spint)((uspint)t & mask);
  t >>= 51;
  tot = (dpint)a[0] * a[4];
  tot += (dpint)a[1] * a[3];
  tot *= 2;
  tot += (dpint)a[2] * a[2];
  t += tot;
  t += (dpint)v0 * p4;
  v4 = (spint)((uspint)t & mask);
  t >>= 51;
  tot = (dpint)a[1] * a[4];
  tot += (dpint)a[2] * a[3];
  tot *= 2;
  t += tot;
  t += (dpint)v1 * p4;
  c[0] = (spint)((uspint)t & mask);
  t >>= 51;
  tot = (dpint)a[2] * a[4];
  tot *= 2;
  tot += (dpint)a[3] * a[3];
  t += tot;
  t += (dpint)v2 * p4;
  c[1] = (spint)((uspint)t & mask);
  t >>= 51;
  tot = (dpint)a[3] * a[4];
  tot *= 2;
  t += tot;
  t += (dpint)v3 * p4;
  c[2] = (spint)((uspint)t & mask);
  t >>= 51;
  tot = (dpint)a[4] * a[4];
  t += tot;
  t += (dpint)v4 * p4;
  c[3] = (spint)((uspint)t & mask);
  t >>= 51;
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
  spint t0[5], t1[5], t2[5], t3[5], t4[5];
  modcpy(w, x);
  modsqr(x, z);
  modmul(x, z, t0);
  modsqr(t0, z);
  modmul(x, z, z);
  modsqr(z, t1);
  modsqr(t1, t3);
  modsqr(t3, t2);
  modcpy(t2, t4);
  for (i = 0; i < 3; i++) {
    modsqr(t4, t4);
  }
  modmul(t2, t4, t2);
  modcpy(t2, t4);
  for (i = 0; i < 6; i++) {
    modsqr(t4, t4);
  }
  modmul(t2, t4, t2);
  modcpy(t2, t4);
  for (i = 0; i < 2; i++) {
    modsqr(t4, t4);
  }
  modmul(t3, t4, t3);
  for (i = 0; i < 13; i++) {
    modsqr(t3, t3);
  }
  modmul(t2, t3, t2);
  modcpy(t2, t3);
  for (i = 0; i < 27; i++) {
    modsqr(t3, t3);
  }
  modmul(t2, t3, t2);
  modmul(z, t2, z);
  modcpy(z, t2);
  for (i = 0; i < 4; i++) {
    modsqr(t2, t2);
  }
  modmul(t1, t2, t1);
  modmul(t0, t1, t0);
  modmul(t1, t0, t1);
  modmul(t0, t1, t0);
  modmul(t1, t0, t2);
  modmul(t0, t2, t0);
  modmul(t1, t0, t1);
  for (i = 0; i < 63; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t1);
  for (i = 0; i < 64; i++) {
    modsqr(t1, t1);
  }
  modmul(t0, t1, t0);
  for (i = 0; i < 57; i++) {
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
  spint c[5] = {0x4cccccccccf5c, 0x1999999999999, 0x3333333333333,
                0x6666666666666, 0xccccccccccc};
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
  return (1 & ((d - 1) >> 51) & (((c0 ^ 1) - 1) >> 51));
}

// is zero?
int modis0(spint *a) {
  sspint d = 0;
  for (int i = 0; i < 5; i++) {
    d |= a[i];
  }
  return (1 & ((d - 1) >> 51));
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
  a[5 - 1] = ((a[5 - 1] << n)) | (a[5 - 2] >> (51 - n));
  for (int i = 5 - 2; i > 0; i--) {
    a[i] = ((a[i] << n) & 0x7ffffffffffff) | (a[i - 1] >> (51 - n));
  }
  a[0] = (a[0] << n) & 0x7ffffffffffff;
}

// shift right by less than a word. Return shifted out part
int modshr(int n, spint *a) {
  spint r = a[0] & (((spint)1 << n) - 1);
  for (int i = 0; i < 5 - 1; i++) {
    a[i] = (a[i] >> n) | ((a[i + 1] << (51 - n)) & 0x7ffffffffffff);
  }
  a[5 - 1] = a[5 - 1] >> n;
  return r;
}

/* API functions calling generated code */
const digit_t p[NWORDS_ORDER] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,0x04ffffffffffffff };

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