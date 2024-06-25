#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#ifdef RADIX_32
typedef int32_t sidigit_t;
#elif defined(RADIX_64)
typedef int64_t sidigit_t;
#endif

#define SELECT(A, B, C)	((-(C) & ((A) ^ (B))) ^ (A))
#define SIGN(D)	(((digit_t)D) >> (RADIX - 1))
#define MASK(B)	((-(digit_t)((B) >= RADIX)) | (((digit_t)1 << ((B) % RADIX)) - 1))
#define LSH(H, L, I)													\
		(H << I) | (L & -(I != 0)) >> ((RADIX - I) & (RADIX - 1))

static digit_t porninstep(sdigit_t m[4], const digit_t f[2], const digit_t g[2], 
		digit_t k, size_t s) {
	digit_t limbx, ai = 1, bi = 0, ci = 0, di = 1;
	digit_t g_lo = g[0], g_hi = g[1], f_lo = f[0], f_hi = f[1];
	digit_t t_lo, t_hi, odd, borrow, xorm;

	/* Unrolling twice gives some small speedup. */
	for (size_t i = 0; i < s; i++) {
		odd = 0 - (g_lo & 1);

		/* g_ -= f_ if g_ is odd */
		t_lo = g_lo, t_hi = g_hi;

		borrow = 0;
		limbx = g_lo - (f_lo & odd);
		borrow = (g_lo < limbx);
		g_lo = limbx;

		limbx = g_hi - (f_hi & odd);
		xorm = limbx - borrow;
		borrow = -((g_hi < limbx) || (borrow && !limbx));
		g_hi = xorm;

		k += ((t_lo & f_lo) >> 1) & borrow;

		/* negate g_-f_ if it borrowed */
		g_lo ^= borrow;
		g_hi ^= borrow;
		limbx = g_lo + (borrow & 1);
		g_hi += (g_lo < limbx);
		g_lo = limbx;

		/* f_=g_ if g_-f_ borrowed */
		f_lo = ((t_lo ^ f_lo) & borrow) ^ f_lo;
		f_hi = ((t_hi ^ f_hi) & borrow) ^ f_hi;

		/* exchange ai and ci if g_-f_ borrowed */
		xorm = (ai ^ ci) & borrow;
		ai ^= xorm;
		ci ^= xorm;

		/* exchange bi and di if g_-f_ borrowed */
		xorm = (bi ^ di) & borrow;
		bi ^= xorm;
		di ^= xorm;

		/* subtract if g_ was odd */
		ai -= ci & odd;
		bi -= di & odd;

		ci <<= 1;
		di <<= 1;
		g_lo >>= 1;
		g_lo |= g_hi << (RADIX - 1);
		g_hi >>= 1;

		k += (f_lo + 2) >> 2;
	}

	m[0] = ai;
	m[1] = bi;
	m[2] = ci;
	m[3] = di;

	return k;
}

bool fp_is_square(const fp_t a) {
	const size_t s = RADIX - 2;
	digit_t f[NWORDS_FIELD + 1], g[NWORDS_FIELD + 1], t[NWORDS_FIELD] = { 0 };
	digit_t t0[NWORDS_FIELD + 1], t1[NWORDS_FIELD + 1];
	digit_t g_[2], f_[2], neg, l, g_hi, g_lo, f_hi, f_lo, mask, k = 0;
	sdigit_t m[4];

	if (fp_is_zero(a)) {
		return 0;
	}

	fp_copy(f, p);
	/* Convert a from Montgomery form. */
	fp_frommont(g, a);

	//int iterations = 2 * NWORDS_FIELD * RADIX;
	//for (size_t i = 0; i < iterations / s; i++) {
	while (!fp_is_zero(g)) {
		f_hi = f[NWORDS_FIELD - 1], f_lo = f[NWORDS_FIELD - 2];
		g_hi = g[NWORDS_FIELD - 1], g_lo = g[NWORDS_FIELD - 2];
		for (int j = NWORDS_FIELD - 2; j >= 0; j--) {
			l = (f_hi | g_hi);
			l = ~l & (l - 1);
			mask = -(l >> (RADIX - 1));
			f_hi = ((f_lo ^ f_hi) & mask) ^ f_hi;
			g_hi = ((g_lo ^ g_hi) & mask) ^ g_hi;
			f_lo = ((f[j] ^ f_lo) & mask) ^ f_lo;
			g_lo = ((g[j] ^ g_lo) & mask) ^ g_lo;
		}
		size_t len = __builtin_clzll(f_hi | g_hi);
		f_[0] = f[0], f_[1] = LSH(f_hi, f_lo, len);
		g_[0] = g[0], g_[1] = LSH(g_hi, g_lo, len);

		k = porninstep(m, f_, g_, k, s);

		t0[NWORDS_FIELD] = bn_muls_low(t0, g, 0, m[0], NWORDS_FIELD);
		t1[NWORDS_FIELD] = bn_muls_low(t1, f, 0, m[1], NWORDS_FIELD);
		bn_addn_low(t0, t0, t1, NWORDS_FIELD + 1);
		neg = SIGN(t0[NWORDS_FIELD]);
		bn_rshs_low(t, t0, NWORDS_FIELD + 1, (RADIX - 2));
		bn_negs_low(t, t, neg, NWORDS_FIELD);

		t0[NWORDS_FIELD] = bn_muls_low(t0, g, 0, m[2], NWORDS_FIELD);
		t1[NWORDS_FIELD] = bn_muls_low(t1, f, 0, m[3], NWORDS_FIELD);
		bn_addn_low(t1, t1, t0, NWORDS_FIELD + 1);
		bn_rshs_low(f, t1, NWORDS_FIELD + 1, (RADIX - 2));
		bn_negs_low(f, f, SIGN(t1[NWORDS_FIELD]), NWORDS_FIELD);

		fp_copy(g, t);
		k += (f[0] >> 1) & neg;
	}

	return (k & 1 ? 0 : 1);
}
