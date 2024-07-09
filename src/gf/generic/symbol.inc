#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#ifdef RADIX_32
typedef int32_t sidigit_t;
#elif defined(RADIX_64)
typedef int64_t sidigit_t;
#endif

#define SIGN(D)	(((digit_t)D) >> (RADIX - 1))
#define SWAP(A, B)	((A) ^= (B), (B) ^= (A), (A) ^= (B))

#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define MIN(A, B) ((A) < (B) ? (A) : (B))

#define MASK(B)	((-(digit_t)((B) >= RADIX)) | (((digit_t)1 << ((B) % RADIX)) - 1))
#define LMASK	(MASK(RADIX >> 1))
#define HMASK	(LMASK << (RADIX >> 1))

bool fp_is_square(const fp_t a) {
	sidigit_t ai, bi, ci, di;
	digit_t n, d, t;
	digit_t t0[NWORDS_FIELD + 1] = { 0 }, t1[NWORDS_FIELD + 1] = { 0 };
	digit_t t2[NWORDS_FIELD + 1] = { 0 }, t3[NWORDS_FIELD + 1] = { 0 };
	unsigned int neg, z, i, j, s = (RADIX >> 1) - 2;

	/* Optimized Pornin's Algorithm by Aleksei Vambol from
	 * https://github.com/privacy-scaling-explorations/halo2curves/pull/95 */

	if (fp_is_zero(a)) {
		return 0;
	}

	/* Convert a from Montgomery form. */
	fp_frommont(t0, a);
	fp_copy(t1, p);
	t = 0;
	j = NWORDS_FIELD;

	while (1) {
		ai = di = 1;
		bi = ci = 0;

		while ((t0[j] == 0) && (t1[j] == 0)) {
			j--;
		}
		if (j == 0) {
			n = t0[0];
			d = t1[0];
			while (n != 0) {
				if (n & 1) {
					if (n < d) {
						SWAP(n, d);
						t ^= (n & d);
					}
					n = (n - d) >> 1;
					t ^= d ^ (d >> 1);
				} else {
					z = __builtin_ctzl(n);
					t ^= (d ^ (d >> 1)) & (z << 1);
					n >>= z;
				}
			}
			return ((d == 1) && !(t & 2));
		}

		z = MIN(__builtin_clzl(t0[j]), __builtin_clzl(t1[j]));
		n = t0[j] << z;
		d = t1[j] << z;
		if (z > (RADIX >> 1)) {
			n |= t0[j - 1] >> z;
			d |= t1[j - 1] >> z;
		}
		n = (n & HMASK) | (t0[0] & LMASK);
		d = (d & HMASK) | (t1[0] & LMASK);

		i = s;
		while (i > 0) {
			if (n & 1) {
				if (n < d) {
					SWAP(ai, ci);
					SWAP(bi, di);
					SWAP(n, d);
					t ^= (n & d);
				}
				n = (n - d) >> 1;
				ai = ai - ci;
				bi = bi - di;
				ci = ci << 1;
				di = di << 1;
				t ^= d ^ (d >> 1);
				i -= 1;
			} else {
				z = MIN(i, __builtin_ctzl(n));
				t ^= (d ^ (d >> 1)) & (z << 1);
				ci = ci << z;
				di = di << z;
				n >>= z;
				i -= z;
			}
		}

#if ARITH_BROADWELL
		// If using the ASM backend, keep precision fixed to use CT functions.
		j = NWORDS_FIELD;
#else
		j = j + 1;
#endif
		t2[j] = bn_muls_low(t2, t0, 0, ai, j);
		t3[j] = bn_muls_low(t3, t1, 0, bi, j);
		bn_addn_low(t3, t3, t2, j + 1);

		t2[j] = bn_muls_low(t2, t0, 0, ci, j);
		t1[j] = bn_muls_low(t1, t1, 0, di, j);
		bn_addn_low(t1, t1, t2, j + 1);

		bn_rshs_low(t1, t1, j + 1, s);
		bn_rshs_low(t0, t3, j + 1, s);

		if (fp_is_zero(t0)) {
			return ((t1[0] == 1) && !(t & 2));
		}

		if (SIGN(t0[j])) {
			t ^= t1[0];
			bn_negs_low(t0, t0, SIGN(t0[j]), j + 1);
		}
		if (SIGN(t1[j])) {
			bn_negs_low(t1, t1, SIGN(t1[j]), j + 1);
		}
	}

	return t;
}
