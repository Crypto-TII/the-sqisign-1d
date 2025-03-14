#include <fp2.h>

extern const digit_t R[NWORDS_FIELD];

/* Arithmetic modulo X^2 + 1 */

void fp2_set(fp2_t* x, const digit_t val)
{
    fp_set(x->re, val);
    fp_set(x->im, 0);
}

bool fp2_is_zero(const fp2_t* a)
{ // Is a GF(p^2) element zero?
  // Returns 1 (true) if a=0, 0 (false) otherwise

    return fp_is_zero(a->re) & fp_is_zero(a->im);
}

bool fp2_is_equal(const fp2_t* a, const fp2_t* b)
{ // Compare two GF(p^2) elements in constant time
  // Returns 1 (true) if a=b, 0 (false) otherwise

    return fp_is_equal(a->re, b->re) & fp_is_equal(a->im, b->im);
}

void fp2_copy(fp2_t* x, const fp2_t* y)
{
    fp_copy(x->re, y->re);
    fp_copy(x->im, y->im);
}

fp2_t fp2_non_residue()
{ // 2 + i is a quadratic non-residue for p1913
    fp_t one = {0};
    fp2_t res;

    one[0] = 1;
    fp_tomont(one, one);
    fp_add(res.re, one, one);
    fp_copy(res.im, one);
    return res;
}

void fp2_add(fp2_t* x, const fp2_t* y, const fp2_t* z)
{
    fp_add(x->re, y->re, z->re);
    fp_add(x->im, y->im, z->im);
}

void fp2_sub(fp2_t* x, const fp2_t* y, const fp2_t* z)
{
    fp_sub(x->re, y->re, z->re);
    fp_sub(x->im, y->im, z->im);
}

void fp2_neg(fp2_t* x, const fp2_t* y)
{
    fp_neg(x->re, y->re);
    fp_neg(x->im, y->im);
}

void fp2_mul(fp2_t* x, const fp2_t* y, const fp2_t* z)
{
    fp_t t0, t1;

    fp_add(t0, y->re, y->im);
    fp_add(t1, z->re, z->im);
    fp_mul(t0, t0, t1);
    fp_mul(t1, y->im, z->im);
    fp_mul(x->re, y->re, z->re);
    fp_sub(x->im, t0, t1);
    fp_sub(x->im, x->im, x->re);
    fp_sub(x->re, x->re, t1);
}

void fp2_sqr(fp2_t* x, const fp2_t* y)
{
    fp_t sum, diff;

    fp_add(sum, y->re, y->im);
    fp_sub(diff, y->re, y->im);
    fp_mul(x->im, y->re, y->im);
    fp_add(x->im, x->im, x->im);
    fp_mul(x->re, sum, diff);
}

void fp2_inv(fp2_t* x)
{
    fp_t t0, t1;

    fp_sqr(t0, x->re);
    fp_sqr(t1, x->im);
    fp_add(t0, t0, t1);
    fp_inv(t0);
    fp_mul(x->re, x->re, t0);
    fp_mul(x->im, x->im, t0);
    fp_neg(x->im, x->im);
}

bool fp2_is_square(const fp2_t* x)
{
    fp_t t0, t1;

    fp_sqr(t0, x->re);
    fp_sqr(t1, x->im);
    fp_add(t0, t0, t1);

    return fp_is_square(t0);
}

void fp2_frob(fp2_t* x, const fp2_t* y)
{
    memcpy((digit_t*)x->re, (digit_t*)y->re, NWORDS_FIELD*RADIX/8);
    fp_neg(x->im, y->im);
}

void fp2_tomont(fp2_t* x, const fp2_t* y)
{
    fp_tomont(x->re, y->re);
    fp_tomont(x->im, y->im);
}

void fp2_frommont(fp2_t* x, const fp2_t* y)
{
    fp_frommont(x->re, y->re);
    fp_frommont(x->im, y->im);
}

void fp2_sqrt(fp2_t* x)
{
    fp_t t, two, sdelta, re, im, tmp;

    if (fp_is_zero(x->im)) {
        fp_copy(re, x->re);
        fp_sqrt(re);
        fp_sqr(t, re);
        fp_frommont(t, t);
        fp_frommont(tmp, x->re);
        if (fp_is_equal(tmp, t)) {
            memcpy((digit_t*)x->re, (digit_t*)re, NWORDS_FIELD*RADIX/8);
            return;
        } else {
            memcpy((digit_t*)x->im, (digit_t*)re, NWORDS_FIELD*RADIX/8);
            fp_set(x->re, 0);
            return;
        }
    }

    // sdelta = sqrt(re^2 + im^2)
    fp_sqr(tmp, x->re);
    fp_sqr(t, x->im);
    fp_add(tmp, tmp, t);
    fp_copy(sdelta, tmp);
    fp_sqrt(sdelta);   // (re^2 + im^2)^((p+1)/4)

    fp_set(two, 2);
    fp_tomont(two, two);

    fp_add(re, x->re, sdelta);
    fp_mul(tmp, re, two);    // 2(x->re+sdelta)
    fp_exp3div4(t, tmp);

    fp_mul(re, t, re);       // (x->re+sdelta) * (2(x->re+sdelta))^((p-3)/4)
    fp_mul(im, t, x->im);    //          x->im * (2(x->re+sdelta))^((p-3)/4)

    fp_mul(t, re, two);      // (2(x->re+sdelta))^((p+1)/4)
    fp_sqr(t, t);
    fp_frommont(t, t);
    fp_frommont(tmp, tmp);
    if (fp_is_equal(tmp, t)) {
        memcpy((digit_t*)x->re, (digit_t*)re, NWORDS_FIELD*RADIX/8);
        memcpy((digit_t*)x->im, (digit_t*)im, NWORDS_FIELD*RADIX/8);
    } else {
#if defined(PRIME_P248)
        memcpy((digit_t*)x->re, (digit_t*)im, NWORDS_FIELD*RADIX/8);
        fp_neg(re, re);
        memcpy((digit_t*)x->im, (digit_t*)re, NWORDS_FIELD*RADIX/8);
#else
        if (fp_is_square(x->im)) {
            memcpy((digit_t*)x->im, (digit_t*)re, NWORDS_FIELD*RADIX/8);
            fp_neg(im, im);
            memcpy((digit_t*)x->re, (digit_t*)im, NWORDS_FIELD*RADIX/8);
        } else {
            memcpy((digit_t*)x->re, (digit_t*)im, NWORDS_FIELD*RADIX/8);
            fp_neg(re, re);
            memcpy((digit_t*)x->im, (digit_t*)re, NWORDS_FIELD*RADIX/8);
        }
#endif
    }
}

// Lexicographic comparison of two field elements. Returns +1 if x > y, -1 if x < y, 0 if x = y
int fp2_cmp(fp2_t* x, fp2_t* y){
    fp2_t a, b;
    fp2_frommont(&a, x);
    fp2_frommont(&b, y);
    for(int i = NWORDS_FIELD-1; i >= 0; i--){
        if(a.re[i] > b.re[i])
            return 1;
        if(a.re[i] < b.re[i])
            return -1;
    }
    for(int i = NWORDS_FIELD-1; i >= 0; i--){
        if(a.im[i] > b.im[i])
            return 1;
        if(a.im[i] < b.im[i])
            return -1;
    }
    return 0;
}