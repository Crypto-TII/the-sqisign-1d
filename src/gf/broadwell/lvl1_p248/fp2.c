#include <fp2.h>

extern void fp2_sq_c0(fp2_t *out, const fp2_t *in);
extern void fp2_sq_c1(fp_t *out, const fp2_t *in);

extern void fp2_mul_c0(fp_t *out, const fp2_t *in0, const fp2_t *in1);
extern void fp2_mul_c1(fp_t *out, const fp2_t *in0, const fp2_t *in1);

/* Arithmetic modulo X^2 + 1 */

void fp2_mul(fp2_t* x, const fp2_t* y, const fp2_t* z)
{
    fp_t t;

    fp2_mul_c0(&t, y, z);              // c0 = a0*b0 - a1*b1
    fp2_mul_c1(&x->im, y, z);          // c1 = a0*b1 + a1*b0 
    x->re[0] = t[0]; x->re[1] = t[1]; x->re[2] = t[2]; x->re[3] = t[3];
}

void fp2_sqr(fp2_t* x, const fp2_t* y) {
    fp2_t t;

    fp2_sq_c0(&t, y);               // c0 = (a0+a1)(a0-a1)
    fp2_sq_c1(&x->im, y);           // c1 = 2a0*a1
    x->re[0] = t.re[0]; x->re[1] = t.re[1]; x->re[2] = t.re[2]; x->re[3] = t.re[3];
}
