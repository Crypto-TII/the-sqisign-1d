#include <assert.h>
#include "include/fp.h"

#define PRECISION 254
#ifdef RADIX_32
const digit_t p[NWORDS_FIELD] =  { 0xffffffff, 0xffffffff, 0x355147FF, 0x252C9E49, 0x87407437, 0x33A6A865, 0x6B95D98C, 0x34E29E28 };
const digit_t R2[NWORDS_FIELD] = { 0x400674D4, 0x233625AE, 0x025A1C2E, 0x20AFD6C1, 0x0920655D, 0x30A841AB, 0x7C30CD3D, 0x0D72E7D6 };
const digit_t pp[NWORDS_FIELD] = { 0x00000001, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
const digit_t pre[NWORDS_FIELD] = { 0xB9713934, 0x14CF4F2D, 0x7482B63C, 0x61C30BED, 0xF05CE6CB, 0x3DA9F365, 0xA8AE3E69, 0x10B5541E };
#elif defined(RADIX_64)
const digit_t p[NWORDS_FIELD] =  { 0xFFFFFFFFFFFFFFFF, 0x252C9E49355147FF, 0x33A6A86587407437, 0x34E29E286B95D98C };
const digit_t R2[NWORDS_FIELD] = { 0x233625AE400674D4, 0x20AFD6C1025A1C2E, 0x30A841AB0920655D, 0x0D72E7D67C30CD3D };
const digit_t pp[NWORDS_FIELD] = { 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 };
const digit_t pre[NWORDS_FIELD] = { 0x2B0A392CB5069140, 0xC5D4CFE16AC471EE, 0xF3D3E4BEC2F7033F, 0x202B6A57AE960B4D };
#endif

void fp_set(digit_t* x, const digit_t val)
{ // Set field element x = val, where val has wordsize

    x[0] = val;
    for (unsigned int i = 1; i < NWORDS_FIELD; i++) {
        x[i] = 0;
    }
}

void fp_mont_setone(digit_t* out1) {
#ifdef RADIX_32
    out1[0] = 0x00000004;
    out1[1] = 0x00000000;
    out1[2] = 0x2abae000;
    out1[3] = 0x6b4d86db;
    out1[4] = 0xe2fe2f23;
    out1[5] = 0x31655e69;
    out1[6] = 0x51a899cf;
    out1[7] = 0x2c75875e;
#elif defined(RADIX_64)
    out1[0] = 0x4;
    out1[1] = UINT64_C(0x6b4d86db2abae000);
    out1[2] = UINT64_C(0x31655e69e2fe2f23);
    out1[3] = UINT64_C(0x2c75875e51a899cf);
#endif
}

bool fp_is_equal(const digit_t* a, const digit_t* b)
{ // Compare two field elements in constant time
  // Returns 1 (true) if a=b, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++)
        r |= a[i] ^ b[i];

    return (bool)is_digit_zero_ct(r);
}

bool fp_is_zero(const digit_t* a)
{ // Is a field element zero?
  // Returns 1 (true) if a=0, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++)
        r |= a[i] ^ 0;

    return (bool)is_digit_zero_ct(r);
}

void fp_copy(digit_t* out, const digit_t* a)
{
    memcpy(out, a, NWORDS_FIELD*RADIX/8);
}

void fp_neg(digit_t* out, const digit_t* a)
{ // Modular negation, out = -a mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(out[i], borrow, ((digit_t*)p)[i], a[i], borrow);
    }
    fp_sub(out, out, (digit_t*)p);
}

void fp_tomont(digit_t* out, const digit_t* a)
{ // Conversion to Montgomery representation
  // out = a*R^2*R^(-1) mod p = a*R mod p, where a in [0, p-1].

    fp_mul(out, a, (digit_t*)&R2);
}

void fp_frommont(digit_t* out, const digit_t* a)
{ // Conversion from Montgomery representation to standard representation
  // out = a*R^(-1) mod p, where a in [0, p-1].
    digit_t one[NWORDS_FIELD] = {0};

    one[0] = 1;
    fp_mul(out, a, one);
}

void MUL(digit_t* out, const digit_t a, const digit_t b)
{ // Digit multiplication, digit*digit -> 2-digit result 
  // Inputs: a, b in [0, 2^w-1], where w is the computer wordsize 
  // Output: 0 < out < 2^(2w)-1    
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), mask_high = (digit_t)(-1) << (sizeof(digit_t)*4);

    al = a & mask_low;                        // Low part
    ah = a >> (sizeof(digit_t)*4);            // High part
    bl = b & mask_low;
    bh = b >> (sizeof(digit_t)*4);

    albl = al * bl;
    albh = al * bh;
    ahbl = ah * bl;
    ahbh = ah * bh;
    out[0] = albl & mask_low;                 // out00

    res1 = albl >> (sizeof(digit_t)*4);
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;
    temp = res1 + res2 + res3;
    carry = temp >> (sizeof(digit_t)*4);
    out[0] ^= temp << (sizeof(digit_t)*4);    // out01   

    res1 = ahbl >> (sizeof(digit_t)*4);
    res2 = albh >> (sizeof(digit_t)*4);
    res3 = ahbh & mask_low;
    temp = res1 + res2 + res3 + carry;
    out[1] = temp & mask_low;                 // out10 
    carry = temp & mask_high;
    out[1] ^= (ahbh & mask_high) + carry;     // out11
}

digit_t mp_shiftr(digit_t* x, const unsigned int shift, const unsigned int nwords)
{ // Multiprecision right shift
    digit_t bit_out = x[0] & 1;

    for (unsigned int i = 0; i < nwords-1; i++) {
        SHIFTR(x[i+1], x[i], shift, x[i], RADIX);
    }
    x[nwords-1] >>= shift;
    return bit_out;
}

void mp_shiftl(digit_t* x, const unsigned int shift, const unsigned int nwords)
{ // Multiprecision left shift

    assert(shift < RADIX*nwords);

    int shift_words = shift / RADIX;

    for (int i = nwords-1; i > shift_words; i--) {
        SHIFTL(x[i-shift_words], x[i-shift_words-1], shift % RADIX, x[i], RADIX);
    }

    x[shift_words] = x[0] << (shift % RADIX);

    for (int i = shift_words - 1; i >= 0; i--) {
        x[i] = 0;
    }
}

static void fp_exp3div4(digit_t* out, const digit_t* a)
{ // Fixed exponentiation out = a^((p-3)/4) mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
  // Requirement: p = 3(mod 4)
    fp_t p_t, acc;
    digit_t bit;

    memcpy((digit_t*)p_t, (digit_t*)p, NWORDS_FIELD*RADIX/8);
    memcpy((digit_t*)acc, (digit_t*)a, NWORDS_FIELD*RADIX/8);
    mp_shiftr(p_t, 1, NWORDS_FIELD);
    mp_shiftr(p_t, 1, NWORDS_FIELD);
    fp_set(out, 1);
    fp_tomont(out, out);

    for (int i = 0; i < NWORDS_FIELD*RADIX-2; i++) {
        bit = p_t[0] & 1;
        mp_shiftr(p_t, 1, NWORDS_FIELD);
        if (bit == 1) {
            fp_mul(out, out, acc);
        }
        fp_sqr(acc, acc);
    }
}

void _fp_inv(digit_t* a)
{ // Modular inversion, out = x^-1*R mod p, where R = 2^(w*nwords), w is the computer wordsize and nwords is the number of words to represent p
  // Input: a=xR in [0, p-1] 
  // Output: out in [0, p-1]. It outputs 0 if the input does not have an inverse  
  // Requirement: Ceiling(Log(p)) < w*nwords
    fp_t t;

    fp_exp3div4(t, a);
    fp_sqr(t, t);
    fp_sqr(t, t);
    fp_mul(a, t, a);    // a^(p-2)
}

#include "../../inversion.inc"

bool _fp_is_square(const digit_t* a)
{ // Is field element a square?
  // Output: out = 0 (false), 1 (true)
    fp_t t, one;

    fp_exp3div4(t, a);
    fp_sqr(t, t);
    fp_mul(t, t, a);    // a^((p-1)/2)
    fp_frommont(t, t);
    fp_set(one, 1);

    return fp_is_equal(t, one);
}

#include "../../symbol.inc"

void fp_sqrt(digit_t* a)
{ // Square root computation, out = a^((p+1)/4) mod p
    fp_t t;

    fp_exp3div4(t, a);
    fp_mul(a, t, a);    // a^((p+1)/4)
}

void fp_to_digit_array(digit_t* out, const digit_t* a) {
    memcpy(out, a, NWORDS_FIELD*RADIX/8);
}

void fp_from_digit_array(digit_t* out, const digit_t* a) {
    memcpy(out, a, NWORDS_FIELD*RADIX/8);
}