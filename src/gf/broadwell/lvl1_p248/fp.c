#include <fp.h>

#ifdef RADIX_32
const digit_t p[NWORDS_FIELD]   = { 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x04ffffff };
const digit_t R2[NWORDS_FIELD]  = { 0x33333d70, 0x33333333, 0x33333333, 0x33333333, 0x33333333, 0x33333333, 0x33333333, 0x03333333 };
const digit_t pre[NWORDS_FIELD] = { 0xd441b778, 0x86b4b5bd, 0x69c1878d, 0xbea87c37, 0x34ef11ef, 0xba9109a9, 0xdd00b817, 0x00efabb8 };
#elif defined(RADIX_64)
const digit_t p[NWORDS_FIELD]   = { 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x04ffffffffffffff };
const digit_t R2[NWORDS_FIELD]  = { 0x3333333333333d70, 0x3333333333333333, 0x3333333333333333, 0x0333333333333333 };
const digit_t pre[NWORDS_FIELD] = { 0x49be0f19ce5d45ab, 0xd4f438a18723cb36, 0x8774530e19f41cf2, 0x01c670115ff9bac1 };
#endif

void fp_mont_setone(digit_t* out1) {
#ifdef RADIX_32
    out1[0] = 0x00000033;
    out1[1] = 0x00000000;
    out1[2] = 0x00000000;
    out1[3] = 0x00000000;
    out1[4] = 0x00000000;
    out1[5] = 0x00000000;
    out1[6] = 0x00000000;
    out1[7] = 0x01000000;
#elif defined(RADIX_64)
    out1[0] = 0x33;
    out1[1] = UINT64_C(0x0);
    out1[2] = UINT64_C(0x0);
    out1[3] = UINT64_C(0x100000000000000);
#endif
}

#include "../../generic/inversion.inc"
#include "../../generic/symbol.inc"
