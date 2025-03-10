#include <fp.h>

#ifdef RADIX_32
const digit_t p[NWORDS_FIELD] =  { 0xffffffff, 0xffffffff, 0x355147FF, 0x252C9E49, 0x87407437, 0x33A6A865, 0x6B95D98C, 0x34E29E28 };
const digit_t R2[NWORDS_FIELD] = { 0x400674D4, 0x233625AE, 0x025A1C2E, 0x20AFD6C1, 0x0920655D, 0x30A841AB, 0x7C30CD3D, 0x0D72E7D6 };
const digit_t pre[NWORDS_FIELD] = { 0xB9713934, 0x14CF4F2D, 0x7482B63C, 0x61C30BED, 0xF05CE6CB, 0x3DA9F365, 0xA8AE3E69, 0x10B5541E };
#elif defined(RADIX_64)
const digit_t p[NWORDS_FIELD] =  { 0xFFFFFFFFFFFFFFFF, 0x252C9E49355147FF, 0x33A6A86587407437, 0x34E29E286B95D98C };
const digit_t R2[NWORDS_FIELD] = { 0x233625AE400674D4, 0x20AFD6C1025A1C2E, 0x30A841AB0920655D, 0x0D72E7D67C30CD3D };
const digit_t pre[NWORDS_FIELD] = { 0x2B0A392CB5069140, 0xC5D4CFE16AC471EE, 0xF3D3E4BEC2F7033F, 0x202B6A57AE960B4D };
#endif

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

#include "../../generic/inversion.inc"
#include "../../generic/symbol.inc"
