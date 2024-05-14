#include <ec_params.h>

const digit_t p_plus_minus_bitlength[P_LEN + M_LEN] =
	{2, 4, 6, 7, 7, 9, 10, 3, 3, 5, 6, 6, 7, 7, 8, 10, 10, 10, 10, 12, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 16, 16, 16, 19};

const digit_t STRATEGY4[] =
	{32, 17, 9, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1};

const digit_t sizeI[] =
	{0, 3, 3, 4, 6, 11, 14, 1, 1, 2, 3, 4, 4, 6, 6, 13, 14, 14, 18, 26, 32, 38, 48, 54, 60, 62, 62, 64, 72, 92, 118, 126, 130, 310};
const digit_t sizeJ[] =
	{0, 1, 3, 4, 4, 10, 14, 1, 1, 2, 3, 3, 4, 5, 6, 12, 13, 13, 14, 24, 32, 33, 37, 49, 56, 60, 62, 62, 64, 63, 102, 125, 126, 256};
const digit_t sizeK[] =
	{1, 0, 2, 1, 3, 10, 21, 0, 1, 0, 0, 2, 4, 3, 3, 9, 2, 5, 0, 21, 28, 21, 11, 6, 75, 21, 82, 59, 75, 21, 21, 123, 0, 396};

#if 0
#elif 8*DIGIT_LEN == 16
const digit_t TWOpFm1[NWORDS_ORDER] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t THREEpE[NWORDS_ORDER] = {0xc6d1, 0x8b91, 0x3e46, 0x215, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t THREEpFdiv2[NWORDS_ORDER] = {0xfb50, 0xe494, 0xcb9a, 0x8f33, 0xb608, 0x3073, 0x2b5e, 0x2, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t p_cofactor_for_2f[NWORDS_ORDER] = {0xd371, 0x3456, 0x4ca6, 0xba67, 0x3c0c, 0x514d, 0xfd78, 0xf532, 0x2506, 0x62c8, 0xff1d, 0x3e80, 0x9328, 0xdf31, 0x15e, 0x73c9, 0xb49b, 0xd801, 0x460a, 0x4de3, 0x5443, 0xaca3, 0x12, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t p_cofactor_for_3g[NWORDS_ORDER] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x75a2, 0x9d68, 0xf46, 0x3e00, 0x3ad4, 0xfa62, 0x10ca, 0x8f78, 0x32c5, 0xe115, 0x27f3, 0xc006, 0xda1d, 0x8880, 0x9ba9, 0x8, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t p_cofactor_for_6fg[NWORDS_ORDER] = {0x3ad1, 0x4eb4, 0x7a3, 0x1f00, 0x1d6a, 0x7d31, 0x865, 0xc7bc, 0x9962, 0xf08a, 0x13f9, 0xe003, 0x6d0e, 0xc440, 0x4dd4, 0x4, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
#elif 8*DIGIT_LEN == 32
const digit_t TWOpFm1[NWORDS_ORDER] = {0x0, 0x0, 0x0, 0x0, 0x10000, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t THREEpE[NWORDS_ORDER] = {0x8b91c6d1, 0x2153e46, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t THREEpFdiv2[NWORDS_ORDER] = {0xe494fb50, 0x8f33cb9a, 0x3073b608, 0x22b5e, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t p_cofactor_for_2f[NWORDS_ORDER] = {0x3456d371, 0xba674ca6, 0x514d3c0c, 0xf532fd78, 0x62c82506, 0x3e80ff1d, 0xdf319328, 0x73c9015e, 0xd801b49b, 0x4de3460a, 0xaca35443, 0x12, 0x0, 0x0, 0x0, 0x0};
const digit_t p_cofactor_for_3g[NWORDS_ORDER] = {0x0, 0x0, 0x0, 0x0, 0x75a20000, 0xf469d68, 0x3ad43e00, 0x10cafa62, 0x32c58f78, 0x27f3e115, 0xda1dc006, 0x9ba98880, 0x8, 0x0, 0x0, 0x0};
const digit_t p_cofactor_for_6fg[NWORDS_ORDER] = {0x4eb43ad1, 0x1f0007a3, 0x7d311d6a, 0xc7bc0865, 0xf08a9962, 0xe00313f9, 0xc4406d0e, 0x44dd4, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
#elif 8*DIGIT_LEN == 64
const digit_t TWOpFm1[NWORDS_ORDER] = {0x0, 0x0, 0x10000, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t THREEpE[NWORDS_ORDER] = {0x2153e468b91c6d1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t THREEpFdiv2[NWORDS_ORDER] = {0x8f33cb9ae494fb50, 0x22b5e3073b608, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
const digit_t p_cofactor_for_2f[NWORDS_ORDER] = {0xba674ca63456d371, 0xf532fd78514d3c0c, 0x3e80ff1d62c82506, 0x73c9015edf319328, 0x4de3460ad801b49b, 0x12aca35443, 0x0, 0x0};
const digit_t p_cofactor_for_3g[NWORDS_ORDER] = {0x0, 0x0, 0xf469d6875a20000, 0x10cafa623ad43e00, 0x27f3e11532c58f78, 0x9ba98880da1dc006, 0x8, 0x0};
const digit_t p_cofactor_for_6fg[NWORDS_ORDER] = {0x1f0007a34eb43ad1, 0xc7bc08657d311d6a, 0xe00313f9f08a9962, 0x44dd4c4406d0e, 0x0, 0x0, 0x0, 0x0};
#endif

// differential addition chains
const digit_t DACS[34] = {0, 0, 42, 72, 48, 3472, 6800, 0, 2, 10, 48, 106, 20, 417, 180, 2692, 1536, 6484, 4618, 10336, 40977, 100898, 68224, 10336, 271530, 279681, 336544, 344148, 853014, 672065, 1049610, 2884160, 2232464, 11536384, };
const int DAC_LEN[34] = {0, 3, 6, 7, 8, 12, 13, 1, 2, 4, 6, 7, 7, 9, 9, 12, 12, 13, 13, 15, 16, 17, 17, 18, 19, 19, 19, 19, 20, 20, 21, 22, 22, 25, };

//quadratic residues
const digit_t NONRES[NONRES_LEN] = { 12, 17, 22, 24, 26, 34, 35, 36, 37, 40, 44, 45, 46, 49, 51, 53, 55, 56, 58, 59, 60, 61, 62, 63, 64, 65, 66, 69, 70, 71, 74, 78, 82, 85, 86, 89, 90, 94, 99, 100, 101, 105, 106, 107, 109, 115, 116, 118, 119, 121, 122, 123, 125, 128, 130, 131, 133, 134, 136, 138, 139, 140, 143, 147, 148, 154, 155, 156, 157, 159, 161, 162, 164, 165, 166, 167, 171, 172, 174, 176, 177, 178, 180, 181, 187, 189, 190, 191, 193, 194, 196, 197, 198, 201, 212, 213, 215, 216, 217, 218, 220, 221, 223, 226, 228, 229, 230, 233, 235, 237, 238, 240, 244, 246, 250, 252, 254, 257, 258, 259, 260, 262, 264, 266, 267, 272, 273, 274 };
