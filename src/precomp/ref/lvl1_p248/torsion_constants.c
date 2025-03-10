#include <stddef.h>
#include <stdint.h>
#include <torsion_constants.h>
#if 0
#elif 8*DIGIT_LEN == 16
const digit_t TORSION_PLUS_EVEN_POWER = 0xf8;
const digit_t TORSION_ODD_PRIMES[2] = {0x5, 0x3};
const digit_t TORSION_ODD_POWERS[2] = {0x1, 0x4};
const size_t TORSION_PLUS_ODD_POWERS[1] = {0x1};
const digit_t TORSION_PLUS_2POWER_DIGITS[NWORDS_ORDER] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x100};
const digit_t TORSION_PLUS_3POWER_DIGITS[NWORDS_ORDER] = {0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
#elif 8*DIGIT_LEN == 32
const digit_t TORSION_PLUS_EVEN_POWER = 0xf8;
const digit_t TORSION_ODD_PRIMES[2] = {0x5, 0x3};
const digit_t TORSION_ODD_POWERS[2] = {0x1, 0x4};
const size_t TORSION_PLUS_ODD_POWERS[1] = {0x1};
const digit_t TORSION_PLUS_2POWER_DIGITS[NWORDS_ORDER] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x1000000};
const digit_t TORSION_PLUS_3POWER_DIGITS[NWORDS_ORDER] = {0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
#elif 8*DIGIT_LEN == 64
const digit_t TORSION_PLUS_EVEN_POWER = 0xf8;
const digit_t TORSION_ODD_PRIMES[2] = {0x5, 0x3};
const digit_t TORSION_ODD_POWERS[2] = {0x1, 0x4};
const size_t TORSION_PLUS_ODD_POWERS[1] = {0x1};
const digit_t TORSION_PLUS_2POWER_DIGITS[NWORDS_ORDER] = {0x0, 0x0, 0x0, 0x100000000000000};
const digit_t TORSION_PLUS_3POWER_DIGITS[NWORDS_ORDER] = {0x1, 0x0, 0x0, 0x0};
#endif
#if defined(ENABLE_SIGN)
#if 0
#elif 8*DIGIT_LEN == 16
const digit_t TORSION_PLUS_ODD_PRIMES[1] = {0x5};
const digit_t TORSION_MINUS_ODD_PRIMES[1] = {0x3};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x4};
const size_t DEGREE_COMMITMENT_POWERS[2] = {0x1, 0x0};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0x4ff}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x195}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[2] = {{{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x51}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x51}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x100}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x100}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x100}}};
#elif 8*DIGIT_LEN == 32
const digit_t TORSION_PLUS_ODD_PRIMES[1] = {0x5};
const digit_t TORSION_MINUS_ODD_PRIMES[1] = {0x3};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x4};
const size_t DEGREE_COMMITMENT_POWERS[2] = {0x1, 0x0};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0x4ffffff}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x195}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[2] = {{{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x51}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x51}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x1000000}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x1000000}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x1000000}}};
#elif 8*DIGIT_LEN == 64
const digit_t TORSION_PLUS_ODD_PRIMES[1] = {0x5};
const digit_t TORSION_MINUS_ODD_PRIMES[1] = {0x3};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x4};
const size_t DEGREE_COMMITMENT_POWERS[2] = {0x1, 0x0};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xffffffffffffffff,0xffffffffffffffff,0xffffffffffffffff,0x4ffffffffffffff}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x195}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[2] = {{{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x51}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x51}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x100000000000000}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x100000000000000}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x5}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x1}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x100000000000000}}};
#endif
#endif
