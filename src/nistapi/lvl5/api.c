// SPDX-License-Identifier: Apache-2.0

#include <api.h>
#include <sig.h>

#if defined(ENABLE_SIGN)
int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {

    return sqisign_keypair(pk, sk);
}

int
crypto_sign(unsigned char *sm, size_t *smlen,
            const unsigned char *m, size_t mlen,
            const unsigned char *sk) {
    return sqisign_sign(sm, smlen, m, mlen, sk);
}
#endif

int
crypto_sign_open(unsigned char *m, size_t *mlen,
                 const unsigned char *sm, size_t smlen,
                 const unsigned char *pk) {
    return sqisign_open(m, mlen, sm, smlen, pk);
}
