#!/usr/bin/env sage
proof.all(False)  # faster

from sage.misc.banner import require_version
if not require_version(9, 8, print_message=True):
    exit('')

################################################################

from parameters import lvl, f, p

################################################################

logp = ceil(log(p, 2))
tors2part = (p+1).p_primary_part(2)
tors3part = (p+1).p_primary_part(3)

# Degree of the challenge isogeny
if(lvl == 1 ): secpar = 128 
elif(lvl == 3 ): secpar = 192 
elif(lvl == 5 ): secpar = 256 
else: raise NotImplementedError(f'security level {lvl} not found')

#XXX first load the constants from klpt_constants.h
import re
klpt_consts = dict()
for l in open('include/klpt_constants.h'):
    m = re.search(r'#define *([^ ]+) *([x0-9]+)$', l)
    if m:
        k,v = m.groups()
        klpt_consts[k] = int(v, 0)

defs = dict()

fp2sz = (logp + 63)//64*8 * 2
defs['FP2_ENCODED_BYTES'] = fp2sz
defs['EC_CURVE_ENCODED_BYTES'] = fp2sz  # just the A
defs['EC_POINT_ENCODED_BYTES'] = fp2sz  # just the x
defs['EC_BASIS_ENCODED_BYTES'] = 3 * defs['EC_POINT_ENCODED_BYTES']

defs['CHAIN_LENGTH'] = klpt_consts['SQISIGN_keygen_length']

defs['QUAT_ALG_ELEM_ENCODED_BITS'] = ceil(((logp/4) + klpt_consts['KLPT_keygen_length'])/2  +55)  #TODO FIXME figure this out XXX XXX
defs['QUAT_ALG_ELEM_ENCODED_BYTES'] = (defs['QUAT_ALG_ELEM_ENCODED_BITS'] + 7)//8
defs['ID2ISO_LONG_TWO_ISOG_ENCODED_BYTES'] = defs['CHAIN_LENGTH'] * (defs['EC_CURVE_ENCODED_BYTES'] + defs['EC_POINT_ENCODED_BYTES'] + 2)
defs['UNCOMPRESSED_LONG_TWO_ISOG_ENCODED_BYTES'] = defs['CHAIN_LENGTH'] * defs['EC_POINT_ENCODED_BYTES']

defs['ZIP_CHAIN_LEN'] = klpt_consts['SQISIGN_signing_length']
defs['ID2ISO_COMPRESSED_LONG_TWO_ISOG_ZIP_CHAIN_BYTES'] = (f + 7) // 8
defs['ID2ISO_COMPRESSED_LONG_TWO_ISOG_BYTES'] = defs['ZIP_CHAIN_LEN'] * defs['ID2ISO_COMPRESSED_LONG_TWO_ISOG_ZIP_CHAIN_BYTES'] + 1

defs['SIGNATURE_LEN'] = defs['ID2ISO_COMPRESSED_LONG_TWO_ISOG_BYTES'] + ((tors2part*tors3part).bit_length()+7)//8 + 1 + (tors2part.bit_length()+7)//8 + (tors3part.bit_length()+7)//8
defs['SMART_SIGNATURE_LEN'] = defs['ID2ISO_COMPRESSED_LONG_TWO_ISOG_BYTES'] + secpar//8 + secpar//8
defs['UNCOMPRESSED_SIGNATURE_LEN'] = defs['UNCOMPRESSED_LONG_TWO_ISOG_ENCODED_BYTES'] + defs['EC_CURVE_ENCODED_BYTES']
defs['PARALLEL_SIGNATURE_LEN'] = defs['UNCOMPRESSED_LONG_TWO_ISOG_ENCODED_BYTES'] + 3*defs['EC_CURVE_ENCODED_BYTES']
defs['COMPRESSED_PARALLEL_SIGNATURE_LEN'] = defs['SMART_SIGNATURE_LEN'] + 2*defs['EC_CURVE_ENCODED_BYTES']
defs['PUBLICKEY_BYTES'] = defs['EC_CURVE_ENCODED_BYTES']
defs['SECRETKEY_BYTES'] = defs['EC_CURVE_ENCODED_BYTES'] + 5*defs['QUAT_ALG_ELEM_ENCODED_BYTES'] + defs['EC_POINT_ENCODED_BYTES'] + defs['EC_BASIS_ENCODED_BYTES'] + defs['EC_BASIS_ENCODED_BYTES']

size_privkey = defs['SECRETKEY_BYTES']
size_pubkey = defs['PUBLICKEY_BYTES']
size_signature = defs['SIGNATURE_LEN']

algname = f'lvl{lvl}'

################################################################

with open('include/encoded_sizes.h','w') as hfile:
    for k,v in defs.items():
        v = ZZ(v)
        print(f'#define {k} {v}', file=hfile)

api = f'''
// SPDX-License-Identifier: Apache-2.0

#ifndef api_h
#define api_h

#include <stddef.h>

#define CRYPTO_SECRETKEYBYTES {size_privkey:4}
#define CRYPTO_PUBLICKEYBYTES {size_pubkey:4}
#define CRYPTO_BYTES          {size_signature:4}

#define CRYPTO_ALGNAME "{algname}"

#if defined(ENABLE_SIGN) || defined(PQM4)
int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk);

int
crypto_sign(unsigned char *sm, size_t *smlen,
            const unsigned char *m, size_t mlen,
            const unsigned char *sk);
#endif

int
crypto_sign_open(unsigned char *m, size_t *mlen,
                 const unsigned char *sm, size_t smlen,
                 const unsigned char *pk);

#endif /* api_h */
'''.strip()

with open(f'../../../nistapi/lvl{lvl}/api.h', 'w') as f:
    print(api, file=f)

