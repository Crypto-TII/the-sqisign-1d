#!/usr/bin/env sage
proof.all(False)  # faster

from sage.misc.banner import require_version
if not require_version(10, 0, print_message=True):
    exit('')
    
################################################################

# The following benchmarks are required for computing optimal strategies.
# The costs can be specified in arbitrary units and can be obtained by
# runing src/ec/ref/lvl1/test/mont.test or left as "None" to use a
# generic estimate.
# ++++++++++ 
# p2 =  None   # cost of xdbl
# q4 = None    # cost of xeval4

##SqiSign lvl 1 costs
p2 = 1766
q4 = 2452

##SqiSign lvl 3 costs
# p2 = 4071
# q4 = 5224

##SqiSign lvl 5 costs
# p2 = 7755
# q4 = 9847
# ++++++++++

################################################################

from parameters import p, f, Tpls, Tmin
g = valuation(Tpls, 3)
T = Tpls * Tmin
Lpls = [x[0] for x in factor(Tpls)]
Lmin = [x[0] for x in factor(Tmin)]
Epls = [x[1] for x in factor(Tpls)]
Emin = [x[1] for x in factor(Tmin)]
L = Lpls + Lmin
E = Epls + Emin
################################################################

if p % 4 != 3:
    raise NotImplementedError('requires p ≡ 3 (mod 4)')
if 3 not in Lpls or 3 in Lmin:
    raise NotImplementedError('power of 3 must be in the positive torsion')

################################################################

from math import log, ceil, floor, sqrt

def strategy(n, p, q):
    S = { 1: [] }
    C = { 1: 0 }
    for i in range(2, n+1):
        b, cost = min(((b, C[i-b] + C[b] + b*p + (i-b)*q) for b in range(1,i)),
                      key=lambda t: t[1])
        S[i] = [b] + S[i-b] + S[b]
        C[i] = cost
    return S[n]

def ijk(l):
    i,j,k = ceil(sqrt((l-1)/4)), floor((l-1)/4/ceil(sqrt((l-1)/4))), int((l-1)/2 - 2*ceil(sqrt((l-1)/4))*floor((l-1)/4/ceil(sqrt((l-1)/4))))
    assert i < 2**16 and j < 2**16 and k < 2**16
    return i,j,k

def dac_search(target,r0,r1,r2,chain,chainlen,best,bestlen):
  if chainlen >= bestlen: return best,bestlen
  if r2 > target: return best,bestlen
  if r2<<(bestlen-1-chainlen) < target: return best,bestlen
  if r2 == target: return chain,chainlen
  chain *= 2
  chainlen += 1
  best,bestlen = dac_search(target,r0,r2,r0+r2,chain+1,chainlen,best,bestlen)
  best,bestlen = dac_search(target,r1,r2,r1+r2,chain,chainlen,best,bestlen)
  return best,bestlen

def dac(target):
  best = None
  bestlen = 0
  while best == None:
    bestlen += 1
    best,bestlen = dac_search(target,1,2,3,0,0,best,bestlen)
  return best,bestlen

def print_dac(dac, name, len, file):
    print(f'static char {name}[{max(len,1)}] = '+'\"', end='', file=file)
    for bit in dac: print(f'{bit}', end='', file=file)
    if len==0: print('0', end='', file=file)
    print('\";', file=file)

bL = '{'+', '.join([str(ceil(log(l,2))) for l in L])+'}'
strategy4 = strategy(f//2-1, 2*p2, q4)
strategy4 = '{'+', '.join([str(s) for s in strategy4])+'}'

## The following primes have IJK that have been optimized externally by an exhaustive dac_search
if p == 0x34e29e286b95d98c33a6a86587407437252c9e49355147ffffffffffffffffff: # Nist lvl1
    sizeI = [0, 2, 4, 6, 6, 7, 12, 14, 28, 1, 2, 3, 3, 5, 6, 6, 6, 6, 9, 8, 10, 12, 12, 12, 16, 16, 18, 22]
    sizeJ = [0, 2, 3, 4, 4, 7, 10, 13, 17, 1, 1, 1, 3, 4, 4, 4, 5, 5, 6, 7, 9, 8, 10, 12, 16, 16, 16, 22]
    sizeK = [1, 3, 5, 2, 6, 0, 5, 7, 4, 1, 1, 0, 0, 4, 0, 5, 5, 8, 3, 7, 11, 2, 9, 15, 4, 12, 20, 18]
elif p == 0x3df6eeeab0871a2c6ae604a45d10ad665bc2e0a90aeb751c722f669356ea4684c6174c1ffffffffffffffffffffffff: # Nist lvl3
    sizeI = [0, 1, 2, 3, 6, 14, 14, 1, 3, 5, 5, 7, 8, 8, 10, 11, 14, 19, 26, 30, 35, 40, 44, 44, 52, 88, 114, 114]
    sizeJ = [0, 1, 1, 3, 6, 9, 13, 1, 1, 4, 5, 6, 7, 7, 6, 10, 10, 16, 23, 28, 32, 32, 32, 33, 45, 76, 87, 104]
    sizeK = [1, 1, 1, 5, 6, 2, 16, 0, 0, 4, 6, 2, 4, 7, 0, 1, 4, 6, 0, 5, 18, 13, 30, 2, 18, 12, 3, 8]
elif p == 0x255946a8869bc68c15b0036936e79202bdbe6326507d01fe3ac5904a0dea65faf0a29a781974ce994c68ada6e1ffffffffffffffffffffffffffffffffffff: # Nist lvl5
    sizeI = [0, 3, 3, 4, 6, 11, 14, 1, 1, 2, 3, 4, 4, 6, 6, 13, 14, 14, 18, 26, 32, 38, 48, 54, 60, 62, 62, 64, 72, 92, 118, 126, 130, 310]
    sizeJ = [0, 1, 3, 4, 4, 10, 14, 1, 1, 2, 3, 3, 4, 5, 6, 12, 13, 13, 14, 24, 32, 33, 37, 49, 56, 60, 62, 62, 64, 63, 102, 125, 126, 256]
    sizeK = [1, 0, 2, 1, 3, 10, 21, 0, 1, 0, 0, 2, 4, 3, 3, 9, 2, 5, 0, 21, 28, 21, 11, 6, 75, 21, 82, 59, 75, 21, 21, 123, 0, 396]

## For other primes we perform a new good-enough search:
else:
    sizeI = []
    sizeJ = []
    sizeK = []
    for l in L:
        I,J,K = ijk(l)
        sizeI.append(I)
        sizeJ.append(J)
        sizeK.append(K)

sImax = max(sizeI)
sJmax = max(sizeJ)
sKmax = max(sizeK)
sizeI = '{'+', '.join([str(s) for s in sizeI])+'}'
sizeJ = '{'+', '.join([str(s) for s in sizeJ])+'}'
sizeK = '{'+', '.join([str(s) for s in sizeK])+'}'

# Small quadratic non-residues
Fp2.<i> = GF((p,2), modulus=[1,0,1])
NONRES_LEN = 128
NONRES = []
k = 0
while len(NONRES) < NONRES_LEN:
    k+=1
    if not (1+k*i).is_square():
        NONRES.append(k)

################################################################

FpEls = {
    'TWOpFm1': 2**(f-1),
    'THREEpE': 3**(g//2),
    'THREEpFdiv2': (3**g)//2,
    'p_cofactor_for_2f': (p+1)//(2**f),
    'p_cofactor_for_3g': (p+1)//(3**g),
    'p_cofactor_for_6fg': (p+1)//(2**f*3**g)
}

DACS = [dac(l) for l in L]

from cformat import FpEl, Object, ObjectFormatter
objs = ObjectFormatter(
    [Object('digit_t', f'{k}[NWORDS_ORDER]', FpEl(v, p))
    for k,v in FpEls.items()]
    )

################################################################

with open('include/ec_params.h', 'w') as hfile:
    with open('ec_params.c','w') as cfile:

        print('#include <ec_params.h>', file=cfile)
        print('',file=cfile)

        print('#ifndef EC_PARAMS_H', file=hfile)
        print('#define EC_PARAMS_H', file=hfile)
        print('',file=hfile)
        print('#include <tutil.h>', file=hfile)
        print('#include <fp_constants.h>', file=hfile)
        print('',file=hfile)
        print(f'#define POWER_OF_2 {f}', file=hfile)
        print(f'#define POWER_OF_3 {g}', file=hfile)
        print('',file=hfile)
        print('#define scaled 1', file=hfile)
        print('#define gap 83', file=hfile)
        print('',file=hfile)
        print(f'#define P_LEN {len(Lpls)}', file=hfile)
        print(f'#define M_LEN {len(Lmin)}', file=hfile)
        print('',file=hfile)

        print('const digit_t p_plus_minus_bitlength[P_LEN + M_LEN] =\n\t'+bL+';', file=cfile)
        print('',file=cfile)
        print(f'const digit_t STRATEGY4[] =\n\t'+strategy4+';', file=cfile)
        print('',file=cfile)
        print(f'const digit_t sizeI[] =\n\t'+sizeI+';', file=cfile)
        print(f'const digit_t sizeJ[] =\n\t'+sizeJ+';', file=cfile)
        print(f'const digit_t sizeK[] =\n\t'+sizeK+';', file=cfile)
        print('',file=cfile)

        print('extern const digit_t p_plus_minus_bitlength[P_LEN + M_LEN];', file=hfile)
        print('',file=hfile)
        print(f'extern const digit_t STRATEGY4[];', file=hfile)
        print('',file=hfile)
        print(f'extern const digit_t sizeI[];', file=hfile)
        print(f'extern const digit_t sizeJ[];', file=hfile)
        print(f'extern const digit_t sizeK[];', file=hfile)
        print('',file=hfile)

        print(f'#define sI_max {sImax}', file=hfile)
        print(f'#define sJ_max {sJmax}', file=hfile)
        print(f'#define sK_max {max(sKmax, 41)}', file=hfile)
        print('',file=hfile)
        print(f'#define ceil_log_sI_max {ceil(log(sImax,2))}', file=hfile)
        print(f'#define ceil_log_sJ_max {ceil(log(sJmax,2))}', file=hfile)
        print('',file=hfile)

        objs.implementation(file=cfile)
        print('',file=cfile)
        objs.header(file=hfile)
        print('',file=hfile)

        print(f'#define P_COFACTOR_FOR_2F_BITLENGTH {ceil(log((p+1)//(2**f),2))}', file=hfile)
        print(f'#define P_COFACTOR_FOR_3G_BITLENGTH {ceil(log((p+1)//(3**g),2))}', file=hfile)
        print(f'#define P_COFACTOR_FOR_6FG_BITLENGTH {ceil(log((p+1)//(2**f*3**g),2))}', file=hfile)
        print('',file=hfile)

        print('// differential addition chains', file= hfile)
        print(f'extern const digit_t DACS[{len(DACS)}];', file=hfile)
        print(f'extern const int DAC_LEN[{len(DACS)}];', file=hfile)
        print('',file=hfile)

        print('// differential addition chains', file= cfile)
        print(f'const digit_t DACS[{len(DACS)}] = '+'{', end= '', file=cfile)
        for dac in DACS:
            print(f'{dac[0]}, ', end= '', file=cfile)
        print('};', file=cfile)
        print(f'const int DAC_LEN[{len(DACS)}] = '+'{', end= '', file=cfile)
        for dac in DACS:
            print(f'{dac[1]}, ', end= '', file=cfile)
        print('};', file=cfile)
        print('',file=cfile)

        print('//quadratic residues',file=hfile)
        print(f'#define NONRES_LEN {NONRES_LEN}',file=hfile)
        print('extern const digit_t NONRES[NONRES_LEN];',file=hfile)
        print('//quadratic residues',file=cfile)
        print('const digit_t NONRES[NONRES_LEN] = { '+', '.join([str(x) for x in NONRES])+' };',file=cfile)

        print('',file=hfile)
        print('//quadratic residues',file=hfile)
        print(f'#define NONRES_LEN {NONRES_LEN}',file=hfile)
        print('extern const digit_t NONRES[NONRES_LEN];',file=hfile)
        print('',file=cfile)
        print('//quadratic residues',file=cfile)
        print('const digit_t NONRES[NONRES_LEN] = { '+', '.join([str(x) for x in NONRES])+' };',file=cfile)

        print('',file=hfile)
        print('#endif', file=hfile)
