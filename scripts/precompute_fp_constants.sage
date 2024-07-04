#!/usr/bin/env sage
proof.all(False)  # faster

from sage.misc.banner import require_version
if not require_version(10, 0, print_message=True):
    exit('')

################################################################

from parameters import p

################################################################

if p % 4 != 3:
    raise NotImplementedError('requires p ≡ 3 (mod 4)')

################################################################

if (p == 0x34e29e286b95d98c33a6a86587407437252c9e49355147ffffffffffffffffff):
    # lvl1
    radix_map = { 32: 29, 64: 52 }
if (p == 0x4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff):
    # lvl1_p248
    radix_map = { 32: 29, 64: 51 }
elif (p == 0x3df6eeeab0871a2c6ae604a45d10ad665bc2e0a90aeb751c722f669356ea4684c6174c1ffffffffffffffffffffffff):
    # lvl3
    radix_map = { 32: 28, 64: 55 }
elif (p == 0x255946a8869bc68c15b0036936e79202bdbe6326507d01fe3ac5904a0dea65faf0a29a781974ce994c68ada6e1ffffffffffffffffffffffffffffffffffff):
    # lvl5
    radix_map = { 32: 28, 64: 56 }

################################################################

from math import ceil, log

with open('include/fp_constants.h', 'w') as hfile:

    print('#ifndef FP_CONSTANTS_H', file=hfile)
    print('#define FP_CONSTANTS_H', file=hfile)
    print('',file=hfile)
    print('#if 0',file=hfile)
    for RADIX in (16,32,64):
        print(f'#elif 8*DIGIT_LEN == {RADIX}',file=hfile)
        print('',file=hfile)
        if RADIX == 16:
            print(f'#define NWORDS_FIELD {ceil(log(p,2**RADIX))}',file=hfile)
        else:
            if RADIX == 32:
                print(f'#if defined(ARITH_REF) || defined (ARITH_M4) || defined(ARITH_BROADWELL)', file=hfile)
            else:
                print(f'#if defined(ARITH_REF) || defined(ARITH_BROADWELL)', file=hfile)
            print(f'#define NWORDS_FIELD {ceil(log(p,2**RADIX))}',file=hfile)
            print(f'#elif defined(ARITH_MIKE)', file=hfile)
            print(f'#define NWORDS_FIELD {ceil(log(p,2**radix_map[RADIX]))}',file=hfile)
            print('#endif', file=hfile)
        print(f'#define NWORDS_ORDER {ceil(log(p+1,2**RADIX))}',file=hfile)
        print('',file=hfile)
    print('#endif', file=hfile)
    print('',file=hfile)
    print(f'#define BITS {RADIX*ceil(log(p,2**RADIX))}',file=hfile)
    print(f'#define LOG2P {ceil(log(p,2))}',file=hfile)
    print('',file=hfile)
    print('#endif', file=hfile)
