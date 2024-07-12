#!/bin/zsh

# This script should be run in the root folder of the repository.

declare -A FIAT_CRYPTO_SOURCES
declare -A FIAT_CRYPTO_FP2_SOURCES
declare -A M4_ARITHMETIC_SOURCES
declare -A PRECOMP_GF_DIRS

FIAT_CRYPTO_SOURCES[1]="fp_p1913_32.c"
FIAT_CRYPTO_SOURCES[1_p248_smart]="fp_p248_32.c"
FIAT_CRYPTO_SOURCES[1_p248_uncompressed]="fp_p248_32.c"
FIAT_CRYPTO_SOURCES[3]="fp_p47441_32.c"
FIAT_CRYPTO_SOURCES[5]="fp_p318233_32.c"

FIAT_CRYPTO_FP2_SOURCES[1]="fp2_p1913.c"
FIAT_CRYPTO_FP2_SOURCES[1_p248_smart]="fp2_p248.c"
FIAT_CRYPTO_FP2_SOURCES[1_p248_uncompressed]="fp2_p248.c"
FIAT_CRYPTO_FP2_SOURCES[3]="fp2_p47441.c"
FIAT_CRYPTO_FP2_SOURCES[5]="fp2_p318233.c"

M4_ARITHMETIC_SOURCES[1]="fp_p1913_32_m4.c"
M4_ARITHMETIC_SOURCES[1_p248_smart]="fp_p248_32_m4.c"
M4_ARITHMETIC_SOURCES[1_p248_uncompressed]="fp_p248_32_m4.c"
M4_ARITHMETIC_SOURCES[3]="fp_p47441_32_m4.c"
M4_ARITHMETIC_SOURCES[5]="fp_p318233_32_m4.c"

PRECOMP_GF_DIRS[1]="lvl1"
PRECOMP_GF_DIRS[1_p248_smart]="lvl1_p248"
PRECOMP_GF_DIRS[1_p248_uncompressed]="lvl1_p248"
PRECOMP_GF_DIRS[3]="lvl3"
PRECOMP_GF_DIRS[5]="lvl5"

for LEVEL in 1 1_p248_smart 1_p248_uncompressed 3 5
do
    for ARITHMETIC in fiat_crypto mikes_arithmetic m4_arithmetic
    do
        LVL=lvl${LEVEL}
        DST_PATH=sqisign${LEVEL}/${ARITHMETIC}

        if [ -d ${DST_PATH} ]; then
            echo Destination folder already exists. Delete it before running this script. Aborting.
            exit 1
        fi

        mkdir -p ${DST_PATH}

        cp include/*.h ${DST_PATH}/
        cp src/sqisign.c ${DST_PATH}/

        cp pqm4/${LVL}/pqm4_api.c ${DST_PATH}/
        if [ ${LEVEL} = "1_p248_smart" ] || [ ${LEVEL} = "1_p248_uncompressed" ]; then
            cp pqm4/${LVL}/api.h ${DST_PATH}/
        else
            cp src/nistapi/${LVL}/api.h ${DST_PATH}/
        fi

        cp src/common/generic/{mem.c,include/tutil.h} ${DST_PATH}/

        for FILE in include/ec_params.h include/encoded_sizes.h include/fp_constants.h include/klpt_constants.h \
                    include/torsion_constants.h ec_params.c torsion_constants.c
        do
            cp src/precomp/ref/${PRECOMP_GF_DIRS[${LEVEL}]}/${FILE} ${DST_PATH}/
        done

        if [ ${ARITHMETIC} = "fiat_crypto" ]; then
            ARITH_FLAG="ARITH_REF"

            cp src/gf/generic/{bn.inc,inversion.inc,symbol.inc} ${DST_PATH}/
            cp src/gf/ref/include/*.h ${DST_PATH}/
            cp src/gf/ref/${PRECOMP_GF_DIRS[${LEVEL}]}/${FIAT_CRYPTO_SOURCES[${LEVEL}]} ${DST_PATH}/
            cp src/gf/ref/${PRECOMP_GF_DIRS[${LEVEL}]}/fp2.c ${DST_PATH}/${FIAT_CRYPTO_FP2_SOURCES[${LEVEL}]}
            cp src/gf/ref/gfx/{fp.c,fp2.c} ${DST_PATH}/
        elif [ ${ARITHMETIC} = "mikes_arithmetic" ]; then
            ARITH_FLAG="ARITH_MIKE"
    
            cp src/gf/generic/{bn.inc,inversion.inc,symbol.inc} ${DST_PATH}/
            cp src/gf/mike/gfx/*.c ${DST_PATH}/
            cp src/gf/mike/include/*.h ${DST_PATH}/
            cp src/gf/mike/${PRECOMP_GF_DIRS[${LEVEL}]}/*.c ${DST_PATH}/
        elif [ ${ARITHMETIC} = "m4_arithmetic" ]; then
            ARITH_FLAG="ARITH_M4"
    
            cp src/gf/mike/gfx/*.c ${DST_PATH}/
            cp src/gf/mike/include/*.h ${DST_PATH}/
            cp pqm4/${LVL}/${M4_ARITHMETIC_SOURCES[${LEVEL}]} ${DST_PATH}/
        fi

        echo "elf/mupq_crypto_sign_sqisign${LEVEL}_${ARITHMETIC}_%.elf: CPPFLAGS+=-DRADIX_32 -D${ARITH_FLAG}" > ${DST_PATH}/config.mk
        echo "obj/libmupq_crypto_sign_sqisign${LEVEL}_${ARITHMETIC}.a: CPPFLAGS+=-DRADIX_32 -D${ARITH_FLAG}" >> ${DST_PATH}/config.mk

        for FILE in ecx/basis.c ecx/ec.c ecx/isog_chains.c ecx/kps.c ecx/xeval.c ecx/xisog.c \
                    include/curve_extras.h include/ec.h include/isog.h include/sdacs.h
        do
            cp src/ec/ref/$FILE ${DST_PATH}/
        done

        cp src/id2iso/ref/id2isox/*.c ${DST_PATH}/
        cp src/id2iso/ref/include/*.h ${DST_PATH}/

        cp src/protocols/ref/protocolsx/{encode.c,verif.c} ${DST_PATH}/
        cp src/protocols/ref/include/*.h ${DST_PATH}/
    done
done
