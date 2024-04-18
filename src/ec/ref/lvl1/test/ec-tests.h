#ifndef EC_TESTS_H
#define EC_TESTS_H

#include "test_extras.h"
#include <curve_extras.h>
#include <stdio.h>
#include <string.h>
#include <bench.h>       //////// NOTE: enable later
#include "test-basis.h"
#include "ec_params.h"
#include <torsion_constants.h>

// Global constants
extern const digit_t p[NWORDS_FIELD];

// Benchmark and test parameters  
static int BENCH_LOOPS = 1000;       // Number of iterations per bench
static int TEST_LOOPS  = 512;       // Number of iterations per test


bool ec_test()
{ // Tests for ecc arithmetic
    bool OK = true;
    int passed;
    ec_point_t P = {0}, Q = {0}, R = {0}, S = {0}, SS = {0}, PQ = {0};
    ec_point_t AC = {0};
    digit_t k[NWORDS_ORDER] = {0}, l[NWORDS_ORDER] = {0};

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Testing ecc functions: \n\n"); 

    // Point doubling
    passed = 1;
#ifdef RADIX_32
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    P.x.re[0] = 0x861BD329; P.x.re[1] = 0xDFD70ED0; P.x.re[2] = 0x8C7F5540; P.x.re[3] = 0x20ACD375;
    P.x.re[4] = 0x7277F80A; P.x.re[5] = 0x3DCCDC00; P.x.re[6] = 0x2981DCE1; P.x.re[7] = 0x18D6D2A2;
    P.x.im[0] = 0x3F08F38C; P.x.im[1] = 0x3C23730A; P.x.im[2] = 0xFD3D954D; P.x.im[3] = 0x98BB973A;
    P.x.im[4] = 0x2829AE8A; P.x.im[5] = 0x8D98ADFC; P.x.im[6] = 0x6369AFBA; P.x.im[7] = 0x21A2464D;
#elif defined(ARITH_MIKE)
    P.x.re[0] = 0x61BD329; P.x.re[1] = 0x1EB87684; P.x.re[2] = 0x1FD55037; P.x.re[3] = 0x19A6EB18; P.x.re[4] = 0x1F80A20A;
    P.x.re[5] = 0xE00393B; P.x.re[6] = 0x1384F733; P.x.re[7] = 0x1445303B; P.x.re[8] = 0x18D6D2;
    P.x.im[0] = 0x1F08F38C; P.x.im[1] = 0x11B9851; P.x.im[2] = 0xF65534F; P.x.im[3] = 0x172E75FA; P.x.im[4] = 0x1AE8A98B;
    P.x.im[5] = 0x16FE1414; P.x.im[6] = 0x1EEA3662; P.x.im[7] = 0x9AC6D35; P.x.im[8] = 0x21A246;
#endif
#elif defined(RADIX_64)
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    P.x.re[0] = 0xDFD70ED0861BD329; P.x.re[1] = 0x20ACD3758C7F5540; P.x.re[2] = 0x3DCCDC007277F80A; P.x.re[3] = 0x18D6D2A22981DCE1;
    P.x.im[0] = 0x3C23730A3F08F38C; P.x.im[1] = 0x98BB973AFD3D954D; P.x.im[2] = 0x8D98ADFC2829AE8A; P.x.im[3] = 0x21A2464D6369AFBA;
#elif defined(ARITH_MIKE)
    P.x.re[0] = 0x70ED0861BD329; P.x.re[1] = 0x758C7F5540DFD; P.x.re[2] = 0x277F80A20ACD3; P.x.re[3] = 0xDCE13DCCDC007; P.x.re[4] = 0x18D6D2A22981;
    P.x.im[0] = 0x3730A3F08F38C; P.x.im[1] = 0x3AFD3D954D3C2; P.x.im[2] = 0x829AE8A98BB97; P.x.im[3] = 0xAFBA8D98ADFC2; P.x.im[4] = 0x21A2464D6369;
#endif
#endif

    P.z.re[0] = 0x01;

    AC.z.re[0] = 0x01;
    fp2_tomont(&AC.z, &AC.z);
        
    fp2_tomont(&R.x, &P.x);
    fp2_tomont(&R.z, &P.z);
    xDBL(&S, &R, &AC);
    fp2_copy(&SS.x, &S.x);    // Copy of S = SS <- 2P 
    fp2_copy(&SS.z, &S.z);
    fp2_inv(&S.z);
    fp2_mul(&S.x, &S.x, &S.z);
    fp2_frommont(&S.x, &S.x);

#ifdef RADIX_32
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    R.x.re[0] = 0x4AF90FC8; R.x.re[1] = 0x5950EE0A; R.x.re[2] = 0xA0A98B08; R.x.re[3] = 0x16488065;
    R.x.re[4] = 0x29DA0FD1; R.x.re[5] = 0xCE653222; R.x.re[6] = 0x781EE204; R.x.re[7] = 0x270A35FF;
    R.x.im[0] = 0x9EC57F6B; R.x.im[1] = 0x564447FD; R.x.im[2] = 0x4294F729; R.x.im[3] = 0x2EE24E98;
    R.x.im[4] = 0x0E972C71; R.x.im[5] = 0x53A6C736; R.x.im[6] = 0x928A7C7E; R.x.im[7] = 0x4FCF4B9;
#elif defined(ARITH_MIKE)
    R.x.re[0] = 0xAF90FC8; R.x.re[1] = 0xA877052; R.x.re[2] = 0xA62C216; R.x.re[3] = 0x1100CB41; R.x.re[4] = 0xFD1164;
    R.x.re[5] = 0x191114ED; R.x.re[6] = 0x8133994; R.x.re[7] = 0x1FEF03DC; R.x.re[8] = 0x270A35;
    R.x.im[0] = 0x1EC57F6B; R.x.im[1] = 0x12223FEC; R.x.im[2] = 0x53DCA55; R.x.im[3] = 0x49D3085; R.x.im[4] = 0x12C712EE;
    R.x.im[5] = 0x39B074B; R.x.im[6] = 0x11F94E9B; R.x.im[7] = 0x1732514F; R.x.im[8] = 0x4FCF4;
#endif
#elif defined(RADIX_64)
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    R.x.re[0] = 0x5950EE0A4AF90FC8; R.x.re[1] = 0x16488065A0A98B08; R.x.re[2] = 0xCE65322229DA0FD1; R.x.re[3] = 0x270A35FF781EE204;
    R.x.im[0] = 0x564447FD9EC57F6B; R.x.im[1] = 0x2EE24E984294F729; R.x.im[2] = 0x53A6C7360E972C71; R.x.im[3] = 0x4FCF4B9928A7C7E;
#elif defined(ARITH_MIKE)
    R.x.re[0] = 0x0EE0A4AF90FC8; R.x.re[1] = 0x65A0A98B08595; R.x.re[2] = 0x9DA0FD1164880; R.x.re[3] = 0xE204CE6532222; R.x.re[4] = 0x270A35FF781E;
    R.x.im[0] = 0x447FD9EC57F6B; R.x.im[1] = 0x984294F729564; R.x.im[2] = 0xE972C712EE24E; R.x.im[3] = 0x7C7E53A6C7360; R.x.im[4] =  0x4FCF4B9928A;
#endif
#endif

    if (compare_words((digit_t*)&R.x, (digit_t*)&S.x, NWORDS_FIELD*2)!=0) { passed=0; goto out0; }
    
#ifdef RADIX_32
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    Q.x.re[0] = 0x70C70053; Q.x.re[1] = 0xC46076A6; Q.x.re[2] = 0x3AB9ED13; Q.x.re[3] = 0x97517AFA;
    Q.x.re[4] = 0x42EDF993; Q.x.re[5] = 0x349644C9; Q.x.re[6] = 0x6F29AF9E; Q.x.re[7] = 0xBB4A4DB;
    Q.x.im[0] = 0xB5A15BB0; Q.x.im[1] = 0x8B47629F; Q.x.im[2] = 0x953C1A10; Q.x.im[3] = 0x4EC6E809;
    Q.x.im[4] = 0x6CBB84D6; Q.x.im[5] = 0x1F83F0EC; Q.x.im[6] = 0xD33265D3; Q.x.im[7] = 0x1D8417C1;
#elif defined(ARITH_MIKE)
    Q.x.re[0] = 0x10C70053; Q.x.re[1] = 0x303B533; Q.x.re[2] = 0xE7B44F1; Q.x.re[3] = 0x2F5F475; Q.x.re[4] = 0x1F993975;
    Q.x.re[5] = 0x264A176; Q.x.re[6] = 0x1E78D259; Q.x.re[7] = 0x1B6DE535; Q.x.re[8] = 0xBB4A4;
    Q.x.im[0] = 0x15A15BB0; Q.x.im[1] = 0x1A3B14FD; Q.x.im[2] = 0xF068422; Q.x.im[3] = 0xDD0132A; Q.x.im[4] = 0x184D64EC;
    Q.x.im[5] = 0x1876365D; Q.x.im[6] = 0x174C7E0F; Q.x.im[7] = 0x183A664C; Q.x.im[8] = 0x1D8417;
#endif
#elif defined(RADIX_64)
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    Q.x.re[0] = 0xC46076A670C70053; Q.x.re[1] = 0x97517AFA3AB9ED13; Q.x.re[2] = 0x349644C942EDF993; Q.x.re[3] = 0xBB4A4DB6F29AF9E;
    Q.x.im[0] = 0x8B47629FB5A15BB0; Q.x.im[1] = 0x4EC6E809953C1A10; Q.x.im[2] = 0x1F83F0EC6CBB84D6; Q.x.im[3] = 0x1D8417C1D33265D3;
#elif defined(ARITH_MIKE)
    Q.x.re[0] = 0x076A670C70053; Q.x.re[1] = 0xFA3AB9ED13C46; Q.x.re[2] = 0x2EDF99397517A; Q.x.re[3] = 0xAF9E349644C94; Q.x.re[4] = 0xBB4A4DB6F29;
    Q.x.im[0] = 0x7629FB5A15BB0; Q.x.im[1] = 0x09953C1A108B4; Q.x.im[2] = 0xCBB84D64EC6E8; Q.x.im[3] = 0x65D31F83F0EC6; Q.x.im[4] = 0x1D8417C1D332;
#endif
#endif
    Q.z.re[0] = 0x01;

#ifdef RADIX_32
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    PQ.x.re[0] = 0x1BE5534F; PQ.x.re[1] = 0x853F66D1; PQ.x.re[2] = 0x52D03D4A; PQ.x.re[3] = 0x27C8FD4E;
    PQ.x.re[4] = 0x0A0C29D2; PQ.x.re[5] = 0xF88EA78D; PQ.x.re[6] = 0xD397A067; PQ.x.re[7] = 0x2F6DFB07;
    PQ.x.im[0] = 0x34434BA1; PQ.x.im[1] = 0xE8DBC4AA; PQ.x.im[2] = 0x2636F8A0; PQ.x.im[3] = 0x7A73AE18;
    PQ.x.im[4] = 0x137868EB; PQ.x.im[5] = 0x419EC260; PQ.x.im[6] = 0x1703D43F; PQ.x.im[7] = 0x129B3E30;
#elif defined(ARITH_MIKE)
    PQ.x.re[0] = 0x1BE5534F; PQ.x.re[1] = 0x9FB3688; PQ.x.re[2] = 0x140F52A1; PQ.x.re[3] = 0x11FA9CA5; PQ.x.re[4] = 0x29D227C;
    PQ.x.re[5] = 0x13C68506; PQ.x.re[6] = 0x19FE23A; PQ.x.re[7] = 0xFA72F4; PQ.x.re[8] = 0x2F6DFB;
    PQ.x.im[0] = 0x14434BA1; PQ.x.im[1] = 0x6DE2551; PQ.x.im[2] = 0xDBE283A; PQ.x.im[3] = 0x75C304C; PQ.x.im[4] = 0x68EB7A7;
    PQ.x.im[5] = 0x13009BC; PQ.x.im[6] = 0x10FD067B; PQ.x.im[7] = 0x602E07A; PQ.x.im[8] = 0x129B3E;
#endif
#elif defined(RADIX_64)
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    PQ.x.re[0] = 0x853F66D11BE5534F; PQ.x.re[1] = 0x27C8FD4E52D03D4A; PQ.x.re[2] = 0xF88EA78D0A0C29D2; PQ.x.re[3] = 0x2F6DFB07D397A067;
    PQ.x.im[0] = 0xE8DBC4AA34434BA1; PQ.x.im[1] = 0x7A73AE182636F8A0; PQ.x.im[2] = 0x419EC260137868EB; PQ.x.im[3] = 0x129B3E301703D43F;
#elif defined(ARITH_MIKE)
    PQ.x.re[0] = 0xF66D11BE5534F; PQ.x.re[1] = 0x4E52D03D4A853; PQ.x.re[2] = 0xA0C29D227C8FD; PQ.x.re[3] = 0xA067F88EA78D0; PQ.x.re[4] = 0x2F6DFB07D397;
    PQ.x.im[0] = 0xBC4AA34434BA1; PQ.x.im[1] = 0x182636F8A0E8D; PQ.x.im[2] = 0x37868EB7A73AE; PQ.x.im[3] = 0xD43F419EC2601; PQ.x.im[4] = 0x129B3E301703;
#endif
#endif
    PQ.z.re[0] = 0x01;

    fp2_tomont(&S.x, &Q.x);
    fp2_tomont(&S.z, &Q.z);
    fp2_tomont(&PQ.x, &PQ.x);
    fp2_tomont(&PQ.z, &PQ.z);
    xADD(&S, &SS, &S, &PQ);
    fp2_inv(&S.z);
    fp2_mul(&S.x, &S.x, &S.z);
    fp2_frommont(&S.x, &S.x);

#ifdef RADIX_32
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    R.x.re[0] = 0x93AB4FF9; R.x.re[1] = 0xED0BEB8F; R.x.re[2] = 0x80CD49BF; R.x.re[3] = 0x27CF508B;
    R.x.re[4] = 0xFA04B2BA; R.x.re[5] = 0x38A6134D; R.x.re[6] = 0xE109EF1F; R.x.re[7] = 0x27B4CB15;
    R.x.im[0] = 0xFD227BDE; R.x.im[1] = 0x6F731BA6; R.x.im[2] = 0x341167F8; R.x.im[3] = 0x14C12335;
    R.x.im[4] = 0x7866E27A; R.x.im[5] = 0xECA7B60F; R.x.im[6] = 0x52880457; R.x.im[7] = 0x2A7A79A1;
#elif defined(ARITH_MIKE)
    R.x.re[0] = 0x13AB4FF9; R.x.re[1] = 0x85F5C7C; R.x.re[2] = 0x13526FFB; R.x.re[3] = 0x1EA11701; R.x.re[4] = 0xB2BA27C;
    R.x.re[5] = 0x9A6FD02; R.x.re[6] = 0x1C7CE298; R.x.re[7] = 0x2BC213D; R.x.re[8] = 0x27B4CB;
    R.x.im[0] = 0x1D227BDE; R.x.im[1] = 0x1B98DD37; R.x.im[2] = 0x459FE1B; R.x.im[3] = 0x2466A68; R.x.im[4] = 0xE27A14C;
    R.x.im[5] = 0x1B07BC33; R.x.im[6] = 0x115FB29E; R.x.im[7] = 0x142A5100; R.x.im[8] = 0x2A7A79;
#endif
#elif defined(RADIX_64)
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    R.x.re[0] = 0xED0BEB8F93AB4FF9; R.x.re[1] = 0x27CF508B80CD49BF; R.x.re[2] = 0x38A6134DFA04B2BA; R.x.re[3] = 0x27B4CB15E109EF1F;
    R.x.im[0] = 0x6F731BA6FD227BDE; R.x.im[1] = 0x14C12335341167F8; R.x.im[2] = 0xECA7B60F7866E27A; R.x.im[3] = 0x2A7A79A152880457;
#elif defined(ARITH_MIKE)
    R.x.re[0] = 0xBEB8F93AB4FF9; R.x.re[1] = 0x8B80CD49BFED0; R.x.re[2] = 0xA04B2BA27CF50; R.x.re[3] = 0xEF1F38A6134DF; R.x.re[4] = 0x27B4CB15E109;
    R.x.im[0] = 0x31BA6FD227BDE; R.x.im[1] = 0x35341167F86F7; R.x.im[2] = 0x866E27A14C123; R.x.im[3] = 0x0457ECA7B60F7; R.x.im[4] = 0x2A7A79A15288;
#endif
#endif

    if (compare_words((digit_t*)&R.x, (digit_t*)&S.x, NWORDS_FIELD*2) != 0) { passed = 0; goto out0; }
    
    fp2_tomont(&R.x, &P.x);
    fp2_tomont(&R.z, &P.z);
    k[0] = 126;
    xMUL(&S, &R, k, (ec_curve_t*)&AC);
    fp2_inv(&S.z);
    fp2_mul(&S.x, &S.x, &S.z);
    fp2_frommont(&S.x, &S.x);

#ifdef RADIX_32
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    R.x.re[0] = 0x1203A147; R.x.re[1] = 0xDE80F87A; R.x.re[2] = 0x928A3B2D; R.x.re[3] = 0xD59E1215;
    R.x.re[4] = 0xA5A8CE46; R.x.re[5] = 0xD5A67F83; R.x.re[6] = 0x488C9CDF; R.x.re[7] = 0xA11E162;
    R.x.im[0] = 0x9A26741B; R.x.im[1] = 0x9417D0D7; R.x.im[2] = 0xF0FE5EEC; R.x.im[3] = 0x8B1F47D6;
    R.x.im[4] = 0xB054CE36; R.x.im[5] = 0xE52188DC; R.x.im[6] = 0xC3148AB3; R.x.im[7] = 0x1A8075A6;
#elif defined(ARITH_MIKE)
    R.x.re[0] = 0x1203A147; R.x.re[1] = 0x1407C3D0; R.x.re[2] = 0x28ECB77; R.x.re[3] = 0x1C242B25; R.x.re[4] = 0xCE46D59;
    R.x.re[5] = 0x1FC1D2D4; R.x.re[6] = 0x137F5699; R.x.re[7] = 0xC491193; R.x.re[8] = 0xA11E1;
    R.x.im[0] = 0x1A26741B; R.x.im[1] = 0xBE86BC; R.x.im[2] = 0x1F97BB25; R.x.im[3] = 0x1E8FADE1; R.x.im[4] = 0xCE368B1;
    R.x.im[5] = 0x46E582A; R.x.im[6] = 0xACF9486; R.x.im[7] = 0x14D86291; R.x.im[8] = 0x1A8075;
#endif
#elif defined(RADIX_64)
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    R.x.re[0] = 0xDE80F87A1203A147; R.x.re[1] = 0xD59E1215928A3B2D; R.x.re[2] = 0xD5A67F83A5A8CE46; R.x.re[3] = 0xA11E162488C9CDF;
    R.x.im[0] = 0x9417D0D79A26741B; R.x.im[1] = 0x8B1F47D6F0FE5EEC; R.x.im[2] = 0xE52188DCB054CE36; R.x.im[3] = 0x1A8075A6C3148AB3;
#elif defined(ARITH_MIKE)
    R.x.re[0] = 0x0F87A1203A147; R.x.re[1] = 0x15928A3B2DDE8; R.x.re[2] = 0x5A8CE46D59E12; R.x.re[3] = 0x9CDFD5A67F83A; R.x.re[4] =  0xA11E162488C;
    R.x.im[0] = 0x7D0D79A26741B; R.x.im[1] = 0xD6F0FE5EEC941; R.x.im[2] = 0x054CE368B1F47; R.x.im[3] = 0x8AB3E52188DCB; R.x.im[4] = 0x1A8075A6C314;
#endif
#endif

    if (compare_words((digit_t*)&R.x, (digit_t*)&S.x, NWORDS_FIELD*2) != 0) { passed = 0; goto out0; }
    
    fp2_tomont(&R.x, &P.x);
    fp2_tomont(&R.z, &P.z);
#ifdef RADIX_32
    k[0] = 0xC6B2D8CD;
    k[1] = 0xE77AD6B6;
    k[2] = 0x00F38D12;
    k[3] = 0xDE43A0B6;
    k[4] = 0x97E17CE2;
    k[5] = 0xA35F4A78;
    k[6] = 0x614D1237;
    k[7] = 0x10ACB62E;
#elif defined(RADIX_64)
    k[0] = 0xE77AD6B6C6B2D8CD;
    k[1] = 0xDE43A0B600F38D12;
    k[2] = 0xA35F4A7897E17CE2;
    k[3] = 0x10ACB62E614D1237;
#endif
    xMUL(&S, &R, k, (ec_curve_t*)&AC);
    fp2_inv(&S.z);
    fp2_mul(&S.x, &S.x, &S.z);
    fp2_frommont(&S.x, &S.x);

#ifdef RADIX_32
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    R.x.re[0] = 0x68A3E7C0; R.x.re[1] = 0xD3938B0A; R.x.re[2] = 0x208A0595; R.x.re[3] = 0xE0667113;
    R.x.re[4] = 0x84E9CB60; R.x.re[5] = 0x258F314C; R.x.re[6] = 0xCA59AB71; R.x.re[7] = 0x14984BA7;
    R.x.im[0] = 0xEE3BFEF4; R.x.im[1] = 0xFE728423; R.x.im[2] = 0xE21AE0E4; R.x.im[3] = 0xBF68C42F;
    R.x.im[4] = 0x528609CA; R.x.im[5] = 0xA8FAF9C9; R.x.im[6] = 0xA1DC0285; R.x.im[7] = 0x1225EC77;
#elif defined(ARITH_MIKE)
    R.x.re[0] = 0x8A3E7C0; R.x.re[1] = 0x1C9C5853; R.x.re[2] = 0x2816574; R.x.re[3] = 0xCE22641; R.x.re[4] = 0x1CB60E06;
    R.x.re[5] = 0x18A64274; R.x.re[6] = 0xDC4963C; R.x.re[7] = 0x14F94B35; R.x.re[8] = 0x14984B;
    R.x.im[0] = 0xE3BFEF4; R.x.im[1] = 0x1394211F; R.x.im[2] = 0x6B8393F; R.x.im[3] = 0x11885FC4; R.x.im[4] = 0x9CABF6;
    R.x.im[5] = 0x1CE4A943; R.x.im[6] = 0xA16A3EB; R.x.im[7] = 0xEF43B80; R.x.im[8] = 0x1225EC;
#endif
#elif defined(RADIX_64)
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    R.x.re[0] = 0xD3938B0A68A3E7C0; R.x.re[1] = 0xE0667113208A0595; R.x.re[2] = 0x258F314C84E9CB60; R.x.re[3] = 0x14984BA7CA59AB71;
    R.x.im[0] = 0xFE728423EE3BFEF4; R.x.im[1] = 0xBF68C42FE21AE0E4; R.x.im[2] = 0xA8FAF9C9528609CA; R.x.im[3] = 0x1225EC77A1DC0285;
#elif defined(ARITH_MIKE)
    R.x.re[0] = 0x38B0A68A3E7C0; R.x.re[1] = 0x13208A0595D39; R.x.re[2] = 0x4E9CB60E06671; R.x.re[3] = 0xAB71258F314C8; R.x.re[4] = 0x14984BA7CA59;
    R.x.im[0] = 0x28423EE3BFEF4; R.x.im[1] = 0x2FE21AE0E4FE7; R.x.im[2] = 0x28609CABF68C4; R.x.im[3] = 0x0285A8FAF9C95; R.x.im[4] = 0x1225EC77A1DC;
#endif
#endif

    if (compare_words((digit_t*)&R.x, (digit_t*)&S.x, NWORDS_FIELD*2) != 0) { passed = 0; goto out0; }
    
    fp2_tomont(&R.x, &Q.x);
    fp2_tomont(&R.z, &Q.z);
#ifdef RADIX_32
    k[0] = 0xC6B2D8CD;
    k[1] = 0xE77AD6B6;
    k[2] = 0x00F38D12;
    k[3] = 0xDE43A0B6;
    k[4] = 0x97E17CE2;
    k[5] = 0xA35F4A78;
    k[6] = 0x614D1237;
    k[7] = 0x10ACB62E;
    l[0] = 0xC6B2D8C0;
    l[1] = 0x34AB78B6;
    l[2] = 0xD00F38D1;
    l[3] = 0xDE6B2D8C;
    l[4] = 0x97E17CE2;
    l[5] = 0xA35F4A78;
    l[6] = 0x89614D13;
    l[7] = 0x20ACF4A7;
#elif defined(RADIX_64)
    k[0] = 0xE77AD6B6C6B2D8CD;
    k[1] = 0xDE43A0B600F38D12;
    k[2] = 0xA35F4A7897E17CE2;
    k[3] = 0x10ACB62E614D1237;
    l[0] = 0x34AB78B6C6B2D8C0;
    l[1] = 0xDE6B2D8CD00F38D1;
    l[2] = 0xA35F4A7897E17CE2;
    l[3] = 0x20ACF4A789614D13;
#endif
    fp2_inv(&SS.z);
    fp2_mul(&SS.x, &SS.x, &SS.z);
    fp2_copy(&SS.z, &R.z);
    xDBLMUL(&S, &R, k, &SS, l, &PQ, (ec_curve_t*)&AC);
    fp2_inv(&S.z);
    fp2_mul(&S.x, &S.x, &S.z);
    fp2_frommont(&S.x, &S.x);

#ifdef RADIX_32
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    R.x.re[0] = 0x609B992F; R.x.re[1] = 0x554E1ADC; R.x.re[2] = 0xF8CC4C42; R.x.re[3] = 0xE407D961;
    R.x.re[4] = 0xED5A68CE; R.x.re[5] = 0x1CF626AF; R.x.re[6] = 0xEE110483; R.x.re[7] = 0x6D02692;
    R.x.im[0] = 0x831C8997; R.x.im[1] = 0x16FB094E; R.x.im[2] = 0x1DC5F702; R.x.im[3] = 0xFDE4ECF3;
    R.x.im[4] = 0x8DFAD7B4; R.x.im[5] = 0x89303D86; R.x.im[6] = 0x1346F22D; R.x.im[7] = 0xC91ACE8;
#elif defined(ARITH_MIKE)
    R.x.re[0] = 0x9B992F; R.x.re[1] = 0xA70D6E3; R.x.re[2] = 0x13131095; R.x.re[3] = 0xFB2C3F1; R.x.re[4] = 0x68CEE40;
    R.x.re[5] = 0x1357F6AD; R.x.re[6] = 0x120C73D8; R.x.re[7] = 0x125DC220; R.x.re[8] = 0x6D026;
    R.x.im[0] = 0x31C8997; R.x.im[1] = 0x17D84A74; R.x.im[2] = 0x117DC085; R.x.im[3] = 0x9D9E63B; R.x.im[4] = 0xD7B4FDE;
    R.x.im[5] = 0x1EC346FD; R.x.im[6] = 0x8B624C0; R.x.im[7] = 0x1D0268DE; R.x.im[8] = 0xC91AC;
#endif
#elif defined(RADIX_64)
#if defined(ARITH_REF) || defined(ARITH_BROADWELL)
    R.x.re[0] = 0x554E1ADC609B992F; R.x.re[1] = 0xE407D961F8CC4C42; R.x.re[2] = 0x1CF626AFED5A68CE; R.x.re[3] = 0x6D02692EE110483;
    R.x.im[0] = 0x16FB094E831C8997; R.x.im[1] = 0xFDE4ECF31DC5F702; R.x.im[2] = 0x89303D868DFAD7B4; R.x.im[3] = 0xC91ACE81346F22D;
#elif defined(ARITH_MIKE)
    R.x.re[0] = 0xE1ADC609B992F; R.x.re[1] = 0x61F8CC4C42554; R.x.re[2] = 0xD5A68CEE407D9; R.x.re[3] = 0x04831CF626AFE; R.x.re[4] = 0x6D02692EE11;
    R.x.im[0] = 0xB094E831C8997; R.x.im[1] = 0xF31DC5F70216F; R.x.im[2] = 0xDFAD7B4FDE4EC; R.x.im[3] = 0xF22D89303D868; R.x.im[4] = 0xC91ACE81346;
#endif
#endif

    if (compare_words((digit_t*)&R.x, (digit_t*)&S.x, NWORDS_FIELD*2) != 0) { passed = 0; goto out0; }
    
out0:
    if (passed==1) printf("  ECC arithmetic tests ............................................ PASSED");
    else { printf("  ECC arithmetic tests... FAILED"); printf("\n"); return false; }
    printf("\n");
 
    return OK;
}

bool dlog_test()
{ // Tests for dlog
    bool OK = true;
    int passed;
    ec_point_t P = {0}, Q = {0}, R = {0}, S = {0}, SS = {0}, PQ = {0};
    ec_curve_t AC = {0};
    ec_basis_t PQ2;
    digit_t scalarP[NWORDS_ORDER], scalarQ[NWORDS_ORDER], k[NWORDS_ORDER] = {0}, l[NWORDS_ORDER] = {0};
    digit_t kt[NWORDS_ORDER], lt[NWORDS_ORDER], f1[NWORDS_ORDER] = {0}, f2[NWORDS_ORDER] = {0}, zero[NWORDS_ORDER] = {0}, tpFdiv2[NWORDS_ORDER] = {0}, tpF[NWORDS_ORDER] = {0};

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Testing dlog functions: \n\n");

    // dlog2 testing
    passed = 1;
    
    fp2_tomont(&P.x, &xP2);
    fp_mont_setone(P.z.re);
    fp_set(P.z.im, 0);
    
    fp2_tomont(&Q.x, &xQ2);
    fp_mont_setone(Q.z.re);
    fp_set(Q.z.im, 0);
    
    fp2_tomont(&PQ.x, &xPQ2);
    fp_mont_setone(PQ.z.re);
    fp_set(PQ.z.im, 0);

    AC.C.re[0] = 0x01;
    memcpy(f1, TWOpFm1, NWORDS_ORDER*RADIX/8);
    memcpy(f2, TORSION_PLUS_2POWER_DIGITS, NWORDS_ORDER*RADIX/8);
    fp2_tomont(&AC.C, &AC.C);

    copy_point(&PQ2.P, &P);
    copy_point(&PQ2.Q, &Q);
    copy_point(&PQ2.PmQ, &PQ);
#ifdef RADIX_32
    k[0] = 0xFFFFFFFF;
    k[1] = 0xFFFFFFFF;
    k[2] = 0x000007FF;
    k[3] = 0x00000000;
    l[0] = 0xFFFFFFFE;
    l[1] = 0xFFFFFFFF;
    l[2] = 0x000007FF;
    l[3] = 0x00000000;
#elif defined(RADIX_64)
    k[0] = 0xFFFFFFFFFFFFFFFF;
    k[1] = 0x00000000000007FF;
    l[0] = 0xFFFFFFFFFFFFFFFE;
    l[1] = 0x00000000000007FF;
#endif

    for (int n = 0; n < TEST_LOOPS; n++)
    {
        k[0] -= 1;
        l[0] -= 2;
        xDBLMUL(&R, &P, k, &Q, l, &PQ, &AC);
        ec_dlog_2(scalarP, scalarQ, &PQ2, &R, &AC);

        memcpy(kt, k, NWORDS_ORDER*RADIX/8);
        memcpy(lt, l, NWORDS_ORDER*RADIX/8);
        if (compare_words(k, f1, NWORDS_ORDER) == 1 ||
           (compare_words(l, f1, NWORDS_ORDER) == 1 && (compare_words(k, zero, NWORDS_ORDER) == 0 || compare_words(k, f1, NWORDS_ORDER) == 0))) {
            if (compare_words(k, zero, NWORDS_ORDER) != 0) {
                sub_test(kt, f2, kt, NWORDS_ORDER);
            }
            if (compare_words(l, zero, NWORDS_ORDER) != 0) {
                sub_test(lt, f2, lt, NWORDS_ORDER);
            }
        }
        if (compare_words((digit_t*)scalarP, (digit_t*)kt, NWORDS_ORDER) != 0 || compare_words((digit_t*)scalarQ, (digit_t*)lt, NWORDS_ORDER) != 0) { passed = 0; break; }
    }

    if (passed == 1) printf("  dlog2 tests ..................................................... PASSED");
    else { printf("  dlog2 tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    // dlog3 testing
    passed = 1;
    
    fp2_tomont(&P.x, &xP3);
    fp_mont_setone(P.z.re);
    fp_set(P.z.im, 0);
    
    fp2_tomont(&Q.x, &xQ3);
    fp_mont_setone(Q.z.re);
    fp_set(Q.z.im, 0);
    
    fp2_tomont(&PQ.x, &xPQ3);
    fp_mont_setone(PQ.z.re);
    fp_set(PQ.z.im, 0);

    AC.C.re[0] = 0x01;
    memcpy(tpFdiv2, THREEpFdiv2, NWORDS_ORDER*RADIX/8);
    memcpy(tpF, TORSION_PLUS_3POWER_DIGITS, NWORDS_ORDER*RADIX/8);
    fp2_tomont(&AC.C, &AC.C);

    copy_point(&PQ2.P, &P);
    copy_point(&PQ2.Q, &Q);
    copy_point(&PQ2.PmQ, &PQ);
#ifdef RADIX_32
    k[3] = 0;
    l[3] = 0;
    k[2] = 0;
    l[2] = 0;
    k[1] = 0x02153E46;
    l[1] = 0x02153E46;
    k[0] = 0x8B91C6D1;
    l[0] = 0x8B91C6D0;
#elif defined(RADIX_64)
    k[1] = 0;
    l[1] = 0;
    k[0] = 0x02153E468B91C6D1;
    l[0] = 0x02153E468B91C6D0;
#endif

    for (int n = 0; n < TEST_LOOPS; n++)
    {
        k[0] -= 1;
        l[0] -= 2;
        xDBLMUL(&R, &P, k, &Q, l, &PQ, &AC);
        ec_dlog_3(scalarP, scalarQ, &PQ2, &R, &AC);

        memcpy(kt, k, NWORDS_ORDER*RADIX/8);
        memcpy(lt, l, NWORDS_ORDER*RADIX/8);
        if (compare_words(k, tpFdiv2, NWORDS_ORDER) == 1 ||
           (compare_words(l, tpFdiv2, NWORDS_ORDER) == 1 && compare_words(k, zero, NWORDS_ORDER) == 0)) {
            if (compare_words(k, zero, NWORDS_ORDER) != 0) {
                sub_test(kt, tpF, kt, NWORDS_ORDER);
            }
            if (compare_words(l, zero, NWORDS_ORDER) != 0) {
                sub_test(lt, tpF, lt, NWORDS_ORDER);
            }
        }
        if (compare_words((digit_t*)scalarP, (digit_t*)kt, NWORDS_ORDER) != 0 || compare_words((digit_t*)scalarQ, (digit_t*)lt, NWORDS_ORDER) != 0) { passed = 0; break; }
    }

    if (passed == 1) printf("  dlog3 tests ..................................................... PASSED");
    else { printf("  dlog3 tests... FAILED"); printf("\n"); return false; }
    printf("\n");

    return OK;
}

bool ec_run()
{
    bool OK = true;
    int n;
    unsigned long long cycles, cycles1, cycles2;
    ec_point_t P, Q, R, PQ, AC;
    digit_t k[NWORDS_ORDER], l[NWORDS_ORDER];
        
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("Benchmarking ecc arithmetic: \n\n"); 

    // Point doubling
    cycles = 0;
    for (n=0; n<BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles(); 
        xDBL(&Q, &P, &AC);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);
    }
    printf("  Montgomery x-only doubling runs in .............................. %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // Point addition
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        xADD(&R, &Q, &P, &PQ);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Montgomery x-only addition runs in .............................. %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // Point multiplication
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        xMUL(&Q, &P, k, (ec_curve_t*)&AC);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Montgomery x-only scalar multiplication runs in ................. %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // Point multiplication
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        cycles1 = cpucycles();
        xDBLMUL(&R, &P, k, &Q, l, &PQ, (ec_curve_t*)&AC);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Montgomery x-only double-scalar multiplication runs in .......... %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    return OK;
}

bool dlog_run()
{
    bool OK = true;
    int n;
    unsigned long long cycles, cycles1, cycles2;
    ec_point_t P = {0}, Q = {0}, R = {0}, S = {0}, SS = {0}, PQ = {0};
    ec_curve_t AC = {0};
    ec_basis_t PQ2;
    digit_t scalarP[NWORDS_ORDER], scalarQ[NWORDS_ORDER], k[NWORDS_ORDER] = {0}, l[NWORDS_ORDER] = {0};

    printf("\n--------------------------------------------------------------------------------------------------------\n\n");
    printf("Benchmarking dlog2: \n\n");

    // dlog2 computation
    
    fp2_tomont(&P.x, &xP2);
    fp_mont_setone(P.z.re);
    fp_set(P.z.im, 0);
    
    fp2_tomont(&Q.x, &xQ2);
    fp_mont_setone(Q.z.re);
    fp_set(Q.z.im, 0);
    
    fp2_tomont(&PQ.x, &xPQ2);
    fp_mont_setone(PQ.z.re);
    fp_set(PQ.z.im, 0);

    AC.C.re[0] = 0x01;
    fp2_tomont(&AC.C, &AC.C);

    copy_point(&PQ2.P, &P);
    copy_point(&PQ2.Q, &Q);
    copy_point(&PQ2.PmQ, &PQ);

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        fprandom_test(k); fprandom_test(l);
        xDBLMUL(&R, &P, k, &Q, l, &PQ, &AC);
        cycles1 = cpucycles();
        ec_dlog_2(scalarP, scalarQ, &PQ2, &R, &AC);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  dlog2 runs in ................................................... %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // dlog3 computation

    fp2_tomont(&P.x, &xP3);
    fp_mont_setone(P.z.re);
    fp_set(P.z.im, 0);
    
    fp2_tomont(&Q.x, &xQ3);
    fp_mont_setone(Q.z.re);
    fp_set(Q.z.im, 0);
    
    fp2_tomont(&PQ.x, &xPQ3);
    fp_mont_setone(PQ.z.re);
    fp_set(PQ.z.im, 0);

    copy_point(&PQ2.P, &P);
    copy_point(&PQ2.Q, &Q);
    copy_point(&PQ2.PmQ, &PQ);

    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        fprandom_test(k); fprandom_test(l);
        xDBLMUL(&R, &P, k, &Q, l, &PQ, &AC);
        cycles1 = cpucycles();
        ec_dlog_3(scalarP, scalarQ, &PQ2, &R, &AC);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  dlog3 runs in ................................................... %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    return OK;
}

#endif
