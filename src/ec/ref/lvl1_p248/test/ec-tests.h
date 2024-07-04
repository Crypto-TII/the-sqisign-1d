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
    k[2] = 0xFFFFFFFF;
    k[3] = 0xFFFFFFFF;
    k[4] = 0xFFFFFFFF;
    k[5] = 0xFFFFFFFF;
    k[6] = 0xFFFFFFFF;
    k[7] = 0x00FFFFFF;
    l[0] = 0xFFFFFFFE;
    l[1] = 0xFFFFFFFF;
    l[2] = 0xFFFFFFFF;
    l[3] = 0xFFFFFFFF;
    l[4] = 0xFFFFFFFF;
    l[5] = 0xFFFFFFFF;
    l[6] = 0xFFFFFFFF;
    l[7] = 0x00FFFFFF;
#elif defined(RADIX_64)
    k[0] = 0xFFFFFFFFFFFFFFFF;
    k[1] = 0xFFFFFFFFFFFFFFFF;
    k[2] = 0xFFFFFFFFFFFFFFFF;
    k[3] = 0x00FFFFFFFFFFFFFF;
    l[0] = 0xFFFFFFFFFFFFFFFE;
    l[1] = 0xFFFFFFFFFFFFFFFF;
    l[2] = 0xFFFFFFFFFFFFFFFF;
    l[3] = 0x00FFFFFFFFFFFFFF;
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

    // Odd cofactor ladder
    ec_point_t A3, A24;
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        fp2random_test(&A3.x);
        fp2random_test(&A3.z);
        fp2random_test(&A24.x);
        fp2random_test(&A24.z);
        fp2random_test(&P.x);
        fp2random_test(&P.z);
        cycles1 = cpucycles();
        xMULv2(&P, &P, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, &A24);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Clearing odd cofactor via Montgomery ladder .................... %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // Odd cofactor DACs
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        fp2random_test(&A3.x);
        fp2random_test(&A3.z);
        fp2random_test(&A24.x);
        fp2random_test(&A24.z);
        fp2random_test(&P.x);
        fp2random_test(&P.z);
        cycles1 = cpucycles();
        for(int i = 0; i < POWER_OF_3; i++){
            xTPL(&P, &P, &A3);
        }
        for(int i = 1; i < P_LEN; i++){
            for(int j = 0; j < TORSION_PLUS_ODD_POWERS[i]; j++){
                xMULdac(&P, &P, DACS[i], DAC_LEN[i], &A24);
            }
        }
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Clearing odd cofactor via DACs ................................ %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // Non-3 cofactor ladder
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        fp2random_test(&A3.x);
        fp2random_test(&A3.z);
        fp2random_test(&A24.x);
        fp2random_test(&A24.z);
        fp2random_test(&P.x);
        fp2random_test(&P.z);
        cycles1 = cpucycles();
        xMULv2(&P, &P, p_cofactor_for_3g, P_COFACTOR_FOR_3G_BITLENGTH, &A24);
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Clearing non-3 cofactor via Montgomery ladder .................. %7lld cycles", cycles/BENCH_LOOPS);
    printf("\n");

    // Non-3 cofactor DACs
    cycles = 0;
    for (n = 0; n < BENCH_LOOPS; n++)
    {
        fp2random_test(&A3.x);
        fp2random_test(&A3.z);
        fp2random_test(&A24.x);
        fp2random_test(&A24.z);
        fp2random_test(&P.x);
        fp2random_test(&P.z);
        cycles1 = cpucycles();
        for(int i = 0; i < POWER_OF_2; i++){
            xDBL(&P, &P, &A24);
        }
        for(int i = 1; i < P_LEN; i++){
            for(int j = 0; j < TORSION_PLUS_ODD_POWERS[i]; j++){
                xMULdac(&P, &P, DACS[i], DAC_LEN[i], &A24);
            }
        }
        cycles2 = cpucycles();
        cycles = cycles + (cycles2 - cycles1);
    }
    printf("  Clearing non-3 cofactor via DACs .............................. %7lld cycles", cycles/BENCH_LOOPS);
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

    return OK;
}

#endif
