#include "isog.h"
#include <assert.h>

static inline void AC_to_A24(ec_point_t *A24, ec_curve_t const *E)
{
    // A24 = (A+2C : 4C)
    fp2_add(&A24->z, &E->C, &E->C);
    fp2_add(&A24->x, &E->A, &A24->z);
    fp2_add(&A24->z, &A24->z, &A24->z);
}

static inline void A24_to_AC(ec_curve_t *E, ec_point_t const *A24)
{
    // (A:C) = ((A+2C)*2-4C : 4C)
    fp2_add(&E->A, &A24->x, &A24->x);
    fp2_sub(&E->A, &E->A, &A24->z);
    fp2_add(&E->A, &E->A, &E->A);
    fp2_copy(&E->C, &A24->z);
}

void ec_eval_even(ec_curve_t* image, const ec_isog_even_t* phi,
    ec_point_t* points, unsigned short length){
        ec_point_t Q4, Q, A24, K[3];
        copy_point(&Q4, &phi->kernel);
        AC_to_A24(&A24, &phi->curve);
        for(int i = 0; i < phi->length - 2; i++)
            xDBLv2(&Q4, &Q4, &A24);
        xDBLv2(&Q, &Q4, &A24);
        if (fp2_is_zero(&Q.x)) {
            xisog_4_singular(&A24, Q4, A24, K);
            xeval_4_singular(points, points, length, Q4, K);
            xeval_4_singular(&Q, &phi->kernel, 1, Q4, K);
        }
        else {
            xisog_4(&A24, Q4, K);
            xeval_4(points, points, length, K);
            xeval_4(&Q, &phi->kernel, 1, K);
        }
        ec_eval_even_strategy(image, points, length, &A24, &Q, phi->length-2);
    }

void ec_eval_even_nonzero(ec_curve_t* image, const ec_isog_even_t* phi,
    ec_point_t* points, unsigned short length){
        ec_point_t Q4, A24, K[3];
        copy_point(&Q4, &phi->kernel);
        AC_to_A24(&A24, &phi->curve);
        for(int i = 0; i < phi->length - 2; i++)
            xDBLv2(&Q4, &Q4, &A24);
        xisog_4(&A24, Q4, K);
        xeval_4(points, points, length, K);
        xeval_4(&Q4, &phi->kernel, 1, K);
        ec_eval_even_strategy(image, points, length, &A24, &Q4, phi->length-2);
    }

static void ec_eval_even_strategy(ec_curve_t* image, ec_point_t* points, unsigned short points_len,
    ec_point_t* A24, const ec_point_t *kernel, const int isog_len){
    
    assert(isog_len == POWER_OF_2-2);
        
    uint8_t log2_of_e, tmp;
    fp2_t t0;
    digit_t e_half = (isog_len)>>1;
    for(tmp = e_half, log2_of_e = 0; tmp > 0; tmp>>=1, ++log2_of_e);
    log2_of_e *= 2; // In order to ensure each splits is at most size log2_of_e

    ec_point_t SPLITTING_POINTS[log2_of_e], K2, K[3];
    copy_point(&SPLITTING_POINTS[0], kernel);

    int strategy = 0,    // Current element of the strategy to be used
    i, j;

    int BLOCK = 0,       // Keeps track of point order
    current = 0;         // Number of points being carried
    int XDBLs[log2_of_e]; // Number of doubles performed

    // If walk length is odd, we start with a 2-isogeny
    if(isog_len & 1){
        copy_point(&SPLITTING_POINTS[1], &SPLITTING_POINTS[0]);
        for(i = 0; i < isog_len-1; i++)
            xDBLv2(&SPLITTING_POINTS[1], &SPLITTING_POINTS[1], A24);
        xisog_2(A24, SPLITTING_POINTS[1], K);
        xeval_2(SPLITTING_POINTS, SPLITTING_POINTS, 1, K);
        xeval_2(points, points, points_len, K);
    }
    
    // Chain of 4-isogenies
    for(j = 0; j < (e_half - 1); j++)
    {   
        // Get the next point of order 4
        while (BLOCK != (e_half -  1 - j) )
        {
            // A new split will be added
            current += 1;
            // We set the seed of the new split to be computed and saved
            copy_point(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current - 1]);
            for(i = 0; i < 2*STRATEGY4[strategy]; i++)
                xDBLv2(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], A24);
            XDBLs[current] = STRATEGY4[strategy];  // The number of doublings performed is saved
            BLOCK += STRATEGY4[strategy];          // BLOCK is increased by the number of doublings performed
            strategy += 1;                  // Next, we move to the next element of the strategy
        }

        // Evaluate 4-isogeny
        xisog_4(A24, SPLITTING_POINTS[current], K);
        xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, K);
        xeval_4(points, points, points_len, K);

        BLOCK -= XDBLs[current];  
        XDBLs[current] = 0;      
        current -= 1;            
    }

    // Final 4-isogeny
    xisog_4(A24, SPLITTING_POINTS[current], K);
    xeval_4(points, points, points_len, K);

    // Output curve in the form (A:C)
    A24_to_AC(image, A24);
}

void ec_eval_even_strategy_smart(ec_point_t *Pa, ec_point_t* A24out,
    const ec_point_t* A24in, const ec_point_t *kernel){
        
    uint8_t log2_of_e, tmp;
    fp2_t t0;
    digit_t e_half = (POWER_OF_2)>>1;
    for(tmp = e_half, log2_of_e = 0; tmp > 0; tmp>>=1, ++log2_of_e);
    log2_of_e *= 2; // In order to ensure each splits is at most size log2_of_e

    ec_point_t SPLITTING_POINTS[log2_of_e], K2, B24, K[3];
    copy_point(&SPLITTING_POINTS[0], kernel);
    copy_point(&B24, A24in);

    int strategy = 0,    // Current element of the strategy to be used
    i, j;

    int BLOCK = 0,       // Keeps track of point order
    current = 0;         // Number of points being carried
    int XDBLs[log2_of_e]; // Number of doubles performed

    // If walk length is odd, we start with a 2-isogeny
    if(POWER_OF_2 & 1){
        copy_point(&SPLITTING_POINTS[1], &SPLITTING_POINTS[0]);
        for(i = 0; i < POWER_OF_2-1; i++)
            xDBLv2(&SPLITTING_POINTS[1], &SPLITTING_POINTS[1], &B24);
        xisog_2(&B24, SPLITTING_POINTS[1], K);
        xeval_2(SPLITTING_POINTS, SPLITTING_POINTS, 1, K);
    }
    
    // Chain of 4-isogenies
    for(j = 0; j < (e_half - 1); j++)
    {   
        // Get the next point of order 4
        while (BLOCK != (e_half -  1 - j) )
        {
            // A new split will be added
            current += 1;
            // We set the seed of the new split to be computed and saved
            copy_point(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current - 1]);
            for(i = 0; i < 2*STRATEGY4_SMART[strategy]; i++)
                xDBLv2(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], &B24);
            XDBLs[current] = STRATEGY4_SMART[strategy];  // The number of doublings performed is saved
            BLOCK += STRATEGY4_SMART[strategy];          // BLOCK is increased by the number of doublings performed
            strategy += 1;                  // Next, we move to the next element of the strategy
        }

        // Evaluate 4-isogeny
        if((j == 0) && !(POWER_OF_2 & 1)){ // if this is the first isogeny we need to check that it is not over (0,0)
            xDBLv2(&K2, &SPLITTING_POINTS[current], &B24);
            if(fp2_is_zero(&K2.x)){
                xisog_4_singular(&B24, SPLITTING_POINTS[current], B24, K);
                xeval_4_singular(SPLITTING_POINTS, SPLITTING_POINTS, current, SPLITTING_POINTS[current], K);
            }
            else{
                xisog_4(&B24, SPLITTING_POINTS[current], K);
                xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, K);
            }
        }
        else{
            xisog_4(&B24, SPLITTING_POINTS[current], K);
            xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, K);
        }

        BLOCK -= XDBLs[current];  
        XDBLs[current] = 0;      
        current -= 1;            
    }

    // Final 4-isogeny replaced by a 2-isogeny only
    xDBLv2(&K2, &SPLITTING_POINTS[current], &B24);
    xisog_2(A24out, K2, K);
    xeval_2(Pa, &SPLITTING_POINTS[current], 1, K);
}

void ec_eval_even_strategy_uncompressed(ec_point_t* A24out, ec_point_t* second_to_last_j_inv, ec_point_t* A24in, const ec_point_t *kernel){
    
    uint8_t log2_of_e, tmp;
    fp2_t t0;
    digit_t e_half = (POWER_OF_2)>>1;
    for(tmp = e_half, log2_of_e = 0; tmp > 0; tmp>>=1, ++log2_of_e);
    log2_of_e *= 2; // In order to ensure each splits is at most size log2_of_e

    ec_point_t SPLITTING_POINTS[log2_of_e], K2, B24, K[3];
    copy_point(&SPLITTING_POINTS[0], kernel);
    copy_point(&B24, A24in);

    int strategy = 0,    // Current element of the strategy to be used
    i, j;

    int BLOCK = 0,       // Keeps track of point order
    current = 0;         // Number of points being carried
    int XDBLs[log2_of_e]; // Number of doubles performed
    
    // Chain of 4-isogenies
    for(j = 0; j < (e_half - 1); j++)
    {   
        // Get the next point of order 4
        while (BLOCK != (e_half -  1 - j) )
        {
            // A new split will be added
            current += 1;
            // We set the seed of the new split to be computed and saved
            copy_point(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current - 1]);
            if( current == 1 ){
            for(i = 0; i < 2*STRATEGY4_SMART[strategy] - 1; i++)
                xDBLv2(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], &B24);
            }
            else{
            for(i = 0; i < 2*STRATEGY4_SMART[strategy]; i++)
                xDBLv2(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], &B24);
            }
            XDBLs[current] = STRATEGY4_SMART[strategy];  // The number of doublings performed is saved
            BLOCK += STRATEGY4_SMART[strategy];          // BLOCK is increased by the number of doublings performed
            strategy += 1;                  // Next, we move to the next element of the strategy
        }

        // Evaluate 4-isogeny
        if( j == 0 ){// For the first isogeny we need to check if kernel is over (0,0)
            xDBLv2(&K2, &SPLITTING_POINTS[current], &B24);
            if(fp2_is_zero(&K2.x)){
                xisog_4_singular(&B24, SPLITTING_POINTS[current], B24, K);
                xeval_4_singular(SPLITTING_POINTS, SPLITTING_POINTS, current, SPLITTING_POINTS[current], K);
            }
            else{
                xisog_4(&B24, SPLITTING_POINTS[current], K);
                xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, K);
            }   
        }
        else{
            xisog_4(&B24, SPLITTING_POINTS[current], K);
            xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, K);
        }

        BLOCK -= XDBLs[current];  
        XDBLs[current] = 0;      
        current -= 1;            
    }

    // Record second-to-last j-invariant
    ec_curve_t E;
    A24_to_AC(&E, &B24);
    ec_j_inv_proj(second_to_last_j_inv, &E);

    // Last 2-isogeny
    xisog_2(A24out, SPLITTING_POINTS[current], K);
}

void ec_eval_even_strategy_parallel(ec_point_t* A24out, ec_point_t *A24mid, ec_point_t *K2, ec_point_t* A24in, const ec_point_t *kernel){
    
    uint8_t log2_of_e, tmp;
    fp2_t t0;
    digit_t e_half = (POWER_OF_2)>>1;
    for(tmp = e_half, log2_of_e = 0; tmp > 0; tmp>>=1, ++log2_of_e);
    log2_of_e *= 2; // In order to ensure each splits is at most size log2_of_e

    ec_point_t SPLITTING_POINTS[log2_of_e], B24, K[3];
    copy_point(&SPLITTING_POINTS[0], kernel);
    copy_point(&B24, A24in);

    int strategy = 0,    // Current element of the strategy to be used
    i, j;

    int BLOCK = 0,       // Keeps track of point order
    current = 0;         // Number of points being carried
    int XDBLs[log2_of_e]; // Number of doubles performed
    
    // Chain of 4-isogenies
    for(j = 0; j < (e_half - 1); j++)
    {   
        // Get the next point of order 4
        while (BLOCK != (e_half -  1 - j) )
        {
            // A new split will be added
            current += 1;
            // We set the seed of the new split to be computed and saved
            copy_point(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current - 1]);
            if( current == 1 ){
            for(i = 0; i < 2*STRATEGY4_SMART[strategy] - 1; i++)
                xDBLv2(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], &B24);
            }
            else{
            for(i = 0; i < 2*STRATEGY4_SMART[strategy]; i++)
                xDBLv2(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], &B24);
            }
            XDBLs[current] = STRATEGY4_SMART[strategy];  // The number of doublings performed is saved
            BLOCK += STRATEGY4_SMART[strategy];          // BLOCK is increased by the number of doublings performed
            strategy += 1;                  // Next, we move to the next element of the strategy
        }

        // Evaluate 4-isogeny
        if( j == 0 ){// For the first isogeny we need to check if kernel is over (0,0)
            xDBLv2(K2, &SPLITTING_POINTS[current], &B24);
            if(fp2_is_zero(&K2->x)){
                xisog_4_singular(&B24, SPLITTING_POINTS[current], B24, K);
                xeval_4_singular(SPLITTING_POINTS, SPLITTING_POINTS, current, SPLITTING_POINTS[current], K);
            }
            else{
                xisog_4(&B24, SPLITTING_POINTS[current], K);
                xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, K);
            }   
        }
        else{
            xisog_4(&B24, SPLITTING_POINTS[current], K);
            xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, K);
        }

        // After lambda/2 steps we record the curve
        if( j == (POWER_OF_2_SECPAR >> 1) - 1){
            copy_point(A24mid, &B24);
        }

        BLOCK -= XDBLs[current];  
        XDBLs[current] = 0;      
        current -= 1;            
    }

    // Last 2-isogeny
    xisog_2(A24out, SPLITTING_POINTS[current], K);
}

void ec_eval_even_strategy_chal(ec_curve_t* image, ec_point_t* push_point,
    ec_point_t* A24, const ec_point_t *kernel){
        
    uint8_t log2_of_e, tmp;
    fp2_t t0;
    digit_t e_half = (POWER_OF_2_SECPAR)>>1;
    for(tmp = e_half, log2_of_e = 0; tmp > 0; tmp>>=1, ++log2_of_e);
    log2_of_e *= 2; // In order to ensure each splits is at most size log2_of_e

    ec_point_t SPLITTING_POINTS[log2_of_e], K2, K[3];
    copy_point(&SPLITTING_POINTS[0], kernel);

    int strategy = 0,    // Current element of the strategy to be used
    i, j;

    int BLOCK = 0,       // Keeps track of point order
    current = 0;         // Number of points being carried
    int XDBLs[log2_of_e]; // Number of doubles performed
    
    // Chain of 4-isogenies
    for(j = 0; j < (e_half - 1); j++)
    {   
        // Get the next point of order 4
        while (BLOCK != (e_half -  1 - j) )
        {
            // A new split will be added
            current += 1;
            // We set the seed of the new split to be computed and saved
            copy_point(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current - 1]);
            for(i = 0; i < 2*STRATEGY4_CHAL[strategy]; i++)
                xDBLv2(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], A24);
            XDBLs[current] = STRATEGY4_CHAL[strategy];  // The number of doublings performed is saved
            BLOCK += STRATEGY4_CHAL[strategy];          // BLOCK is increased by the number of doublings performed
            strategy += 1;                  // Next, we move to the next element of the strategy
        }

        // Evaluate 4-isogeny
        xisog_4(A24, SPLITTING_POINTS[current], K);
        xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, K);
        xeval_4(push_point, push_point, 1, K);

        BLOCK -= XDBLs[current];  
        XDBLs[current] = 0;      
        current -= 1;            
    }

    // Final 4-isogeny
    xisog_4(A24, SPLITTING_POINTS[current], K);
    xeval_4(push_point, push_point, 1, K);

    // Output curve in the form (A:C)
    A24_to_AC(image, A24);
}

void ec_eval_even_strategy_chal_uncompressed(ec_curve_t* image, ec_point_t* second_to_last_j_inv,
    const ec_point_t* A24, const ec_point_t *kernel){
        
    uint8_t log2_of_e, tmp;
    fp2_t t0;
    digit_t e_half = (POWER_OF_2_SECPAR)>>1;
    for(tmp = e_half, log2_of_e = 0; tmp > 0; tmp>>=1, ++log2_of_e);
    log2_of_e *= 2; // In order to ensure each splits is at most size log2_of_e

    ec_point_t SPLITTING_POINTS[log2_of_e], K2, B24, K[3];
    copy_point(&SPLITTING_POINTS[0], kernel);
    copy_point(&B24, A24);

    int strategy = 0,    // Current element of the strategy to be used
    i, j;

    int BLOCK = 0,       // Keeps track of point order
    current = 0;         // Number of points being carried
    int XDBLs[log2_of_e]; // Number of doubles performed
    
    // Chain of 4-isogenies
    for(j = 0; j < (e_half - 1); j++)
    {   
        // Get the next point of order 4
        while (BLOCK != (e_half -  1 - j) )
        {
            // A new split will be added
            current += 1;
            // We set the seed of the new split to be computed and saved
            copy_point(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current - 1]);
            for(i = 0; i < 2*STRATEGY4_CHAL[strategy]; i++)
                xDBLv2(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], &B24);
            XDBLs[current] = STRATEGY4_CHAL[strategy];  // The number of doublings performed is saved
            BLOCK += STRATEGY4_CHAL[strategy];          // BLOCK is increased by the number of doublings performed
            strategy += 1;                  // Next, we move to the next element of the strategy
        }

        // Evaluate 4-isogeny
        xisog_4(&B24, SPLITTING_POINTS[current], K);
        xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, K);

        BLOCK -= XDBLs[current];  
        XDBLs[current] = 0;      
        current -= 1;            
    }

    // Second-to-last 2-isogeny
    xDBLv2(&K2, &SPLITTING_POINTS[current], &B24);
    xisog_2(&B24, K2, K);
    xeval_2(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], 1, K);

    // Record second-to-last j-invariant
    A24_to_AC(image, &B24);
    ec_j_inv_proj(second_to_last_j_inv, image);

    // Last 2-isogeny
    xisog_2(&B24, SPLITTING_POINTS[current], K);
    A24_to_AC(image, &B24);
}

void ec_eval_odd(ec_curve_t* image, const ec_isog_odd_t* phi,
        ec_point_t* points, unsigned short length){
        
    ec_point_t ker_plus, ker_minus, P, Kp, A24, B24;
    int i,j,k;

#if defined(ENABLE_SIGN) && !defined(PQM4)
    ec_point_t K[sK_max];
#else
    ec_point_t K[3];
#endif

    AC_to_A24(&A24, &phi->curve);

    // Isogenies with kernel in E[p+1]
    copy_point(&ker_plus, &phi->ker_plus);
    copy_point(&ker_minus, &phi->ker_minus);
    for(i = 0; i < P_LEN; i++){
        copy_point(&P, &ker_plus);
        for(j = i+1; j < P_LEN; j++){
            for(k = 0; k < phi->degree[j]; k++)
                xMULdac(&P, &P, DACS[j], DAC_LEN[j], &A24);
        }
        for(k = 0; k < phi->degree[i]; k++){
            copy_point(&Kp, &P);
            for(j = 0; j < phi->degree[i]-k-1; j++)
                xMULdac(&Kp, &Kp, DACS[i], DAC_LEN[i], &A24);
            kps(i, Kp, A24, K);
            xisog(&B24, i, A24, K);
            xeval(&P, i, P, A24, K);
            xeval(&ker_plus, i, ker_plus, A24, K);
            xeval(&ker_minus, i, ker_minus, A24, K);
            for(j = 0; j < length; j++)
                xeval(&points[j], i, points[j], A24, K);
            copy_point(&A24, &B24);
        }
    }

    // Isogenies with kernel in E[p-1]
    for(i = P_LEN; i < P_LEN+M_LEN; i++){
        copy_point(&P, &ker_minus);
        for(j = i+1; j < P_LEN+M_LEN; j++){
            for(k = 0; k < phi->degree[j]; k++)
                xMULdac(&P, &P, DACS[j], DAC_LEN[j], &A24);
        }
        for(k = 0; k < phi->degree[i]; k++){
            copy_point(&Kp, &P);
            for(j = 0; j < phi->degree[i]-k-1; j++)
                xMULdac(&Kp, &Kp, DACS[i], DAC_LEN[i], &A24);
            kps(i, Kp, A24, K);
            xisog(&B24, i, A24, K);
            xeval(&P, i, P, A24, K);
            xeval(&ker_minus, i, ker_minus, A24, K);
            for(j = 0; j < length; j++)
                xeval(&points[j], i, points[j], A24, K);
            copy_point(&A24, &B24);
        }
    }

    A24_to_AC(image, &A24);
}

void ec_curve_normalize(ec_curve_t *new, ec_isom_t *isom, const ec_curve_t *old){
    fp2_t t0, t1, t2, t3, t4, t5;
    // Compute the other solutions:
    // A'^2 = [ sqrt(A^2-4C^2)*(9C^2-A^2) +- (A^3-3AC^2) ] / [ 2C^2*sqrt(A^2-4C^2) ]
    fp2_sqr(&t0, &old->C);      //C^2
    fp2_add(&t1, &t0, &t0);     //2C^2
    fp2_add(&t2, &t1, &t1);     //4C^2
    fp2_sqr(&t3, &old->A);      //A^2
    fp2_sub(&t2, &t3, &t2);     //A^2-4C^2
    fp2_sqrt(&t2);              //sqrt(A^2-4C^2)
    fp2_add(&t0, &t0, &t1);     //3C^2
    fp2_mul(&t1, &t2, &t1);     //2C^2*sqrt(A^2-4C^2)
    fp2_sub(&t5, &t3, &t0);     //A^2-3C^2
    fp2_mul(&t5, &t5, &old->A);     //A^3-3AC^2
    fp2_add(&t4, &t0, &t0);     //6C^2
    fp2_add(&t0, &t4, &t0);     //9C^2
    fp2_sub(&t0, &t0, &t3);     //9C^2-A^2
    fp2_add(&t3, &t3, &t3);     //2A^2
    fp2_mul(&t3, &t3, &t2);     //2A^2*sqrt(A^2-4C^2)
    fp2_mul(&t2, &t2, &t0);     //sqrt(A^2-4C^2)*(9C^2-A^2)
    fp2_add(&t0, &t2, &t5);     //sqrt(A^2-4C^2)*(9C^2-A^2) + (A^3-3AC^2)
    fp2_sub(&t2, &t2, &t5);     //sqrt(A^2-4C^2)*(9C^2-A^2) - (A^3-3AC^2)
    fp2_inv(&t1);               //1/2C^2*sqrt(A^2-4C^2)
    fp2_mul(&t0, &t0, &t1);     // First solution
    fp2_mul(&t2, &t2, &t1);     // Second solution
    fp2_mul(&t1, &t3, &t1);     // Original solution

    // Chose the lexicographically first solution
    if(fp2_cmp(&t0, &t1)==1)
        fp2_copy(&t0, &t1);
    if(fp2_cmp(&t0, &t2)==1)
        fp2_copy(&t0, &t2);

    // Copy the solution
    fp2_sqrt(&t0);
    ec_curve_t E;
    fp2_copy(&E.A, &t0);
    fp_mont_setone(E.C.re);
    fp_set(E.C.im, 0);
    ec_isomorphism(isom, old, &E);
    fp2_copy(&new->A, &E.A);
    fp2_copy(&new->C, &E.C);
}

void ec_isomorphism(ec_isom_t* isom, const ec_curve_t* from, const ec_curve_t* to){
    fp2_t t0, t1, t2, t3, t4;
    fp2_mul(&t0, &from->A, &to->C);
    fp2_sqr(&t0, &t0);                  //fromA^2toC^2
    fp2_mul(&t1, &to->A, &from->C);
    fp2_sqr(&t1, &t1);                  //toA^2fromC^2
    fp2_mul(&t2, &to->C, &from->C);
    fp2_sqr(&t2, &t2);                  //toC^2fromC^2
    fp2_add(&t3, &t2, &t2);
    fp2_add(&t2, &t3, &t2);             //3toC^2fromC^2
    fp2_sub(&t3, &t2, &t0);             //3toC^2fromC^2-fromA^2toC^2
    fp2_sub(&t4, &t2, &t1);             //3toC^2fromC^2-toA^2fromC^2
    fp2_inv(&t3);
    fp2_mul(&t4, &t4, &t3);
    fp2_sqrt(&t4);                      //lambda^2 constant for SW isomorphism
    fp2_sqr(&t3, &t4);
    fp2_mul(&t3, &t3, &t4);             //lambda^6

    // Check sign of lambda^2, such that lambda^6 has the right sign
    fp2_sqr(&t0, &from->C);
    fp2_add(&t1, &t0, &t0);
    fp2_add(&t1, &t1, &t1);
    fp2_add(&t1, &t1, &t1);
    fp2_add(&t0, &t0, &t1); // 9fromC^2
    fp2_sqr(&t2, &from->A);
    fp2_add(&t2, &t2, &t2); // 2fromA^2
    fp2_sub(&t2, &t2, &t0);
    fp2_mul(&t2, &t2, &from->A); // -9fromC^2fromA+2fromA^3
    fp2_sqr(&t0, &to->C);
    fp2_mul(&t0, &t0, &to->C);
    fp2_mul(&t2, &t2, &t0);     //toC^3* [-9fromC^2fromA+2fromA^3]
    fp2_mul(&t3, &t3, &t2);             //lambda^6*(-9fromA+2fromA^3)*toC^3
    fp2_sqr(&t0, &to->C);
    fp2_add(&t1, &t0, &t0);
    fp2_add(&t1, &t1, &t1);
    fp2_add(&t1, &t1, &t1);
    fp2_add(&t0, &t0, &t1); // 9toC^2
    fp2_sqr(&t2, &to->A);
    fp2_add(&t2, &t2, &t2); // 2toA^2
    fp2_sub(&t2, &t2, &t0);
    fp2_mul(&t2, &t2, &to->A); // -9toC^2toA+2toA^3
    fp2_sqr(&t0, &from->C);
    fp2_mul(&t0, &t0, &from->C);
    fp2_mul(&t2, &t2, &t0);     //fromC^3* [-9toC^2toA+2toA^3]
    if(!fp2_is_equal(&t2, &t3))
        fp2_neg(&t4, &t4);

    // Mont -> SW -> SW -> Mont
    fp_mont_setone(t0.re);
    fp_set(t0.im, 0);
    fp2_add(&isom->D, &t0, &t0);
    fp2_add(&isom->D, &isom->D, &t0);
    fp2_mul(&isom->D, &isom->D, &from->C);
    fp2_mul(&isom->D, &isom->D, &to->C);
    fp2_mul(&isom->Nx, &isom->D, &t4);
    fp2_mul(&t4, &t4, &from->A);
    fp2_mul(&t4, &t4, &to->C);
    fp2_mul(&t0, &to->A, &from->C);
    fp2_sub(&isom->Nz, &t0, &t4);
}

void ec_iso_inv(ec_isom_t* isom){
    fp2_t tmp;
    fp2_copy(&tmp, &isom->D);
    fp2_copy(&isom->D, &isom->Nx);
    fp2_copy(&isom->Nx, &tmp);
    fp2_neg(&isom->Nz, &isom->Nz);
}

void ec_iso_eval(ec_point_t *P, ec_isom_t* isom){
    fp2_t tmp;
    fp2_mul(&P->x, &P->x, &isom->Nx);
    fp2_mul(&tmp, &P->z, &isom->Nz);
    fp2_sub(&P->x, &P->x, &tmp);
    fp2_mul(&P->z, &P->z, &isom->D);
}
