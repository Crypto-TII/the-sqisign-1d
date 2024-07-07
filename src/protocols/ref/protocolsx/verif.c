#include <protocols.h>
#include <inttypes.h>
#include <assert.h>

static inline void copy_point(ec_point_t* A, const ec_point_t* B){
    fp2_copy(&A->x, &B->x);
    fp2_copy(&A->z, &B->z);
}

//XXX FIXME stolen from src/ec/opt/generic/test/isog-test.c
static void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf(HEX_FS, b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf(HEX_FS, b.im[i]);
    printf("\n");
}

static inline void A24_to_AC(ec_curve_t *E, ec_point_t const *A24)
{
    // (A:C) = ((A+2C)*2-4C : 4C)
    fp2_add(&E->A, &A24->x, &A24->x);
    fp2_sub(&E->A, &E->A, &A24->z);
    fp2_add(&E->A, &E->A, &E->A);
    fp2_copy(&E->C, &A24->z);
}

static void curve_print(char *name, ec_curve_t E){
    fp2_t a;
    fp2_copy(&a, &E.C);
    fp2_inv(&a);
    fp2_mul(&a, &a, &E.A);
    fp2_print(name, a);
}
//XXX

static void point_print(char *name, ec_point_t P){
    fp2_t a;
    if(fp2_is_zero(&P.z)){
        printf("%s = INF\n", name);
    }
    else{
    fp2_copy(&a, &P.z);
    fp2_inv(&a);
    fp2_mul(&a, &a, &P.x);
    fp2_print(name, a);
    }
}


void protocols_verif_unpack_chall(ec_curve_t *E2, ec_point_t *dual, const signature_t *sig, const public_key_t *pk)
{
    // Generate 2^f-basis
    ec_basis_t B2;
    ec_curve_to_basis_2(&B2, &pk->E);
    if(sig->zip.bit_first_step){
        // Swap P and Q
        ec_point_t tmp;
        copy_point(&tmp, &B2.P);
        copy_point(&B2.P, &B2.Q);
        copy_point(&B2.Q, &tmp);
    }

    // 2^f-isogeny chains
    ec_point_t K;
    ec_curve_t E;
    ec_isog_even_t isog2;
    fp2_copy(&E.A, &pk->E.A);
    fp2_copy(&E.C, &pk->E.C);
    isog2.length = POWER_OF_2;
    for(size_t i = 0; i < sig->zip.length-1; i++){        
        ec_ladder3pt(&K, sig->zip.zip_chain[i], &B2.Q, &B2.P, &B2.PmQ, &E);
        copy_point(&isog2.kernel, &K);
        fp2_copy(&isog2.curve.A, &E.A);
        fp2_copy(&isog2.curve.C, &E.C);
        ec_eval_even(&E, &isog2, &B2.P, 1);
        ec_complete_basis_2_singularP(&B2, &E, &B2.P);
    }   
    ec_ladder3pt(&K, sig->zip.zip_chain[sig->zip.length-1], &B2.Q, &B2.P, &B2.PmQ, &E);
    copy_point(&isog2.kernel, &K);
    fp2_copy(&isog2.curve.A, &E.A);
    fp2_copy(&isog2.curve.C, &E.C);
    ec_eval_even(&E, &isog2, &B2.P, 1);

    copy_point(dual,&B2.P);
    // Normalize curve
    ec_isom_t isom;
    ec_curve_normalize(E2, &isom, &E);
    ec_iso_eval(dual,&isom);

}

void protocols_verif_unpack_chall_smart(ec_point_t *Pa, ec_point_t *A24, const signature_t *sig, const public_key_smart_t *pk)
{
    ec_point_t K;

    copy_point(Pa, &pk->Pa);
    ec_A24_from_mont_root(A24, Pa);

    for(size_t i = 0; i < sig->zip.length; i++){
        ec_scalar_to_kernel_smart(&K, Pa, sig->zip.zip_chain[i], (i == 0) & (1 - sig->zip.bit_first_step));
        ec_eval_even_strategy_smart(Pa, A24, A24, &K);
    }
}

void protocols_verif_unpack_chall_uncompressed(ec_point_t *jinv, ec_point_t *jprev, const signature_uncompressed_t *sig, const public_key_t *pk)
{
    // Import the public key curve in (A+2C:4C) form
    ec_point_t A24;
    fp2_add(&A24.z, &pk->E.C, &pk->E.C);
    fp2_add(&A24.x, &pk->E.A, &A24.z);
    fp2_add(&A24.z, &A24.z, &A24.z);
    
    for (int i = 0; i < ZIP_CHAIN_LEN; i++)
        ec_eval_even_strategy_uncompressed(&A24, jprev, &A24, &sig->kernel_points[i]);

    // Challenge curve in (A:C) form
    ec_curve_t E;
    fp2_add(&E.A, &A24.x, &A24.x);
    fp2_sub(&E.A, &E.A, &A24.z);
    fp2_add(&E.C, &A24.z, &A24.z);

    // j invariant
    ec_j_inv_proj(jinv, &E);
}

int protocols_verif_from_chall(const signature_t *sig, ec_curve_t const *E2, const ec_point_t *dual, const unsigned char* m, size_t l)
{   

    // Generate 2^f and 3^g bases
    ec_basis_t B2, B3, B6;
    ec_curve_to_basis_6(&B6, E2);
    copy_point(&B3.P, &B6.P);
    copy_point(&B3.Q, &B6.Q);
    copy_point(&B3.PmQ, &B6.PmQ);
    for(size_t i = 0; i < POWER_OF_2; i++){
        ec_dbl(&B3.P, E2, &B3.P);
        ec_dbl(&B3.Q, E2, &B3.Q);
        ec_dbl(&B3.PmQ, E2, &B3.PmQ);
    }
    copy_point(&B2.P, &B6.P);
    copy_point(&B2.Q, &B6.Q);
    copy_point(&B2.PmQ, &B6.PmQ);
    for(size_t i = 0; i < POWER_OF_3; i++){
        ec_point_t tmp;
        ec_dbl(&tmp, E2, &B2.P);
        ec_add(&B2.P, &tmp, &B2.P, &B2.P);
        ec_dbl(&tmp, E2, &B2.Q);
        ec_add(&B2.Q, &tmp, &B2.Q, &B2.Q);
        ec_dbl(&tmp, E2, &B2.PmQ);
        ec_add(&B2.PmQ, &tmp, &B2.PmQ, &B2.PmQ);
    }

    bool const bit2 = sig->s.select23 & 1;
    bool const bit3 = sig->s.select23 & 2;

    if(!bit2){
        // Swap P2 and Q2
        ec_point_t tmp;
        copy_point(&tmp, &B2.P);
        copy_point(&B2.P, &B2.Q);
        copy_point(&B2.Q, &tmp);
    }
    if(!bit3){
        // Swap P3 and Q3
        ec_point_t tmp;
        copy_point(&tmp, &B3.P);
        copy_point(&B3.P, &B3.Q);
        copy_point(&B3.Q, &tmp);
    }

    // Kernels for the dual of the challenge isogeny
    ec_point_t K2, K3;
    ec_ladder3pt(&K2, sig->s.scalar2, &B2.P, &B2.Q, &B2.PmQ, E2);
    ec_ladder3pt(&K3, sig->s.scalar3, &B3.P, &B3.Q, &B3.PmQ, E2);

    // Evaluating the dual isogeny
    ec_isog_odd_t isog3={0};
    ec_point_t push_points[2];
    if (bit2 == bit3) {
        push_points[0] = bit2 ? B6.Q : B6.P;
    }
    else {
        const digit_t *scalars[2];
        scalars[bit3] = TORSION_PLUS_2POWER_DIGITS;
        scalars[bit2] = TORSION_PLUS_3POWER_DIGITS;
        ec_biscalar_mul(&push_points[0], E2, scalars[0], scalars[1], &B6);
    }
    copy_point(&push_points[1], &K3);
    ec_isog_even_t isog2;
    isog2.length = POWER_OF_2;
    copy_point(&isog2.kernel, &K2);


    // entering the cyclicity test zone
    ec_point_t test1,test2;
    test1 = *dual;
    test2 = K2;
    for (int i =0 ; i < TORSION_PLUS_EVEN_POWER-1; i++) {
        ec_dbl(&test1,E2,&test1);
        ec_dbl(&test2,E2,&test2);
    }
    assert(!fp2_is_zero(&test1.z));
    assert(!fp2_is_zero(&test2.z));
    if (ec_is_equal(&test1,&test2)) {
        return 0;
    }      

    fp2_copy(&isog2.curve.A, &E2->A);
    fp2_copy(&isog2.curve.C, &E2->C);
    ec_curve_t E;
    ec_eval_even(&E, &isog2, push_points, 2);
    isog3.degree[0] = POWER_OF_3;
    copy_point(&isog3.ker_plus, &push_points[1]);
    fp2_copy(&isog3.curve.A, &E.A);
    fp2_copy(&isog3.curve.C, &E.C);
    ec_eval_odd(&E, &isog3, push_points, 1);
assert(!fp2_is_zero(&E.C));

    // Normalize curve
    ec_curve_t E1;
    ec_isom_t isom;
    ec_curve_normalize(&E1, &isom, &E);
    ec_iso_eval(&push_points[0], &isom);
assert(!fp2_is_zero(&E1.C));

    // Generate 2^f*3^g basis
    ec_curve_to_basis_6(&B6, &E1);

    // Recover original challenge kernel point from the hash
    ec_point_t K;
    digit_vec_2_t scalars;
    hash_to_challenge(&scalars, &E1, m, l);
    ec_biscalar_mul(&K, &E1, scalars[0], scalars[1], &B6);
assert(!fp2_is_zero(&K.z));

    // Multiply Q by r
    ec_mul(&push_points[0], &E1, sig->r, &push_points[0]);
assert(!fp2_is_zero(&push_points[0].z));

    return ec_is_equal(&K, &push_points[0]);
}


int protocols_verif_from_chall_smart(const signature_t *sig, ec_point_t const *Pa, ec_point_t *A24, const unsigned char* m, size_t l)
{   
    ec_point_t K, push_point;
    ec_curve_t E;

    ec_scalar_to_kernel_secpar_smart(&K, &push_point, Pa, A24, sig->s.scalar2);

    ec_eval_even_strategy_chal(&E, &push_point, A24, &K);

assert(!fp2_is_zero(&E.C));

    // Normalize curve
    ec_curve_t E1;
    ec_isom_t isom;
    ec_curve_normalize(&E1, &isom, &E);
    ec_iso_eval(&push_point, &isom);
assert(!fp2_is_zero(&E1.C));

    // Recover original challenge kernel point from the hash
    hash_to_challenge_smart(&K, &E1, m, l);

assert(!fp2_is_zero(&K.z));

    // Multiply Q by r
    ec_mul(&push_point, &E1, sig->r, &push_point);
assert(!fp2_is_zero(&push_point.z));

    return ec_is_equal(&K, &push_point);
}


int protocols_verif_from_chall_uncompressed(const signature_uncompressed_t *sig, ec_point_t const *jinv, ec_point_t const *jprev, const unsigned char* m, size_t l)
{   
    // Compute the challenge isogeny
    ec_point_t Kchal, B24, j2, j2prev;
    ec_curve_t B;
    hash_to_challenge_smart(&Kchal, &sig->E_COM, m, l);
    fp2_add(&B24.z, &sig->E_COM.C, &sig->E_COM.C);
    fp2_add(&B24.x, &sig->E_COM.A, &B24.z);
    fp2_add(&B24.z, &B24.z, &B24.z);
    ec_eval_even_strategy_chal_uncompressed(&B, &j2prev, &B24, &Kchal);
    ec_j_inv_proj(&j2, &B);

    // Check that second-to-last j-invariants are different but last j-invariants are the same
    return(!ec_is_equal(jprev, &j2prev) && ec_is_equal(jinv, &j2));
}

/** @defgroup verification The verification protocol
 * @{
*/

/**
 * @brief Verifying a signature
 *
 * @param sig: the signature
 * @param pk the public key 
 * @param m the message
 * @returns a bit indicating if the verification succeeded  
    */
int protocols_verif(const signature_t *sig, const public_key_t *pk, const unsigned char* m, size_t l)
{
    //TODO general input checks?

    ec_curve_t E2;
    ec_point_t dual;
    protocols_verif_unpack_chall(&E2, &dual, sig, pk);

    return protocols_verif_from_chall(sig, &E2, &dual, m, l);
}

/**
 * @brief Verifying a signature with smart sampling and hashing
 *
 * @param sig: the signature
 * @param pk the public key 
 * @param m the message
 * @returns a bit indicating if the verification succeeded  
    */
int protocols_verif_smart(const signature_t *sig, const public_key_smart_t *pk, const unsigned char* m, size_t l)
{
    ec_point_t Pa, A24;
    protocols_verif_unpack_chall_smart(&Pa, &A24, sig, pk);
    return protocols_verif_from_chall_smart(sig, &Pa, &A24, m, l);
}

int protocols_verif_uncompressed(const signature_uncompressed_t *sig, const public_key_t *pk, const unsigned char* m, size_t l)
{
    ec_point_t jinv, jprev;
    protocols_verif_unpack_chall_uncompressed(&jinv, &jprev, sig, pk);
    return protocols_verif_from_chall_uncompressed(sig, &jinv, &jprev, m, l);
}