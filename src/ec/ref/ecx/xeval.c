#include "isog.h"
#include "ec.h"
#include <assert.h>

// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
// Traditional isogeny evaluation (xEVAL)

// Degree-2 isogeny evaluation with kenerl generated by P != (0, 0)
void xeval_2(ec_point_t* R, ec_point_t* const Q, const int lenQ)
{
	fp2_t t0, t1, t2;
	for(int j = 0; j < lenQ; j++){
		fp2_add(&t0, &Q[j].x, &Q[j].z);
		fp2_sub(&t1, &Q[j].x, &Q[j].z);
		fp2_mul(&t2, &K[0].x, &t1);
		fp2_mul(&t1, &K[0].z, &t0);
		fp2_add(&t0, &t2, &t1);
		fp2_sub(&t1, &t2, &t1);
		fp2_mul(&R[j].x, &Q[j].x, &t0);
		fp2_mul(&R[j].z, &Q[j].z, &t1);
	}
}

// Degree-4 isogeny evaluation with kenerl generated by P such that [2]P != (0, 0)
void xeval_4(ec_point_t* R, const ec_point_t* Q, const int lenQ)
{
	fp2_t t0, t1;

	for(int i = 0; i < lenQ; i++){
		fp2_add(&t0, &Q[i].x, &Q[i].z);
		fp2_sub(&t1, &Q[i].x, &Q[i].z);
		fp2_mul(&(R[i].x), &t0, &K[1].x);
		fp2_mul(&(R[i].z), &t1, &K[2].x);
		fp2_mul(&t0, &t0, &t1);
		fp2_mul(&t0, &t0, &K[0].x); 
		fp2_add(&t1, &(R[i].x), &(R[i].z));
		fp2_sub(&(R[i].z), &(R[i].x), &(R[i].z));
		fp2_sqr(&t1, &t1);
		fp2_sqr(&(R[i].z), &(R[i].z));
		fp2_add(&(R[i].x), &t0, &t1);
		fp2_sub(&t0, &t0, &(R[i].z));
		fp2_mul(&(R[i].x), &(R[i].x), &t1);
		fp2_mul(&(R[i].z), &(R[i].z), &t0);
	}
}

// Degree-4 isogeny evaluation with kenerl generated by P such that [2]P = (0, 0)
// Must call after xisog_4_singular
void xeval_4_singular(ec_point_t* R, const ec_point_t* Q, const int lenQ, const ec_point_t P)
{
	fp2_t t0, t1, t2;
	for(int i = 0; i < lenQ; i++){
		fp2_add(&t0, &Q[i].x, &Q[i].z);
		fp2_sub(&t2, &Q[i].x, &Q[i].z);
		fp2_sqr(&t0, &t0);
		fp2_sqr(&t2, &t2);
		fp2_sub(&R[i].z, &t0, &t2);
		if(fp2_is_equal(&P.x, &P.z)){
			// Branch for P = (+1,_)
			fp2_copy(&t1, &t2);
		}
		else{
			// Branch for P = (-1,_)
			fp2_copy(&t1, &t0);
			fp2_copy(&t0, &t2);
		}
		fp2_mul(&R[i].x, &R[i].z, &K[0].x);
		fp2_mul(&R[i].z, &R[i].z, &K[1].x);
		fp2_mul(&R[i].z, &R[i].z, &t1);
		fp2_mul(&t1, &t1, &K[0].z);
		fp2_add(&R[i].x, &R[i].x, &t1);
		fp2_mul(&R[i].x, &R[i].x, &t0);
	}
}

// Degree-3 isogeny evaluation with kenerl generated by P.
// Must be called after kps_3(P).
void xeval_3(ec_point_t* R, const ec_point_t Q)
{
	fp2_t t0, t1, t2;
	fp2_add(&t0, &Q.x, &Q.z);
	fp2_sub(&t1, &Q.x, &Q.z);
	fp2_mul(&t0, &K[0].x, &t0);
	fp2_mul(&t1, &K[1].x, &t1);
	fp2_add(&t2, &t0, &t1);
	fp2_sub(&t0, &t1, &t0);
	fp2_sqr(&t2, &t2);
	fp2_sqr(&t0, &t0);
	fp2_mul(&R->x, &Q.x, &t2);
	fp2_mul(&R->z, &Q.z, &t0);
}

#if defined(ENABLE_SIGN)

// CrissCross procedure as described in Hisil and Costello paper
void CrissCross(fp2_t *r0, fp2_t *r1, fp2_t const alpha, fp2_t const beta, fp2_t const gamma, fp2_t const delta)
{
	fp2_t t_1, t_2;

	fp2_mul(&t_1, &alpha, &delta);
    fp2_mul(&t_2, &beta, &gamma);
	fp2_add(&*r0, &t_1, &t_2);
	fp2_sub(&*r1, &t_1, &t_2);
}

// Isogeny evaluation on Montgomery curves
// Recall: K has been computed in Twisted Edwards model and none extra additions are required.
void xeval_t(ec_point_t* Q, int i, ec_point_t const P)
{
	int j;
	int d = ((int)TORSION_ODD_PRIMES[i] - 1) / 2;	// Here, l = 2d + 1

	fp2_t R0, R1, S0, S1, T0, T1;
	fp2_add(&S0, &P.x, &P.z);
	fp2_sub(&S1, &P.x, &P.z);

	CrissCross(&R0, &R1, K[0].z, K[0].x, S0, S1);
	for (j = 1; j < d; j++)
	{
		CrissCross(&T0, &T1, K[j].z, K[j].x, S0, S1);
		fp2_mul(&R0, &T0, &R0);
		fp2_mul(&R1, &T1, &R1);
	};

	fp2_sqr(&R0, &R0);
	fp2_sqr(&R1, &R1);

	fp2_mul(&(Q->x), &P.x, &R0);
	fp2_mul(&(Q->z), &P.z, &R1);
}

// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
// Isogeny evaluation (xEVAL) used in velu SQRT

void xeval_s(ec_point_t* Q, int i, ec_point_t const P, ec_point_t const A)
{
	// =================================================================================
	assert(TORSION_ODD_PRIMES[i] > gap);     // Ensuring velusqrt is used for l_i > gap
	sI = sizeI[i];          // size of I
	sJ = sizeJ[i];          // size of J
	sK = sizeK[i];          // size of K

	assert(sI >= sJ);       // Ensuring #I >= #J
	assert(sK >= 0);        // Recall, it must be that #K >= 0
	assert(sJ > 1);         // ensuring sI >= sJ > 1
	// =================================================================================

	// We require the curve coefficient A = A'/C ... well, a multiple of these ones
	fp2_t Ap;
	fp2_add(&Ap, &A.x, &A.x); // 2A' + 4C
	fp2_sub(&Ap, &Ap, &A.z);   // 2A'
	fp2_add(&Ap, &Ap, &Ap);     // 4A'

	//  --------------------------------------------------------------------------------------------------
	//                   ~~~~~~~~
	//                    |    | 
	// Computing E_J(W) = |    | [ F0(W, x([j]P)) * alpha^2 + F1(W, x([j]P)) * alpha + F2(W, x([j]P)) ]
	//                    j in J 
	// In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates
	// In particular, for a degree-l isogeny construction, we need alpha = X/Z and alpha = Z/X (i.e., 1/alpha)

	//fp2_t EJ_0[sJ][3]; // EJ_0[j][2] factors of one polynomial to be used in a resultant 

	fp2_t XZ_add, XZj_add,
	   XZ_sub, XZj_sub,
	   AXZ2,
	   CXZ2,
	   CX2Z2,
	   t1, t2;

	fp2_add(&XZ_add, &P.x, &P.z);	// X + Z
	fp2_sub(&XZ_sub, &P.x, &P.z);	// X - Z

	fp2_mul(&AXZ2, &P.x, &P.z);	// X * Z
	fp2_sqr(&t1, &P.x);		// X ^ 2
	fp2_sqr(&t2, &P.z);		// Z ^ 2

	fp2_add(&CX2Z2, &t1, &t2);		//      X^2 + Z^2
	fp2_mul(&CX2Z2, &CX2Z2, &A.z);	// C * (X^2 + Z^2)

	fp2_add(&AXZ2, &AXZ2, &AXZ2);	//       2 * (X * Z)
	fp2_mul(&CXZ2, &AXZ2, &A.z);	// C  * [2 * (X * Z)]
	fp2_mul(&AXZ2, &AXZ2, &Ap);		// A' * [2 * (X * Z)]

	int j;
	for (j = 0; j < sJ; j++)
	{
		fp2_add(&XZj_add, &J[j].x, &J[j].z);		// Xj + Zj
		fp2_sub(&XZj_sub, &J[j].x, &J[j].z);		// Xj - Zj

		fp2_mul(&t1, &XZ_sub, &XZj_add);			// (X - Z) * (Xj + Zj)
		fp2_mul(&t2, &XZ_add, &XZj_sub);			// (X + Z) * (Xj - Zj)

		// ...................................
		// Computing the quadratic coefficient
		fp2_sub(&EJ_0[j][2], &t1, &t2);			//       2 * [(X*Zj) - (Z*Xj)]
		fp2_sqr(&EJ_0[j][2], &EJ_0[j][2]);			//     ( 2 * [(X*Zj) - (Z*Xj)] )^2
		fp2_mul(&EJ_0[j][2], &A.z, &EJ_0[j][2]);		// C * ( 2 * [(X*Zj) - (Z*Xj)] )^2

		// ..................................
		// Computing the constant coefficient
		fp2_add(&EJ_0[j][0], &t1, &t2);			//       2 * [(X*Xj) - (Z*Zj)]
		fp2_sqr(&EJ_0[j][0], &EJ_0[j][0]);			//     ( 2 * [(X*Xj) - (Z*Zj)] )^2
		fp2_mul(&EJ_0[j][0], &A.z, &EJ_0[j][0]);		// C * ( 2 * [(X*Xj) - (Z*Zj)] )^2

		// ................................
		// Computing the linear coefficient
	
		// C * [ (-2*Xj*Zj)*(alpha^2 + 1) + (-2*alpha)*(Xj^2 + Zj^2)] + [A' * (-2*Xj*Zj) * (2*X*Z)] where alpha = X/Z
		fp2_add(&t1, &J[j].x, &J[j].z);			//      (Xj + Zj)
		fp2_sqr(&t1, &t1);					//      (Xj + Zj)^2
		fp2_add(&t1, &t1, &t1);				//  2 * (Xj + Zj)^2
		fp2_add(&t1, &t1, &XZJ4[j]);			//  2 * (Xj + Zj)^2 - (4*Xj*Zj) := 2 * (Xj^2 + Zj^2)
		fp2_mul(&t1, &t1, &CXZ2);				// [2 * (Xj^2 + Zj^2)] * (2 * [ C * (X * Z)])

		fp2_mul(&t2, &CX2Z2, &XZJ4[j]);			// [C * (X^2 + Z^2)] * (-4 * Xj * Zj)
		fp2_sub(&t1, &t2, &t1);				// [C * (X^2 + Z^2)] * (-4 * Xj * Zj) - [2 * (Xj^2 + Zj^2)] * (2 * [ C * (X * Z)])

		fp2_mul(&t2, &AXZ2, &XZJ4[j]);			// (2 * [A' * (X * Z)]) * (-4 * Xj * Zj)
		fp2_add(&EJ_0[j][1], &t1, &t2);			// This is our desired equation but multiplied by 2
		fp2_add(&EJ_0[j][1], &EJ_0[j][1], &EJ_0[j][1]);	// This is our desired equation but multiplied by 4
	};

        // ---------------------------------------------------------------------
        // The faster way for multiplying is using a divide-and-conquer approach

	// product tree of EJ_0 (we only require the root)
	product_tree_LENFeq3(ptree_EJ, deg_ptree_EJ, 0, EJ_0, sJ);
	assert( deg_ptree_EJ[0] == (2*sJ) );
	if (!scaled)
	{
		// unscaled remainder tree approach
		multieval_unscaled(leaves, ptree_EJ[0], 2*sJ + 1, rtree_hI, (const fp2_t*)rtree_A, ptree_hI, deg_ptree_hI, 0, sI);
	}
	else
	{
		// scaled remainder tree approach
		fp2_t G[sI_max], G_rev[sI_max];
		poly_redc(G, ptree_EJ[0], 2*sJ + 1, ptree_hI[0], sI + 1, R0, A0);
		for (j = 0; j < sI; j++)
			fp2_copy(&G_rev[j], &G[sI - 1 - j]);

		poly_mul_middle(G_rev, G_rev, sI, R0, sI);
		for (j = 0; j < sI; j++)
			fp2_copy(&G[j], &G_rev[sI - 1 - j]);

		multieval_scaled(leaves, G, ptree_hI, deg_ptree_hI, 0, sI);
        };

	// Finally, we must multiply the leaves of the outpur of remainders
	fp2_t r0;
	product(&r0, (const fp2_t*)leaves, sI);
	// EJ_1 is just reverting the ordering in the coefficients of EJ_0
	for (j = 0; j < sJ; j++){
		fp2_copy(&t1, &ptree_EJ[0][j]);
		fp2_copy(&ptree_EJ[0][j], &ptree_EJ[0][2*sJ - j]);
		fp2_copy(&ptree_EJ[0][2*sJ - j], &t1);
	}

	if (!scaled)
	{
		// unscaled remainder tree approach
		multieval_unscaled(leaves, ptree_EJ[0], 2*sJ + 1, rtree_hI, (const fp2_t*)rtree_A, ptree_hI, deg_ptree_hI, 0, sI);
	}
	else
	{
		// scaled remainder tree approach
		fp2_t G[sI_max], G_rev[sI_max];
		poly_redc(G, ptree_EJ[0], 2*sJ + 1, ptree_hI[0], sI + 1, R0, A0);
		for (j = 0; j < sI; j++)
			fp2_copy(&G_rev[j], &G[sI - 1 - j]);

		poly_mul_middle(G_rev, G_rev, sI, R0, sI);
		for (j = 0; j < sI; j++)
			fp2_copy(&G[j], &G_rev[sI - 1 - j]);

		multieval_scaled(leaves, G, ptree_hI, deg_ptree_hI, 0, sI);
        };
	clear_tree(ptree_EJ, 0, sJ);
	// Finally, we must multiply the leaves of the outpur of remainders
	fp2_t r1;
	product(&r1, (const fp2_t*)leaves, sI);

	// -------------------------------
	// Sometimes the public value sK is equal to zero,
	// Thus for avoing runtime error we add one when sK =0
	fp2_t hK_0[sK_max + 1], hK_1[sK_max + 1], hk_0, hk_1;
	for (j = 0; j < sK; j++)
	{
		fp2_add(&XZj_add, &K[j].x, &K[j].z);	// Xk + Zk
		fp2_sub(&XZj_sub, &K[j].x, &K[j].z);	// Xk - Zk
		fp2_mul(&t1, &XZ_sub, &XZj_add);		// (X - Z) * (Xk + Zk)
		fp2_mul(&t2, &XZ_add, &XZj_sub);		// (X + Z) * (Xk - Zk)

		// Case alpha = X/Z
		fp2_sub(&hK_0[j], &t1, &t2);		// 2 * [(X*Zk) - (Z*Xk)]

		// Case 1/alpha = Z/X
		fp2_add(&hK_1[j], &t1, &t2);		// 2 * [(X*Xk) - (Z*Zk)]
	};

	// hk_0 <- use product to mulitiply all the elements in hK_0
	product(&hk_0, (const fp2_t*)hK_0, sK);
	// hk_1 <- use product to mulitiply all the elements in hK_1
	product(&hk_1, (const fp2_t*)hK_1, sK);

	// ---------------------------------------------------------------------------------
	// Now, unifying all the computations
	fp2_mul(&t1, &hk_1, &r1);				// output of algorithm 2 with 1/alpha = Z/X and without the demoninator
	fp2_sqr(&t1, &t1);
	fp2_mul(&(Q->x), &t1, &P.x);

	fp2_mul(&t2, &hk_0, &r0);				// output of algorithm 2 with alpha = X/Z and without the demoninator
	fp2_sqr(&t2, &t2);
	fp2_mul(&(Q->z), &t2, &P.z);
}

#endif
