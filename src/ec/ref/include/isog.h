#ifndef _ISOG_H_
#define _ISOG_H_

#include "curve_extras.h"
#include "assert.h"

#if defined(ENABLE_SIGN) && !defined(PQM4)
#include "poly.h"

extern int sI, sJ, sK;	// Sizes of each current I, J, and K	

extern fp2_t I[sI_max][2],		// I plays also as the linear factors of the polynomial h_I(X)
			EJ_0[sJ_max][3], EJ_1[sJ_max][3];	// To be used in xisog y xeval

extern ec_point_t J[sJ_max], K[sK_max];		// Finite subsets of the kernel
extern fp2_t XZJ4[sJ_max],		// -4* (Xj * Zj) for each j in J, and x([j]P) = (Xj : Zj)
    rtree_A[(1 << (ceil_log_sI_max+2)) - 1],		// constant multiple of the reciprocal tree computation
    A0;			// constant multiple of the reciprocal R0

extern fp2_t ptree_hI[],		// product tree of h_I(X)
     rtree_hI[],		// reciprocal tree of h_I(X)
     ptree_EJ[];		// product tree of E_J(X)
     
extern fp2_t R0[2*sJ_max + 1];		// Reciprocal of h_I(X) required in the scaled remainder tree approach

extern int deg_ptree_hI[(1 << (ceil_log_sI_max+1)) - 1],	// degree of each noed in the product tree of h_I(X)
    deg_ptree_EJ[(1 << (ceil_log_sJ_max+1)) - 1];	// degree of each node in the product tree of E_J(X)

extern fp2_t leaves[sI_max];		// leaves of the remainder tree, which are required in the Resultant computation

#endif

void kps_4(ec_point_t const P, ec_point_t K[3]);
void kps_3(ec_point_t const P, ec_point_t K[3]);

void xisog_4(ec_point_t* B, ec_point_t const P, ec_point_t K[3]);			// degree-4 isogeny construction
void xisog_4_singular(ec_point_t* B24, ec_point_t const P, ec_point_t A24, ec_point_t K[3]);
void xisog_2(ec_point_t* B, ec_point_t const P, ec_point_t K[3]);			// degree-2 isogeny construction
void xisog_3(ec_point_t* B, ec_point_t K[3]);			// degree-3 isogeny construction

void xeval_4(ec_point_t* R, const ec_point_t* Q, const int lenQ, ec_point_t K[3]);					// degree-4 isogeny evaluation
void xeval_4_singular(ec_point_t* R, const ec_point_t* Q, const int lenQ, const ec_point_t P, ec_point_t K[3]);
void xeval_2(ec_point_t* R, ec_point_t* const Q, const int lenQ, ec_point_t K[3]);	// degree-2 isogeny evaluation
void xeval_3(ec_point_t* R, ec_point_t const Q, ec_point_t K[3]);	// degree-3 isogeny evaluation

// Strategy-based 4-isogeny chain
static void ec_eval_even_strategy(ec_curve_t* image, ec_point_t* points, unsigned short points_len,
    ec_point_t* A24, const ec_point_t *kernel, const int isog_len);

#ifndef ENABLE_SIGN
// Odd isogenies are always degree 3
static inline void kps(int i, ec_point_t const P, ec_point_t const A, ec_point_t* K)	
{
	assert(i==0);
	kps_3(P, K);
}

static inline void xisog(ec_point_t* B, int i, ec_point_t const A, ec_point_t* K)
{
	xisog_3(B, K);
}

static inline void xeval(ec_point_t* Q, int i, ec_point_t const P, ec_point_t const A, ec_point_t* K)
{
	xeval_3(Q, P, K);
}

#else

void eds2mont(ec_point_t* P);						// mapping from Twisted edwards into Montogmery
void yadd(ec_point_t* R, ec_point_t* const P, ec_point_t* const Q, ec_point_t* const PQ);	// differential addition on Twisted edwards model
void CrissCross(fp2_t *r0, fp2_t *r1, fp2_t const alpha, fp2_t const beta, fp2_t const gamma, fp2_t const delta);

void kps_t(int i, ec_point_t const P, ec_point_t const A, ec_point_t* K);	// tvelu formulae
void kps_s(int i, ec_point_t const P, ec_point_t const A, ec_point_t* K);	// svelu formulae

void xisog_t(ec_point_t* B, int i, ec_point_t const A, ec_point_t* K);	// tvelu formulae
void xisog_s(ec_point_t* B, int i, ec_point_t const A, ec_point_t* K);	// svelu formulae

void xeval_t(ec_point_t* Q, int i, ec_point_t const P, ec_point_t* K);			// tvelu formulae
void xeval_s(ec_point_t* Q, int i, ec_point_t const P, ec_point_t const A, ec_point_t* K);	// svelu formulae


// hybrid velu formulae
static inline void kps(int i, ec_point_t const P, ec_point_t const A, ec_point_t* K)	
{
	// Next branch only depends on a fixed public bound (named gap)
	if (TORSION_ODD_PRIMES[i] <= gap)
		kps_t(i, P, A, K);
	else
		kps_s(i, P, A, K);
}

static inline void xisog(ec_point_t* B, int i, ec_point_t const A, ec_point_t* K)
{
	// Next branch only depends on a fixed public bound (named gap)
	if (TORSION_ODD_PRIMES[i] <= gap)
		xisog_t(B, i, A, K);
	else
		xisog_s(B, i, A, K);
}

static inline void xeval(ec_point_t* Q, int i, ec_point_t const P, ec_point_t const A, ec_point_t* K)
{
	// Next branch only depends on a fixed public bound (named gap)
	if (TORSION_ODD_PRIMES[i] <= gap)
		xeval_t(Q, i, P, K);
	else
		xeval_s(Q, i, P, A, K);
}

#endif

void kps_clear(int i);	// Clear memory assigned by KPS

#endif
