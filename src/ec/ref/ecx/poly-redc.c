#define _POLY_MUL_REDC_H_
#include "poly.h"
#include <assert.h>

void reciprocal(poly h, fp2_t *c, const poly f, const int lenf, const int n){
  
  // Writes a polynomial to h and a field element to c such that f*h = c mod x^n
  // REQUIRES h to have space for n terms
  // NOT responsible for terms in h beyond h[n-1]

  int i;

  // Case when f needs to be padded with zeroes
  if(n > lenf)
  {
    fp2_t fpad[n];
    for(i = 0; i < lenf; i++)
      fp2_copy(&fpad[i], &f[i]);
    for(i = lenf; i < n; i++)
      fp2_set(&fpad[i], 0);
    reciprocal(h, c, fpad, n, n);
    return;
  }

  // Trivial case
  if(n == 0)
  {
    fp2_set(&*c, 0);
    return;
  }

  // Case n = 1
  if(n == 1)
  {
    fp2_copy(&*c, &f[0]);
    fp_mont_setone(h[0].re);fp_set(h[0].im,0);
    return;
  }

  // Case n = 2
  if(n == 2)
  {
    fp2_sqr(&*c, &f[0]);
    fp2_copy(&h[0], &f[0]);
    fp2_neg(&h[1], &f[1]);
    return;
  }

  // Case n = 3
  if(n == 3)
  {
    fp2_t t0, t1;

    fp2_sqr(&t0, &f[1]);
    fp2_mul(&t1, &f[0], &f[2]);
    fp2_sub(&t1, &t1, &t0);
    fp2_mul(&t1, &t1, &f[0]);

    reciprocal(h, c, f, 2, 2);
    fp2_mul(&h[0], &h[0], &*c);
    fp2_mul(&h[1], &h[1], &*c);
    fp2_neg(&h[2], &t1);
    fp2_sqr(&*c, &*c);
    return;
  }

  // Case n = 4
  if(n == 4)
  {
    fp2_t t0, t1, t2, t3, g[2];

    reciprocal(g, &t3, f, 2, 2);
    fp2_sqr(&t0, &f[1]);
    fp2_mul(&t1, &g[0], &f[2]);
    fp2_mul(&t2, &g[0], &f[3]);
    fp2_mul(&h[1], &g[1], &f[2]);
    fp2_sub(&t0, &t1, &t0);
    fp2_add(&t1, &t2, &h[1]);
    fp2_mul(&t2, &t0, &g[0]);
    fp2_mul(&h[1], &t0, &g[1]);
    fp2_mul(&h[3], &t1, &g[0]);
    fp2_add(&h[3], &h[1], &h[3]);
    
    fp2_mul(&h[0], &g[0], &t3);
    fp2_mul(&h[1], &g[1], &t3);
    fp2_neg(&h[2], &t2);
    fp2_neg(&h[3], &h[3]);
    fp2_sqr(&*c, &t3);
    return;
  }


  // General case
  // Compute the reciprocal g mod x^m for m = ceil(n/2)
  // Then f*g-c is multiple of x^m so we only care about terms from m to n-1
  const int m = n - (n>>1);
  fp2_t g[m], t[m], t0;

  reciprocal(g, &t0, f, lenf, m);
  poly_mul_middle(t, g, m, f, n);
  poly_mul_low(t, n-m, g, m, &(t[2*m-n]), n-m);
  for(i = 0; i < m; i++)
    fp2_mul(&h[i], &g[i], &t0);
  for(i = m; i < n; i++)
    fp2_neg(&h[i], &t[i-m]);
  fp2_sqr(&*c, &t0);
  return;
}


void poly_redc(poly h, const poly g, const int leng, const poly f, const int lenf,//
	       const poly f_rev_inv, const fp2_t c)
{
  // Computes h(x) =  a * g(x) mod f(x) for some scalar a, writting lenf-1 terms to h.
  // REQUIRES an inverse f_rev_inv such that f_rev*f_rev_inv = c mod x^(leng-lenf+1),
  // where f_rev is the polynomial with the coefficients of f listed in reverse order.
  // The scalar a is equal to c, except for special cases:
  //    - If leng<lenf (no reduction needed) then a = 1
  //    - If lenf = leng = 2, then a = f[1] 
  //    - If lenf = leng = 3, then a = f[2] 
  //    - If lenf=2, leng=3 then a = 2*f[1]^2
  //
  // REQUIRES h to have space for lenf-1 terms
  // NOT responsible for terms in h beyond h[lenf-2]

  int i;
  
  // Case without reduction
  if(leng < lenf)
  {
    for(i = 0; i < leng; i++)
      fp2_copy(&h[i], &g[i]);
    for(i = leng; i < lenf-1; i++)
      fp2_set(&h[i], 0);
    return;
  }

  // Small cases for f linear
  if(lenf == 2)
  {
    if(leng == 2)
    {
      fp2_t t0;
      fp2_mul(&t0, &g[0], &f[1]);
      fp2_mul(&h[0], &g[1], &f[0]);
      fp2_sub(&h[0], &t0, &h[0]);
      return;
    }
    
    if(leng == 3)
    {
      fp2_t f0f1, f02, f12;
      fp2_sqr(&f02, &f[0]);
      fp2_sqr(&f12, &f[1]);
      fp2_sub(&f0f1, &f[0], &f[1]);
      fp2_sqr(&f0f1, &f0f1);
      fp2_sub(&f0f1, &f0f1, &f02);
      fp2_sub(&f0f1, &f0f1, &f12);
      fp2_add(&f02, &f02, &f02);
      fp2_add(&f12, &f12, &f12);
      fp2_mul(&f02, &f02, &g[2]);
      fp2_mul(&f12, &f12, &g[0]);
      fp2_mul(&f0f1, &f0f1, &g[1]);
      fp2_add(&h[0], &f02, &f12);
      fp2_add(&h[0], &h[0], &f0f1);
      return;
    }
  }

  // Small case for f cuadratic
  if(lenf == 3 && leng == 3)
  {
    fp2_t f2g1, f2g0, f1g2;
    fp2_mul(&f2g1, &g[1], &f[2]);
    fp2_mul(&f2g0, &g[0], &f[2]);
    fp2_mul(&f1g2, &g[2], &f[1]);
    fp2_mul(&h[0], &g[2], &f[0]);
    fp2_sub(&h[0], &f2g0, &h[0]);
    fp2_sub(&h[1], &f2g1, &f1g2);
    return;
  }

  // General case
  fp2_t g_reversed[leng], Q[leng - lenf + 1], Q_reversed[leng - lenf + 1];
  
  for(i = 0; i < leng; i++)
    fp2_copy(&g_reversed[i], &g[leng-1-i]);

  poly_mul_low(Q, leng-lenf+1, f_rev_inv, leng-lenf+1, g_reversed, leng-lenf+1);

  for(i = 0; i < leng - lenf + 1; i++)
    fp2_copy(&Q_reversed[i], &Q[leng - lenf - i]);

  poly_mul_low(g_reversed, lenf-1, Q_reversed, leng-lenf+1, f, lenf);

  for(i = 0; i < lenf-1; i++)
  {
    fp2_mul(&h[i], &g[i], &c);
    fp2_sub(&h[i], &h[i], &g_reversed[i]);
  }
  return;
}


#define POLY_TREE_INDEX(i,j,LENF,h) (i)*(1 << (h))*((LENF)-1) + (1 << (i)) + (j)*((LENF)-1)*(1 << ((h)-(i))) + (j) - 1
#define SCALAR_TREE_INDEX(i,j) (1 << (i)) - 1 + (j)
void reciprocal_tree(fp2_t R[], fp2_t A[], const int leng, const fp2_t H[], const int DEG[],//
		     const int i, const int j, const int n, const int h)
{
  // Given the product tree H and degree tree DEG generated by n LINEAR polynomials, 
  // writes the reverse-reciprocal tree to R and field elements to A such that
  // Rev(H_i)*R_i = A_i mod x^(N). The modulus is N = deg(parent) - deg(self)
  // for inner nodes, and N = leng - deg(root) for the root.
  //
  // The structures of R (resp. A) is the same tree structure as H (resp. DEG); see product_tree()
	//
	// This is a recursive function; the outer call should always be made with i=j=0
  //
  // REQUIRES 1 < n < leng <= 2^h+n+1.
	// REQUIRES R to have enough space for 2^h*(LENF + 1 + (LENF-1)*h) elements.
	// REQUIRES A to have enough space for 2^(i+1)-1 elements.

  if(n == 0)
    return;

  const int root = POLY_TREE_INDEX(i,j,2,h);
  const int parent = POLY_TREE_INDEX(i-1,j>>1,2,h);
  const int brother = POLY_TREE_INDEX(i,j^1,2,h);
  const int root_deg = SCALAR_TREE_INDEX(i,j);
  const int parent_deg = SCALAR_TREE_INDEX(i-1,j>>1);
  const int brother_deg = SCALAR_TREE_INDEX(i,j^1);

  int lenr;

  if(root > 0)
    lenr = DEG[parent_deg] - DEG[root_deg];
  else
    lenr = leng - DEG[root_deg];
  
  // ----------------------------------
  // base cases determined by poly_redc
  if(n == 1)
    return;


  // case for computing  g mod f when len(f), len(g) = 3
  if (DEG[root_deg] == 2 && lenr == 1)
  {
    reciprocal_tree(R, A, lenr-1, H, DEG, i+1, 2*j, n-(n>>1), h);
    reciprocal_tree(R, A, lenr-1, H, DEG, i+1, 2*j+1, n>>1, h);
    return;
  }
  
  // ----------------------------------

  int l;
  
  // When the parent's inverse was calculated to a smaller modulus, need to invert from scratch
  if(i == 0 || leng < lenr)
  {
    for(l = 0; l < lenr && l < DEG[root_deg]+1; l++)
      fp2_copy(&R[root+l], &H[root+DEG[root_deg]-l]);
    for(l = DEG[root_deg]+1; l < lenr; l++){
      fp2_set(&R[root+l], 0);
    }
    reciprocal(&R[root], &(A[root_deg]), &R[root], lenr, lenr);
  }
  else
  {
  // When parent's inverse was to a greater/equal modulus, this inverse can be obtained from it
    for(l = 0; l < lenr; l++)
      fp2_copy(&R[root+l], &H[brother+DEG[brother_deg]-l]);
    poly_mul_low(&R[root], lenr, &R[parent], leng, &R[root], lenr);
    fp2_copy(&A[root_deg], &A[parent_deg]);
  }

  // Now move on to the children
  reciprocal_tree(R, A, lenr, H, DEG, i+1, 2*j, n-(n>>1), h);
  reciprocal_tree(R, A, lenr, H, DEG, i+1, 2*j+1, n>>1, h);
  return;
}


void multieval_unscaled(fp2_t REM[], const poly g, const int leng, const fp2_t R[], const fp2_t A[],//
		const fp2_t H[], const int DEG[], const int i, const int j, const int n, const int h)
{
  // Given the product tree H, degree tree DEG, reciprocal tree R, and constant
  // tree A generated by LINEAR polynomials f_0, ... , f_{n-1}, writes the constant term 
  // of c_i*g mod f_i to REM[i]. The constants c_i are unspecified, but are a function
  // only of leng and f_0,...,f_{n-1} so they cancel out when taking the ratios of
  // remainders of different g's of the same length.
	//
	// This is a recursive function; the outer call should always be made with i=j=0
  //
  // REQUIRES REM to have space for n terms and n > 1

  if(n == 0)
    return;

  const int root = POLY_TREE_INDEX(i,j,2,h);
  const int root_deg = SCALAR_TREE_INDEX(i,j);
  
  fp2_t g_mod[DEG[root_deg]];
  poly_redc(g_mod, g, leng, (const poly)&H[root], DEG[root_deg]+1, (const poly)&R[root], A[root_deg]);

  if(n == 1)
  {
    fp2_copy(&REM[0], &g_mod[0]);
    return;
  }
  
  multieval_unscaled(REM, g_mod, DEG[root_deg], R, A, H, DEG, i+1, 2*j, n-(n>>1), h);
  multieval_unscaled(&(REM[n-(n>>1)]), g_mod, DEG[root_deg], R, A, H, DEG, i+1, 2*j+1, n>>1, h);
  return;
}


void multieval_scaled(fp2_t REM[], const poly G, const fp2_t H[], //
			   const int DEG[], const int i, const int j, const int n, const int h)
{
  // Given the product tree H and degree tree DEG generated by LINEAR f_0,...,f_{n-1},
  // writes the constant term of c_i * g mod f_i(x) to REM[i]
  // The constants c_i are unspecified but are only a function of leng and f_0,...,f_{n-1},
  // so they cancel out when taking the ratios of remainders of different g's of the same length.
	//
	// This is a recursive function; the outer call should always be made with i=j=0
  //
  // REQUIRES REM to have space for n terms and n > 1
  // Also REQUIRES a precomputed G := rev((rev(g mod F)) * F_rev_inv mod x^deg(F)-1) where F is
  // the root of the product tree and F_rev_inv is its reverse's reciprocal mod x^deg(F)

  const int brother = POLY_TREE_INDEX(i,j^1,2,h);
  const int root_deg = SCALAR_TREE_INDEX(i,j);
  const int brother_deg = SCALAR_TREE_INDEX(i,j^1);
  const int uncle_deg = SCALAR_TREE_INDEX(i-1,(j>>1)^1);
  fp2_t fg[DEG[brother_deg]+1];

  if(i == 0)
  {
    if(n == 1)
    {
      fp2_copy(&REM[0], &G[DEG[root_deg]-1]);
      return;
    }
    else
    {
      multieval_scaled(REM, G, H, DEG, i+1, 2*j, n-(n>>1), h);
      multieval_scaled(&(REM[n-(n>>1)]), G, H, DEG, i+1, 2*j+1, n>>1, h);
      return;
    }
  }

  if(i > 1)
    poly_mul_middle(fg, (const poly)&H[brother], DEG[brother_deg]+1, G, DEG[uncle_deg]+1);
  else
    poly_mul_middle(fg, (const poly)&H[brother], DEG[brother_deg]+1, G, DEG[0]);
    
  
  if(n == 1)
  {
    fp2_copy(&REM[0], &fg[DEG[brother_deg]]);
    return;
  }

  multieval_scaled(REM, fg, H, DEG, i+1, 2*j, n-(n>>1), h);
  multieval_scaled(&(REM[n-(n>>1)]), fg, H, DEG, i+1, 2*j+1, n>>1, h);
  return;
}