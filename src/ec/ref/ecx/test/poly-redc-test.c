#include "poly.h"
#include <assert.h>
#include <stdio.h>

bool fp2_isequal(fp2_t a, fp2_t b)
{
  return fp_is_equal(a.re, b.re) && fp_is_equal(a.im, b.im);
}

// VERY NOT SECURE (testing only)
void fp2_random(fp2_t *a)
{
  for (int i = 0; i < NWORDS_FIELD; i++)
  {
    a->re[i] = rand();
    a->im[i] = rand();
  }
  // Normalize
  fp2_t one;
  fp_mont_setone(one.re);
  fp_set(one.im, 0);
  fp2_mul(&*a, &*a, &one);
  // Update seed
  srand((unsigned)a->re[0]);
}

int main()
{
  fp2_t fp2_0, fp2_1;
  fp2_set(&fp2_0, 0);
  fp_mont_setone(fp2_1.re);
  fp_set(fp2_1.im, 0);

// Product tree
#define nmax 16
#define logn 4
#define TREE_SIZE(LENF, h) (1 << (h)) * ((LENF) + 1 + ((LENF)-1) * h)
#define POLY_TREE_INDEX(i, j, LENF, h) (i) * (1 << (h)) * ((LENF)-1) + (1 << (i)) + (j) * ((LENF)-1) * (1 << ((h) - (i))) + (j)-1
#define SCALAR_TREE_INDEX(i, j) (1 << (i)) - 1 + (j)
  int DEG[(1 << (logn + 1)) - 1];                             // Degree tree
  fp2_t H[TREE_SIZE(nmax + 1, logn)],                       // Product tree
      F[nmax * (nmax + 1)],                                 // Product tree input
      R[TREE_SIZE(nmax + 1, logn)], A[(1 << (logn + 1)) - 1]; // Residue tree

  fp2_t h[nmax * nmax + 1], f[nmax * nmax + 1], g[nmax * nmax + 1], f_rev[nmax * nmax + 1],
      f_rev_inv[nmax * nmax + 1], g1[nmax * nmax + 1], g2[nmax * nmax + 1], REM1[nmax],
      REM2[nmax], C[nmax], G1[nmax * nmax + 1], G2[nmax * nmax + 1], G1_rev[nmax * nmax + 1],
      G2_rev[nmax * nmax + 1], R0[nmax * nmax + 1];

  int lenf, leng, n, e, iteration, array_size, tree_size, i, j, root, brother, LENF;
  fp2_t c, ratio, A0;

  // TEST FOR RECIPROCAL
  for (lenf = 1; lenf < nmax; lenf++)
  {
    printf("[%3d%%] Testing reciprocals", 100 * lenf / nmax);
    fflush(stdout);
    printf("\r\x1b[K");

    // Get random poly
    for (e = 0; e < lenf; e++)
      fp2_random(&f[e]);

    for (n = 1; n < nmax; n++)
    {
      // Get the reciprocal and multiply them
      memset(h, 0, sizeof(fp2_t) * n);
      reciprocal(h, &c, f, lenf, n);
      poly_mul_low(h, n, f, lenf, h, n);

      // Compare with expected
      assert(fp2_isequal(h[0], c));
      for (e = 1; e < n; e++)
        assert(fp2_is_zero(&h[e]));
    }
  }
  printf("[%3d%%] Tested reciprocals:\t\tNo errors!\n", 100 * lenf / nmax);

  // TEST FOR REDUCTION
  for (lenf = 2; lenf < nmax; lenf++)
  {
    printf("[%3d%%] Testing polynomial reduction", 100 * lenf / nmax);
    fflush(stdout);
    printf("\r\x1b[K");

    // Get random poly for the mod
    for (e = 0; e < lenf; e++)
    {
      fp2_random(&f[e]);
      fp2_copy(&f_rev[lenf - 1 - e], &f[e]);
    }

    for (leng = 1; leng < nmax; leng++)
    {
      // Get random poly to reduce
      for (e = 0; e < leng; e++)
      {
        fp2_random(&g[e]);
      }

      // Get reverse-inverse mod x^(leng-lenf+1)
      if (leng >= lenf)
      {
        reciprocal(f_rev_inv, &c, f_rev, lenf, leng - lenf + 1);
      }
      else
      {
        fp_mont_setone(c.re);
        fp_set(c.im, 0);
      }

      // Compute the reduction
      poly_redc(h, g, leng, f, lenf, f_rev_inv, c);

      // Reduce manually
      int leng_red = leng;
      fp2_t scale, f_e;
      while (leng_red >= lenf)
      {
        fp2_copy(&scale, &f[lenf - 1]);
        fp2_inv(&scale);
        fp2_mul(&scale, &scale, &g[leng_red - 1]);
        for (e = 0; e < lenf; e++)
        {
          fp2_mul(&f_e, &f[e], &scale);
          fp2_sub(&g[e + leng_red - lenf], &g[e + leng_red - lenf], &f_e);
        }
        leng_red--;
      }

      // Rescale manual result
      if (leng < lenf)
      {
        fp_mont_setone(scale.re);
        fp_set(scale.im, 0);
      }
      else if (lenf == 2 && leng == 3)
      {
        fp2_sqr(&scale, &f[1]);
        fp2_add(&scale, &scale, &scale);
      }
      else
        fp2_copy(&scale, &c);
      for (e = 0; e < leng_red; e++)
        fp2_mul(&g[e], &g[e], &scale);

      // Comapre results
      for (e = leng_red - 1; e >= 0; e--)
        assert(fp2_isequal(h[e], g[e]));
      for (e = leng_red; e < lenf - 1; e++)
        assert(fp2_is_zero(&h[e]));
    }
  }
  printf("[%3d%%] Tested polynomial reduction:\tNo errors!\n", 100 * lenf / nmax);

  // TEST FOR RECIPROCAL TREES

  for (tree_size = 2; tree_size < nmax; tree_size++)
  {
    printf("[%3d%%] Testing reciprocal tree:\t\tTree size %d out of %d", 100 * tree_size / nmax, tree_size, nmax);
    fflush(stdout);
    printf("\r\x1b[K");

    // Get random polys
    LENF = 2;
    for (i = 0; i < tree_size; i++)
    {
      for (e = 0; e < LENF; e++)
      {
        fp2_random(&F[LENF * i + e]);
      }
    }

    // Get product tree then reciprocal tree
    product_tree(H, DEG, 0, 0, F, LENF, tree_size, logn);
    leng = DEG[0] + 1 + (rand() % (1 << logn));
    reciprocal_tree(R, A, leng, H, DEG, 0, 0, tree_size, logn);

    // Check the root
    root = 0;
    lenf = leng - DEG[root];
    for (e = 0; e < DEG[root] + 1 && e < lenf; e++)
    {
      fp2_copy(&f[e], &H[root + DEG[root] - e]);
    }
    for (e = DEG[root] + 1; e < lenf; e++)
    {
      fp2_set(&f[e], 0);
    }
    poly_mul_low(f, lenf, f, lenf, &R[root], lenf);
    assert(fp2_isequal(f[0], A[root]));
    for (e = 1; e < lenf; e++)
    {
      assert(fp2_is_zero(&f[e]));
    }

    // Perform random walks
    for (iteration = 0; iteration < nmax - tree_size; iteration++)
    {
      i = 0;
      j = 0;
      n = tree_size;
      while (n > 1)
      {
        i += 1;
        if (rand() & 1)
        {
          j = 2*j;
          n = n - (n >> 1);
        }
        else
        {
          j = 2*j + 1;
          n = n >> 1;
        }

        int root = POLY_TREE_INDEX(i,j,2,logn);
        int root_deg = SCALAR_TREE_INDEX(i,j);
        int brother_deg = SCALAR_TREE_INDEX(i,j^1);

        // Check current node
        lenf = DEG[brother_deg];
        for (e = 0; e < DEG[root_deg] + 1 && e < lenf; e++)
        {
          fp2_copy(&f[e], &H[root + DEG[root_deg] - e]);
        }
        for (e = DEG[root_deg] + 1; e < lenf; e++)
        {
          fp2_set(&f[e], 0);
        }
        poly_mul_low(f, lenf, f, lenf, &R[root], lenf);
        assert(fp2_isequal(f[0], A[root_deg]));
        for (e = 1; e < lenf; e++)
        {
          assert(fp2_is_zero(&f[e]));
        }
      }
    }
  }
  printf("[%3d%%] Tested reciprocal tree:\t\tNo errors!\n", 100 * tree_size / nmax);

  // TEST FOR REMAINDERS
  for (tree_size = 2; tree_size < nmax; tree_size++)
  {
    printf("[%3d%%] Testing batched remainders:\t\tTree size %d out of %d", 100 * tree_size / nmax, tree_size, nmax);
    fflush(stdout);
    printf("\r\x1b[K");

    // Get random polys
    LENF = 2;
    for (i = 0; i < tree_size; i++)
    {
      for (e = 0; e < LENF; e++)
        fp2_random(&F[LENF * i + e]);
    }

    // Get product tree, reciprocal tree, and remainders
    product_tree(H, DEG, 0, 0, F, LENF, tree_size, logn);
    leng = DEG[0] + 1 + (rand() % (1 << logn));
    for (e = 0; e < leng; e++)
    {
      fp2_random(&g1[e]);
      fp2_random(&g2[e]);
    }
    reciprocal_tree(R, A, leng, H, DEG, 0, 0, tree_size, logn);
    multieval_unscaled(REM1, g1, leng, R, A, H, DEG, 0, 0, tree_size, logn);
    multieval_unscaled(REM2, g2, leng, R, A, H, DEG, 0, 0, tree_size, logn);

    for (i = 0; i < tree_size; i++)
    {
      // Get ratio of the remainder
      fp2_inv(&REM1[i]);
      fp2_mul(&ratio, &REM1[i], &REM2[i]);

      // Compute remainders manually
      for (e = 0; e < LENF; e++)
        fp2_copy(&f_rev[e], &F[LENF * i + LENF - 1 - e]);
      reciprocal(f_rev_inv, &c, f_rev, LENF, leng - LENF + 1);
      poly_redc(h, g1, leng, &F[LENF * i], LENF, f_rev_inv, c);
      fp2_copy(&REM1[i], &h[0]);
      poly_redc(h, g2, leng, &F[LENF * i], LENF, f_rev_inv, c);
      fp2_copy(&REM2[i], &h[0]);

      // Compare results
      fp2_inv(&REM1[i]);
      fp2_mul(&REM1[i], &REM1[i], &REM2[i]);
      assert(fp2_is_equal(&(REM1[i]), &ratio));
    }
  }
  printf("[%3d%%] Tested batched remainders:\tNo errors!\n", 100 * tree_size / nmax);

  // TEST FOR SCALED REMAINDER TREE
  for (tree_size = 1; tree_size < nmax; tree_size++)
  {
    printf("[%3d%%] Testing scaled remainder tree:\tTree size %d out of %d", 100 * tree_size / nmax, tree_size, nmax);
    fflush(stdout);
    printf("\r\x1b[K");

    // Get random polys
    LENF = 2;
    for (i = 0; i < tree_size; i++)
    {
      for (e = 0; e < LENF; e++)
        fp2_random(&F[LENF * i + e]);
    }

    // Get random polys to reduce
    product_tree(H, DEG, 0, 0, F, LENF, tree_size, logn);
    leng = DEG[0] + 1 + (rand() % nmax);
    for (e = 0; e < leng; e++)
    {
      fp2_random(&g1[e]);
      fp2_random(&g2[e]);
    }

    // Get the required initial nodes
    for (e = 0; e < DEG[0] + 1; e++)
      fp2_copy(&f_rev[e], &H[DEG[0] - e]);
    if (DEG[0] > leng - DEG[0])
      reciprocal(R0, &A0, f_rev, DEG[0] + 1, DEG[0]);
    else
      reciprocal(R0, &A0, f_rev, DEG[0] + 1, leng - DEG[0]);
    poly_redc(G1, g1, leng, H, DEG[0] + 1, R0, A0);
    poly_redc(G2, g2, leng, H, DEG[0] + 1, R0, A0);
    for (e = 0; e < DEG[0]; e++)
    {
      fp2_copy(&G1_rev[e], &G1[DEG[0] - 1 - e]);
      fp2_copy(&G2_rev[e], &G2[DEG[0] - 1 - e]);
    }
    poly_mul_middle(G1_rev, G1_rev, DEG[0], R0, DEG[0]);
    poly_mul_middle(G2_rev, G2_rev, DEG[0], R0, DEG[0]);
    for (e = 0; e < DEG[0]; e++)
    {
      fp2_copy(&G1[e], &G1_rev[DEG[0] - 1 - e]);
      fp2_copy(&G2[e], &G2_rev[DEG[0] - 1 - e]);
    }

    // Compute the scaled remainder trees
    multieval_scaled(REM1, G1, H, DEG, 0, 0, tree_size, logn);
    multieval_scaled(REM2, G2, H, DEG, 0, 0, tree_size, logn);

    for (i = 0; i < tree_size; i++)
    {
      // Get ratio of the remainder
      fp2_inv(&REM1[i]);
      fp2_mul(&ratio, &REM1[i], &REM2[i]);

      // Compute remainders manually
      for (e = 0; e < LENF; e++)
        fp2_copy(&f_rev[e], &F[LENF * i + LENF - 1 - e]);
      reciprocal(f_rev_inv, &c, f_rev, LENF, leng - LENF + 1);
      poly_redc(h, g1, leng, &F[LENF * i], LENF, f_rev_inv, c);
      fp2_copy(&REM1[i], &h[0]);
      poly_redc(h, g2, leng, &F[LENF * i], LENF, f_rev_inv, c);
      fp2_copy(&REM2[i], &h[0]);

      // Compare results
      fp2_inv(&REM1[i]);
      fp2_mul(&REM1[i], &REM1[i], &REM2[i]);
      assert(fp2_is_equal(&(REM1[i]), &ratio));
    }
  }
  printf("[%3d%%] Tested scaled remainder tree:\tNo errors!\n", 100 * tree_size / nmax);

  printf("-- All tests passed.\n");
}
