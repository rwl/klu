/**
 * KLU: a sparse LU factorization algorithm.
 * Copyright (C) 2004-2009, Timothy A. Davis.
 * Copyright (C) 2011-2012, Richard W. Lincoln.
 * http://www.cise.ufl.edu/research/sparse/klu
 *
 * -------------------------------------------------------------------------
 *
 * KLU is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * KLU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */
part of edu.ufl.cise.klu.tdouble;

/**
 * Allocate Symbolic object, and check input matrix. Not user callable.
 */
KLU_symbolic _alloc_symbolic(final int n, final Int32List Ap,
                             final Int32List Ai, final KLU_common Common) {
  if (Common == null) {
    return (null);
  }
  Common.status = KLU_OK;

  /* A is n-by-n, with n > 0.  Ap [0] = 0 and nz = Ap [n] >= 0 required.
   * Ap [j] <= Ap [j+1] must hold for all j = 0 to n-1.  Row indices in Ai
   * must be in the range 0 to n-1, and no duplicate entries can be present.
   * The list of row indices in each column of A need not be sorted.
   */

  if (n <= 0 || Ap == null || Ai == null) {
    /* Ap and Ai must be present, and n must be > 0 */
    Common.status = KLU_INVALID;
    return (null);
  }

  final nz = Ap[n];
  if (Ap[0] != 0 || nz < 0) {
    /* nz must be >= 0 and Ap [0] must equal zero */
    Common.status = KLU_INVALID;
    return (null);
  }

  for (int j = 0; j < n; j++) {
    if (Ap[j] > Ap[j + 1]) {
      /* column pointers must be non-decreasing */
      Common.status = KLU_INVALID;
      return (null);
    }
  }
  final P = malloc_int(n, Common);
  if (Common.status < KLU_OK) {
    /* out of memory */
    Common.status = KLU_OUT_OF_MEMORY;
    return (null);
  }
  for (int i = 0; i < n; i++) {
    P[i] = EMPTY;
  }
  for (int j = 0; j < n; j++) {
    final pend = Ap[j + 1];
    for (int p = Ap[j]; p < pend; p++) {
      final i = Ai[p];
      if (i < 0 || i >= n || P[i] == j) {
        /* row index out of range, or duplicate entry */
        //klu_free (P, n, Integer.class, Common) ;
        Common.status = KLU_INVALID;
        return (null);
      }
      /* flag row i as appearing in column j */
      P[i] = j;
    }
  }

  /* ---------------------------------------------------------------------- */
  /* allocate the Symbolic object */
  /* ---------------------------------------------------------------------- */

  KLU_symbolic Symbolic;
  try {
    //Symbolic = klu_malloc (KLU_symbolic.class, 1, Common) ;
    Symbolic = new KLU_symbolic();
  } on OutOfMemoryError catch (e) {
    /* out of memory */
    //klu_free (P, n, sizeof (int), Common) ;
    Common.status = KLU_OUT_OF_MEMORY;
    return (null);
  }

  final Q = malloc_int(n, Common);
  final R = malloc_int(n + 1, Common);
  final Lnz = malloc_dbl(n, Common);

  Symbolic.n = n;
  Symbolic.nz = nz;
  Symbolic.P = P;
  Symbolic.Q = Q;
  Symbolic.R = R;
  Symbolic.Lnz = Lnz;

  if (Common.status < KLU_OK) {
    /* out of memory */
    //klu_free_symbolic (Symbolic, Common) ;
    Symbolic = null;
    Common.status = KLU_OUT_OF_MEMORY;
    return (null);
  }

  return (Symbolic);
}

/**
 * Analyzes a matrix using given P and Q.
 *
 * Given an input permutation P and Q, create the Symbolic object. BTF can
 * be done to modify the user's P and Q (does not perform the max transversal;
 * just finds the strongly-connected components).
 *
 * Order the matrix with BTF (or not), then use natural or given ordering
 * P and Q on the blocks. P and Q are interpreted as identity if `null`.
 *
 * User's row permutation [Puser] is size [n] and may be `null`. User's
 * column permutation [Quser] is size [n] and may be `null`.
 */
KLU_symbolic analyze_given(final int n, final Int32List Ap,
                           final Int32List Ai, final Int32List Puser,
                           final Int32List Quser, final KLU_common Common) {

  /* ---------------------------------------------------------------------- */
  /* determine if input matrix is valid, and get # of nonzeros */
  /* ---------------------------------------------------------------------- */

  final Symbolic = _alloc_symbolic(n, Ap, Ai, Common);
  if (Symbolic == null) {
    return (null);
  }
  final P = Symbolic.P;
  final Q = Symbolic.Q;
  final R = Symbolic.R;
  final Lnz = Symbolic.Lnz;
  final nz = Symbolic.nz;

  /* ---------------------------------------------------------------------- */
  /* Q = Quser, or identity if Quser is null */
  /* ---------------------------------------------------------------------- */

  if (Quser == null) {
    for (int k = 0; k < n; k++) {
      Q[k] = k;
    }
  } else {
    for (int k = 0; k < n; k++) {
      Q[k] = Quser[k];
    }
  }

  /* ---------------------------------------------------------------------- */
  /* get the control parameters for BTF and ordering method */
  /* ---------------------------------------------------------------------- */

  int do_btf = Common.btf;
  do_btf = (do_btf != 0) ? TRUE : FALSE;
  Symbolic.ordering = 2;
  Symbolic.do_btf = do_btf;

  /* ---------------------------------------------------------------------- */
  /* find the block triangular form, if requested */
  /* ---------------------------------------------------------------------- */

  int maxblock, nzoff, nblocks;
  if (do_btf != 0) {

    /* ------------------------------------------------------------------ */
    /* get workspace for BTF_strongcomp */
    /* ------------------------------------------------------------------ */

    Int32List Pinv, Work, Bi;
    int k1, k2, nk, oldcol;

    Work = malloc_int(4 * n, Common);
    Pinv = malloc_int(n, Common);
    if (Puser != null) {
      Bi = malloc_int(nz + 1, Common);
    } else {
      Bi = Ai;
    }

    if (Common.status < KLU_OK) {
      /* out of memory */
      //klu_free (Work, 4*n, sizeof (int), Common) ;
      Work = null;
      //klu_free (Pinv, n, sizeof (int), Common) ;
      Pinv = null;
      if (Puser != null) {
        //klu_free (Bi, nz+1, sizeof (int), Common) ;
        Bi = null;
      }
      //klu_free_symbolic (Symbolic, Common) ;
      Common.status = KLU_OUT_OF_MEMORY;
      return (null);
    }

    /* ------------------------------------------------------------------ */
    /* B = Puser * A */
    /* ------------------------------------------------------------------ */

    if (Puser != null) {
      for (int k = 0; k < n; k++) {
        Pinv[Puser[k]] = k;
      }
      for (int p = 0; p < nz; p++) {
        Bi[p] = Pinv[Ai[p]];
      }
    }

    /* ------------------------------------------------------------------ */
    /* find the strongly-connected components */
    /* ------------------------------------------------------------------ */

    /* modifies Q, and determines P and R */
    nblocks = btf.strongcomp(n, Ap, Bi, Q, P, R);

    /* ------------------------------------------------------------------ */
    /* P = P * Puser */
    /* ------------------------------------------------------------------ */

    if (Puser != null) {
      for (int k = 0; k < n; k++) {
        Work[k] = Puser[P[k]];
      }
      for (int k = 0; k < n; k++) {
        P[k] = Work[k];
      }
    }

    /* ------------------------------------------------------------------ */
    /* Pinv = inverse of P */
    /* ------------------------------------------------------------------ */

    for (int k = 0; k < n; k++) {
      Pinv[P[k]] = k;
    }

    /* ------------------------------------------------------------------ */
    /* analyze each block */
    /* ------------------------------------------------------------------ */

    nzoff = 0;
    /* nz in off-diagonal part */
    maxblock = 1;
    /* size of the largest block */

    for (int block = 0; block < nblocks; block++) {

      /* -------------------------------------------------------------- */
      /* the block is from rows/columns k1 to k2-1 */
      /* -------------------------------------------------------------- */

      k1 = R[block];
      k2 = R[block + 1];
      nk = k2 - k1;
      PRINTF("BLOCK $block, k1 $k1 k2-1 $k2-1 nk $nk\n");
      maxblock = MAX(maxblock, nk);

      /* -------------------------------------------------------------- */
      /* scan the kth block, C */
      /* -------------------------------------------------------------- */

      for (int k = k1; k < k2; k++) {
        oldcol = Q[k];
        final pend = Ap[oldcol + 1];
        for (int p = Ap[oldcol]; p < pend; p++) {
          if (Pinv[Ai[p]] < k1) {
            nzoff++;
          }
        }
      }

      /* fill-in not estimated */
      Lnz[block] = EMPTY_D;
    }

    /* ------------------------------------------------------------------ */
    /* free all workspace */
    /* ------------------------------------------------------------------ */

    //klu_free (Work, 4*n, sizeof (int), Common) ;
    Work = null;
    //klu_free (Pinv, n, sizeof (int), Common) ;
    Pinv = null;
    if (Puser != null) {
      //klu_free (Bi, nz+1, sizeof (int), Common) ;
      Bi = null;
    }

  } else {

    /* ------------------------------------------------------------------ */
    /* BTF not requested */
    /* ------------------------------------------------------------------ */

    nzoff = 0;
    nblocks = 1;
    maxblock = n;
    R[0] = 0;
    R[1] = n;
    Lnz[0] = EMPTY_D;

    /* ------------------------------------------------------------------ */
    /* P = Puser, or identity if Puser is null */
    /* ------------------------------------------------------------------ */

    for (int k = 0; k < n; k++) {
      P[k] = (Puser == null) ? k : Puser[k];
    }
  }

  /* ---------------------------------------------------------------------- */
  /* return the symbolic object */
  /* ---------------------------------------------------------------------- */

  Symbolic.nblocks = nblocks;
  Symbolic.maxblock = maxblock;
  Symbolic.lnz = EMPTY_D;
  Symbolic.unz = EMPTY_D;
  Symbolic.nzoff = nzoff;

  return (Symbolic);
}
