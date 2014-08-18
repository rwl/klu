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

void factor2(final Int32List Ap, final Int32List Ai, final Float64List Ax,
             final KLU_symbolic Symbolic, final KLU_numeric Numeric,
             final KLU_common Common) {
  Int32List lnz_block = new Int32List(1);
  Int32List unz_block = new Int32List(1);

  /* ---------------------------------------------------------------------- */
  /* initializations */
  /* ---------------------------------------------------------------------- */

  /* get the contents of the Symbolic object */
  final n = Symbolic.n;
  final P = Symbolic.P;
  final Q = Symbolic.Q;
  final R = Symbolic.R;
  final Lnz = Symbolic.Lnz;
  final nblocks = Symbolic.nblocks;
  final nzoff = Symbolic.nzoff;

  final Pnum = Numeric.Pnum;
  final Offp = Numeric.Offp;
  final Offi = Numeric.Offi;
  final Offx = Numeric.Offx;

  final Lip = Numeric.Lip;
  final Uip = Numeric.Uip;
  final Llen = Numeric.Llen;
  final Ulen = Numeric.Ulen;
  final LUbx = Numeric.LUbx;
  final Udiag = Numeric.Udiag;

  final Rs = Numeric.Rs;
  final Pinv = Numeric.Pinv;
  final X = Numeric.Xwork; // X is of size n
  final Iwork = Numeric.Iwork; // 5*maxblock for KLU_factor
  //Pblock = Iwork + 5*((int) Symbolic.maxblock) ;
  final Pblock = new Int32List(Symbolic.maxblock);// 1*maxblock for Pblock
  Common.nrealloc = 0;
  final scale = Common.scale;
  int max_lnz_block = 1;
  int max_unz_block = 1;

  /* compute the inverse of P from symbolic analysis. Will be updated to
   * become the inverse of the numerical factorization when the factorization
   * is done, for use in KLU_refactor */
  if (!NDEBUG) {
    for (int k = 0; k < n; k++) {
      Pinv[k] = EMPTY;
    }
  }
  for (int k = 0; k < n; k++) {
    ASSERT(P[k] >= 0 && P[k] < n);
    Pinv[P[k]] = k;
  }
  if (!NDEBUG) {
    for (int k = 0; k < n; k++) ASSERT(Pinv[k] != EMPTY);
  }

  int lnz = 0;
  int unz = 0;
  Common.noffdiag = 0;
  Offp[0] = 0;

  /* ---------------------------------------------------------------------- */
  /* optionally check input matrix and compute scale factors */
  /* ---------------------------------------------------------------------- */

  if (scale >= 0) {
    /* use Pnum as workspace. NOTE: scale factors are not yet permuted
     * according to the final pivot row ordering, so Rs [oldrow] is the
     * scale factor for A (oldrow,:), for the user's matrix A.  Pnum is
     * used as workspace in KLU_scale.  When the factorization is done,
     * the scale factors are permuted according to the final pivot row
     * permutation, so that Rs [k] is the scale factor for the kth row of
     * A(p,q) where p and q are the final row and column permutations. */
    klu_scale(scale, n, Ap, Ai, Ax, Rs, Pnum, Common);
    if (Common.status < KLU_OK) {
      /* matrix is invalid */
      return;
    }
  }

  if (!NDEBUG) {
    if (scale > 0) {
      for (int k = 0; k < n; k++) PRINTF("Rs [$k] ${Rs [k]}\n");
    }
  }

  /* ---------------------------------------------------------------------- */
  /* factor each block using klu */
  /* ---------------------------------------------------------------------- */

  for (int block = 0; block < nblocks; block++) {

    /* ------------------------------------------------------------------ */
    /* the block is from rows/columns k1 to k2-1 */
    /* ------------------------------------------------------------------ */

    final k1 = R[block];
    final k2 = R[block + 1];
    final nk = k2 - k1;
    PRINTF("FACTOR BLOCK $block, k1 $k1 k2-1 ${k2-1} nk $nk\n");

    if (nk == 1) {

      /* -------------------------------------------------------------- */
      /* singleton case */
      /* -------------------------------------------------------------- */

      int poff = Offp[k1];
      final oldcol = Q[k1];
      final pend = Ap[oldcol + 1];
      //CLEAR (s) ;
      double s = 0.0;

      if (scale <= 0) {
        /* no scaling */
        for (int p = Ap[oldcol]; p < pend; p++) {
          final oldrow = Ai[p];
          final newrow = Pinv[oldrow];
          if (newrow < k1) {
            Offi[poff] = oldrow;
            Offx[poff] = Ax[p];
            poff++;
          } else {
            ASSERT(newrow == k1);
            PRINTF("singleton block $block");
            PRINT_ENTRY(Ax[p]);
            s = Ax[p];
          }
        }
      } else {
        /* row scaling.  NOTE: scale factors are not yet permuted
         * according to the pivot row permutation, so Rs [oldrow] is
         * used below.  When the factorization is done, the scale
         * factors are permuted, so that Rs [newrow] will be used in
         * klu_solve, klu_tsolve, and klu_rgrowth */
        for (int p = Ap[oldcol]; p < pend; p++) {
          final oldrow = Ai[p];
          final newrow = Pinv[oldrow];
          if (newrow < k1) {
            Offi[poff] = oldrow;
            //SCALE_DIV_ASSIGN (Offx [poff], Ax [p], Rs [oldrow]) ;
            Offx[poff] = Ax[p] / Rs[oldrow];
            poff++;
          } else {
            ASSERT(newrow == k1);
            PRINTF("singleton block $block ");
            PRINT_ENTRY(Ax[p]);
            s = Ax[p] / Rs[oldrow];
            //SCALE_DIV_ASSIGN (s, Ax [p], Rs [oldrow]) ;
          }
        }
      }

      Udiag[k1] = s;

      if (IS_ZERO(s)) {
        /* singular singleton */
        Common.status = KLU_SINGULAR;
        Common.numerical_rank = k1;
        Common.singular_col = oldcol;
        if (Common.halt_if_singular == 1) {
          return;
        }
      }

      Offp[k1 + 1] = poff;
      Pnum[k1] = P[k1];
      lnz++;
      unz++;

    } else {

      /* -------------------------------------------------------------- */
      /* construct and factorize the kth block */
      /* -------------------------------------------------------------- */

      double lsize;
      if (Lnz[block] < 0) {
        /* COLAMD was used - no estimate of fill-in */
        /* use 10 times the nnz in A, plus n */
        lsize = -(Common.initmem);
      } else {
        lsize = Common.initmem_amd * Lnz[block] + nk;
      }

      /* allocates 1 arrays: LUbx [block] */
      Numeric.LUsize[block] = kernel_factor(nk, Ap, Ai, Ax, Q, lsize, LUbx, block, Udiag, k1, Llen, k1, Ulen, k1, Lip, k1, Uip, k1, Pblock, lnz_block, unz_block, X, Iwork, k1, Pinv, Rs, Offp, Offi, Offx, Common);

      if (Common.status < KLU_OK || (Common.status == KLU_SINGULAR && Common.halt_if_singular == 1)) {
        /* out of memory, invalid inputs, or singular */
        return;
      }

      PRINTF("\n----------------------- L $block:\n");
      if (!NDEBUG) ASSERT(_valid_LU(nk, TRUE, Lip, k1, Llen, k1, LUbx[block]));
      PRINTF("\n----------------------- U $block:\n");
      if (!NDEBUG) ASSERT(_valid_LU(nk, FALSE, Uip, k1, Ulen, k1, LUbx[block]));

      /* -------------------------------------------------------------- */
      /* get statistics */
      /* -------------------------------------------------------------- */

      lnz += lnz_block[0];
      unz += unz_block[0];
      max_lnz_block = MAX(max_lnz_block, lnz_block[0]);
      max_unz_block = MAX(max_unz_block, unz_block[0]);

      if (Lnz[block] == EMPTY) {
        /* revise estimate for subsequent factorization */
        Lnz[block] = MAX(lnz_block[0], unz_block[0]);
      }

      /* -------------------------------------------------------------- */
      /* combine the klu row ordering with the symbolic pre-ordering */
      /* -------------------------------------------------------------- */

      PRINTF("Pnum, 1-based:\n");
      for (int k = 0; k < nk; k++) {
        ASSERT(k + k1 < n);
        ASSERT(Pblock[k] + k1 < n);
        Pnum[k + k1] = P[Pblock[k] + k1];
        PRINTF("Pnum ($k + $k1 + 1 = ${k+k1+1}) = ${Pnum [k+k1]} + 1 = ${Pnum [k+k1]+1}\n");
      }

      /* the local pivot row permutation Pblock is no longer needed */
    }
  }
  ASSERT(nzoff == Offp[n]);
  PRINTF("\n------------------- Off diagonal entries:\n");
  if (!NDEBUG) ASSERT_INT(_valid(n, Offp, Offi, Offx));

  Numeric.lnz = lnz;
  Numeric.unz = unz;
  Numeric.max_lnz_block = max_lnz_block;
  Numeric.max_unz_block = max_unz_block;

  /* compute the inverse of Pnum */
  if (!NDEBUG) {
    for (int k = 0; k < n; k++) {
      Pinv[k] = EMPTY;
    }
  }
  for (int k = 0; k < n; k++) {
    ASSERT(Pnum[k] >= 0 && Pnum[k] < n);
    Pinv[Pnum[k]] = k;
  }
  if (!NDEBUG) {
    for (int k = 0; k < n; k++) ASSERT(Pinv[k] != EMPTY);
  }

  /* permute scale factors Rs according to pivotal row order */
  if (scale > 0) {
    for (int k = 0; k < n; k++) {
      //REAL (X [k]) = Rs [Pnum [k]] ;
      X[k] = Rs[Pnum[k]];
    }
    for (int k = 0; k < n; k++) {
      Rs[k] = X[k]; //REAL (X [k]) ;
    }
  }

  PRINTF("\n------------------- Off diagonal entries, old:\n");
  if (!NDEBUG) ASSERT_INT(_valid(n, Offp, Offi, Offx));

  /* apply the pivot row permutations to the off-diagonal entries */
  for (int p = 0; p < nzoff; p++) {
    ASSERT(Offi[p] >= 0 && Offi[p] < n);
    Offi[p] = Pinv[Offi[p]];
  }

  PRINTF("\n------------------- Off diagonal entries, new:\n");
  if (!NDEBUG) ASSERT_INT(_valid(n, Offp, Offi, Offx));

  if (!NDEBUG) {
    PRINTF("\n ############# KLU_BTF_FACTOR done, nblocks $nblocks\n");
    double ss;
    final Udiag = Numeric.Udiag;
    for (int block = 0; block < nblocks && Common.status == KLU_OK; block++) {
      final k1 = R[block];
      final k2 = R[block + 1];
      final nk = k2 - k1;
      PRINTF("\n======================KLU_factor output: k1 $k1 k2 $k2 nk $nk\n");
      if (nk == 1) {
        PRINTF("singleton  ");
        /* ENTRY_PRINT (singleton [block]) ; */
        ss = Udiag[k1];
        PRINT_ENTRY(ss);
      } else {
        final Lip = Numeric.Lip;
        final Lip_offset = k1;
        final Llen = Numeric.Llen;
        final Llen_offset = k1;
        final LU = Numeric.LUbx[block];// as Float64List ;
        PRINTF("\n---- L block $block\n");
        if (!NDEBUG) ASSERT(_valid_LU(nk, TRUE, Lip, Lip_offset, Llen, Llen_offset, LU));
        final Uip = Numeric.Uip;
        final Uip_offset = k1;
        final Ulen = Numeric.Ulen;
        final Ulen_offset = k1;
        PRINTF("\n---- U block $block\n");
        if (!NDEBUG) ASSERT(_valid_LU(nk, FALSE, Uip, Uip_offset, Ulen, Ulen_offset, LU));
      }
    }
  }
}

/**
 * Factor the matrix, after ordering and analyzing it with KLU_analyze
 * or KLU_analyze_given.
 *
 * Returns `null` on error or a valid [KLU_numeric] object if successful.
 */
KLU_numeric factor(final Int32List Ap, final Int32List Ai,
                   final Float64List Ax, final KLU_symbolic Symbolic,
                   final KLU_common Common) {
  if (Common == null) {
    return (null);
  }
  Common.status = KLU_OK;
  Common.numerical_rank = EMPTY;
  Common.singular_col = EMPTY;

  /* ---------------------------------------------------------------------- */
  /* get the contents of the Symbolic object */
  /* ---------------------------------------------------------------------- */

  /* check for a valid Symbolic object */
  if (Symbolic == null) {
    Common.status = KLU_INVALID;
    return (null);
  }

  final n = Symbolic.n;
  final nzoff = Symbolic.nzoff;
  final nblocks = Symbolic.nblocks;
  final maxblock = Symbolic.maxblock;
  PRINTF("KLU_factor:  n $n nzoff $nzoff nblocks $nblocks maxblock $maxblock\n");

  /* ---------------------------------------------------------------------- */
  /* get control parameters and make sure they are in the proper range */
  /* ---------------------------------------------------------------------- */

  Common.initmem_amd = MAX(1.0, Common.initmem_amd);
  Common.initmem = MAX(1.0, Common.initmem);
  Common.tol = MIN(Common.tol, 1.0);
  Common.tol = MAX(0.0, Common.tol);
  Common.memgrow = MAX(1.0, Common.memgrow);

  /* ---------------------------------------------------------------------- */
  /* allocate the Numeric object  */
  /* ---------------------------------------------------------------------- */

  /* this will not cause int overflow (already checked by KLU_symbolic) */
  final n1 = n + 1;
  final nzoff1 = nzoff + 1;

  //Numeric = klu_malloc (sizeof (KLU_numeric), 1, Common) ;
  KLU_numeric Numeric;
  try {
    Numeric = new KLU_numeric();
  } on OutOfMemoryError catch (e) {
    Common.status = KLU_OUT_OF_MEMORY;
    return (null);
  }
  Numeric.n = n;
  Numeric.nblocks = nblocks;
  Numeric.nzoff = nzoff;
  Numeric.Pnum = malloc_int(n, Common);
  Numeric.Offp = malloc_int(n1, Common);
  Numeric.Offi = malloc_int(nzoff1, Common);
  Numeric.Offx = malloc_dbl(nzoff1, Common);

  Numeric.Lip = malloc_int(n, Common);
  Numeric.Uip = malloc_int(n, Common);
  Numeric.Llen = malloc_int(n, Common);
  Numeric.Ulen = malloc_int(n, Common);

  Numeric.LUsize = malloc_int(nblocks, Common);

  //Numeric.LUbx = klu_malloc (nblocks, sizeof (Float64List), Common) ;
  Numeric.LUbx = new List<Float64List>(nblocks);
  if (Numeric.LUbx != null) {
    for (int k = 0; k < nblocks; k++) {
      Numeric.LUbx[k] = null;
    }
  }

  Numeric.Udiag = malloc_dbl(n, Common);

  if (Common.scale > 0) {
    Numeric.Rs = malloc_dbl(n, Common);
  } else {
    /* no scaling */
    Numeric.Rs = null;
  }

  Numeric.Pinv = malloc_int(n, Common);

  /* allocate permanent workspace for factorization and solve.  Note that the
   * solver will use an Xwork of size 4n, whereas the factorization codes use
   * an Xwork of size n and integer space (Iwork) of size 6n. KLU_condest
   * uses an Xwork of size 2n.  Total size is:
   *
   *    n*sizeof(double) + max (6*maxblock*sizeof(Int), 3*n*sizeof(double))
   */
  //s = klu_mult_size_t (n, sizeof (double), ok) ;
  final s = n;
  //n3 = klu_mult_size_t (n, 3 * sizeof (double), ok) ;
  final n3 = 3 * n;
  //b6 = klu_mult_size_t (maxblock, 6 * sizeof (Int), ok) ;
  final b6 = 6 * maxblock;
  final ok = new Int32List.fromList([TRUE]);
  Numeric.worksize = add_size_t(s, MAX(n3, b6), ok);
  try {
    if (ok[0] == 0) throw new OutOfMemoryError();

    //Numeric.Work = klu_malloc (Numeric.worksize, 1, Common) ;
    Numeric.Work = new Float64List(Numeric.worksize);
    Numeric.Xwork = Numeric.Work;
    //Numeric.Iwork = (Int[]) ((Float64List) Numeric.Xwork + n) ;
    Numeric.Iwork = new Int32List(b6);
  } on OutOfMemoryError catch (e) {
    /* out of memory or problem too large */
    Common.status = ok[0] == 1 ? KLU_OUT_OF_MEMORY : KLU_TOO_LARGE;
    //klu_free_numeric (Numeric, Common) ;
    Numeric = null;
    return (null);
  }

  /* ---------------------------------------------------------------------- */
  /* factorize the blocks */
  /* ---------------------------------------------------------------------- */

  factor2(Ap, Ai, Ax, Symbolic, Numeric, Common);

  /* ---------------------------------------------------------------------- */
  /* return or free the Numeric object */
  /* ---------------------------------------------------------------------- */

  if (Common.status < KLU_OK) {
    /* out of memory or inputs invalid */
    //klu_free_numeric (Numeric, Common) ;
    Numeric = null;
  } else if (Common.status == KLU_SINGULAR) {
    if (Common.halt_if_singular == 1) {
      /* Matrix is singular, and the Numeric object is only partially
       * defined because we halted early.  This is the default case for
       * a singular matrix. */
      //klu_free_numeric (Numeric, Common) ;
      Numeric = null;

    }
  } else if (Common.status == KLU_OK) {
    /* successful non-singular factorization */
    Common.numerical_rank = n;
    Common.singular_col = n;
  }
  return (Numeric);
}
