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
 * Orders and analyzes a matrix.
 *
 * Order the matrix using BTF (or not), and then AMD, COLAMD, the natural
 * ordering, or the user-provided-function on the blocks.  Does not support
 * using a given ordering (use klu_analyze_given for that case).
 *
 * A is [n]-by-[n]. Column pointers [Ap] is size `n+1`. Row indices [Ai] is
 * size `nz`. [nblocks] is the number of blocks. [Pbtf] is the BTF row
 * permutation. [Qbtf] is the BTF column permutation. R is size `n+1`, but
 * only `Rbtf[0..nblocks]` is used. [ordering] may be 0, 1, or 3 for this
 * routine. [P] is size [n]. [Q] is size [n]. [Lnz] is size [n], but only
 * `Lnz[0..nblocks-1]` is used. [Pblk] is size `maxblock`. [Cp] is size
 * `maxblock+1`. Ci is size `MAX(nz+1, Cilen)`. [Cilen] is size `nz+1`, or
 * `COLAMD_recommend(nz,n,n)` for COLAMD. [Pinv] is size `maxblock`. Returns
 * [KLU_OK] or < 0 if an error occurrs.
 */
int analyze_worker(final int n, final Int32List Ap, final Int32List Ai,
                   final int nblocks, final Int32List Pbtf,
                   final Int32List Qbtf, final Int32List R,
                   final int ordering, final Int32List P, final Int32List Q,
                   final Float64List Lnz, final Int32List Pblk,
                   final Int32List Cp, final Int32List Ci, final int Cilen,
                   final Int32List Pinv, final KLU_symbolic Symbolic,
                   final KLU_common Common) {
  final amd_Info = new List<num>(amd.AMD_INFO);
  double lnz1, flops1;
  int ok,
      err = KLU_INVALID;
  final cstats = new Int32List(colamd.COLAMD_STATS);

  /* ---------------------------------------------------------------------- */
  /* initializations */
  /* ---------------------------------------------------------------------- */

  /* compute the inverse of Pbtf */
  if (!NDEBUG) {
    for (int k = 0; k < n; k++) {
      P[k] = EMPTY;
      Q[k] = EMPTY;
      Pinv[k] = EMPTY;
    }
  }
  for (int k = 0; k < n; k++) {
    ASSERT(Pbtf[k] >= 0 && Pbtf[k] < n);
    Pinv[Pbtf[k]] = k;
  }
  if (!NDEBUG) {
    for (int k = 0; k < n; k++) ASSERT(Pinv[k] != EMPTY);
  }
  int nzoff = 0;
  double lnz = 0.0;
  int maxnz = 0;
  double flops = 0.0;
  /* only computed by AMD */
  Symbolic.symmetry = EMPTY_D;

  /* ---------------------------------------------------------------------- */
  /* order each block */
  /* ---------------------------------------------------------------------- */

  for (int block = 0; block < nblocks; block++) {

    /* ------------------------------------------------------------------ */
    /* the block is from rows/columns k1 to k2-1 */
    /* ------------------------------------------------------------------ */

    final k1 = R[block];
    final k2 = R[block + 1];
    final nk = k2 - k1;
    PRINTF("BLOCK $block, k1 $k1 k2-1 ${k2-1} nk $nk\n");

    /* ------------------------------------------------------------------ */
    /* construct the kth block, C */
    /* ------------------------------------------------------------------ */

    Lnz[block] = EMPTY_D;
    int pc = 0;
    for (int k = k1; k < k2; k++) {
      final newcol = k - k1;
      Cp[newcol] = pc;
      final oldcol = Qbtf[k];
      final pend = Ap[oldcol + 1];
      for (int p = Ap[oldcol]; p < pend; p++) {
        int newrow = Pinv[Ai[p]];
        if (newrow < k1) {
          nzoff++;
        } else {
          /* (newrow, newcol) is an entry in the block */
          ASSERT(newrow < k2);
          newrow -= k1;
          Ci[pc++] = newrow;
        }
      }
    }
    Cp[nk] = pc;
    maxnz = MAX(maxnz, pc);
    if (!NDEBUG) ASSERT_INT(_valid(nk, Cp, Ci, null));

    /* ------------------------------------------------------------------ */
    /* order the block C */
    /* ------------------------------------------------------------------ */

    if (nk <= 3) {

      /* -------------------------------------------------------------- */
      /* use natural ordering for tiny blocks (3-by-3 or less) */
      /* -------------------------------------------------------------- */

      for (int k = 0; k < nk; k++) {
        Pblk[k] = k;
      }
      lnz1 = nk * (nk + 1) / 2;
      flops1 = nk * (nk - 1) / 2 + (nk - 1) * nk * (2 * nk - 1) / 6;
      ok = TRUE;

    } else if (ordering == 0) {

      /* -------------------------------------------------------------- */
      /* order the block with AMD (C+C') */
      /* -------------------------------------------------------------- */

      final result = amd.order(nk, Cp, Ci, Pblk, null, amd_Info);
      ok = (result >= amd.AMD_OK) ? 1 : 0;
      if (result == amd.AMD_OUT_OF_MEMORY) {
        err = KLU_OUT_OF_MEMORY;
      }

      /* account for memory usage in AMD */
      Common.mempeak = MAX(Common.mempeak, Common.memusage +
          amd_Info[amd.AMD_MEMORY].toInt());

      /* get the ordering statistics from AMD */
      lnz1 = (amd_Info[amd.AMD_LNZ]) + nk;
      flops1 = 2 * amd_Info[amd.AMD_NMULTSUBS_LU] + amd_Info[amd.AMD_NDIV];
      if (pc == maxnz) {
        /* get the symmetry of the biggest block */
        Symbolic.symmetry = amd_Info[amd.AMD_SYMMETRY];
      }

    } else if (ordering == 1) {

      /* -------------------------------------------------------------- */
      /* order the block with COLAMD (C) */
      /* -------------------------------------------------------------- */

      /* order (and destroy) Ci, returning column permutation in Cp.
       * COLAMD "cannot" fail since the matrix has already been checked,
       * and Ci allocated. */

      ok = colamd.colamd(nk, nk, Cilen, Ci, Cp, null, cstats);
      lnz1 = EMPTY_D;
      flops1 = EMPTY_D;

      /* copy the permutation from Cp to Pblk */
      for (int k = 0; k < nk; k++) {
        Pblk[k] = Cp[k];
      }

    } else {

      /* -------------------------------------------------------------- */
      /* pass the block to the user-provided ordering function */
      /* -------------------------------------------------------------- */

      lnz1 = Common.user_order.order(nk, Cp, Ci, Pblk, Common);
      flops1 = EMPTY_D;
      ok = (lnz1 != 0) ? 1 : 0;
    }

    if (ok != 1) {
      return (err); // ordering method failed
    }

    /* ------------------------------------------------------------------ */
    /* keep track of nnz(L) and flops statistics */
    /* ------------------------------------------------------------------ */

    Lnz[block] = lnz1;
    lnz = (lnz == EMPTY || lnz1 == EMPTY) ? EMPTY : (lnz + lnz1);
    flops = (flops == EMPTY || flops1 == EMPTY) ? EMPTY : (flops + flops1);

    /* ------------------------------------------------------------------ */
    /* combine the preordering with the BTF ordering */
    /* ------------------------------------------------------------------ */

    PRINTF("Pblk, 1-based:\n");
    for (int k = 0; k < nk; k++) {
      ASSERT(k + k1 < n);
      ASSERT(Pblk[k] + k1 < n);
      Q[k + k1] = Qbtf[Pblk[k] + k1];
    }
    for (int k = 0; k < nk; k++) {
      ASSERT(k + k1 < n);
      ASSERT(Pblk[k] + k1 < n);
      P[k + k1] = Pbtf[Pblk[k] + k1];
    }
  }

  PRINTF("nzoff $nzoff  Ap[n] ${Ap [n]}\n");
  ASSERT(nzoff >= 0 && nzoff <= Ap[n]);

  /* return estimates of # of nonzeros in L including diagonal */
  Symbolic.lnz = lnz; // EMPTY if COLAMD used
  Symbolic.unz = lnz;
  Symbolic.nzoff = nzoff;
  Symbolic.est_flops = flops; // EMPTY if COLAMD or user-ordering used
  return (KLU_OK);
}

/**
 * Orders the matrix with or without BTF, then orders each block with AMD,
 * COLAMD, or the user ordering function. Does not handle the natural or
 * given ordering cases.
 *
 * A is [n]-by-[n]. Column pointers [Ap] is size `n+1`. Row indices [Ai] is
 * size [nz]. Returns `null` if error, or a valid [KLU_symbolic] object if
 * successful.
 */
KLU_symbolic order_and_analyze(final int n, final Int32List Ap,
                               final Int32List Ai, final KLU_common Common) {
  final work = new Float64List(1);
  final structural_rank = new Int32List(1);

  /* ---------------------------------------------------------------------- */
  /* allocate the Symbolic object, and check input matrix */
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

  final ordering = Common.ordering;
  int Cilen;
  if (ordering == 1) {
    /* COLAMD */
    Cilen = colamd.COLAMD_recommended(nz, n, n);
  } else if (ordering == 0 || (ordering == 3 && Common.user_order != null)) {
    /* AMD or user ordering function */
    Cilen = nz + 1;
  } else {
    /* invalid ordering */
    Common.status = KLU_INVALID;
    //klu_free_symbolic (Symbolic, Common) ;
    return (null);
  }

  /* AMD memory management routines */
  //amd_malloc  = Common.malloc_memory ;
  //amd_free    = Common.free_memory ;
  //amd_calloc  = Common.calloc_memory ;
  //amd_realloc = Common.realloc_memory ;

  /* ---------------------------------------------------------------------- */
  /* allocate workspace for BTF permutation */
  /* ---------------------------------------------------------------------- */

  final Pbtf = malloc_int(n, Common);
  final Qbtf = malloc_int(n, Common);
  if (Common.status < KLU_OK) {
    //KLU_free (Pbtf, n, sizeof (int), Common) ;
    //KLU_free (Qbtf, n, sizeof (int), Common) ;
    //klu_free_symbolic (Symbolic, Common) ;
    return (null);
  }

  /* ---------------------------------------------------------------------- */
  /* get the common parameters for BTF and ordering method */
  /* ---------------------------------------------------------------------- */

  int do_btf = Common.btf;
  do_btf = (do_btf != 0) ? TRUE : FALSE;
  Symbolic.ordering = ordering;
  Symbolic.do_btf = do_btf;
  Symbolic.structural_rank = EMPTY;

  /* ---------------------------------------------------------------------- */
  /* find the block triangular form (if requested) */
  /* ---------------------------------------------------------------------- */

  Common.work = 0.0;

  int nblocks, maxblock;
  if (do_btf != 0) {
    //Work = klu_malloc_int (5*n, Common) ;
    if (Common.status < KLU_OK) {
      /* out of memory */
      //klu_free (Pbtf, n, sizeof (int), Common) ;
      //klu_free (Qbtf, n, sizeof (int), Common) ;
      //klu_free_symbolic (Symbolic, Common) ;
      return (null);
    }

    nblocks = btf.order(n, Ap, Ai, Common.maxwork, work, Pbtf, Qbtf, R,
        structural_rank);
    Symbolic.structural_rank = structural_rank[0];
    Common.structural_rank = Symbolic.structural_rank;
    Common.work += work[0];

    //klu_free (Work, 5*n, sizeof (int), Common) ;

    /* unflip Qbtf if the matrix does not have full structural rank */
    if (Symbolic.structural_rank < n) {
      for (int k = 0; k < n; k++) {
        Qbtf[k] = btf.BTF_UNFLIP(Qbtf[k]);
      }
    }

    /* find the size of the largest block */
    maxblock = 1;
    for (int block = 0; block < nblocks; block++) {
      final k1 = R[block];
      final k2 = R[block + 1];
      final nk = k2 - k1;
      PRINTF("block $block size $nk\n");
      maxblock = MAX(maxblock, nk);
    }
  } else {
    /* BTF not requested */
    nblocks = 1;
    maxblock = n;
    R[0] = 0;
    R[1] = n;
    for (int k = 0; k < n; k++) {
      Pbtf[k] = k;
      Qbtf[k] = k;
    }
  }

  Symbolic.nblocks = nblocks;

  PRINTF("maxblock size $maxblock\n");
  Symbolic.maxblock = maxblock;

  /* ---------------------------------------------------------------------- */
  /* allocate more workspace, for analyze_worker */
  /* ---------------------------------------------------------------------- */

  final Pblk = malloc_int(maxblock, Common);
  final Cp = malloc_int(maxblock + 1, Common);
  final Ci = malloc_int(MAX(Cilen, nz + 1), Common);
  final Pinv = malloc_int(n, Common);

  /* ---------------------------------------------------------------------- */
  /* order each block of the BTF ordering, and a fill-reducing ordering */
  /* ---------------------------------------------------------------------- */

  if (Common.status == KLU_OK) {
    PRINTF(("calling analyze_worker\n"));
    Common.status = analyze_worker(n, Ap, Ai, nblocks, Pbtf, Qbtf, R,
        ordering, P, Q, Lnz, Pblk, Cp, Ci, Cilen, Pinv, Symbolic, Common);
    PRINTF("analyze_worker done\n");
  }

  /* ---------------------------------------------------------------------- */
  /* free all workspace */
  /* ---------------------------------------------------------------------- */

  //klu_free (Pblk, maxblock, sizeof (int), Common) ;
  //klu_free (Cp, maxblock+1, sizeof (int), Common) ;
  //klu_free (Ci, MAX (Cilen, nz+1), sizeof (int), Common) ;
  //klu_free (Pinv, n, sizeof (int), Common) ;
  //klu_free (Pbtf, n, sizeof (int), Common) ;
  //klu_free (Qbtf, n, sizeof (int), Common) ;

  /* ---------------------------------------------------------------------- */
  /* return the symbolic object */
  /* ---------------------------------------------------------------------- */

  if (Common.status < KLU_OK) {
    //klu_free_symbolic (Symbolic, Common) ;
  }
  return (Symbolic);
}

/**
 * Order the matrix with BTF (or not), then order each block with AMD,
 * COLAMD, a natural ordering, or with a user-provided ordering function.
 *
 * A is [n]-by-[n]. Column pointers [Ap] is size `n+1`. Row indices [Ai] is
 * size `nz`. Returns `null` if error, or a valid [KLU_symbolic] object if
 * successful.
 */
KLU_symbolic analyze(final int n, final Int32List Ap, final Int32List Ai,
                     final KLU_common Common) {
  /* ---------------------------------------------------------------------- */
  /* get the control parameters for BTF and ordering method */
  /* ---------------------------------------------------------------------- */

  if (Common == null) {
    return (null);
  }
  Common.status = KLU_OK;
  Common.structural_rank = EMPTY;

  /* ---------------------------------------------------------------------- */
  /* order and analyze */
  /* ---------------------------------------------------------------------- */

  if (Common.ordering == 2) {
    /* natural ordering */
    return (analyze_given(n, Ap, Ai, null, null, Common));
  } else {
    /* order with P and Q */
    return (order_and_analyze(n, Ap, Ai, Common));
  }
}
