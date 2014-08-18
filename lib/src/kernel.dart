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

/**
 * Sparse left-looking LU factorization, with partial pivoting.  Based on
 * Gilbert & Peierl's method, with a non-recursive DFS and with Eisenstat &
 * Liu's symmetric pruning.  No user-callable routines are in this file.
 */
part of edu.ufl.cise.klu.tdouble;

/**
 * Does a depth-first-search, starting at node j.
 *
 * [j] is the node at which to start the DFS. [k] is the mark value, for the
 * Flag array. `Pinv[i] = k` if row `i` is `k`th pivot row, or [EMPTY] if
 * row `i` is not yet pivotal. [Llen] is size `n` and `Llen[k] =` # nonzeros
 * in column `k` of `L`. [Lip] is size `n` and Lip[k] is the position in LU
 * of column [k] of `L`. [Stack] is size `n`. `Flag[i] == k` means `i` is
 * marked. [Lpend] is for symmetric pruning. [top] is the top of the stack on
 * input. `Li` row index array of the `k`th column is [Lik]. [Ap_pos] keeps
 * track of position in adj list during DFS.
 */
int _dfs(final int j, final int k, final Int32List Pinv, final Int32List Llen,
         final int Llen_offset, final Int32List Lip, final int Lip_offset,
         final Int32List Stack, final Int32List Flag, final Int32List Lpend,
         int top, final Float64List LU, final Float64List Lik,
         final int Lik_offset, final Int32List plength,
         final Int32List Ap_pos) {
  int pos, l_length = plength[0];

  int head = 0;
  Stack[0] = j;
  ASSERT(Flag[j] != k);

  while (head >= 0) {
    final j = Stack[head];
    final jnew = Pinv[j];
    ASSERT(jnew >= 0 && jnew < k); // j is pivotal

    if (Flag[j] != k) { // a node is not yet visited
      /* first time that j has been visited */
      Flag[j] = k;
      PRINTF("[ start dfs at $j : new $jnew\n");
      /* set Ap_pos [head] to one past the last entry in col j to scan */
      Ap_pos[head] = (Lpend[jnew] == EMPTY) ? Llen[Llen_offset + jnew] : Lpend[jnew];
    }

    /* add the adjacent nodes to the recursive stack by iterating through
     * until finding another non-visited pivotal node */
    final Li = LU;
    int Li_offset = Lip[Lip_offset + jnew];
    for (pos = --Ap_pos[head]; pos >= 0; --pos) {
      final i = Li[Li_offset + pos].toInt();
      if (Flag[i] != k) {
        /* node i is not yet visited */
        if (Pinv[i] >= 0) {
          /* keep track of where we left off in the scan of the
           * adjacency list of node j so we can restart j where we
           * left off. */
          Ap_pos[head] = pos;

          /* node i is pivotal; push it onto the recursive stack
           * and immediately break so we can recurse on node i. */
          Stack[++head] = i;
          break;
        } else {
          /* node i is not pivotal (no outgoing edges). */
          /* Flag as visited and store directly into L,
           * and continue with current node j. */
          Flag[i] = k;
          Lik[Lik_offset + l_length] = i.toDouble();
          l_length++;
        }
      }
    }

    if (pos == -1) {
      /* if all adjacent nodes of j are already visited, pop j from
       * recursive stack and push j onto output stack */
      head--;
      Stack[--top] = j;
      PRINTF("  end   dfs at $j ] head : $head\n");
    }
  }

  plength[0] = l_length;
  return (top);
}

/**
 * Finds the pattern of x, for the solution of Lx=b.
 *
 * L is [n]-by-[n], where `n >= 0`. [k] is also used as the mark value for
 * the [Flag] array. `Pinv[i] = k` if `i` is `k`th pivot row, or EMPTY if row
 * `i` is not yet pivotal. [Stack] is size [n]. [Flag] is size [n]. Initially,
 * all of `Flag[0..n-1] < k`. After [_lsolve_symbolic] is done, `Flag[i] == k`
 * if `i` is in the pattern of the output, and `Flag[0..n-1] <= k`. [Lpend] is
 * for symmetric pruning. [Ap_pos] is a workspace used in dfs. [LU] factors
 * (pattern and values). [lup] is a pointer to free space in LU. [Llen] size
 * [n] and `Llen[k] =` number of nonzeros in column `k` of `L`. [Lip] is size
 * [n] and `Lip[k]` is position in [LU] of column [k] of `L`. [k1] is the
 * block of A and is from `k1` to `k2-1`. [PSinv] inverse of P from symbolic
 * factorization.
 */
int _lsolve_symbolic(final int n, final int k, final Int32List Ap,
                     final Int32List Ai, final Int32List Q,
                     final Int32List Pinv, final Int32List Stack,
                     final Int32List Flag, final Int32List Lpend,
                     final Int32List Ap_pos, final Float64List LU,
                     final int lup, final Int32List Llen,
                     final int Llen_offset, final Int32List Lip,
                     final int Lip_offset, final int k1,
                     final Int32List PSinv) {
  int top = n;
  final l_length = new Int32List.fromList([0]);
  final Lik = LU;
  final Lik_offset = lup;

  /* ---------------------------------------------------------------------- */
  /* BTF factorization of A (k1:k2-1, k1:k2-1) */
  /* ---------------------------------------------------------------------- */

  final kglobal = k + k1; // column k of the block is col kglobal of A
  final oldcol = Q[kglobal]; // Q must be present for BTF case
  final pend = Ap[oldcol + 1];
  for (int p = Ap[oldcol]; p < pend; p++) {
    final i = PSinv[Ai[p]] - k1;
    if (i < 0) continue; // skip entry outside the block

    /* (i,k) is an entry in the block.  start a DFS at node i */
    PRINTF("\n ===== DFS at node $i in b, inew: ${Pinv [i]}\n");
    if (Flag[i] != k) {
      if (Pinv[i] >= 0) {
        top = _dfs(i, k, Pinv, Llen, Llen_offset, Lip, Lip_offset, Stack,
            Flag, Lpend, top, LU, Lik, Lik_offset, l_length, Ap_pos);
      } else {
        /* i is not pivotal, and not flagged. Flag and put in L */
        Flag[i] = k;
        Lik[Lik_offset + l_length[0]] = i.toDouble();
        l_length[0]++;
      }
    }
  }

  /* If Llen [k] is zero, the matrix is structurally singular */
  Llen[Llen_offset + k] = l_length[0];
  return (top);
}

/**
 * Construct the kth column of A, and the off-diagonal part, if requested.
 * Scatter the numerical values into the workspace X, and construct the
 * corresponding column of the off-diagonal matrix.
 *
 * [k] the column of `A` (or the column of the block) to get.
 * [Q] is the column pre-ordering.
 * [k1] the block of `A` is from `k1` to `k2-1`. [PSinv] is the inverse of
 * [P] from symbolic factorization. [Rs] are the scale factors for `A`.
 * [scale] 0: no scaling, nonzero: scale the rows with [Rs]. [Offp]
 * off-diagonal matrix (modified by this routine).
 */
void _construct_column(final int k, final Int32List Ap, final Int32List Ai,
                       final Float64List Ax, final Int32List Q,
                       final Float64List X, final int k1,
                       final Int32List PSinv, final Float64List Rs,
                       final int scale, final Int32List Offp,
                       final Int32List Offi, final Float64List Offx) {
//  double aik;
//  int i, p, pend, oldcol, kglobal, poff, oldrow;

  /* ---------------------------------------------------------------------- */
  /* Scale and scatter the column into X. */
  /* ---------------------------------------------------------------------- */

  final kglobal = k + k1; // column k of the block is col kglobal of A
  int poff = Offp[kglobal]; // start of off-diagonal column
  final oldcol = Q[kglobal];
  final pend = Ap[oldcol + 1];

  if (scale <= 0) {
    /* no scaling */
    for (int p = Ap[oldcol]; p < pend; p++) {
      final oldrow = Ai[p];
      final i = PSinv[oldrow] - k1;
      final aik = Ax[p];
      if (i < 0) {
        /* this is an entry in the off-diagonal part */
        Offi[poff] = oldrow;
        Offx[poff] = aik;
        poff++;
      } else {
        /* (i,k) is an entry in the block. scatter into X */
        X[i] = aik;
      }
    }
  } else {
    /* row scaling */
    for (int p = Ap[oldcol]; p < pend; p++) {
      final oldrow = Ai[p];
      final i = PSinv[oldrow] - k1;
      double aik = Ax[p];
      aik = SCALE_DIV(aik, Rs[oldrow]);
      if (i < 0) {
        /* this is an entry in the off-diagonal part */
        Offi[poff] = oldrow;
        Offx[poff] = aik;
        poff++;
      } else {
        /* (i,k) is an entry in the block. scatter into X */
        X[i] = aik;
      }
    }
  }

  Offp[kglobal + 1] = poff;
  /* start of the next col of off-diag part */
}

/**
 * Computes the numerical values of `x`, for the solution of `Lx=b`. Note that
 * `x` may include explicit zeros if numerical cancelation occurs. L is
 * assumed to be unit-diagonal, with possibly unsorted columns (but the first
 * entry in the column must always be the diagonal entry).
 *
 * `Pinv[i] = k` if `i` is `k`th pivot row, or EMPTY if row `i` is not yet
 * pivotal. [LU] factors (pattern and values). [Stack] is the stack for dfs.
 * [Lip] is size [n] and `Lip[k]` is position in LU of column `k` of `L`.
 * [top] is the top of the stack on input. A is [n]-by-[n]. [Llen] is size
 * [n] and `Llen[k] =` number of nonzeros in column `k` of `L`. X is size [n],
 * initially zero. On output, `X[Ui[up1..up-1]]` and `X[Li[lp1..lp-1]]`
 * contains the solution.
 */
void _lsolve_numeric(final Int32List Pinv, final Float64List LU,
                     final Int32List Stack, final Int32List Lip,
                     final int Lip_offset, final int top, final int n,
                     final Int32List Llen, final int Llen_offset,
                     final Float64List X) {
  final len = new Int32List(1);
  final Li_offset = new Int32List(1);
  final Lx_offset = new Int32List(1);

  /* solve Lx=b */
  for (int s = top; s < n; s++) {
    /* forward solve with column j of L */
    final j = Stack[s];
    final jnew = Pinv[j];
    ASSERT(jnew >= 0);
    final xj = X[j];
    final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, jnew, len);
    final Lx = Li;
    ASSERT(Lip[Lip_offset + jnew] <= Lip[Lip_offset + jnew + 1]);
    for (int p = 0; p < len[0]; p++) {
      //MULT_SUB (X [Li [p]], Lx [p], xj) ;
      X[Li[Li_offset[0] + p].toInt()] -= Lx[Lx_offset[0] + p] * xj;
    }
  }
}

/**
 * Find a pivot via partial pivoting, and scale the column of L.
 */
int _lpivot(final int diagrow, final Int32List p_pivrow,
            final Float64List p_pivot, final Float64List p_abs_pivot,
            final double tol, final Float64List X, final Float64List LU,
            final Int32List Lip, final int Lip_offset, final Int32List Llen,
            final int Llen_offset, final int k, final int n,
            final Int32List Pinv, final Int32List p_firstrow,
            final KLU_common Common) {
  double pivot;
  Float64List Lx;
//  double abs_pivot, xabs;
//  int p, i, ppivrow, pdiag, pivrow, last_row_index, firstrow;
  /*Int32List*/Float64List Li;
  final len = new Int32List(1);
  final Li_offset = new Int32List(1);
  final Lx_offset = new Int32List(1);

  int pivrow = EMPTY;
  if (Llen[Llen_offset + k] == 0) {
    /* matrix is structurally singular */
    if (Common.halt_if_singular != 0) {
      return (FALSE);
    }
    int firstrow;
    for (firstrow = p_firstrow[0]; firstrow < n; firstrow++) {
      PRINTF("check $firstrow\n");
      if (Pinv[firstrow] < 0) {
        /* found the lowest-numbered non-pivotal row. Pick it. */
        pivrow = firstrow;
        PRINTF("Got pivotal row: $pivrow\n");
        break;
      }
    }
    ASSERT(pivrow >= 0 && pivrow < n);
    pivot = 0.0; //CLEAR (pivot) ;
    p_pivrow[0] = pivrow;
    p_pivot[0] = pivot;
    p_abs_pivot[0] = 0.0;
    p_firstrow[0] = firstrow;
    return (FALSE);
  }

  int pdiag = EMPTY;
  int ppivrow = EMPTY;
  double abs_pivot = EMPTY.toDouble();
  final i = Llen[Llen_offset + k] - 1;
  Li = Lx = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);
  final last_row_index = Li[Li_offset[0] + i].toInt();

  /* decrement the length by 1 */
  Llen[Llen_offset + k] = i;
  Li = Lx = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);

  /* look in Li [0 ..Llen [k] - 1 ] for a pivot row */
  for (int p = 0; p < len[0]; p++) {
    /* gather the entry from X and store in L */
    final i = Li[Li_offset[0] + p].toInt();
    final x = X[i];
    CLEAR(X, i);

    Lx[Lx_offset[0] + p] = x;
    //ABS (xabs, x) ;
    final xabs = ABS(x);

    /* find the diagonal */
    if (i == diagrow) {
      pdiag = p;
    }

    /* find the partial-pivoting choice */
    if (xabs > abs_pivot) {
      abs_pivot = xabs;
      ppivrow = p;
    }
  }

  //ABS (xabs, X [last_row_index]) ;
  final xabs = ABS(X[last_row_index]);
  if (xabs > abs_pivot) {
    abs_pivot = xabs;
    ppivrow = EMPTY;
  }

  /* compare the diagonal with the largest entry */
  if (last_row_index == diagrow) {
    if (xabs >= tol * abs_pivot) {
      abs_pivot = xabs;
      ppivrow = EMPTY;
    }
  } else if (pdiag != EMPTY) {
    //ABS (xabs, Lx [pdiag]) ;
    final xabs = ABS(Lx[Lx_offset[0] + pdiag]);
    if (xabs >= tol * abs_pivot) {
      /* the diagonal is large enough */
      abs_pivot = xabs;
      ppivrow = pdiag;
    }
  }

  if (ppivrow != EMPTY) {
    pivrow = Li[Li_offset[0] + ppivrow].toInt();
    pivot = Lx[Lx_offset[0] + ppivrow];
    /* overwrite the ppivrow values with last index values */
    Li[Li_offset[0] + ppivrow] = last_row_index.toDouble();
    Lx[Lx_offset[0] + ppivrow] = X[last_row_index];
  } else {
    pivrow = last_row_index;
    pivot = X[last_row_index];
  }
  CLEAR(X, last_row_index);

  p_pivrow[0] = pivrow;
  p_pivot[0] = pivot;
  p_abs_pivot[0] = abs_pivot;
  ASSERT(pivrow >= 0 && pivrow < n);

  if (IS_ZERO(pivot) && Common.halt_if_singular != 0) {
    /* numerically singular case */
    return (FALSE);
  }

  /* divide L by the pivot value */
  for (int p = 0; p < Llen[Llen_offset + k]; p++) {
    //DIV (Lx [p], Lx [p], pivot) ;
    Lx[Lx_offset[0] + p] /= pivot;
  }

  return (TRUE);
}

/**
 * Prune the columns of L to reduce work in subsequent depth-first searches.
 *
 * `Lpend[j]` marks symmetric pruning point for `L(:,j)`. `Pinv[i] = k` if
 * row `i` is `k`th pivot row, or [EMPTY] if row `i` is not yet pivotal.
 * Prune using column [k] of `U`. [pivrow] is the current pivot row. [LU]
 * factors (pattern and values). Column pointers [Uip] for `U` is size `n`.
 * Column pointers [Lip] for `L` is size n. Column length [Ulen] of `U` is
 * size `n`. Column length [Llen] of `L` is size `n`.
 */
void _prune(final Int32List Lpend, final Int32List Pinv, final int k,
            final int pivrow, final Float64List LU, final Int32List Uip,
            final int Uip_offset, final Int32List Lip, final int Lip_offset,
            final Int32List Ulen, final int Ulen_offset, final Int32List Llen,
            final int Llen_offset) {
  final llen = new Int32List(1);
  final ulen = new Int32List(1);
  final Li_offset = new Int32List(1);
  final Lx_offset = new Int32List(1);
  final Ui_offset = new Int32List(1);
  final Ux_offset = new Int32List(1);

  /* check to see if any column of L can be pruned */
  final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, ulen);
  for (int p = 0; p < ulen[0]; p++) {
    final j = Ui[Ui_offset[0] + p].toInt();
    ASSERT(j < k);
    PRINTF("$j is pruned: ${Lpend [j] != EMPTY ? 1 : 0}. Lpend[j] ${Lpend [j]} Lip[j+1] ${Lip [Lip_offset + j+1]}\n");
    if (Lpend[j] == EMPTY) {
      /* scan column j of L for the pivot row */
      final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, j, llen);
      final Lx = Li;
      for (int p2 = 0; p2 < llen[0]; p2++) {
        if (pivrow == Li[Li_offset[0] + p2]) {
          /* found it!  This column can be pruned */
          if (!NDEBUG) {
            PRINTF("==== PRUNE: col j $j of L\n");
            {
              int p3;
              for (p3 = 0; p3 < Llen[Llen_offset + j]; p3++) {
                PRINTF("before: ${Li [Li_offset[0] + p3].toInt()}  pivotal: ${Pinv [Li [Li_offset[0] + p3].toInt()] >= 0 ? 1 : 0}\n");
              }
            }
          }

          /* partition column j of L. The unit diagonal of L
           * is not stored in the column of L. */
          int phead = 0;
          int ptail = Llen[Llen_offset + j];
          while (phead < ptail) {
            final i = Li[Li_offset[0] + phead].toInt();
            if (Pinv[i] >= 0) {
              /* leave at the head */
              phead++;
            } else {
              /* swap with the tail */
              ptail--;
              Li[Li_offset[0] + phead] = Li[Li_offset[0] + ptail];
              Li[Li_offset[0] + ptail] = i.toDouble();
              final x = Lx[Lx_offset[0] + phead];
              Lx[Lx_offset[0] + phead] = Lx[Lx_offset[0] + ptail];
              Lx[Lx_offset[0] + ptail] = x;
            }
          }

          /* set Lpend to one past the last entry in the
           * first part of the column of L.  Entries in
           * Li [0 ... Lpend [j]-1] are the only part of
           * column j of L that needs to be scanned in the DFS.
           * Lpend [j] was EMPTY; setting it >= 0 also flags
           * column j as pruned. */
          Lpend[j] = ptail;

          if (!NDEBUG) {
            int p3;
            for (p3 = 0; p3 < Llen[Llen_offset + j]; p3++) {
              if (p3 == Lpend[j]) PRINTF(("----\n"));
              PRINTF("after: ${Li [Li_offset[0] + p3].toInt()}  pivotal: ${Pinv [Li [Li_offset[0] + p3].toInt()] >= 0 ? 1 : 0}\n");
            }
          }

          break;
        }
      }
    }
  }
}

/**
 * A is [n]-by-[n]. Column pointers [Ap] for A is size `n+1`. Row indices
 * [Ai] for A is size `nz = Ap[n]`. Values [Ax] of A is size `nz`.
 * Optional input permutation [Q] is size `n`. [lusize] is the initial size
 * of LU on input. Inverse row permutation [Pinv] is size [n] and
 * `Pinv[i] = k` if row `i` is the `k`th pivot row.
 * Row permutation [P] is size [n] where `P[k] = i` if row `i` is the `k`th
 * pivot row. LU array [p_LU] is size [lusize] on input.
 * Diagonal of U [Udiag] is size [n].
 * Column length [Llen] of L is size [n].
 * Column length [Ulen] of U is size [n].
 * Column pointers [Lip] for L is size [n].
 * Column pointers [Uip] for U is size [n].
 * [lnz] is size 1, size of L. [unz] is size 1, size of U.
 * [X] is size [n], undefined on input, zero on output.
 * [Lpend] is a size n workspace, for pruning only.
 * [k1] is the block of `A` is from `k1` to `k2-1`.
 * [PSinv] is the inverse of P from symbolic factorization.
 * [Rs] is scale factors for `A`. [Offp] is the off-diagonal matrix (modified
 * by this routine). Returns the final size of LU on output.
 */
int _kernel(final int n, final Int32List Ap, final Int32List Ai,
            final Float64List Ax, final Int32List Q, int lusize,
            final Int32List Pinv, final Int32List P,
            final List<Float64List> p_LU, final Float64List Udiag,
            final int Udiag_offset, final Int32List Llen,
            final int Llen_offset, final Int32List Ulen,
            final int Ulen_offset, final Int32List Lip,
            final int Lip_offset, final Int32List Uip,
            final int Uip_offset, final Int32List lnz,
            final Int32List unz, final Float64List X,
            final Int32List Stack, final Int32List Flag,
            final Int32List Ap_pos, final Int32List Lpend,
            final int k1, final Int32List PSinv, final Float64List Rs,
            final Int32List Offp, final Int32List Offi,
            final Float64List Offx, final KLU_common Common) {
  final pivot = new Float64List(1);
  final abs_pivot = new Float64List(1);
  Float64List Ux;
  /*Int32List*/Float64List Li, Ui;
  Float64List LU; // LU factors (pattern and values)
  final len = new Int32List(1);
  final firstrow = new Int32List(1);
  final pivrow = new Int32List.fromList([0]);
  int newlusize;
  final Ui_offset = new Int32List(1);
  final Ux_offset = new Int32List(1);
  final Li_offset = new Int32List(1);
  final Lx_offset = new Int32List(1);

  Float64List Lx; // only used when debugging

  ASSERT(Common != null);
  final scale = Common.scale;
  final tol = Common.tol;
  final memgrow = Common.memgrow;
  lnz[0] = 0;
  unz[0] = 0;
  pivot[0] = 0.0; //CLEAR (pivot) ;

  /* ---------------------------------------------------------------------- */
  /* get initial Li, Lx, Ui, and Ux */
  /* ---------------------------------------------------------------------- */

  PRINTF("input: lusize $lusize \n");
  ASSERT(lusize > 0);
  LU = p_LU[0];

  /* ---------------------------------------------------------------------- */
  /* initializations */
  /* ---------------------------------------------------------------------- */

  firstrow[0] = 0;
  int lup = 0;

  for (int k = 0; k < n; k++) {
    /* X [k] = 0 ; */
    CLEAR(X, k);
    Flag[k] = EMPTY;
    Lpend[k] = EMPTY; // flag k as not pruned
  }

  /* ---------------------------------------------------------------------- */
  /* mark all rows as non-pivotal and determine initial diagonal mapping */
  /* ---------------------------------------------------------------------- */

  /* PSinv does the symmetric permutation, so don't do it here */
  for (int k = 0; k < n; k++) {
    P[k] = k;
    Pinv[k] = FLIP(k); // mark all rows as non-pivotal
  }
  /* initialize the construction of the off-diagonal matrix */
  Offp[0] = 0;

  /* P [k] = row means that UNFLIP (Pinv [row]) = k, and visa versa.
   * If row is pivotal, then Pinv [row] >= 0.  A row is initially "flipped"
   * (Pinv [k] < EMPTY), and then marked "unflipped" when it becomes
   * pivotal. */

  if (!NDEBUG) {
    for (int k = 0; k < n; k++) {
      PRINTF("Initial P [$k] = ${P [k]}\n");
    }
  }

  /* ---------------------------------------------------------------------- */
  /* factorize */
  /* ---------------------------------------------------------------------- */

  for (int k = 0; k < n; k++) {

    PRINTF("\n\n==================================== k: $k\n");

    /* ------------------------------------------------------------------ */
    /* determine if LU factors have grown too big */
    /* ------------------------------------------------------------------ */

    /* (n - k) entries for L and k entries for U */
    //nunits = DUNITS (Integer, n - k) + DUNITS (Integer, k) +
    //	DUNITS (Double, n - k) + DUNITS (Double, k) ;
    final nunits = ((n - k) + (k) + (n - k) + (k)).toDouble();

    /* LU can grow by at most 'nunits' entries if the column is dense */
    PRINTF("lup $lup lusize $lusize lup+nunits: ${lup+nunits}\n");
    double xsize = lup + nunits;
    if (xsize > lusize) {
      /* check here how much to grow */
      xsize = memgrow * lusize + 4 * n + 1;
      if (INT_OVERFLOW(xsize)) {
        PRINTF("Matrix is too large (int overflow)\n");
        Common.status = KLU_TOO_LARGE;
        return (lusize);
      }
      newlusize = (memgrow * lusize + 2 * n + 1).toInt();
      /* Future work: retry mechanism in case of malloc failure */
      LU = realloc_dbl(newlusize, lusize, LU, Common);
      Common.nrealloc++;
      p_LU[0] = LU;
      if (Common.status == KLU_OUT_OF_MEMORY) {
        PRINTF("Matrix is too large (LU)\n");
        return (lusize);
      }
      lusize = newlusize;
      PRINTF("inc LU to $lusize done\n");
    }

    /* ------------------------------------------------------------------ */
    /* start the kth column of L and U */
    /* ------------------------------------------------------------------ */

    Lip[Lip_offset + k] = lup;

    /* ------------------------------------------------------------------ */
    /* compute the nonzero pattern of the kth column of L and U */
    /* ------------------------------------------------------------------ */

    if (!NDEBUG) {
      for (int i = 0; i < n; i++) {
        ASSERT(Flag[i] < k);
        /* ASSERT (X [i] == 0) ; */
        ASSERT(IS_ZERO(X[i]));
      }
    }

    final top = _lsolve_symbolic(n, k, Ap, Ai, Q, Pinv, Stack, Flag, Lpend, Ap_pos, LU, lup, Llen, Llen_offset, Lip, Lip_offset, k1, PSinv);

    if (!NDEBUG) {
      PRINTF("--- in U:\n");
      for (int p = top; p < n; p++) {
        PRINTF("pattern of X for U: $p : ${Stack [p]} pivot row: ${Pinv [Stack [p]]}\n");
        ASSERT(Flag[Stack[p]] == k);
      }
      PRINTF("--- in L:\n");
      Li = LU;
      Li_offset[0] = Lip[Lip_offset + k];
      for (int p = 0; p < Llen[Llen_offset + k]; p++) {
        PRINTF("pattern of X in L: $p : ${Li [Li_offset[0] + p]} pivot row: ${Pinv [Li [Li_offset[0] + p].toInt()]}\n");
        ASSERT(Flag[Li[Li_offset[0] + p].toInt()] == k);
      }
      int p = 0;
      for (int i = 0; i < n; i++) {
        ASSERT(Flag[i] <= k);
        if (Flag[i] == k) p++;
      }
    }

    /* ------------------------------------------------------------------ */
    /* get the column of the matrix to factorize and scatter into X */
    /* ------------------------------------------------------------------ */

    _construct_column(k, Ap, Ai, Ax, Q, X, k1, PSinv, Rs, scale, Offp, Offi, Offx);

    /* ------------------------------------------------------------------ */
    /* compute the numerical values of the kth column (s = L \ A (:,k)) */
    /* ------------------------------------------------------------------ */

    _lsolve_numeric(Pinv, LU, Stack, Lip, Lip_offset, top, n, Llen, Llen_offset, X);

    if (!NDEBUG) {
      for (int p = top; p < n; p++) {
        PRINTF("X for U ${Stack [p]} : ");
        PRINT_ENTRY(X[Stack[p]]);
      }
      Li = LU;
      Li_offset[0] = Lip[Lip_offset + k];
      for (int p = 0; p < Llen[Llen_offset + k]; p++) {
        PRINTF("X for L ${Li [Li_offset[0] + p].toInt()} : ");
        PRINT_ENTRY(X[Li[Li_offset[0] + p].toInt()]);
      }
    }

    /* ------------------------------------------------------------------ */
    /* partial pivoting with diagonal preference */
    /* ------------------------------------------------------------------ */

    /* determine what the "diagonal" is */
    final diagrow = P[k]; // might already be pivotal
    PRINTF("k $k, diagrow = $diagrow, UNFLIP (diagrow) = ${UNFLIP (diagrow)}\n");

    /* find a pivot and scale the pivot column */
    if (_lpivot(diagrow, pivrow, pivot, abs_pivot, tol, X, LU, Lip, Lip_offset, Llen, Llen_offset, k, n, Pinv, firstrow, Common) == 0) {
      /* matrix is structurally or numerically singular */
      Common.status = KLU_SINGULAR;
      if (Common.numerical_rank == EMPTY) {
        Common.numerical_rank = k + k1;
        Common.singular_col = Q[k + k1];
      }
      if (Common.halt_if_singular != 0) {
        /* do not continue the factorization */
        return (lusize);
      }
    }

    /* we now have a valid pivot row, even if the column has NaN's or
     * has no entries on or below the diagonal at all. */
    PRINTF("\nk $k : Pivot row ${pivrow[0]} : ");
    PRINT_ENTRY(pivot[0]);
    ASSERT(pivrow[0] >= 0 && pivrow[0] < n);
    ASSERT(Pinv[pivrow[0]] < 0);

    /* set the Uip pointer */
    //Uip [Uip_offset + k] = Lip [Lip_offset + k] +
    //		UNITS (Integer, Llen [Llen_offset + k]) +
    //		UNITS (Double, Llen [Llen_offset + k]) ;
    Uip[Uip_offset + k] = Lip[Lip_offset + k] + Llen[Llen_offset + k] + Llen[Llen_offset + k];

    /* move the lup pointer to the position where indices of U
     * should be stored */
    //lup += UNITS (Integer, Llen [Llen_offset + k]) +
    //		UNITS (Double, Llen [Llen_offset + k]) ;
    lup += Llen[Llen_offset + k] + Llen[Llen_offset + k];

    Ulen[Ulen_offset + k] = n - top;

    /* extract Stack [top..n-1] to Ui and the values to Ux and clear X */
    Ui = Ux = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len);

    int p = top;
    for (int i = 0; p < n; p++, i++) {
      final j = Stack[p];
      Ui[Ui_offset[0] + i] = Pinv[j].toDouble();
      Ux[Ux_offset[0] + i] = X[j];
      //CLEAR (X [j]) ;
      X[j] = 0.0;
    }

    /* position the lu index at the starting point for next column */
    //lup += UNITS (Integer, Ulen [Ulen_offset + k]) +
    //		UNITS (Double, Ulen [Ulen_offset + k]) ;
    lup += Ulen[Ulen_offset + k] + Ulen[Ulen_offset + k];

    /* U(k,k) = pivot */
    Udiag[Udiag_offset + k] = pivot[0];

    /* ------------------------------------------------------------------ */
    /* log the pivot permutation */
    /* ------------------------------------------------------------------ */

    ASSERT(UNFLIP(Pinv[diagrow]) < n);
    ASSERT(P[UNFLIP(Pinv[diagrow])] == diagrow);

    if (pivrow[0] != diagrow) {
      /* an off-diagonal pivot has been chosen */
      Common.noffdiag++;
      PRINTF(">>>>>>>>>>>>>>>>> pivrow ${pivrow[0]} k $k off-diagonal\n");
      if (Pinv[diagrow] < 0) {
        /* the former diagonal row index, diagrow, has not yet been
         * chosen as a pivot row.  Log this diagrow as the "diagonal"
         * entry in the column kbar for which the chosen pivot row,
         * pivrow, was originally logged as the "diagonal" */
        final kbar = FLIP(Pinv[pivrow[0]]);
        P[kbar] = diagrow;
        Pinv[diagrow] = FLIP(kbar);
      }
    }
    P[k] = pivrow[0];
    Pinv[pivrow[0]] = k;

    if (!NDEBUG) {
      for (int i = 0; i < n; i++) {
        ASSERT(IS_ZERO(X[i]));
      }
      Ui = Ux = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len);
      for (p = 0; p < len[0]; p++) {
        PRINTF("Column %k of U: ${Ui [Ui_offset[0] + p].toInt()} : ");
        PRINT_ENTRY(Ux[Ux_offset[0] + p]);
      }

      Li = Lx = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);
      for (p = 0; p < len[0]; p++) {
        PRINTF("Column $k of L: ${Li [Li_offset[0] + p].toInt()} : ");
        PRINT_ENTRY(Lx[Lx_offset[0] + p]);
      }
    }

    /* ------------------------------------------------------------------ */
    /* symmetric pruning */
    /* ------------------------------------------------------------------ */

    _prune(Lpend, Pinv, k, pivrow[0], LU, Uip, Uip_offset, Lip, Lip_offset, Ulen, Ulen_offset, Llen, Llen_offset);

    lnz[0] += Llen[Llen_offset + k] + 1; // 1 added to lnz for diagonal
    unz[0] += Ulen[Ulen_offset + k] + 1; // 1 added to unz for diagonal
  }

  /* ---------------------------------------------------------------------- */
  /* finalize column pointers for L and U, and put L in the pivotal order */
  /* ---------------------------------------------------------------------- */

  for (int p = 0; p < n; p++) {
    Li = LU;
    Li_offset[0] = Lip[Lip_offset + p];
    for (int i = 0; i < Llen[Llen_offset + p]; i++) {
      Li[Li_offset[0] + i] = Pinv[Li[Li_offset[0] + i].toInt()].toDouble();
    }
  }

  if (!NDEBUG) {
    for (int i = 0; i < n; i++) {
      PRINTF("P [$i] = ${P [i]}   Pinv [$i] = ${Pinv [i]}\n");
    }
    for (int i = 0; i < n; i++) {
      ASSERT(Pinv[i] >= 0 && Pinv[i] < n);
      ASSERT(P[i] >= 0 && P[i] < n);
      ASSERT(P[Pinv[i]] == i);
      ASSERT(IS_ZERO(X[i]));
    }
  }

  /* ---------------------------------------------------------------------- */
  /* shrink the LU factors to just the required size */
  /* ---------------------------------------------------------------------- */

  newlusize = lup;
  ASSERT(newlusize <= lusize);

  /* this cannot fail, since the block is descreasing in size */
  LU = realloc_dbl(newlusize, lusize, LU, Common);
  p_LU[0] = LU;
  return (newlusize);
}
