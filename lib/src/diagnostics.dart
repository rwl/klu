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
 * Linear algebraic diagnostics.
 */
part of edu.ufl.cise.klu.tdouble;

/**
 * Compute the reciprocal pivot growth factor. In MATLAB notation:
 *
 *   rgrowth = min (max (abs ((R \ A (p,q)) - F))) ./ max (abs (U)))
 *
 * Takes O(|A|+|U|) time.
 *
 * Returns [TRUE] if successful, [FALSE] otherwise.
 */
int rgrowth(final Int32List Ap, final Int32List Ai, final Float64List Ax,
            final KLU_symbolic Symbolic, final KLU_numeric Numeric,
            final KLU_common Common) {
  final len = new Int32List(1);
  final Ui_offset = new Int32List(1);
  final Ux_offset = new Int32List(1);

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  if (Common == null) {
    return (FALSE);
  }

  if (Symbolic == null || Ap == null || Ai == null || Ax == null) {
    Common.status = KLU_INVALID;
    return (FALSE);
  }

  if (Numeric == null) {
    /* treat this as a singular matrix */
    Common.rgrowth = 0.0;
    Common.status = KLU_SINGULAR;
    return (TRUE);
  }
  Common.status = KLU_OK;

  /* ---------------------------------------------------------------------- */
  /* compute the reciprocal pivot growth */
  /* ---------------------------------------------------------------------- */

  final Aentry = Ax;// as Float64List;
  final Pinv = Numeric.Pinv;
  final Rs = Numeric.Rs;
  final Q = Symbolic.Q;
  Common.rgrowth = 1.0;

  for (int i = 0; i < Symbolic.nblocks; i++) {
    final k1 = Symbolic.R[i];
    final k2 = Symbolic.R[i + 1];
    final nk = k2 - k1;
    if (nk == 1) {
      continue; // skip singleton blocks
    }
    final LU = Numeric.LUbx[i];
    final Uip = Numeric.Uip;
    final Uip_offset = k1;
    final Ulen = Numeric.Ulen;
    final Ulen_offset = k1;
    final Ukk = Numeric.Udiag;
    final Ukk_offset = k1;
    double min_block_rgrowth = 1.0;
    for (int j = 0; j < nk; j++) {
      double max_ai = 0.0;
      double max_ui = 0.0;
      final oldcol = Q[j + k1];
      final pend = Ap[oldcol + 1];
      for (int k = Ap[oldcol]; k < pend; k++) {
        final oldrow = Ai[k];
        final newrow = Pinv[oldrow];
        if (newrow < k1) {
          continue; // skip entry outside the block
        }
        ASSERT(newrow < k2);
        double aik;
        if (Rs != null) {
          //SCALE_DIV_ASSIGN (aik, Aentry [k], Rs [newrow]) ;
          aik = Aentry[k] / Rs[newrow];
        } else {
          aik = Aentry[k];
        }
        //ABS (temp, aik) ;
        final temp = ABS(aik);
        if (temp > max_ai) {
          max_ai = temp;
        }
      }

      final Ux = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset,
          Ui_offset, Ux_offset, j, len);
      for (int k = 0; k < len[0]; k++) {
        //ABS (temp, Ux [k]) ;
        final temp = ABS(Ux[Ux_offset[0] + k]);
        if (temp > max_ui) {
          max_ui = temp;
        }
      }
      /* consider the diagonal element */
      //ABS (temp, Ukk [j]) ;
      final temp = ABS(Ukk[Ukk_offset + j]);
      if (temp > max_ui) {
        max_ui = temp;
      }

      /* if max_ui is 0, skip the column */
      if (SCALAR_IS_ZERO(max_ui)) {
        continue;
      }
      final temp2 = max_ai / max_ui;
      if (temp2 < min_block_rgrowth) {
        min_block_rgrowth = temp2;
      }
    }

    if (min_block_rgrowth < Common.rgrowth) {
      Common.rgrowth = min_block_rgrowth;
    }
  }
  return (TRUE);
}

/**
 * Estimate the condition number. Uses Higham and Tisseur's algorithm
 * (A block algorithm for matrix 1-norm estimation, with applications to
 * 1-norm pseudospectra, SIAM J. Matrix Anal. Appl., 21(4):1185-1201, 2000.
 *
 * Takes about `O(|A|+5*(|L|+|U|))` time.
 *
 * Returns [TRUE] if successful, [FALSE] otherwise.
 */
int condest(final Int32List Ap, final Float64List Ax,
            final KLU_symbolic Symbolic, final KLU_numeric Numeric,
            final KLU_common Common) {

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  if (Common == null) {
    return (FALSE);
  }
  if (Symbolic == null || Ap == null || Ax == null) {
    Common.status = KLU_INVALID;
    return (FALSE);
  }
  double abs_value = 0.0;
  if (Numeric == null) {
    /* treat this as a singular matrix */
    Common.condest = 1 / abs_value;
    Common.status = KLU_SINGULAR;
    return (TRUE);
  }
  Common.status = KLU_OK;

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */

  final n = Symbolic.n;
  final Udiag = Numeric.Udiag;

  /* ---------------------------------------------------------------------- */
  /* check if diagonal of U has a zero on it */
  /* ---------------------------------------------------------------------- */

  for (int i = 0; i < n; i++) {
    //ABS (abs_value, Udiag [i]) ;
    abs_value = ABS(Udiag[i]);
    if (SCALAR_IS_ZERO(abs_value)) {
      Common.condest = 1 / abs_value;
      Common.status = KLU_SINGULAR;
      return (TRUE);
    }
  }

  /* ---------------------------------------------------------------------- */
  /* compute 1-norm (maximum column sum) of the matrix */
  /* ---------------------------------------------------------------------- */

  double anorm = 0.0;
  final Aentry = Ax;// as Float64List ;
  for (int i = 0; i < n; i++) {
    final pend = Ap[i + 1];
    double csum = 0.0;
    for (int j = Ap[i]; j < pend; j++) {
      //ABS (abs_value, Aentry [j]) ;
      abs_value = ABS(Aentry[j]);
      csum += abs_value;
    }
    if (csum > anorm) {
      anorm = csum;
    }
  }

  /* ---------------------------------------------------------------------- */
  /* compute estimate of 1-norm of inv (A) */
  /* ---------------------------------------------------------------------- */

  /* get workspace (size 2*n double's) */
  final X = Numeric.Xwork;
  /* size n space used in KLU_solve, tsolve */
  //X += n ;                       /* X is size n */
  int X_offset = n;
  //S = X + n ;                    /* S is size n */
  final S = X;
  int S_offset = 2 * n;

  for (int i = 0; i < n; i++) {
    CLEAR(S, S_offset + i);
    CLEAR(X, X_offset + i);
    //REAL (X [i]) = 1.0 / ((double) n) ;
    X[X_offset + i] = 1.0 / n;
  }
  int jmax = 0;

  double ainv_norm = 0.0;
  for (int i = 0; i < 5; i++) {
    if (i > 0) {
      /* X [jmax] is the largest entry in X */
      for (int j = 0; j < n; j++) {
        CLEAR(X, X_offset + j);
      }
      //REAL (X [jmax]) = 1 ;
      X[X_offset + jmax] = 1.0;
    }

    solve(Symbolic, Numeric, n, 1, X/* as Float64List*/, X_offset, Common);
    final est_old = ainv_norm;
    ainv_norm = 0.0;

    for (int j = 0; j < n; j++) {
      /* ainv_norm += ABS (X [j]) ;*/
      //ABS (abs_value, X [j]) ;
      abs_value = ABS(X[X_offset + j]);
      ainv_norm += abs_value;
    }

    int unchanged = TRUE;

    for (int j = 0; j < n; j++) {
      double s = (X[X_offset + j] >= 0) ? 1.0 : -1.0;
      if (s != S[S_offset + j]) // s != REAL (S [j])
      {
        S[S_offset + j] = s;
        unchanged = FALSE;
      }
    }

    if (i > 0 && (ainv_norm <= est_old || unchanged == 1)) {
      break;
    }

    for (int j = 0; j < n; j++) {
      X[j] = S[S_offset + j];
    }

    /* do a transpose solve */
    tsolve(Symbolic, Numeric, n, 1, X, X_offset, Common);

    /* jnew = the position of the largest entry in X */
    int jnew = 0;
    double Xmax = 0.0;
    for (int j = 0; j < n; j++) {
      //ABS (xj, X [j]) ;
      final xj = ABS(X[X_offset + j]);
      if (xj > Xmax) {
        Xmax = xj;
        jnew = j;
      }
    }
    if (i > 0 && jnew == jmax) {
      /* the position of the largest entry did not change
       * from the previous iteration */
      break;
    }
    jmax = jnew;
  }

  /* ---------------------------------------------------------------------- */
  /* compute another estimate of norm(inv(A),1), and take the largest one */
  /* ---------------------------------------------------------------------- */

  for (int j = 0; j < n; j++) {
    CLEAR(X, X_offset + j);
    if (j % 2 != 0) {
      //REAL (X [j]) = 1 + ((double) j) / ((double) (n-1)) ;
      X[X_offset + j] = 1 + j / n - 1;
    } else {
      //REAL (X [j]) = -1 - ((double) j) / ((double) (n-1)) ;
      X[X_offset + j] = -1 - j / n - 1;
    }
  }

  solve(Symbolic, Numeric, n, 1, X/* as Float64List*/, X_offset, Common);

  double est_new = 0.0;
  for (int j = 0; j < n; j++) {
    /* est_new += ABS (X [j]) ;*/
    //ABS (abs_value, X [j]) ;
    abs_value = ABS(X[X_offset + j]);
    est_new += abs_value;
  }
  est_new = 2 * est_new / (3 * n);
  ainv_norm = MAX(est_new, ainv_norm);

  /* ---------------------------------------------------------------------- */
  /* compute estimate of condition number */
  /* ---------------------------------------------------------------------- */

  Common.condest = ainv_norm * anorm;
  return (TRUE);
}

/**
 * Compute the flop count for the LU factorization (in Common.flops)
 *
 * Returns [TRUE] if successful, [FALSE] otherwise.
 */
int flops(final KLU_symbolic Symbolic, final KLU_numeric Numeric,
          final KLU_common Common) {

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  if (Common == null) {
    return (FALSE);
  }
  Common.flops = EMPTY_D;
  if (Numeric == null || Symbolic == null) {
    Common.status = KLU_INVALID;
    return (FALSE);
  }
  Common.status = KLU_OK;

  /* ---------------------------------------------------------------------- */
  /* get the contents of the Symbolic object */
  /* ---------------------------------------------------------------------- */

  final R = Symbolic.R;
  final nblocks = Symbolic.nblocks;

  /* ---------------------------------------------------------------------- */
  /* get the contents of the Numeric object */
  /* ---------------------------------------------------------------------- */

  final LUbx = Numeric.LUbx;// as List<Float64List> ;

  /* ---------------------------------------------------------------------- */
  /* compute the flop count */
  /* ---------------------------------------------------------------------- */

  double flops = 0.0;
  for (int block = 0; block < nblocks; block++) {
    final k1 = R[block];
    final nk = R[block + 1] - k1;
    if (nk > 1) {
      final Llen = Numeric.Llen;
      final Llen_offset = k1;
      final Uip = Numeric.Uip;
      final Uip_offset = k1;
      final Ulen = Numeric.Ulen;
      final Ulen_offset = k1;
      final LU = LUbx[block];
      final Ui_offset = new Int32List(1);
      for (int k = 0; k < nk; k++) {
        /* compute kth column of U, and update kth column of A */
        final Ui = GET_I_POINTER(LU, Uip, Uip_offset, Ui_offset, k);
        final ulen = Ulen[Ulen_offset + k];
        for (int p = 0; p < ulen; p++) {
          flops += 2 * Llen[Llen_offset + Ui[Ui_offset[0] + p].toInt().toInt()];
        }
        /* gather and divide by pivot to get kth column of L */
        flops += Llen[Llen_offset + k];
      }
    }
  }
  Common.flops = flops;
  return (TRUE);
}

/**
 * Compute a really cheap estimate of the reciprocal of the condition number,
 * condition number, min(abs(diag(U))) / max(abs(diag(U))).  If U has a zero
 * pivot, or a NaN pivot, rcond will be zero.  Takes O(n) time.
 *
 * Result is in `Common.rcond`. Returns [TRUE] if successful, [FALSE] otherwise.
 */
int rcond(final KLU_symbolic Symbolic, final KLU_numeric Numeric,
          final KLU_common Common) {

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  if (Common == null) {
    return (FALSE);
  }
  if (Symbolic == null) {
    Common.status = KLU_INVALID;
    return (FALSE);
  }
  if (Numeric == null) {
    Common.rcond = 0.0;
    Common.status = KLU_SINGULAR;
    return (TRUE);
  }
  Common.status = KLU_OK;

  /* ---------------------------------------------------------------------- */
  /* compute rcond */
  /* ---------------------------------------------------------------------- */

  final n = Symbolic.n;
  final Udiag = Numeric.Udiag;
  double umin = 0.0;
  double umax = 0.0;
  for (int j = 0; j < n; j++) {
    /* get the magnitude of the pivot */
    //ABS (ukk, Udiag [j]) ;
    final ukk = ABS(Udiag[j]);
    if (SCALAR_IS_NAN(ukk) || SCALAR_IS_ZERO(ukk)) {
      /* if NaN, or zero, the rcond is zero */
      Common.rcond = 0.0;
      Common.status = KLU_SINGULAR;
      return (TRUE);
    }
    if (j == 0) {
      /* first pivot entry */
      umin = ukk;
      umax = ukk;
    } else {
      /* subsequent pivots */
      umin = MIN(umin, ukk);
      umax = MAX(umax, ukk);
    }
  }

  Common.rcond = umin / umax;
  if (SCALAR_IS_NAN(Common.rcond) || SCALAR_IS_ZERO(Common.rcond)) {
    /* this can occur if umin or umax are Inf or NaN */
    Common.rcond = 0.0;
    Common.status = KLU_SINGULAR;
  }
  return (TRUE);
}
