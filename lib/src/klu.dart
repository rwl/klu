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
 * KLU factorizes P*A into L*U, using the Gilbert-Peierls algorithm[1], with
 * optional symmetric pruning by Eisenstat and Liu[2].  The code is by Tim
 * Davis.  This algorithm is what appears as the default sparse LU routine in
 * MATLAB version 6.0, and still appears in MATLAB 6.5 as[L,U,P] = lu (A).
 * Note that no column ordering is provided (see COLAMD or AMD for suitable
 * orderings).  SuperLU is based on this algorithm, except that it adds the
 * use of dense matrix operations on "supernodes" (adjacent columns with
 * identical).  This code doesn't use supernodes, thus its name ("Kent" LU,
 * as in "Clark Kent", in contrast with Super-LU...).  This algorithm is slower
 * than SuperLU and UMFPACK for large matrices with lots of nonzeros in their
 * factors (such as for most finite-element problems).  However, for matrices
 * with very sparse LU factors, this algorithm is typically faster than both
 * SuperLU and UMFPACK, since in this case there is little chance to exploit
 * dense matrix kernels (the BLAS).
 *
 * Only one block of A is factorized, in the BTF form.  The input n is the
 * size of the block; k1 is the first row and column in the block.
 *
 * NOTE: no error checking is done on the inputs.  This version is not meant to
 * be called directly by the user.  Use klu_factor instead.
 *
 * No fill-reducing ordering is provided.  The ordering quality of
 * klu_kernel_factor is the responsibility of the caller.  The input A must
 * pre-permuted to reduce fill-in, or fill-reducing input permutation Q must
 * be provided.
 *
 * The input matrix A must be in compressed-column form, with either sorted
 * or unsorted row indices.  Row indices for column j of A is in
 * Ai[Ap[j] ... Ap[j+1]-1] and the same range of indices in Ax holds the
 * numerical values.  No duplicate entries are allowed.
 *
 * Copyright 2004-2009, Tim Davis.  All rights reserved.  See the README
 * file for details on permitted use.  Note that no code from The MathWorks,
 * Inc, or from SuperLU, or from any other source appears here.  The code is
 * written from scratch, from the algorithmic description in Gilbert & Peierls'
 * and Eisenstat & Liu's journal papers[1,2].
 *
 * If an input permutation Q is provided, the factorization L*U = A (P,Q)
 * is computed, where P is determined by partial pivoting, and Q is the input
 * ordering.  If the pivot tolerance is less than 1, the "diagonal" entry that
 * KLU attempts to choose is the diagonal of A (Q,Q).  In other words, the
 * input permutation is applied symmetrically to the input matrix.  The output
 * permutation P includes both the partial pivoting ordering and the input
 * permutation.  If Q is null, then it is assumed to be the identity
 * permutation.  Q is not modified.
 *
 * 1. Gilbert, J. R. and Peierls, T., "Sparse Partial Pivoting in Time
 *    Proportional to Arithmetic Operations," SIAM J. Sci. Stat. Comp.,
 *    vol 9, pp.  862-874, 1988.
 * 2. Eisenstat, S. C. and Liu, J. W. H., "Exploiting Structural Symmetry in
 *    Unsymmetric Sparse Symbolic Factorization," SIAM J. Matrix Analysis &
 *    Applic., vol 13, pp.  202-211, 1992.
 */
library edu.ufl.cise.klu.tdouble;

import 'dart:typed_data';

import './common/common.dart';

import 'package:btf/btf.dart' as btf;
import 'package:colamd/colamd.dart' as colamd;
import 'package:amd/amd.dart' as amd;

part 'analyze.dart';
part 'analyze_given.dart';
part 'defaults.dart';
part 'diagnostics.dart';
part 'dump.dart';
part 'extract.dart';
part 'factor.dart';
part 'internal.dart';
part 'kernel.dart';
part 'memory.dart';
part 'refactor.dart';
part 'scale.dart';
part 'solve.dart';
part 'sort.dart';
part 'tsolve.dart';
part 'version.dart';

/**
 * A is [n]-by-[n]. [n] must be > 0. Column pointers [Ap] for A is size n+1.
 * Row indices [Ai] for A is size `nz = Ap[n]`. Values of A are [Ax] and
 * size `nz`. Optional column permutation [Q] is size [n]. [Lsize] is
 * estimate of number of nonzeros in L. [p_LU] is the row indices and values
 * of L and U. The diagonal of U [Udiag] is size [n]. [Llen] is the column
 * length of L and size [n]. [Ulen] is the column length of U and size [n].
 * [Lip] is the column pointers for L and size [n]. [Uip] column pointers
 * for U and size [n]. [P] is the row permutation and size [n]. [lnz] is
 * size of L. [unz] is the size of U. [X] is size [n] and zero on output.
 * [Work] is size [5n]. [k1] is the block of A from [k1] to [k2-1]. [PSinv]
 * is the inverse of P from symbolic factorization. [Rs] are the scale
 * factors for A. [Offp] is the off-diagonal matrix (modified by this
 * routine).
 */
int kernel_factor(int n, final Int32List Ap, final Int32List Ai,
                  final Float64List Ax, final Int32List Q, double Lsize,
                  final List<Float64List> p_LU, final int block,
                  final Float64List Udiag, final int Udiag_offset,
                  final Int32List Llen, final int Llen_offset,
                  final Int32List Ulen, final int Ulen_offset,
                  final Int32List Lip, final int Lip_offset,
                  final Int32List Uip, final int Uip_offset,
                  final Int32List P, final Int32List lnz, final Int32List unz,
                  final Float64List X, final Int32List Work, final int k1,
                  final Int32List PSinv, final Float64List Rs,
                  final Int32List Offp, final Int32List Offi,
                  final Float64List Offx, final KLU_common Common) {
  List<Float64List> LU = new List<Float64List>(1);//[] ;
  ASSERT(Common != null);

  /* ---------------------------------------------------------------------- */
  /* get control parameters, or use defaults */
  /* ---------------------------------------------------------------------- */

  n = MAX(1, n);
  final anz = Ap[n + k1] - Ap[k1];

  int lsize;
  if (Lsize <= 0) {
    Lsize = -Lsize;
    Lsize = MAX(Lsize, 1.0);
    lsize = (Lsize * anz + n).toInt();
  } else {
    lsize = Lsize.toInt();
  }

  int usize = lsize;

  lsize = MAX(n + 1, lsize);
  usize = MAX(n + 1, usize);

  double maxlnz = ((n) * (n) + (n)) / 2.0;
  maxlnz = MIN(maxlnz, (INT_MAX));
  lsize = MIN(maxlnz.toInt(), lsize);
  usize = MIN(maxlnz.toInt(), usize);

  PRINTF("Welcome to klu: n $n anz $anz k1 $k1 lsize $lsize usize $usize maxlnz $maxlnz\n");

  /* ---------------------------------------------------------------------- */
  /* allocate workspace and outputs */
  /* ---------------------------------------------------------------------- */

  /* return arguments are not yet assigned */
  p_LU[block] = null;

  /* these computations are safe from int overflow */
  //W = Work ;
  //Pinv = W ;      //W += n ;
  //int Pinv_offset = n ;
  final Pinv = new Int32List(n);
  //Stack = W ;     //W += n ;
  //int Stack_offset = 2*n ;
  final Stack = new Int32List(n);
  //Flag = W ;      //W += n ;
  //int Flag_offset = 3*n ;
  final Flag = new Int32List(n);
  //Lpend = W ;     //W += n ;
  //int Lpend_offset = 4*n ;
  final Lpend = new Int32List(n);
  //Ap_pos = W ;    //W += n ;
  //int Ap_pos_offset = 5*n ;
  final Ap_pos = new Int32List(n);

  //dunits = DUNITS (Integer, lsize) + DUNITS (Double, lsize) +
  //		 DUNITS (Integer, usize) + DUNITS (Double, usize) ;
  final dunits = (lsize + lsize + usize + usize).toDouble();
  int lusize = dunits.toInt();
  final ok = INT_OVERFLOW(dunits) ? FALSE : TRUE;
  LU[0] = ok != 0 ? malloc_dbl(lusize, Common) : null;
  if (LU[0] == null) {
    /* out of memory, or problem too large */
    Common.status = KLU_OUT_OF_MEMORY;
    lusize = 0;
    return (lusize);
  }

  /* ---------------------------------------------------------------------- */
  /* factorize */
  /* ---------------------------------------------------------------------- */

  /* with pruning, and non-recursive depth-first-search */
  lusize = _kernel(n, Ap, Ai, Ax, Q, lusize, Pinv, P, LU, Udiag, Udiag_offset,
      Llen, Llen_offset, Ulen, Ulen_offset, Lip, Lip_offset, Uip, Uip_offset,
      lnz, unz, X, Stack, Flag, Ap_pos, Lpend, k1, PSinv, Rs, Offp, Offi, Offx,
      Common);

  /* ---------------------------------------------------------------------- */
  /* return LU factors, or return nothing if an error occurred */
  /* ---------------------------------------------------------------------- */

  if (Common.status < KLU_OK) {
    //LU = KLU_free (LU, lusize, sizeof (double), Common) ;
    LU[0] = null;
    lusize = 0;
  }
  p_LU[block] = LU[0];
  PRINTF(" in klu noffdiag ${Common.noffdiag}\n");
  return (lusize);
}

/**
 * Solve `Lx=b`. Assumes L is unit lower triangular and where the unit
 * diagonal entry is *not* stored. Overwrites `B` with the solution [X]. B is
 * [n]-by-[nrhs] and is stored in *row* form with row dimension [nrhs]. [nrhs]
 * must be in the range 1 to 4.
 *
 * [X] right-hand-side on input, solution to `Lx=b` on output.
 */
void lsolve(final int n, final Int32List Lip, final int Lip_offset,
            final Int32List Llen, final int Llen_offset, final Float64List LU,
            final int nrhs, final Float64List X, final int X_offset) {
  final x = new Float64List(4);
  final len = new Int32List(1);
  final Li_offset = new Int32List(1);
  final Lx_offset = new Int32List(1);

  switch (nrhs) {

    case 1:
      for (int k = 0; k < n; k++) {
        x[0] = X[X_offset + k];
        final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);
        final Lx = Li;
        /* unit diagonal of L is not stored*/
        for (int p = 0; p < len[0]; p++) {
          //MULT_SUB (X [Li [p]], Lx [p], x [0]) ;
          X[X_offset + (Li[Li_offset[0] + p]).toInt()] -= Lx[Lx_offset[0] + p] * x[0];
        }
      }
      break;

    case 2:

      for (int k = 0; k < n; k++) {
        x[0] = X[X_offset + 2 * k];
        x[1] = X[X_offset + 2 * k + 1];
        final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);
        final Lx = Li;
        for (int p = 0; p < len[0]; p++) {
          final i = (Li[Li_offset[0] + p]).toInt();
          final lik = Lx[Lx_offset[0] + p];
          //MULT_SUB (X [2*i], lik, x [0]) ;
          X[X_offset + 2 * i] -= lik * x[0];
          //MULT_SUB (X [2*i + 1], lik, x [1]) ;
          X[X_offset + 2 * i + 1] -= lik * x[1];
        }
      }
      break;

    case 3:

      for (int k = 0; k < n; k++) {
        x[0] = X[X_offset + 3 * k];
        x[1] = X[X_offset + 3 * k + 1];
        x[2] = X[X_offset + 3 * k + 2];
        final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);
        final Lx = Li;
        for (int p = 0; p < len[0]; p++) {
          final i = (Li[Li_offset[0] + p]).toInt();
          final lik = Lx[Lx_offset[0] + p];
          //MULT_SUB (X [3*i], lik, x [0]) ;
          X[X_offset + 3 * i] -= lik * x[0];
          //MULT_SUB (X [3*i + 1], lik, x [1]) ;
          X[X_offset + 3 * i + 1] -= lik * x[1];
          //MULT_SUB (X [3*i + 2], lik, x [2]) ;
          X[X_offset + 3 * i + 2] -= lik * x[2];
        }
      }
      break;

    case 4:

      for (int k = 0; k < n; k++) {
        x[0] = X[X_offset + 4 * k];
        x[1] = X[X_offset + 4 * k + 1];
        x[2] = X[X_offset + 4 * k + 2];
        x[3] = X[X_offset + 4 * k + 3];
        final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);
        final Lx = Li;
        for (int p = 0; p < len[0]; p++) {
          final i = (Li[Li_offset[0] + p]).toInt();
          final lik = Lx[Lx_offset[0] + p];
          //MULT_SUB (X [4*i], lik, x [0]) ;
          X[X_offset + 4 * i] -= lik * x[0];
          //MULT_SUB (X [4*i + 1], lik, x [1]) ;
          X[X_offset + 4 * i + 1] -= lik * x[1];
          //MULT_SUB (X [4*i + 2], lik, x [2]) ;
          X[X_offset + 4 * i + 2] -= lik * x[2];
          //MULT_SUB (X [4*i + 3], lik, x [3]) ;
          X[X_offset + 4 * i + 3] -= lik * x[3];
        }
      }
      break;

  }
}

/**
 * Solve `Ux=b`. Assumes `U` is non-unit upper triangular and where the
 * diagonal entry is *not* stored. Overwrites `B` with the solution [X]. B is
 * [n]-by-[nrhs] and is stored in *row* form with row dimension [nrhs]. [nrhs]
 * must be in the range 1 to 4.
 *
 * [X] is the right-hand-side on input, solution to `Ux=b` on output.
 */
void usolve(final int n, final Int32List Uip, final int Uip_offset,
            final Int32List Ulen, final int Ulen_offset, final Float64List LU,
            final Float64List Udiag, final int Udiag_offset, final int nrhs,
            final Float64List X, final int X_offset) {
  final x = new Float64List(4);
  final len = new Int32List(1);
  final Ui_offset = new Int32List(1);
  final Ux_offset = new Int32List(1);

  switch (nrhs) {

    case 1:

      for (int k = n - 1; k >= 0; k--) {
        final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len);
        final Ux = Ui;
        //DIV (x [0], X [k], Udiag [k]) ;
        x[0] = X[X_offset + k] / Udiag[Udiag_offset + k];
        X[X_offset + k] = x[0];
        for (int p = 0; p < len[0]; p++) {
          //MULT_SUB (X [Ui [p]], Ux [p], x [0]) ;
          X[X_offset + (Ui[Ui_offset[0] + p]).toInt()] -= Ux[Ux_offset[0] + p] * x[0];

        }
      }

      break;

    case 2:

      for (int k = n - 1; k >= 0; k--) {
        final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len);
        final Ux = Ui;
        final ukk = Udiag[Udiag_offset + k];
        //DIV (x [0], X [2*k], ukk) ;
        x[0] = X[X_offset + 2 * k] / ukk;
        //DIV (x [1], X [2*k + 1], ukk) ;
        x[1] = X[X_offset + 2 * k + 1] / ukk;

        X[X_offset + 2 * k] = x[0];
        X[X_offset + 2 * k + 1] = x[1];
        for (int p = 0; p < len[0]; p++) {
          final i = Ui[Ui_offset[0] + p].toInt();
          final uik = Ux[Ux_offset[0] + p];
          //MULT_SUB (X [2*i], uik, x [0]) ;
          X[X_offset + 2 * i] -= uik * x[0];
          //MULT_SUB (X [2*i + 1], uik, x [1]) ;
          X[X_offset + 2 * i + 1] -= uik * x[1];
        }
      }

      break;

    case 3:

      for (int k = n - 1; k >= 0; k--) {
        final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len);
        final Ux = Ui;
        final ukk = Udiag[Udiag_offset + k];

        //DIV (x [0], X [3*k], ukk) ;
        x[0] = X[X_offset + 3 * k] / ukk;
        //DIV (x [1], X [3*k + 1], ukk) ;
        x[1] = X[X_offset + 3 * k + 1] / ukk;
        //DIV (x [2], X [3*k + 2], ukk) ;
        x[2] = X[X_offset + 3 * k + 2] / ukk;

        X[3 * k] = x[0];
        X[3 * k + 1] = x[1];
        X[3 * k + 2] = x[2];
        for (int p = 0; p < len[0]; p++) {
          final i = Ui[Ui_offset[0] + p].toInt();
          final uik = Ux[Ux_offset[0] + p];
          //MULT_SUB (X [3*i], uik, x [0]) ;
          X[X_offset + 3 * i] -= uik * x[0];
          //MULT_SUB (X [3*i + 1], uik, x [1]) ;
          X[X_offset + 3 * i + 1] -= uik * x[1];
          //MULT_SUB (X [3*i + 2], uik, x [2]) ;
          X[X_offset + 3 * i + 2] -= uik * x[2];
        }
      }

      break;

    case 4:

      for (int k = n - 1; k >= 0; k--) {
        final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len);
        final Ux = Ui;
        final ukk = Udiag[Udiag_offset + k];

        //DIV (x [0], X [4*k], ukk) ;
        x[0] = X[X_offset + 4 * k] / ukk;
        //DIV (x [1], X [4*k + 1], ukk) ;
        x[1] = X[X_offset + 4 * k + 1] / ukk;
        //DIV (x [2], X [4*k + 2], ukk) ;
        x[2] = X[X_offset + 4 * k + 2] / ukk;
        //DIV (x [3], X [4*k + 3], ukk) ;
        x[3] = X[X_offset + 4 * k + 3] / ukk;

        X[X_offset + 4 * k] = x[0];
        X[X_offset + 4 * k + 1] = x[1];
        X[X_offset + 4 * k + 2] = x[2];
        X[X_offset + 4 * k + 3] = x[3];
        for (int p = 0; p < len[0]; p++) {
          final i = Ui[Ui_offset[0] + p].toInt();
          final uik = Ux[Ux_offset[0] + p];

          //MULT_SUB (X [4*i], uik, x [0]) ;
          X[X_offset + 4 * i] -= uik * x[0];
          //MULT_SUB (X [4*i + 1], uik, x [1]) ;
          X[X_offset + 4 * i + 1] -= uik * x[1];
          //MULT_SUB (X [4*i + 2], uik, x [2]) ;
          X[X_offset + 4 * i + 2] -= uik * x[2];
          //MULT_SUB (X [4*i + 3], uik, x [3]) ;
          X[X_offset + 4 * i + 3] -= uik * x[3];
        }
      }

      break;

  }
}

/**
 * Solve `L'x=b`. Assumes `L` is unit lower triangular and where the unit
 * diagonal entry is *not* stored. Overwrites `B` with the solution [X]. B is
 * [n]-by-[nrhs] and is stored in *row* form with row dimension [nrhs]. [nrhs]
 * must in the range 1 to 4.
 *
 * [X] right-hand-side on input, solution to `L'x=b` on output.
 */
void ltsolve(final int n, final Int32List Lip, final int Lip_offset,
             final Int32List Llen, final int Llen_offset, final Float64List LU,
             final int nrhs, final Float64List X, final int X_offset) {
  final x = new Float64List(4);
  final len = new Int32List(1);
  final Li_offset = new Int32List(1);
  final Lx_offset = new Int32List(1);

  switch (nrhs) {

    case 1:

      for (int k = n - 1; k >= 0; k--) {
        final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);
        final Lx = Li;
        x[0] = X[X_offset + k];
        for (int p = 0; p < len[0]; p++) {
          {
            //MULT_SUB (x [0], Lx [p], X [Li [p]]) ;
            x[0] -= Lx[Lx_offset[0] + p] * X[X_offset + Li[Li_offset[0] + p].toInt()];
          }
        }
        X[X_offset + k] = x[0];
      }
      break;

    case 2:

      for (int k = n - 1; k >= 0; k--) {
        x[0] = X[X_offset + 2 * k];
        x[1] = X[X_offset + 2 * k + 1];
        final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);
        final Lx = Li;
        for (int p = 0; p < len[0]; p++) {
          final i = Li[Li_offset[0] + p].toInt();
          //{
            final lik = Lx[Lx_offset[0] + p];
          //}
          //MULT_SUB (x [0], lik, X [2*i]) ;
          x[0] -= lik * X[X_offset + 2 * i];
          //MULT_SUB (x [1], lik, X [2*i + 1]) ;
          x[1] -= lik * X[X_offset + 2 * i + 1];
        }
        X[X_offset + 2 * k] = x[0];
        X[X_offset + 2 * k + 1] = x[1];
      }
      break;

    case 3:

      for (int k = n - 1; k >= 0; k--) {
        x[0] = X[X_offset + 3 * k];
        x[1] = X[X_offset + 3 * k + 1];
        x[2] = X[X_offset + 3 * k + 2];
        final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);
        final Lx = Li;
        for (int p = 0; p < len[0]; p++) {
          final i = Li[Li_offset[0] + p].toInt();
          //{
            final lik = Lx[Lx_offset[0] + p];
          //}
          //MULT_SUB (x [0], lik, X [3*i]) ;
          x[0] -= lik * X[X_offset + 3 * i];
          //MULT_SUB (x [1], lik, X [3*i + 1]) ;
          x[1] -= lik * X[X_offset + 3 * i + 1];
          //MULT_SUB (x [2], lik, X [3*i + 2]) ;
          x[2] -= lik * X[X_offset + 3 * i + 2];
        }
        X[X_offset + 3 * k] = x[0];
        X[X_offset + 3 * k + 1] = x[1];
        X[X_offset + 3 * k + 2] = x[2];
      }
      break;

    case 4:

      for (int k = n - 1; k >= 0; k--) {
        x[0] = X[X_offset + 4 * k];
        x[1] = X[X_offset + 4 * k + 1];
        x[2] = X[X_offset + 4 * k + 2];
        x[3] = X[X_offset + 4 * k + 3];
        final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len);
        final Lx = Li;
        for (int p = 0; p < len[0]; p++) {
          final i = Li[Li_offset[0] + p].toInt();
          //{
            final lik = Lx[Lx_offset[0] + p];
          //}
          //MULT_SUB (x [0], lik, X [4*i]) ;
          x[0] -= lik * X[X_offset + 4 * i];
          //MULT_SUB (x [1], lik, X [4*i + 1]) ;
          x[1] -= lik * X[X_offset + 4 * i + 1];
          //MULT_SUB (x [2], lik, X [4*i + 2]) ;
          x[2] -= lik * X[X_offset + 4 * i + 2];
          //MULT_SUB (x [3], lik, X [4*i + 3]) ;
          x[3] -= lik * X[X_offset + 4 * i + 3];
        }
        X[X_offset + 4 * k] = x[0];
        X[X_offset + 4 * k + 1] = x[1];
        X[X_offset + 4 * k + 2] = x[2];
        X[X_offset + 4 * k + 3] = x[3];
      }
      break;
  }
}

/**
 * Solve `U'x=b`. Assumes `U` is non-unit upper triangular and where the
 * diagonal entry is stored (and appears last in each column of `U`).
 * Overwrites `B` with the solution [X]. `B` is [n]-by-[nrhs] and is stored
 * in *row* form with row dimension [nrhs]. [nrhs] must be in the range 1 to 4.
 *
 * [X] right-hand-side on input, solution to `Ux=b` on output.
 */
void utsolve(final int n, final Int32List Uip, final int Uip_offset,
             final Int32List Ulen, final int Ulen_offset, final Float64List LU,
             final Float64List Udiag, final int Udiag_offset, final int nrhs,
             final Float64List X, final int X_offset) {
  final x = new Float64List(4);
  final len = new Int32List(1);
  final Ui_offset = new Int32List(1);
  final Ux_offset = new Int32List(1);

  switch (nrhs) {

    case 1:

      for (int k = 0; k < n; k++) {
        final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len);
        final Ux = Ui;
        x[0] = X[X_offset + k];
        for (int p = 0; p < len[0]; p++) {
          {
            //MULT_SUB (x [0], Ux [p], X [Ui [p]]) ;
            x[0] -= Ux[Ux_offset[0] + p] * X[X_offset + Ui[Ui_offset[0] + p].toInt()];
          }
        }
        //{
          final ukk = Udiag[k];
        //}
        //DIV (X [k], x [0], ukk) ;
        X[X_offset + k] = x[0] / ukk;
      }
      break;

    case 2:

      for (int k = 0; k < n; k++) {
        final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len);
        final Ux = Ui;
        x[0] = X[X_offset + 2 * k];
        x[1] = X[X_offset + 2 * k + 1];
        for (int p = 0; p < len[0]; p++) {
          final i = Ui[Ui_offset[0] + p].toInt();
          //{
          final uik = Ux[Ux_offset[0] + p];
          //}
          //MULT_SUB (x [0], uik, X [2*i]) ;
          x[0] -= uik * X[X_offset + 2 * i];
          //MULT_SUB (x [1], uik, X [2*i + 1]) ;
          x[1] -= uik * X[X_offset + 2 * i + 1];
        }
        //{
          final ukk = Udiag[k];
        //}
        //DIV (X [2*k], x [0], ukk) ;
        X[X_offset + 2 * k] = x[0] / ukk;
        //DIV (X [2*k + 1], x [1], ukk) ;
        X[X_offset + 2 * k + 1] = x[1] / ukk;
      }
      break;

    case 3:

      for (int k = 0; k < n; k++) {
        final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len);
        final Ux = Ui;
        x[0] = X[X_offset + 3 * k];
        x[1] = X[X_offset + 3 * k + 1];
        x[2] = X[X_offset + 3 * k + 2];
        for (int p = 0; p < len[0]; p++) {
          final i = Ui[Ui_offset[0] + p].toInt();
          //{
            final uik = Ux[Ux_offset[0] + p];
          //}
          //MULT_SUB (x [0], uik, X [3*i]) ;
          x[0] -= uik * X[X_offset + 3 * i];
          //MULT_SUB (x [1], uik, X [3*i + 1]) ;
          x[1] -= uik * X[X_offset + 3 * i + 1];
          //MULT_SUB (x [2], uik, X [3*i + 2]) ;
          x[2] -= uik * X[X_offset + 3 * i + 2];
        }
        //{
          final ukk = Udiag[k];
        //}
        //DIV (X [3*k], x [0], ukk) ;
        X[X_offset + 3 * k] = x[0] / ukk;
        //DIV (X [3*k + 1], x [1], ukk) ;
        X[X_offset + 3 * k + 1] = x[1] / ukk;
        //DIV (X [3*k + 2], x [2], ukk) ;
        X[X_offset + 3 * k + 2] = x[2] / ukk;
      }
      break;

    case 4:

      for (int k = 0; k < n; k++) {
        final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len);
        final Ux = Ui;
        x[0] = X[X_offset + 4 * k];
        x[1] = X[X_offset + 4 * k + 1];
        x[2] = X[X_offset + 4 * k + 2];
        x[3] = X[X_offset + 4 * k + 3];
        for (int p = 0; p < len[0]; p++) {
          final i = Ui[Ui_offset[0] + p].toInt();
          //{
            final uik = Ux[Ux_offset[0] + p];
          //}
          //MULT_SUB (x [0], uik, X [4*i]) ;
          x[0] -= uik * X[X_offset + 4 * i];
          //MULT_SUB (x [1], uik, X [4*i + 1]) ;
          x[1] -= uik * X[X_offset + 4 * i + 1];
          //MULT_SUB (x [2], uik, X [4*i + 2]) ;
          x[2] -= uik * X[X_offset + 4 * i + 2];
          //MULT_SUB (x [3], uik, X [4*i + 3]) ;
          x[3] -= uik * X[X_offset + 4 * i + 3];
        }
        //{
          final ukk = Udiag[k];
        //}
        //DIV (X [4*k], x [0], ukk) ;
        X[X_offset + 4 * k] = x[0] / ukk;
        //DIV (X [4*k + 1], x [1], ukk) ;
        X[X_offset + 4 * k + 1] = x[1] / ukk;
        //DIV (X [4*k + 2], x [2], ukk) ;
        X[X_offset + 4 * k + 2] = x[2] / ukk;
        //DIV (X [4*k + 3], x [3], ukk) ;
        X[X_offset + 4 * k + 3] = x[3] / ukk;
      }
      break;
  }
}
