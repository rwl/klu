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
 * Solve `Ax=b` using the symbolic and numeric objects from [analyze]
 * (or [analyze_given]) and [factor]. Note that no iterative refinement is
 * performed. Uses [Numeric.Xwork] as workspace (undefined on input and
 * output), of size `4n` double's (note that columns 2 to 4 of Xwork overlap
 * with [Numeric.Iwork]).
 *
 * [d] is the leading dimension of [B]. [nrhs] s the number of
 * right-hand-sides. [B] is the right-hand-side on input, overwritten with
 * solution to `Ax=b` on output. Size `n*nrhs`, in column-oriented form, with
 * leading dimension d.
 */
int solve(final KLU_symbolic Symbolic, final KLU_numeric Numeric, final int d,
          final int nrhs, final Float64List B, int B_offset,
          final KLU_common Common) {
  final x = new Float64List(4);

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  if (Common == null) {
    return (FALSE);
  }
  if (Numeric == null || Symbolic == null || d < Symbolic.n || nrhs < 0 || B == null) {
    Common.status = KLU_INVALID;
    return (FALSE);
  }
  Common.status = KLU_OK;

  /* ---------------------------------------------------------------------- */
  /* get the contents of the Symbolic object */
  /* ---------------------------------------------------------------------- */

  final Bz = B;
  final n = Symbolic.n;
  final nblocks = Symbolic.nblocks;
  final Q = Symbolic.Q;
  final R = Symbolic.R;

  /* ---------------------------------------------------------------------- */
  /* get the contents of the Numeric object */
  /* ---------------------------------------------------------------------- */

  ASSERT(nblocks == Numeric.nblocks);
  final Pnum = Numeric.Pnum;
  final Offp = Numeric.Offp;
  final Offi = Numeric.Offi;
  final Offx = Numeric.Offx;

  final Lip = Numeric.Lip;
  final Llen = Numeric.Llen;
  final Uip = Numeric.Uip;
  final Ulen = Numeric.Ulen;
  final LUbx = Numeric.LUbx;
  final Udiag = Numeric.Udiag;

  final Rs = Numeric.Rs;
  final X = Numeric.Xwork;

  if (!NDEBUG) ASSERT_INT(_valid(n, Offp, Offi, Offx));

  /* ---------------------------------------------------------------------- */
  /* solve in chunks of 4 columns at a time */
  /* ---------------------------------------------------------------------- */

  for (int chunk = 0; chunk < nrhs; chunk += 4) {

    /* ------------------------------------------------------------------ */
    /* get the size of the current chunk */
    /* ------------------------------------------------------------------ */

    final nr = MIN(nrhs - chunk, 4);

    /* ------------------------------------------------------------------ */
    /* scale and permute the right hand side, X = P*(R\B) */
    /* ------------------------------------------------------------------ */

    if (Rs == null) {

      /* no scaling */
      switch (nr) {

        case 1:

          for (int k = 0; k < n; k++) {
            X[k] = Bz[B_offset + Pnum[k]];
          }
          break;

        case 2:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            X[2 * k] = Bz[B_offset + i];
            X[2 * k + 1] = Bz[B_offset + i + d];
          }
          break;

        case 3:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            X[3 * k] = Bz[B_offset + i];
            X[3 * k + 1] = Bz[B_offset + i + d];
            X[3 * k + 2] = Bz[B_offset + i + d * 2];
          }
          break;

        case 4:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            X[4 * k] = Bz[B_offset + i];
            X[4 * k + 1] = Bz[B_offset + i + d];
            X[4 * k + 2] = Bz[B_offset + i + d * 2];
            X[4 * k + 3] = Bz[B_offset + i + d * 3];
          }
          break;
      }

    } else {

      switch (nr) {

        case 1:

          for (int k = 0; k < n; k++) {
            //SCALE_DIV_ASSIGN (X [k], Bz  [Pnum [k]], Rs [k]) ;
            X[k] = Bz[B_offset + Pnum[k]] / Rs[k];
          }
          break;

        case 2:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            final rs = Rs[k];
            //SCALE_DIV_ASSIGN (X [2*k], Bz [i], rs) ;
            X[2 * k] = Bz[B_offset + i] / rs;
            //SCALE_DIV_ASSIGN (X [2*k + 1], Bz [i + d], rs) ;
            X[2 * k + 1] = Bz[B_offset + i + d] / rs;
          }
          break;

        case 3:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            final rs = Rs[k];
            //SCALE_DIV_ASSIGN (X [3*k], Bz [i], rs) ;
            X[3 * k] = Bz[B_offset + i] / rs;
            //SCALE_DIV_ASSIGN (X [3*k + 1], Bz [i + d], rs) ;
            X[3 * k + 1] = Bz[B_offset + i + d] / rs;
            //SCALE_DIV_ASSIGN (X [3*k + 2], Bz [i + d*2], rs) ;
            X[3 * k + 2] = Bz[B_offset + i + d * 2] / rs;
          }
          break;

        case 4:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            final rs = Rs[k];
            //SCALE_DIV_ASSIGN (X [4*k], Bz [i], rs) ;
            X[4 * k] = Bz[B_offset + i] / rs;
            //SCALE_DIV_ASSIGN (X [4*k + 1], Bz [i + d], rs) ;
            X[4 * k + 1] = Bz[B_offset + i + d] / rs;
            //SCALE_DIV_ASSIGN (X [4*k + 2], Bz [i + d*2], rs) ;
            X[4 * k + 2] = Bz[B_offset + i + d * 2] / rs;
            //SCALE_DIV_ASSIGN (X [4*k + 3], Bz [i + d*3], rs) ;
            X[4 * k + 3] = Bz[B_offset + i + d * 3] / rs;
          }
          break;
      }
    }

    /* ------------------------------------------------------------------ */
    /* solve X = (L*U + Off)\X */
    /* ------------------------------------------------------------------ */

    for (int block = nblocks - 1; block >= 0; block--) {

      /* -------------------------------------------------------------- */
      /* the block of size nk is from rows/columns k1 to k2-1 */
      /* -------------------------------------------------------------- */

      final k1 = R[block];
      final k2 = R[block + 1];
      final nk = k2 - k1;
      PRINTF("solve $block, k1 $k1 k2-1 ${k2-1} nk $nk\n");

      /* solve the block system */
      if (nk == 1) {
        final s = Udiag[k1];
        switch (nr) {

          case 1:
            //DIV (X [k1], X [k1], s) ;
            X[k1] = X[k1] / s;
            break;

          case 2:
            //DIV (X [2*k1], X [2*k1], s) ;
            X[2 * k1] = X[2 * k1] / s;
            //DIV (X [2*k1 + 1], X [2*k1 + 1], s) ;
            X[2 * k1 + 1] = X[2 * k1 + 1] / s;
            break;

          case 3:
            //DIV (X [3*k1], X [3*k1], s) ;
            X[3 * k1] = X[3 * k1] / s;
            //DIV (X [3*k1 + 1], X [3*k1 + 1], s) ;
            X[3 * k1 + 1] = X[3 * k1 + 1] / s;
            //DIV (X [3*k1 + 2], X [3*k1 + 2], s) ;
            X[3 * k1 + 2] = X[3 * k1 + 2] / s;
            break;

          case 4:
            //DIV (X [4*k1], X [4*k1], s) ;
            X[4 * k1] = X[4 * k1] / s;
            //DIV (X [4*k1 + 1], X [4*k1 + 1], s) ;
            X[4 * k1 + 1] = X[4 * k1 + 1] / s;
            //DIV (X [4*k1 + 2], X [4*k1 + 2], s) ;
            X[4 * k1 + 2] = X[4 * k1 + 2] / s;
            //DIV (X [4*k1 + 3], X [4*k1 + 3], s) ;
            X[4 * k1 + 3] = X[4 * k1 + 3] / s;
            break;

        }
      } else {
        lsolve(nk, Lip, k1, Llen, k1, LUbx[block], nr, X, nr * k1);
        usolve(nk, Uip, k1, Ulen, k1, LUbx[block], Udiag, k1, nr, X, nr * k1);
      }

      /* -------------------------------------------------------------- */
      /* block back-substitution for the off-diagonal-block entries */
      /* -------------------------------------------------------------- */

      if (block > 0) {
        switch (nr) {

          case 1:

            for (int k = k1; k < k2; k++) {
              final pend = Offp[k + 1];
              x[0] = X[k];
              for (int p = Offp[k]; p < pend; p++) {
                //MULT_SUB (X [Offi [p]], Offx [p], x [0]) ;
                X[Offi[p]] -= Offx[p] * x[0];
              }
            }
            break;

          case 2:

            for (int k = k1; k < k2; k++) {
              final pend = Offp[k + 1];
              x[0] = X[2 * k];
              x[1] = X[2 * k + 1];
              for (int p = Offp[k]; p < pend; p++) {
                final i = Offi[p];
                final offik = Offx[p];
                //MULT_SUB (X [2*i], offik, x [0]) ;
                X[2 * i] -= offik * x[0];
                //MULT_SUB (X [2*i + 1], offik, x [1]) ;
                X[2 * i + 1] -= offik * x[1];
              }
            }
            break;

          case 3:

            for (int k = k1; k < k2; k++) {
              final pend = Offp[k + 1];
              x[0] = X[3 * k];
              x[1] = X[3 * k + 1];
              x[2] = X[3 * k + 2];
              for (int p = Offp[k]; p < pend; p++) {
                final i = Offi[p];
                final offik = Offx[p];
                //MULT_SUB (X [3*i], offik, x [0]) ;
                X[3 * i] -= offik * x[0];
                //MULT_SUB (X [3*i + 1], offik, x [1]) ;
                X[3 * i + 1] -= offik * x[1];
                //MULT_SUB (X [3*i + 2], offik, x [2]) ;
                X[3 * i + 2] -= offik * x[2];
              }
            }
            break;

          case 4:

            for (int k = k1; k < k2; k++) {
              final pend = Offp[k + 1];
              x[0] = X[4 * k];
              x[1] = X[4 * k + 1];
              x[2] = X[4 * k + 2];
              x[3] = X[4 * k + 3];
              for (int p = Offp[k]; p < pend; p++) {
                final i = Offi[p];
                final offik = Offx[p];
                //MULT_SUB (X [4*i], offik, x [0]) ;
                X[4 * i] -= offik * x[0];
                //MULT_SUB (X [4*i + 1], offik, x [1]) ;
                X[4 * i + 1] -= offik * x[1];
                //MULT_SUB (X [4*i + 2], offik, x [2]) ;
                X[4 * i + 2] -= offik * x[2];
                //MULT_SUB (X [4*i + 3], offik, x [3]) ;
                X[4 * i + 3] -= offik * x[3];
              }
            }
            break;
        }
      }
    }

    /* ------------------------------------------------------------------ */
    /* permute the result, Bz  = Q*X */
    /* ------------------------------------------------------------------ */

    switch (nr) {

      case 1:

        for (int k = 0; k < n; k++) {
          Bz[B_offset + Q[k]] = X[k];
        }
        break;

      case 2:

        for (int k = 0; k < n; k++) {
          final i = Q[k];
          Bz[B_offset + i] = X[2 * k];
          Bz[B_offset + i + d] = X[2 * k + 1];
        }
        break;

      case 3:

        for (int k = 0; k < n; k++) {
          final i = Q[k];
          Bz[B_offset + i] = X[3 * k];
          Bz[B_offset + i + d] = X[3 * k + 1];
          Bz[B_offset + i + d * 2] = X[3 * k + 2];
        }
        break;

      case 4:

        for (int k = 0; k < n; k++) {
          final i = Q[k];
          Bz[B_offset + i] = X[4 * k];
          Bz[B_offset + i + d] = X[4 * k + 1];
          Bz[B_offset + i + d * 2] = X[4 * k + 2];
          Bz[B_offset + i + d * 3] = X[4 * k + 3];
        }
        break;
    }

    /* ------------------------------------------------------------------ */
    /* go to the next chunk of B */
    /* ------------------------------------------------------------------ */

    B_offset += d * 4;
  }
  return (TRUE);
}
