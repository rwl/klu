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
 * Solve `A'x=b` using the symbolic and numeric objects from [analyze]
 * (or [analyze_given]) and [factor]. Note that no iterative refinement is
 * performed. Uses [Numeric.Xwork] as workspace (undefined on input and
 * output), of size `4n` double's (note that columns 2 to 4 of `Xwork`
 * overlap with [Numeric.Iwork]).
 *
 * [d] is the leading dimension of [B]. [nrhs] is the number of
 * right-hand-sides. [B] is the right-hand-side on input, overwritten with
 * solution to `Ax=b` on output. Size `n*nrhs`, in column-oriented form,
 * with leading dimension d.
 */
int tsolve(final KLU_symbolic Symbolic, final KLU_numeric Numeric,
           final int d, final int nrhs, final Float64List B,
           int B_offset, final KLU_common Common) {
  final x = new Float64List(4);
//  double offik, s;
//  double rs;
//  Float64List Rs;
//  Float64List Offx, X, Bz, Udiag;
//  Int32List Q, R, Pnum, Offp, Offi, Lip, Uip, Llen, Ulen;
//  List<Float64List> LUbx;
//  int k1, k2, nk, k, block, pend, n, p, nblocks, chunk, nr, i;

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
    /* permute the right hand side, X = Q'*B */
    /* ------------------------------------------------------------------ */

    switch (nr) {

      case 1:

        for (int k = 0; k < n; k++) {
          X[k] = Bz[B_offset + Q[k]];
        }
        break;

      case 2:

        for (int k = 0; k < n; k++) {
          final i = Q[k];
          X[2 * k] = Bz[B_offset + i];
          X[2 * k + 1] = Bz[B_offset + i + d];
        }
        break;

      case 3:

        for (int k = 0; k < n; k++) {
          final i = Q[k];
          X[3 * k] = Bz[B_offset + i];
          X[3 * k + 1] = Bz[B_offset + i + d];
          X[3 * k + 2] = Bz[B_offset + i + d * 2];
        }
        break;

      case 4:

        for (int k = 0; k < n; k++) {
          final i = Q[k];
          X[4 * k] = Bz[B_offset + i];
          X[4 * k + 1] = Bz[B_offset + i + d];
          X[4 * k + 2] = Bz[B_offset + i + d * 2];
          X[4 * k + 3] = Bz[B_offset + i + d * 3];
        }
        break;

    }

    /* ------------------------------------------------------------------ */
    /* solve X = (L*U + Off)'\X */
    /* ------------------------------------------------------------------ */

    for (int block = 0; block < nblocks; block++) {

      /* -------------------------------------------------------------- */
      /* the block of size nk is from rows/columns k1 to k2-1 */
      /* -------------------------------------------------------------- */

      final k1 = R[block];
      final k2 = R[block + 1];
      final nk = k2 - k1;
      PRINTF("tsolve $block, k1 $k1 k2-1 ${k2-1} nk $nk\n");

      /* -------------------------------------------------------------- */
      /* block back-substitution for the off-diagonal-block entries */
      /* -------------------------------------------------------------- */

      if (block > 0) {
        switch (nr) {

          case 1:

            for (int k = k1; k < k2; k++) {
              final pend = Offp[k + 1];
              for (int p = Offp[k]; p < pend; p++) {
                {
                  //MULT_SUB (X [k], Offx [p], X [Offi [p]]) ;
                  X[k] -= Offx[p] * X[Offi[p]];
                }
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
                //{
                  final offik = Offx[p];
                //}
                //MULT_SUB (x [0], offik, X [2*i]) ;
                x[0] -= offik * X[2 * i];
                //MULT_SUB (x [1], offik, X [2*i + 1]) ;
                x[1] -= offik * X[2 * i + 1];
              }
              X[2 * k] = x[0];
              X[2 * k + 1] = x[1];
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
                //{
                  final offik = Offx[p];
                //}
                //MULT_SUB (x [0], offik, X [3*i]) ;
                x[0] -= offik * X[3 * i];
                //MULT_SUB (x [1], offik, X [3*i + 1]) ;
                x[1] -= offik * X[3 * i + 1];
                //MULT_SUB (x [2], offik, X [3*i + 2]) ;
                x[2] -= offik * X[3 * i + 2];
              }
              X[3 * k] = x[0];
              X[3 * k + 1] = x[1];
              X[3 * k + 2] = x[2];
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
                //{
                  final offik = Offx[p];
                //}
                //MULT_SUB (x [0], offik, X [4*i]) ;
                x[0] -= offik * X[4 * i];
                //MULT_SUB (x [1], offik, X [4*i + 1]) ;
                x[1] -= offik * X[4 * i + 1];
                //MULT_SUB (x [2], offik, X [4*i + 2]) ;
                x[2] -= offik * X[4 * i + 2];
                //MULT_SUB (x [3], offik, X [4*i + 3]) ;
                x[3] -= offik * X[4 * i + 3];
              }
              X[4 * k] = x[0];
              X[4 * k + 1] = x[1];
              X[4 * k + 2] = x[2];
              X[4 * k + 3] = x[3];
            }
            break;
        }
      }

      /* -------------------------------------------------------------- */
      /* solve the block system */
      /* -------------------------------------------------------------- */

      if (nk == 1) {
        //{
          final s = Udiag[k1];
        //}
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
        utsolve(nk, Uip, k1, Ulen, k1, LUbx[block], Udiag, k1, nr, X, nr * k1);
        ltsolve(nk, Lip, k1, Llen, k1, LUbx[block], nr, X, nr * k1);
      }
    }

    /* ------------------------------------------------------------------ */
    /* scale and permute the result, Bz  = P'(R\X) */
    /* ------------------------------------------------------------------ */

    if (Rs == null) {

      /* no scaling */
      switch (nr) {

        case 1:

          for (int k = 0; k < n; k++) {
            Bz[B_offset + Pnum[k]] = X[k];
          }
          break;

        case 2:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            Bz[B_offset + i] = X[2 * k];
            Bz[B_offset + i + d] = X[2 * k + 1];
          }
          break;

        case 3:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            Bz[B_offset + i] = X[3 * k];
            Bz[B_offset + i + d] = X[3 * k + 1];
            Bz[B_offset + i + d * 2] = X[3 * k + 2];
          }
          break;

        case 4:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            Bz[B_offset + i] = X[4 * k];
            Bz[B_offset + i + d] = X[4 * k + 1];
            Bz[B_offset + i + d * 2] = X[4 * k + 2];
            Bz[B_offset + i + d * 3] = X[4 * k + 3];
          }
          break;
      }

    } else {

      switch (nr) {

        case 1:

          for (int k = 0; k < n; k++) {
            //SCALE_DIV_ASSIGN (Bz [Pnum [k]], X [k], Rs [k]) ;
            Bz[B_offset + Pnum[k]] = X[k] / Rs[k];
          }
          break;

        case 2:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            final rs = Rs[k];
            //SCALE_DIV_ASSIGN (Bz [i], X [2*k], rs) ;
            Bz[B_offset + i] = X[2 * k] / rs;
            //SCALE_DIV_ASSIGN (Bz [i + d], X [2*k + 1], rs) ;
            Bz[B_offset + i + d] = X[2 * k + 1] / rs;
          }
          break;

        case 3:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            final rs = Rs[k];
            //SCALE_DIV_ASSIGN (Bz [i], X [3*k], rs) ;
            Bz[B_offset + i] = X[3 * k] / rs;
            //SCALE_DIV_ASSIGN (Bz [i + d], X [3*k + 1], rs) ;
            Bz[B_offset + i + d] = X[3 * k + 1] / rs;
            //SCALE_DIV_ASSIGN (Bz [i + d*2], X [3*k + 2], rs) ;
            Bz[B_offset + i + d * 2] = X[3 * k + 2] / rs;
          }
          break;

        case 4:

          for (int k = 0; k < n; k++) {
            final i = Pnum[k];
            final rs = Rs[k];
            //SCALE_DIV_ASSIGN (Bz [i], X [4*k], rs) ;
            Bz[B_offset + i] = X[4 * k] / rs;
            //SCALE_DIV_ASSIGN (Bz [i + d], X [4*k + 1], rs) ;
            Bz[B_offset + i + d] = X[4 * k + 1] / rs;
            //SCALE_DIV_ASSIGN (Bz [i + d*2], X [4*k + 2], rs) ;
            Bz[B_offset + i + d * 2] = X[4 * k + 2] / rs;
            //SCALE_DIV_ASSIGN (Bz [i + d*3], X [4*k + 3], rs) ;
            Bz[B_offset + i + d * 3] = X[4 * k + 3] / rs;
          }
          break;
      }
    }

    /* ------------------------------------------------------------------ */
    /* go to the next chunk of B */
    /* ------------------------------------------------------------------ */

    B_offset += d * 4;
  }
  return (TRUE);
}
