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
 * Factor the matrix, after ordering and analyzing it with KLU_analyze, and
 * factoring it once with KLU_factor. This routine cannot do any numerical
 * pivoting. The pattern of the input matrix (Ap, Ai) must be identical to
 * the pattern given to KLU_factor.
 */
part of edu.ufl.cise.klu.tdouble;

/**
 * Returns `true` if successful, `false` otherwise.
 */
int refactor(final Int32List Ap, final Int32List Ai, final Float64List Ax,
             final KLU_symbolic Symbolic, final KLU_numeric Numeric,
             final KLU_common Common) {
  final ulen = new Int32List(1);
  final Ui_offset = new Int32List(1);
  final Ux_offset = new Int32List(1);
  final llen = new Int32List(1);
  final Li_offset = new Int32List(1);
  final Lx_offset = new Int32List(1);

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  if (Common == null) {
    return (FALSE);
  }
  Common.status = KLU_OK;

  if (Numeric == null) {
    /* invalid Numeric object */
    Common.status = KLU_INVALID;
    return (FALSE);
  }

  Common.numerical_rank = EMPTY;
  Common.singular_col = EMPTY;

  final Az = new Float64List.fromList(Ax);

  /* ---------------------------------------------------------------------- */
  /* get the contents of the Symbolic object */
  /* ---------------------------------------------------------------------- */

  final n = Symbolic.n;
  final Q = Symbolic.Q;
  final R = Symbolic.R;
  final nblocks = Symbolic.nblocks;
  final maxblock = Symbolic.maxblock;

  /* ---------------------------------------------------------------------- */
  /* get the contents of the Numeric object */
  /* ---------------------------------------------------------------------- */

  final Pnum = Numeric.Pnum;
  final Offp = Numeric.Offp;
  final Offi = Numeric.Offi;
  final Offx = Numeric.Offx;

  final LUbx = Numeric.LUbx;

  final scale = Common.scale;
  if (scale > 0) {
    /* factorization was not scaled, but refactorization is scaled */
    if (Numeric.Rs == null) {
      Numeric.Rs = malloc_dbl(n, Common);
      if (Common.status < KLU_OK) {
        Common.status = KLU_OUT_OF_MEMORY;
        return (FALSE);
      }
    }
  } else {
    /* no scaling for refactorization; ensure Numeric.Rs is freed.  This
     * does nothing if Numeric.Rs is already null. */
    //Numeric.Rs = KLU_free (Numeric.Rs, n, sizeof (double), Common) ;
    Numeric.Rs = null;
  }
  final Rs = Numeric.Rs;

  final Pinv = Numeric.Pinv;
  final X = Numeric.Xwork;
  Common.nrealloc = 0;
  final Udiag = Numeric.Udiag;
  final nzoff = Symbolic.nzoff;

  /* ---------------------------------------------------------------------- */
  /* check the input matrix compute the row scale factors, Rs */
  /* ---------------------------------------------------------------------- */

  /* do no scale, or check the input matrix, if scale < 0 */
  if (scale >= 0) {
    /* check for out-of-range indices, but do not check for duplicates */
    if (klu_scale(scale, n, Ap, Ai, Ax, Rs, null, Common) == 0) {
      return (FALSE);
    }
  }

  /* ---------------------------------------------------------------------- */
  /* clear workspace X */
  /* ---------------------------------------------------------------------- */

  for (int k = 0; k < maxblock; k++) {
    /* X [k] = 0 ; */
    CLEAR(X, k);
  }

  int poff = 0;

  /* ---------------------------------------------------------------------- */
  /* factor each block */
  /* ---------------------------------------------------------------------- */

  if (scale <= 0) {

    /* ------------------------------------------------------------------ */
    /* no scaling */
    /* ------------------------------------------------------------------ */

    for (int block = 0; block < nblocks; block++) {

      /* -------------------------------------------------------------- */
      /* the block is from rows/columns k1 to k2-1 */
      /* -------------------------------------------------------------- */

      final k1 = R[block];
      final k2 = R[block + 1];
      final nk = k2 - k1;

      if (nk == 1) {

        /* ---------------------------------------------------------- */
        /* singleton case */
        /* ---------------------------------------------------------- */

        final oldcol = Q[k1];
        final pend = Ap[oldcol + 1];
        double s = 0.0; //CLEAR (s) ;
        for (int p = Ap[oldcol]; p < pend; p++) {
          final newrow = Pinv[Ai[p]] - k1;
          if (newrow < 0 && poff < nzoff) {
            /* entry in off-diagonal block */
            Offx[poff] = Az[p];
            poff++;
          } else {
            /* singleton */
            s = Az[p];
          }
        }
        Udiag[k1] = s;

      } else {

        /* ---------------------------------------------------------- */
        /* construct and factor the kth block */
        /* ---------------------------------------------------------- */

        final Lip = Numeric.Lip;
        final Lip_offset = k1;
        final Llen = Numeric.Llen;
        final Llen_offset = k1;
        final Uip = Numeric.Uip;
        final Uip_offset = k1;
        final Ulen = Numeric.Ulen;
        final Ulen_offset = k1;
        final LU = LUbx[block];

        for (int k = 0; k < nk; k++) {

          /* ------------------------------------------------------ */
          /* scatter kth column of the block into workspace X */
          /* ------------------------------------------------------ */

          final oldcol = Q[k + k1];
          final pend = Ap[oldcol + 1];
          for (int p = Ap[oldcol]; p < pend; p++) {
            final newrow = Pinv[Ai[p]] - k1;
            if (newrow < 0 && poff < nzoff) {
              /* entry in off-diagonal block */
              Offx[poff] = Az[p];
              poff++;
            } else {
              /* (newrow,k) is an entry in the block */
              X[newrow] = Az[p];
            }
          }

          /* ------------------------------------------------------ */
          /* compute kth column of U, and update kth column of A */
          /* ------------------------------------------------------ */

          final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, ulen);
          final Ux = Ui;
          for (int up = 0; up < ulen[0]; up++) {
            final j = Ui[Ui_offset[0] + up].toInt();
            final ujk = X[j];
            /* X [j] = 0 ; */
            CLEAR(X, j);
            Ux[Ux_offset[0] + up] = ujk;
            final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, j, llen);
            final Lx = Li;
            for (int p = 0; p < llen[0]; p++) {
              //MULT_SUB (X [Li [p]], Lx [p], ujk) ;
              X[Li[Li_offset[0] + p].toInt()] -= Lx[Lx_offset[0] + p] * ujk;
            }
          }
          /* get the diagonal entry of U */
          final ukk = X[k];
          /* X [k] = 0 ; */
          CLEAR(X, k);
          /* singular case */
          if (IS_ZERO(ukk)) {
            /* matrix is numerically singular */
            Common.status = KLU_SINGULAR;
            if (Common.numerical_rank == EMPTY) {
              Common.numerical_rank = k + k1;
              Common.singular_col = Q[k + k1];
            }
            if (Common.halt_if_singular != 0) {
              /* do not continue the factorization */
              return (FALSE);
            }
          }
          Udiag[k + k1] = ukk;
          /* gather and divide by pivot to get kth column of L */
          final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, llen);
          final Lx = Li;
          for (int p = 0; p < llen[0]; p++) {
            final i = Li[Li_offset[0] + p].toInt();
            //DIV (Lx [p], X [i], ukk) ;
            Lx[Lx_offset[0] + p] = X[i] / ukk;
            CLEAR(X, i);
          }

        }
      }
    }

  } else {

    /* ------------------------------------------------------------------ */
    /* scaling */
    /* ------------------------------------------------------------------ */

    for (int block = 0; block < nblocks; block++) {

      /* -------------------------------------------------------------- */
      /* the block is from rows/columns k1 to k2-1 */
      /* -------------------------------------------------------------- */

      final k1 = R[block];
      final k2 = R[block + 1];
      final nk = k2 - k1;

      if (nk == 1) {

        /* ---------------------------------------------------------- */
        /* singleton case */
        /* ---------------------------------------------------------- */

        final oldcol = Q[k1];
        final pend = Ap[oldcol + 1];
        double s = 0.0; //CLEAR (s) ;
        for (int p = Ap[oldcol]; p < pend; p++) {
          final oldrow = Ai[p];
          final newrow = Pinv[oldrow] - k1;
          if (newrow < 0 && poff < nzoff) {
            /* entry in off-diagonal block */
            Offx[poff] = Az[p] / Rs[oldrow];
            //SCALE_DIV_ASSIGN (Offx [poff], Az [p], Rs [oldrow]) ;
            poff++;
          } else {
            /* singleton */
            s = Az[p] / Rs[oldrow];
            //SCALE_DIV_ASSIGN (s, Az [p], Rs [oldrow]) ;
          }
        }
        Udiag[k1] = s;

      } else {

        /* ---------------------------------------------------------- */
        /* construct and factor the kth block */
        /* ---------------------------------------------------------- */

        final Lip = Numeric.Lip;
        final Lip_offset = k1;
        final Llen = Numeric.Llen;
        final Llen_offset = k1;
        final Uip = Numeric.Uip;
        final Uip_offset = k1;
        final Ulen = Numeric.Ulen;
        final Ulen_offset = k1;
        final LU = LUbx[block];

        for (int k = 0; k < nk; k++) {

          /* ------------------------------------------------------ */
          /* scatter kth column of the block into workspace X */
          /* ------------------------------------------------------ */

          final oldcol = Q[k + k1];
          final pend = Ap[oldcol + 1];
          for (int p = Ap[oldcol]; p < pend; p++) {
            final oldrow = Ai[p];
            final newrow = Pinv[oldrow] - k1;
            if (newrow < 0 && poff < nzoff) {
              /* entry in off-diagonal part */
              //SCALE_DIV_ASSIGN (Offx [poff], Az [p], Rs [oldrow]);
              Offx[poff] = Az[p] / Rs[oldrow];
              poff++;
            } else {
              /* (newrow,k) is an entry in the block */
              //SCALE_DIV_ASSIGN (X [newrow], Az [p], Rs [oldrow]) ;
              X[newrow] = Az[p] / Rs[oldrow];
            }
          }

          /* ------------------------------------------------------ */
          /* compute kth column of U, and update kth column of A */
          /* ------------------------------------------------------ */

          final Ui = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, ulen);
          final Ux = Ui;
          for (int up = 0; up < ulen[0]; up++) {
            final j = Ui[Ui_offset[0] + up].toInt();
            final ujk = X[j];
            /* X [j] = 0 ; */
            CLEAR(X, j);
            Ux[Ux_offset[0] + up] = ujk;
            final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, j, llen);
            final Lx = Li;
            for (int p = 0; p < llen[0]; p++) {
              //MULT_SUB (X [Li [p]], Lx [p], ujk) ;
              X[Li[Li_offset[0] + p].toInt()] -= Lx[Lx_offset[0] + p] * ujk;
            }
          }
          /* get the diagonal entry of U */
          final ukk = X[k];
          /* X [k] = 0 ; */
          CLEAR(X, k);
          /* singular case */
          if (IS_ZERO(ukk)) {
            /* matrix is numerically singular */
            Common.status = KLU_SINGULAR;
            if (Common.numerical_rank == EMPTY) {
              Common.numerical_rank = k + k1;
              Common.singular_col = Q[k + k1];
            }
            if (Common.halt_if_singular != 0) {
              /* do not continue the factorization */
              return (FALSE);
            }
          }
          Udiag[k + k1] = ukk;
          /* gather and divide by pivot to get kth column of L */
          final Li = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, llen);
          final Lx = Li;
          for (int p = 0; p < llen[0]; p++) {
            final i = Li[Li_offset[0] + p].toInt();
            //DIV (Lx [p], X [i], ukk) ;
            Lx[Lx_offset[0] + p] = X[i] / ukk;
            CLEAR(X, i);
          }
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /* permute scale factors Rs according to pivotal row order */
  /* ---------------------------------------------------------------------- */

  if (scale > 0) {
    for (int k = 0; k < n; k++) {
      X[k] = Rs[Pnum[k]];
      //REAL (X [k]) = Rs [Pnum [k]] ;
    }
    for (int k = 0; k < n; k++) {
      Rs[k] = X[k];
      //Rs [k] = REAL (X [k]) ;
    }
  }

  if (!NDEBUG) {
    ASSERT(Offp[n] == poff);
    ASSERT(Symbolic.nzoff == poff);
    PRINTF(("\n------------------- Off diagonal entries, new:\n"));
    if (!NDEBUG) ASSERT_INT(_valid(n, Offp, Offi, Offx));
    if (Common.status == KLU_OK) {
      PRINTF("\n ########### KLU_BTF_REFACTOR done, nblocks $nblocks\n");
      for (int block = 0; block < nblocks; block++) {
        final k1 = R[block];
        final k2 = R[block + 1];
        final nk = k2 - k1;
        PRINTF("\n================KLU_refactor output: k1 $k1 k2 $k2 nk $nk\n");
        if (nk == 1) {
          PRINTF("singleton  ");
          PRINT_ENTRY(Udiag[k1]);
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

  return (TRUE);
}
