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
 * Extract KLU factorization into conventional compressed-column matrices.
 * If any output array is `null`, that part of the LU factorization is not
 * extracted (this is not an error condition).
 *
 * `nnz(L) = Numeric.lnz, nnz(U) = Numeric.unz, and nnz(F) = Numeric.Offp[n]`
 *
 * [Lp] is size `n+1`. [Li] is size `nnz(L)`. Lx is size `nnz(L)`. Up is
 * size `n+1`. [Ui] is size `nnz(U)`. [Ux] is size `nnz(U)`. Fp is size `n+1`.
 * [Fi] is size `nnz(F)`. [Fx] is size `nnz(F)`. Row permutation [P] is size
 * `n`. Column permutation [Q] is size `n`. Scale factors [Rs] is size `n`.
 * Block boundaries [R] is size `nblocks+1`.
 */
int extract(final KLU_numeric Numeric, final KLU_symbolic Symbolic,
            final Int32List Lp, final Int32List Li, final Float64List Lx,
            final Int32List Up, final Int32List Ui, final Float64List Ux,
            final Int32List Fp, final Int32List Fi, final Float64List Fx,
            final Int32List P, final Int32List Q, final Float64List Rs,
            final Int32List R, final KLU_common Common) {
  final len = new Int32List(1);
  final Li2_offset = new Int32List(1);
  final Lx2_offset = new Int32List(1);
  final Ui2_offset = new Int32List(1);
  final Ux2_offset = new Int32List(1);

  if (Common == null) {
    return (FALSE);
  }

  if (Symbolic == null || Numeric == null) {
    Common.status = KLU_INVALID;
    return (FALSE);
  }

  Common.status = KLU_OK;
  final n = Symbolic.n;
  final nblocks = Symbolic.nblocks;

  /* ---------------------------------------------------------------------- */
  /* extract scale factors */
  /* ---------------------------------------------------------------------- */

  if (Rs != null) {
    if (Numeric.Rs != null) {
      for (int i = 0; i < n; i++) {
        Rs[i] = Numeric.Rs[i];
      }
    } else {
      /* no scaling */
      for (int i = 0; i < n; i++) {
        Rs[i] = 1.0;
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /* extract block boundaries */
  /* ---------------------------------------------------------------------- */

  if (R != null) {
    for (int block = 0; block <= nblocks; block++) {
      R[block] = Symbolic.R[block];
    }
  }

  /* ---------------------------------------------------------------------- */
  /* extract final row permutation */
  /* ---------------------------------------------------------------------- */

  if (P != null) {
    for (int k = 0; k < n; k++) {
      P[k] = Numeric.Pnum[k];
    }
  }

  /* ---------------------------------------------------------------------- */
  /* extract column permutation */
  /* ---------------------------------------------------------------------- */

  if (Q != null) {
    for (int k = 0; k < n; k++) {
      Q[k] = Symbolic.Q[k];
    }
  }

  /* ---------------------------------------------------------------------- */
  /* extract each block of L */
  /* ---------------------------------------------------------------------- */

  int nz;
  if (Lp != null && Li != null && Lx != null) {
    nz = 0;
    for (int block = 0; block < nblocks; block++) {
      final k1 = Symbolic.R[block];
      final k2 = Symbolic.R[block + 1];
      final nk = k2 - k1;
      if (nk == 1) {
        /* singleton block */
        Lp[k1] = nz;
        Li[nz] = k1;
        Lx[nz] = 1.0;
        nz++;
      } else {
        /* non-singleton block */
        final LU = Numeric.LUbx[block];
        final Lip = Numeric.Lip;
        final Lip_offset = k1;
        final Llen = Numeric.Llen;
        int Llen_offset = k1;
        for (int kk = 0; kk < nk; kk++) {
          Lp[k1 + kk] = nz;
          /* add the unit diagonal entry */
          Li[nz] = k1 + kk;
          Lx[nz] = 1.0;
          nz++;
          final Li2 = GET_POINTER(LU, Lip, Lip_offset, Llen, Llen_offset,
              Li2_offset, Lx2_offset, kk, len);
          final Lx2 = Li2;
          for (int p = 0; p < len[0]; p++) {
            Li[nz] = k1 + Li2[Li2_offset[0] + p].toInt();
            Lx[nz] = Lx2[Lx2_offset[0] + p]; //REAL (Lx2 [p]) ;
            nz++;
          }
        }
      }
    }
    Lp[n] = nz;
    ASSERT(nz == Numeric.lnz);
  }

  /* ---------------------------------------------------------------------- */
  /* extract each block of U */
  /* ---------------------------------------------------------------------- */

  if (Up != null && Ui != null && Ux != null) {
    nz = 0;
    for (int block = 0; block < nblocks; block++) {
      final k1 = Symbolic.R[block];
      final k2 = Symbolic.R[block + 1];
      final nk = k2 - k1;
      final Ukk = Numeric.Udiag;
      final Ukk_offset = k1;
      if (nk == 1) {
        /* singleton block */
        Up[k1] = nz;
        Ui[nz] = k1;
        Ux[nz] = Ukk[Ukk_offset + 0]; //REAL (Ukk [0]) ;
        nz++;
      } else {
        /* non-singleton block */
        final LU = Numeric.LUbx[block];
        final Uip = Numeric.Uip;
        final Uip_offset = k1;
        final Ulen = Numeric.Ulen;
        final Ulen_offset = k1;
        for (int kk = 0; kk < nk; kk++) {
          Up[k1 + kk] = nz;
          final Ui2 = GET_POINTER(LU, Uip, Uip_offset, Ulen, Ulen_offset,
              Ui2_offset, Ux2_offset, kk, len);
          final Ux2 = Ui2;
          for (int p = 0; p < len[0]; p++) {
            Ui[nz] = k1 + Ui2[Ui2_offset[0] + p].toInt();
            Ux[nz] = Ux2[Ux2_offset[0] + p]; //REAL (Ux2 [p]) ;
            nz++;
          }
          /* add the diagonal entry */
          Ui[nz] = k1 + kk;
          Ux[nz] = Ukk[Ukk_offset + kk]; //REAL (Ukk [kk]) ;
          nz++;
        }
      }
    }
    Up[n] = nz;
    ASSERT(nz == Numeric.unz);
  }

  /* ---------------------------------------------------------------------- */
  /* extract the off-diagonal blocks, F */
  /* ---------------------------------------------------------------------- */

  if (Fp != null && Fi != null && Fx != null) {
    for (int k = 0; k <= n; k++) {
      Fp[k] = Numeric.Offp[k];
    }
    nz = Fp[n];
    for (int k = 0; k < nz; k++) {
      Fi[k] = Numeric.Offi[k];
    }
    for (int k = 0; k < nz; k++) {
      Fx[k] = Numeric.Offx[k];
      //Fx [k] = REAL (((Float64List) Numeric.Offx) [k]) ;
    }
  }

  return (TRUE);
}
