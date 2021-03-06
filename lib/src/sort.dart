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
 * Sort L or U using a double-transpose.
 */
void _sort(final int n, final Int32List Xip, final int Xip_offset,
           final Int32List Xlen, final int Xlen_offset, final Float64List LU,
           final Int32List Tp, final Int32List Tj, final Float64List Tx,
           final Int32List W) {
  final len = new Int32List(1);
  final Xi_offset = new Int32List(1);
  final Xx_offset = new Int32List(1);

  ASSERT(_valid_LU(n, FALSE, Xip, Xip_offset, Xlen, Xlen_offset, LU));

  /* count the number of entries in each row of L or U */
  for (int i = 0; i < n; i++) {
    W[i] = 0;
  }
  for (int j = 0; j < n; j++) {
    final Xi = GET_POINTER(LU, Xip, Xip_offset, Xlen, Xlen_offset, Xi_offset, Xx_offset, j, len);
    final Xx = Xi;
    for (int p = 0; p < len[0]; p++) {
      W[Xi[Xi_offset[0] + p].toInt()]++;
    }
  }

  /* construct the row pointers for T */
  int nz = 0;
  for (int i = 0; i < n; i++) {
    Tp[i] = nz;
    nz += W[i];
  }
  Tp[n] = nz;
  for (int i = 0; i < n; i++) {
    W[i] = Tp[i];
  }

  /* transpose the matrix into Tp, Ti, Tx */
  for (int j = 0; j < n; j++) {
    final Xi = GET_POINTER(LU, Xip, Xip_offset, Xlen, Xlen_offset, Xi_offset, Xx_offset, j, len);
    final Xx = Xi;
    for (int p = 0; p < len[0]; p++) {
      final tp = W[Xi[Xi_offset[0] + p].toInt()]++;
      Tj[tp] = j;
      Tx[tp] = Xx[Xx_offset[0] + p];
    }
  }

  /* transpose the matrix back into Xip, Xlen, Xi, Xx */
  for (int j = 0; j < n; j++) {
    W[j] = 0;
  }
  for (int i = 0; i < n; i++) {
    final pend = Tp[i + 1];
    for (int p = Tp[i]; p < pend; p++) {
      final j = Tj[p];
      final Xi = GET_POINTER(LU, Xip, Xip_offset, Xlen, Xlen_offset, Xi_offset, Xx_offset, j, len);
      final Xx = Xi;
      final xlen = W[j]++;
      Xi[Xi_offset[0] + xlen] = i.toDouble();
      Xx[Xx_offset[0] + xlen] = Tx[p];
    }
  }

  ASSERT(_valid_LU(n, FALSE, Xip, Xip_offset, Xlen, Xlen_offset, LU));
}

/**
 * Sorts the columns of `L` and `U` so that the row indices appear in strictly
 * increasing order.
 */
int sort(final KLU_symbolic Symbolic, final KLU_numeric Numeric,
         final KLU_common Common) {
  if (Common == null) {
    return (FALSE);
  }
  Common.status = KLU_OK;

  final R = Symbolic.R;
  final nblocks = Symbolic.nblocks;
  final maxblock = Symbolic.maxblock;

  final Lip = Numeric.Lip;
  final Llen = Numeric.Llen;
  final Uip = Numeric.Uip;
  final Ulen = Numeric.Ulen;
  final LUbx = new List<Float64List>.from(Numeric.LUbx);

  final m1 = (maxblock.toInt()) + 1;

  /* allocate workspace */
  final nz = MAX(Numeric.max_lnz_block, Numeric.max_unz_block);
  final W = malloc_int(maxblock, Common);
  final Tp = malloc_int(m1, Common);
  final Ti = malloc_int(nz, Common);
  final Tx = malloc_dbl(nz, Common);

  PRINTF("\n======================= Start sort:\n");

  if (Common.status == KLU_OK) {
    /* sort each block of L and U */
    for (int block = 0; block < nblocks; block++) {
      final k1 = R[block];
      final nk = R[block + 1] - k1;
      if (nk > 1) {
        PRINTF("\n-------------------block: $block nk $nk\n");
        _sort(nk, Lip, k1, Llen, k1, LUbx[block], Tp, Ti, Tx, W);
        _sort(nk, Uip, k1, Ulen, k1, LUbx[block], Tp, Ti, Tx, W);
      }
    }
  }

  PRINTF("\n======================= sort done.\n");

  /* free workspace */
  //KLU_free (W, maxblock, sizeof (Int), Common) ;
  //KLU_free (Tp, m1, sizeof (Int), Common) ;
  //KLU_free (Ti, nz, sizeof (Int), Common) ;
  //KLU_free (Tx, nz, sizeof (double), Common) ;

  return (Common.status == KLU_OK ? 1 : 0);
}
