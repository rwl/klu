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
 * Debug routines for klu. Only used when NDEBUG is true.
 */
part of edu.ufl.cise.klu.tdouble;

/**
 * Check if a column-form matrix is valid or not. The matrix A is
 * [n]-by-[n].  The row indices of entries in column j are in
 * `Ai[Ap [j] ... Ap [j+1]-1]`. Required conditions are:
 *
 * 1. `n >= 0`
 * 2. `nz = Ap[n_col] >= 0` number of entries in the matrix
 * 3. `Ap[0] == 0`
 * 4. `Ap[j] <= Ap[j+1]` for all j in the range 0 to n_col.
 * 5. row indices in `Ai[Ap [j] ... Ap [j+1]-1]` must be in the
 *    range 0 to n_row-1, and no duplicate entries can exist
 *    (duplicates not checked here).
 *
 * Not user-callable.  Only used when debugging.
 */
int _valid(final int n, final Int32List Ap, final Int32List Ai,
           final Float64List Ax) {
  PRINTF("\ncolumn oriented matrix, n = $n\n");
  if (n <= 0) {
    PRINTF("n must be >= 0: $n\n");
    return (FALSE);
  }
  final nz = Ap[n];
  if (Ap[0] != 0 || nz < 0) {
    /* column pointers must start at Ap [0] = 0, and Ap [n] must be >= 0 */
    PRINTF("column 0 pointer bad or nz < 0\n");
    return (FALSE);
  }
  for (int j = 0; j < n; j++) {
    final p1 = Ap[j];
    final p2 = Ap[j + 1];
    PRINTF("\nColumn: $j p1: $p1 p2: $p2\n");
    if (p1 > p2) {
      /* column pointers must be ascending */
      PRINTF("column $j pointer bad\n");
      return (FALSE);
    }
    for (int p = p1; p < p2; p++) {
      final i = Ai[p];
      PRINTF("row: $i");
      if (i < 0 || i >= n) {
        /* row index out of range */
        PRINTF("index out of range, col $j row $i\n");
        return (FALSE);
      }
      if (Ax != null) {
        PRINT_ENTRY(Ax[p]);
      }
      PRINTF("\n");
    }
  }
  return (TRUE);
}

/**
 * This function does the same validity tests as KLU_valid but for the
 * LU factor storage format. The flag flag_test_start_ptr is used to
 * test if Xip [0] = 0. This is not applicable for U. So when calling this
 * function for U, the flag should be set to false. Only used when debugging.
 */
_valid_LU(final int n, final int flag_test_start_ptr, final Int32List Xip,
          final int Xip_offset, final Int32List Xlen, final int Xlen_offset,
          final Float64List LU) {
  final len = new Int32List(1);
  final Xi_offset = new Int32List(1);
  final Xx_offset = new Int32List(1);

  PRINTF("\ncolumn oriented matrix, n = $n\n");
  if (n <= 0) {
    PRINTF("n must be >= 0: $n\n");
    return (FALSE);
  }
  if (flag_test_start_ptr != 0 && Xip[Xip_offset + 0] != 0) {
    /* column pointers must start at Xip [0] = 0*/
    PRINTF("column 0 pointer bad\n");
    return (FALSE);
  }

  for (int j = 0; j < n; j++) {
    final p1 = Xip[Xip_offset + j];
    final p2 = Xip[Xip_offset + j + 1];
    PRINTF("\nColumn: $j p1: $p1 p2: $p2\n");
    if (p1 > p2) {
      /* column pointers must be ascending */
      PRINTF("column $j pointer bad\n");
      return (FALSE);
    }
    final Xi = GET_POINTER(LU, Xip, Xip_offset, Xlen, Xlen_offset, Xi_offset, Xx_offset, j, len);
    final Xx = Xi;
    for (int p = 0; p < len[0]; p++) {
      final i = Xi[Xi_offset[0] + p].toInt();
      PRINTF("row: $i");
      if (i < 0 || i >= n) {
        /* row index out of range */
        PRINTF("index out of range, col $j row $i\n");
        return (FALSE);
      }
      if (Xx != null) {
        PRINT_ENTRY(Xx[Xx_offset[0] + p]);
      }
      PRINTF("\n");
    }
  }

  return (TRUE);
}
