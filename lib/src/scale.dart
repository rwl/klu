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
 * Scale a matrix and check to see if it is valid. Can be called by the user.
 * This is called by KLU_factor and KLU_refactor. Returns `true` if the input
 * matrix is valid, `false` otherwise. If the [W] input argument is non-null,
 * then the input matrix is checked for duplicate entries.
 *
 * [scale] methods:
 *  * `<0`: no scaling, do not compute Rs, and do not check input matrix.
 *  * `0`: no scaling
 *  * `1`: the scale factor for row i is sum (abs (A (i,:)))
 *  * `2`: or more: the scale factor for row i is max(abs (A (i,:)))
 */
int klu_scale(final int scale, final int n, final Int32List Ap,
              final Int32List Ai, final Float64List Ax, final Float64List Rs,
              final Int32List W, final KLU_common Common) {
//  double a;
//  Float64List Az;
//  int row, col, p, pend;
//  bool check_duplicates;

  /* ---------------------------------------------------------------------- */
  /* check inputs */
  /* ---------------------------------------------------------------------- */

  if (Common == null) {
    return (FALSE);
  }
  Common.status = KLU_OK;

  if (scale < 0) {
    /* return without checking anything and without computing the
     * scale factors */
    return (TRUE);
  }

  final Az = new Float64List.fromList(Ax);

  if (n <= 0 || Ap == null || Ai == null || Az == null || (scale > 0 && Rs == null)) {
    /* Ap, Ai, Ax and Rs must be present, and n must be > 0 */
    Common.status = KLU_INVALID;
    return (FALSE);
  }
  if (Ap[0] != 0 || Ap[n] < 0) {
    /* nz = Ap [n] must be >= 0 and Ap [0] must equal zero */
    Common.status = KLU_INVALID;
    return (FALSE);
  }
  for (int col = 0; col < n; col++) {
    if (Ap[col] > Ap[col + 1]) {
      /* column pointers must be non-decreasing */
      Common.status = KLU_INVALID;
      return (FALSE);
    }
  }

  /* ---------------------------------------------------------------------- */
  /* scale */
  /* ---------------------------------------------------------------------- */

  if (scale > 0) {
    /* initialize row sum or row max */
    for (int row = 0; row < n; row++) {
      Rs[row] = 0.0;
    }
  }

  /* check for duplicates only if W is present */
  final check_duplicates = (W != null);
  if (check_duplicates) {
    for (int row = 0; row < n; row++) {
      W[row] = EMPTY;
    }
  }

  for (int col = 0; col < n; col++) {
    final pend = Ap[col + 1];
    for (int p = Ap[col]; p < pend; p++) {
      final row = Ai[p];
      if (row < 0 || row >= n) {
        /* row index out of range, or duplicate entry */
        Common.status = KLU_INVALID;
        return (FALSE);
      }
      if (check_duplicates) {
        if (W[row] == col) {
          /* duplicate entry */
          Common.status = KLU_INVALID;
          return (FALSE);
        }
        /* flag row i as appearing in column col */
        W[row] = col;
      }
      //ABS (a, Az [p]) ;
      final a = ABS(Az[p]);
      if (scale == 1) {
        /* accumulate the abs. row sum */
        Rs[row] += a;
      } else if (scale > 1) {
        /* find the max abs. value in the row */
        Rs[row] = MAX(Rs[row], a);
      }
    }
  }

  if (scale > 0) {
    /* do not scale empty rows */
    for (int row = 0; row < n; row++) {
      /* matrix is singular */
      PRINTF("Rs [$row] = ${Rs [row]}\n");

      if (Rs[row] == 0.0) {
        PRINTF("Row $row of A is all zero\n");
        Rs[row] = 1.0;
      }
    }
  }

  return (TRUE);
}
