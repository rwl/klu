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

part of edu.ufl.cise.klu.common;

/**
 * Symbolic object - contains the pre-ordering computed by klu_analyze.
 *
 * 64-bit version.
 */
class KLU_l_symbolic {
  /* A (P,Q) is in upper block triangular form.  The kth block goes from
     * row/col index R [k] to R [k+1]-1. The estimated number of nonzeros
     * in the L factor of the kth block is Lnz [k].
     */

  /* only computed if the AMD ordering is chosen: */
  /** symmetry of largest block */
  double symmetry;
  /** est. factorization flop count */
  double est_flops;
  /** estimated nz in L and U, including diagonals */
  double lnz, unz;
  /** size n, but only Lnz [0..nblocks-1] is used */
  Float64List Lnz;

  /* computed for all orderings: */
  /** input matrix A is n-by-n */
  int n;
  /** # entries in input matrix */
  int nz;
  /** nz in off-diagonal blocks */
  int nzoff;
  /** number of blocks */
  int nblocks;
  /** size of largest block */
  int maxblock;
  /** ordering used (AMD, COLAMD, or GIVEN) */
  int ordering;
  /** whether or not BTF preordering was requested */
  int do_btf;

  Int32List P; // size n
  Int32List Q; // size n
  Int32List R; // size n+1, but only R [0..nblocks] is used

  /* only computed if BTF preordering requested */
  /** 0 to n-1 if the matrix is structurally rank
   * deficient.  -1 if not computed.  n if the matrix has
   * full structural rank */
  int structural_rank;

}
