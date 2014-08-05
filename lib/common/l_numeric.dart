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
 * License aint with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

part of edu.ufl.cise.klu.common;

/**
 * Numeric object - contains the factors computed by klu_factor.
 *
 * 64-bit version.
 */
class KLU_l_numeric {

  /* LU factors of each block, the pivot row permutation, and the
   * entries in the off-diagonal blocks */

  /** A is n-by-n */
  int n;
  /** number of diagonal blocks */
  int nblocks;
  /** actual nz in L, including diagonal */
  int lnz;
  /** actual nz in U, including diagonal */
  int unz;
  /** max actual nz in L in any one block, incl. diag */
  int max_lnz_block;
  /** max actual nz in U in any one block, incl. diag */
  int max_unz_block;
  /** size n. final pivot permutation */
  Int32List Pnum;
  /** size n. inverse of final pivot permutation */
  Int32List Pinv;

  /* LU factors of each block */
  /** size n. pointers into LUbx[block] for L */
  Int32List Lip;
  /** size n. pointgers into LUbx[block] for U */
  Int32List Uip;
  /** size n. Llen [k] = # of entries in kth column of L */
  Int32List Llen;
  /** size n. Ulen [k] = # of entries in kth column of U */
  Int32List Ulen;
  /** L and U indices and entries (excl. diagonal of U) */
  List<Float64List> LUbx;
  /** size of each LUbx [block], in sizeof (Unit) */
  Int32List LUsize;
  /** diagonal of U */
  Float64List Udiag;

  /* scale factors; can be NULL if no scaling */
  /** size n. Rs [i] is scale factor for row i */
  Float64List Rs;

  /* permanent workspace for factorization and solve */
  /** size (in bytes) of Work */
  Int32List worksize;
  /** workspace */
  Float64List Work;
  /** alias into Numeric.Work */
  Float64List Xwork;
  /** alias into Numeric->Work */
  Int32List Iwork;

  /* off-diagonal entries in a conventional compressed-column sparse matrix */
  /** size n+1, column pointers */
  Int32List Offp;
  /** size nzoff, row indices */
  Int32List Offi;
  /** size nzoff, numerical values */
  Float64List Offx;
  int nzoff;

}
