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
 * Numeric object - contains the factors computed by klu_factor.
 */
class KLU_numeric
{

	/* LU factors of each block, the pivot row permutation, and the
	 * entries in the off-diagonal blocks */

	int n;             /* A is n-by-n */
	int nblocks;       /* number of diagonal blocks */
	int lnz;           /* actual nz in L, including diagonal */
	int unz;           /* actual nz in U, including diagonal */
	int max_lnz_block; /* max actual nz in L in any one block, incl. diag */
	int max_unz_block; /* max actual nz in U in any one block, incl. diag */
	List<int> Pnum;        /* size n. final pivot permutation */
	List<int> Pinv;        /* size n. inverse of final pivot permutation */

	/* LU factors of each block */
	List<int> Lip;         /* size n. pointers into LUbx[block] for L */
	List<int> Uip;         /* size n. pointers into LUbx[block] for U */
	List<int> Llen;        /* size n. Llen [k] = # of entries in kth column of L */
	List<int> Ulen;        /* size n. Ulen [k] = # of entries in kth column of U */
  List<List<double>> LUbx;     /* L and U indices and entries (excl. diagonal of U) */  // FIXME: indices are doubles and cast to integers
	List<int> LUsize;   /* size of each LUbx [block], in sizeof (Unit) */
	List<double> Udiag;      /* diagonal of U */

	/* scale factors; can be NULL if no scaling */
	List<double> Rs;       /* size n. Rs [i] is scale factor for row i */

	/* permanent workspace for factorization and solve */
	int worksize; /* size (in bytes) of Work */
	List<double> Work;       /* workspace */
	List<double> Xwork;      /* alias into Numeric->Work */
	List<int> Iwork;       /* alias into Numeric->Work */

	/* off-diagonal entries in a conventional compressed-column sparse matrix */
	List<int> Offp;        /* size n+1, column pointers */
	List<int> Offi;        /* size nzoff, row indices */
	List<double> Offx;       /* size nzoff, numerical values */
	int nzoff;

}
