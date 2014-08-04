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
class KLU_l_numeric
{

	/* LU factors of each block, the pivot row permutation, and the
     * entries in the off-diagonal blocks */

    int n;             /* A is n-by-n */
    int nblocks;       /* number of diagonal blocks */
    int lnz;           /* actual nz in L, including diagonal */
    int unz;           /* actual nz in U, including diagonal */
    int max_lnz_block; /* max actual nz in L in any one block, incl. diag */
    int max_unz_block; /* max actual nz in U in any one block, incl. diag */
    Int32List Pnum;        /* size n. final pivot permutation */
    Int32List Pinv;        /* size n. inverse of final pivot permutation */

    /* LU factors of each block */
    Int32List Lip;         /* size n. pointers into LUbx[block] for L */
    Int32List Uip;         /* size n. pointgers into LUbx[block] for U */
    Int32List Llen;        /* size n. Llen [k] = # of entries in kth column of L */
    Int32List Ulen;        /* size n. Ulen [k] = # of entries in kth column of U */
    List<Float64List> LUbx;      /* L and U indices and entries (excl. diagonal of U) */
    Int32List LUsize;    /* size of each LUbx [block], in sizeof (Unit) */
    Float64List Udiag;       /* diagonal of U */

    /* scale factors; can be NULL if no scaling */
    Float64List Rs;        /* size n. Rs [i] is scale factor for row i */

    /* permanent workspace for factorization and solve */
    Int32List worksize;  /* size (in bytes) of Work */
    Float64List Work;        /* workspace */
    Float64List Xwork;       /* alias into Numeric->Work */
    Int32List Iwork;       /* alias into Numeric->Work */

    /* off-diagonal entries in a conventional compressed-column sparse matrix */
    Int32List Offp;        /* size n+1, column pointers */
    Int32List Offi;        /* size nzoff, row indices */
    Float64List Offx;        /* size nzoff, numerical values */
    int nzoff;

}
