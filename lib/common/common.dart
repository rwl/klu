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

library edu.ufl.cise.klu.common;

import 'dart:typed_data';

part 'defaults.dart';
part 'l_common.dart';
part 'l_defaults.dart';
part 'l_numeric.dart';
part 'l_symbolic.dart';
part 'numeric.dart';
part 'symbolic.dart';
part 'user_order.dart';
part 'version.dart';

/**
 * KLU control parameters and statistics.
 */
class KLU_common
{

	/* ---------------------------------------------------------------------- */
	/* parameters */
	/* ---------------------------------------------------------------------- */

	double tol;             /* pivot tolerance for diagonal preference */
	double memgrow;         /* realloc memory growth size for LU factors */
	double initmem_amd;     /* init. memory size with AMD: c*nnz(L) + n */
	double initmem;         /* init. memory size: c*nnz(A) + n */
	double maxwork;         /* maxwork for BTF, <= 0 if no limit */

	int btf;                /* use BTF pre-ordering, or not */
	int ordering;           /* 0: AMD, 1: COLAMD, 2: user P and Q,
                                    * 3: user function */
	int scale;              /* row scaling: -1: none (and no error check),
                                    * 0: none, 1: sum, 2: max */

	/* memory management routines */
//	void *(*malloc_memory) (size_t);            /* pointer to malloc */
//	void *(*realloc_memory) (void *, size_t);   /* pointer to realloc */
//	void (*free_memory) (void *);               /* pointer to free */
//	void *(*calloc_memory) (size_t, size_t);    /* pointer to calloc */

	/* pointer to user ordering function */
	KLU_user_order user_order;


	/* pointer to user data, passed unchanged as the last parameter to the
	 * user ordering function (optional, the user function need not use this
	 * information). */
	Object user_data;

	int halt_if_singular;       /* how to handle a singular matrix:
	    * FALSE: keep going.  Return a Numeric object with a zero U(k,k).  A
	    *   divide-by-zero may occur when computing L(:,k).  The Numeric object
	    *   can be passed to klu_solve (a divide-by-zero will occur).  It can
	    *   also be safely passed to klu_refactor.
	    * TRUE: stop quickly.  klu_factor will free the partially-constructed
	    *   Numeric object.  klu_refactor will not free it, but will leave the
	    *   numerical values only partially defined.  This is the default. */

	/* ---------------------------------------------------------------------- */
	/* statistics */
	/* ---------------------------------------------------------------------- */

	int status;                 /* KLU_OK if OK, < 0 if error */
	int nrealloc;               /* # of reallocations of L and U */

	int structural_rank;        /* 0 to n-1 if the matrix is structurally rank
	    * deficient (as determined by maxtrans).  -1 if not computed.  n if the
	    * matrix has full structural rank.  This is computed by klu_analyze
	    * if a BTF preordering is requested. */

	int numerical_rank;         /* First k for which a zero U(k,k) was found,
	    * if the matrix was singular (in the range 0 to n-1).  n if the matrix
	    * has full rank. This is not a true rank-estimation.  It just reports
	    * where the first zero pivot was found.  -1 if not computed.
	    * Computed by klu_factor and klu_refactor. */

	int singular_col;           /* n if the matrix is not singular.  If in the
	    * range 0 to n-1, this is the column index of the original matrix A that
	    * corresponds to the column of U that contains a zero diagonal entry.
	    * -1 if not computed.  Computed by klu_factor and klu_refactor. */

	int noffdiag;       /* # of off-diagonal pivots, -1 if not computed */

	double flops;       /* actual factorization flop count, from klu_flops */
	double rcond;       /* crude reciprocal condition est., from klu_rcond */
	double condest;     /* accurate condition est., from klu_condest */
	double rgrowth;     /* reciprocal pivot rgrowth, from klu_rgrowth */
	double work;        /* actual work done in BTF, in klu_analyze */

	int memusage;    /* current memory usage, in bytes */
	int mempeak;     /* peak memory usage, in bytes */

}
