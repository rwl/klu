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
 * Sets default parameters for KLU.
 */
int defaults(KLU_common Common) {
  if (Common == null) {
    return (FALSE);
  }

  /* parameters */
  /** pivot tolerance for diagonal */
  Common.tol = 0.001;
  /** realloc size ratio increase for LU factors */
  Common.memgrow = 1.2;
  /** init. mem with AMD:  c*nnz(L) + n */
  Common.initmem_amd = 1.2;
  /** init. mem otherwise: c*nnz(A) + n */
  Common.initmem = 10.0;
  /** use BTF pre-ordering, or not */
  Common.btf = TRUE;
  /** no limit to work done by btf_order */
  Common.maxwork = 0.0;
  /** 0: AMD, 1: COLAMD, 2: user-provided P and Q,
   * 3: user-provided function */
  Common.ordering = 0;
  /** scale: -1: none, and do not check for errors
  * in the input matrix in KLU_refactor.
  * 0: none, but check for errors,
  * 1: sum, 2: max */
  Common.scale = 2;
  /** quick halt if matrix is singular */
  Common.halt_if_singular = TRUE;

  /* memory management routines */
  //Common.malloc_memory  = malloc ;
  //Common.calloc_memory  = calloc ;
  //Common.free_memory    = free ;
  //Common.realloc_memory = realloc ;

  /* user ordering function and optional argument */
  Common.user_order = null;
  Common.user_data = null;

  /* statistics */
  Common.status = KLU_OK;
  Common.nrealloc = 0;
  Common.structural_rank = EMPTY;
  Common.numerical_rank = EMPTY;
  Common.noffdiag = EMPTY;
  Common.flops = EMPTY_D;
  Common.rcond = EMPTY_D;
  Common.condest = EMPTY_D;
  Common.rgrowth = EMPTY_D;
  Common.work = 0.0;
  /* work done by btf_order */

  Common.memusage = 0;
  Common.mempeak = 0;

  return (TRUE);
}
