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

library edu.ufl.cise.klu.tdouble;

//import edu.ufl.cise.klu.common.KLU_common;
import './common/common.dart';

import 'package:btf/btf.dart' as btf;
import 'package:colamd/colamd.dart' as colamd;
import 'package:amd/amd.dart' as amd;

//import static edu.ufl.cise.klu.tdouble.Dklu_memory.klu_malloc_dbl;
//import static edu.ufl.cise.klu.tdouble.Dklu_kernel.klu_kernel;

part 'analyze.dart';
part 'analyze_given.dart';
part 'defaults.dart';
part 'diagnostics.dart';
part 'dump.dart';
part 'extract.dart';
part 'factor.dart';
part 'internal.dart';
part 'kernel.dart';
part 'memory.dart';
part 'refactor.dart';
part 'scale.dart';
part 'solve.dart';
part 'sort.dart';
part 'tsolve.dart';
part 'version.dart';

/**
 * KLU: factorizes P*A into L*U, using the Gilbert-Peierls algorithm[1], with
 * optional symmetric pruning by Eisenstat and Liu[2].  The code is by Tim
 * Davis.  This algorithm is what appears as the default sparse LU routine in
 * MATLAB version 6.0, and still appears in MATLAB 6.5 as[L,U,P] = lu (A).
 * Note that no column ordering is provided (see COLAMD or AMD for suitable
 * orderings).  SuperLU is based on this algorithm, except that it adds the
 * use of dense matrix operations on "supernodes" (adjacent columns with
 * identical).  This code doesn't use supernodes, thus its name ("Kent" LU,
 * as in "Clark Kent", in contrast with Super-LU...).  This algorithm is slower
 * than SuperLU and UMFPACK for large matrices with lots of nonzeros in their
 * factors (such as for most finite-element problems).  However, for matrices
 * with very sparse LU factors, this algorithm is typically faster than both
 * SuperLU and UMFPACK, since in this case there is little chance to exploit
 * dense matrix kernels (the BLAS).
 *
 * Only one block of A is factorized, in the BTF form.  The input n is the
 * size of the block; k1 is the first row and column in the block.
 *
 * NOTE: no error checking is done on the inputs.  This version is not meant to
 * be called directly by the user.  Use klu_factor instead.
 *
 * No fill-reducing ordering is provided.  The ordering quality of
 * klu_kernel_factor is the responsibility of the caller.  The input A must
 * pre-permuted to reduce fill-in, or fill-reducing input permutation Q must
 * be provided.
 *
 * The input matrix A must be in compressed-column form, with either sorted
 * or unsorted row indices.  Row indices for column j of A is in
 * Ai[Ap[j] ... Ap[j+1]-1] and the same range of indices in Ax holds the
 * numerical values.  No duplicate entries are allowed.
 *
 * Copyright 2004-2009, Tim Davis.  All rights reserved.  See the README
 * file for details on permitted use.  Note that no code from The MathWorks,
 * Inc, or from SuperLU, or from any other source appears here.  The code is
 * written from scratch, from the algorithmic description in Gilbert & Peierls'
 * and Eisenstat & Liu's journal papers[1,2].
 *
 * If an input permutation Q is provided, the factorization L*U = A (P,Q)
 * is computed, where P is determined by partial pivoting, and Q is the input
 * ordering.  If the pivot tolerance is less than 1, the "diagonal" entry that
 * KLU attempts to choose is the diagonal of A (Q,Q).  In other words, the
 * input permutation is applied symmetrically to the input matrix.  The output
 * permutation P includes both the partial pivoting ordering and the input
 * permutation.  If Q is null, then it is assumed to be the identity
 * permutation.  Q is not modified.
 *
 *[1] Gilbert, J. R. and Peierls, T., "Sparse Partial Pivoting in Time
 *      Proportional to Arithmetic Operations," SIAM J. Sci. Stat. Comp.,
 *      vol 9, pp.  862-874, 1988.
 *[2] Eisenstat, S. C. and Liu, J. W. H., "Exploiting Structural Symmetry in
 *      Unsymmetric Sparse Symbolic Factorization," SIAM J. Matrix Analysis &
 *      Applic., vol 13, pp.  202-211, 1992.
 */
//public class Dklu extends Dklu_internal {

/**
 *
 * @param n A is n-by-n. n must be > 0.
 * @param Ap size n+1, column pointers for A
 * @param Ai size nz = Ap[n], row indices for A
 * @param Ax size nz, values of A
 * @param Q size n, optional column permutation
 * @param Lsize estimate of number of nonzeros in L
 * @param p_LU row indices and values of L and U
 * @param Udiag size n, diagonal of U
 * @param Llen size n, column length of L
 * @param Ulen size n, column length of U
 * @param Lip size n, column pointers for L
 * @param Uip size n, column pointers for U
 * @param P row permutation, size n
 * @param lnz size of L
 * @param unz size of U
 * @param X size n double's, zero on output
 * @param Work size 5n int's
 * @param k1 the block of A is from k1 to k2-1
 * @param PSinv inverse of P from symbolic factorization
 * @param Rs scale factors for A
 * @param Offp off-diagonal matrix (modified by this routine)
 * @param Offi
 * @param Offx
 * @param Common
 * @return
 */
int kernel_factor(int n, List<int> Ap, List<int> Ai,
		List<double> Ax, List<int> Q, double Lsize,
		List<List<double>> p_LU, int block,
		List<double> Udiag, int Udiag_offset, List<int> Llen, int Llen_offset,
		List<int> Ulen, int Ulen_offset, List<int> Lip, int Lip_offset,
		List<int> Uip, int Uip_offset, List<int> P, List<int> lnz, List<int> unz,
		List<double> X, List<int> Work, int k1, List<int> PSinv, List<double> Rs,
		List<int> Offp, List<int> Offi, List<double> Offx, KLU_common Common)
{
	double maxlnz, dunits ;
	List<List<double>> LU = new List<List<double>>(1);//[] ;
	List<int> Pinv, Lpend, Stack, Flag, Ap_pos ;
	int lsize, usize, anz, ok ;
	int lusize ;
	ASSERT (Common != null) ;

	/* ---------------------------------------------------------------------- */
	/* get control parameters, or use defaults */
	/* ---------------------------------------------------------------------- */

	n = MAX (1, n) ;
	anz = Ap [n+k1] - Ap [k1] ;

	if (Lsize <= 0)
	{
		Lsize = -Lsize ;
		Lsize = MAX (Lsize, 1.0) ;
		lsize = (Lsize * anz + n).toInt() ;
	}
	else
	{
		lsize = Lsize.toInt() ;
	}

	usize = lsize ;

	lsize  = MAX (n+1, lsize) ;
	usize  = MAX (n+1, usize) ;

	maxlnz = ((n) * (n) + (n)) / 2.0 ;
	maxlnz = MIN (maxlnz, (INT_MAX)) ;
	lsize  = MIN (maxlnz.toInt(), lsize) ;
	usize  = MIN (maxlnz.toInt(), usize) ;

	PRINTF ("Welcome to klu: n $n anz $anz k1 $k1 lsize $lsize usize $usize maxlnz $maxlnz\n") ;

	/* ---------------------------------------------------------------------- */
	/* allocate workspace and outputs */
	/* ---------------------------------------------------------------------- */

	/* return arguments are not yet assigned */
	p_LU [block] = null ;

	/* these computations are safe from int overflow */
	//W = Work ;
	//Pinv = W ;      //W += n ;
	//int Pinv_offset = n ;
	Pinv = new List<int>(n) ;
	//Stack = W ;     //W += n ;
	//int Stack_offset = 2*n ;
	Stack = new List<int>(n) ;
	//Flag = W ;      //W += n ;
	//int Flag_offset = 3*n ;
	Flag = new List<int>(n) ;
	//Lpend = W ;     //W += n ;
	//int Lpend_offset = 4*n ;
	Lpend = new List<int>(n) ;
	//Ap_pos = W ;    //W += n ;
	//int Ap_pos_offset = 5*n ;
	Ap_pos = new List<int>(n) ;

	//dunits = DUNITS (Integer, lsize) + DUNITS (Double, lsize) +
	//		 DUNITS (Integer, usize) + DUNITS (Double, usize) ;
	dunits = (lsize + lsize + usize + usize).toDouble() ;
	lusize = dunits.toInt() ;
	ok = INT_OVERFLOW (dunits) ? FALSE : TRUE ;
	LU [0] = ok != 0 ? malloc_dbl (lusize, Common) : null ;
	if (LU [0] == null)
	{
		/* out of memory, or problem too large */
		Common.status = KLU_OUT_OF_MEMORY ;
		lusize = 0 ;
		return (lusize) ;
	}

	/* ---------------------------------------------------------------------- */
	/* factorize */
	/* ---------------------------------------------------------------------- */

	/* with pruning, and non-recursive depth-first-search */
	lusize = kernel (n, Ap, Ai, Ax, Q, lusize,
			Pinv, P, LU, Udiag, Udiag_offset, Llen, Llen_offset,
			Ulen, Ulen_offset, Lip, Lip_offset, Uip, Uip_offset,
			lnz, unz, X, Stack, Flag, Ap_pos, Lpend,
			k1, PSinv, Rs, Offp, Offi, Offx, Common) ;

	/* ---------------------------------------------------------------------- */
	/* return LU factors, or return nothing if an error occurred */
	/* ---------------------------------------------------------------------- */

	if (Common.status < KLU_OK)
	{
		//LU = KLU_free (LU, lusize, sizeof (double), Common) ;
		LU [0] = null ;
		lusize = 0 ;
	}
	p_LU [block] = LU [0] ;
	PRINTF (" in klu noffdiag ${Common.noffdiag}\n") ;
	return (lusize) ;
}

/**
 * Solve Lx=b.  Assumes L is unit lower triangular and where the unit diagonal
 * entry is NOT stored.  Overwrites B  with the solution X.  B is n-by-nrhs
 * and is stored in ROW form with row dimension nrhs.  nrhs must be in the
 * range 1 to 4.
 *
 * @param n
 * @param Lip
 * @param Llen
 * @param LU
 * @param nrhs
 * @param X right-hand-side on input, solution to Lx=b on output
 */
void lsolve(int n, List<int> Lip, int Lip_offset,
		List<int> Llen, int Llen_offset, List<double> LU, int nrhs,
		List<double> X, int X_offset)
{
	List<double> x = new List<double>(4) ;
	double lik ;
	/*List<int>*/List<double> Li ;
	List<double> Lx ;
	int k, p, i ;
	List<int> len = new List<int>(1) ;
	List<int> Li_offset = new List<int>(1) ;
	List<int> Lx_offset = new List<int>(1) ;

	switch (nrhs)
	{

		case 1:
			for (k = 0 ; k < n ; k++)
			{
				x [0] = X [X_offset + k] ;
				Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len) ;
				/* unit diagonal of L is not stored*/
				for (p = 0 ; p < len[0] ; p++)
				{
					//MULT_SUB (X [Li [p]], Lx [p], x [0]) ;
					X [X_offset + (Li [Li_offset[0] + p]).toInt()] -= Lx [Lx_offset[0] + p] * x [0] ;
				}
			}
			break ;

		case 2:

			for (k = 0 ; k < n ; k++)
			{
				x [0] = X [X_offset + 2*k    ] ;
				x [1] = X [X_offset + 2*k + 1] ;
				Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len) ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = (Li [Li_offset[0] + p]).toInt() ;
					lik = Lx [Lx_offset[0] + p] ;
					//MULT_SUB (X [2*i], lik, x [0]) ;
					X [X_offset + 2*i] -= lik * x [0] ;
					//MULT_SUB (X [2*i + 1], lik, x [1]) ;
					X [X_offset + 2*i + 1] -= lik * x [1] ;
				}
			}
			break ;

		case 3:

			for (k = 0 ; k < n ; k++)
			{
				x [0] = X [X_offset + 3*k    ] ;
				x [1] = X [X_offset + 3*k + 1] ;
				x [2] = X [X_offset + 3*k + 2] ;
				Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len) ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = (Li [Li_offset[0] + p]).toInt() ;
					lik = Lx [Lx_offset[0] + p] ;
					//MULT_SUB (X [3*i], lik, x [0]) ;
					X [X_offset + 3*i] -= lik * x [0] ;
					//MULT_SUB (X [3*i + 1], lik, x [1]) ;
					X [X_offset + 3*i + 1] -= lik * x [1] ;
					//MULT_SUB (X [3*i + 2], lik, x [2]) ;
					X [X_offset + 3*i + 2] -= lik * x [2] ;
				}
			}
			break ;

		case 4:

			for (k = 0 ; k < n ; k++)
			{
				x [0] = X [X_offset + 4*k    ] ;
				x [1] = X [X_offset + 4*k + 1] ;
				x [2] = X [X_offset + 4*k + 2] ;
				x [3] = X [X_offset + 4*k + 3] ;
				Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len) ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = (Li [Li_offset[0] + p]).toInt() ;
					lik = Lx [Lx_offset[0] + p] ;
					//MULT_SUB (X [4*i], lik, x [0]) ;
					X [X_offset + 4*i] -= lik * x [0] ;
					//MULT_SUB (X [4*i + 1], lik, x [1]) ;
					X [X_offset + 4*i + 1] -= lik * x [1] ;
					//MULT_SUB (X [4*i + 2], lik, x [2]) ;
					X [X_offset + 4*i + 2] -= lik * x [2] ;
					//MULT_SUB (X [4*i + 3], lik, x [3]) ;
					X [X_offset + 4*i + 3] -= lik * x [3] ;
				}
			}
			break ;

	}
}

/**
 * Solve Ux=b.  Assumes U is non-unit upper triangular and where the diagonal
 * entry is NOT stored.  Overwrites B with the solution X.  B is n-by-nrhs
 * and is stored in ROW form with row dimension nrhs.  nrhs must be in the
 * range 1 to 4.
 *
 * @param n
 * @param Uip
 * @param Ulen
 * @param LU
 * @param Udiag
 * @param nrhs
 * @param X right-hand-side on input, solution to Ux=b on output
 */
void usolve(int n, List<int> Uip, int Uip_offset,
		List<int> Ulen, int Ulen_offset, List<double> LU,
		List<double> Udiag, int Udiag_offset, int nrhs,
		List<double> X, int X_offset)
{
	List<double> x = new List<double>(4) ;
	double uik, ukk ;
	/*List<int>*/List<double> Ui ;
	List<double> Ux ;
	int k, p, i ;
	List<int> len = new List<int>(1) ;
	List<int> Ui_offset = new List<int>(1) ;
	List<int> Ux_offset = new List<int>(1) ;

	switch (nrhs)
	{

		case 1:

			for (k = n-1 ; k >= 0 ; k--)
			{
				Ui = Ux = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset,
						Ui_offset, Ux_offset, k, len) ;
				//DIV (x [0], X [k], Udiag [k]) ;
				x [0] = X [X_offset + k] / Udiag [Udiag_offset + k] ;
				X [X_offset + k] = x [0] ;
				for (p = 0 ; p < len[0] ; p++)
				{
					//MULT_SUB (X [Ui [p]], Ux [p], x [0]) ;
					X [X_offset + (Ui [Ui_offset[0] + p]).toInt()] -= Ux [Ux_offset[0] + p] * x [0] ;

				}
			}

			break ;

		case 2:

			for (k = n-1 ; k >= 0 ; k--)
			{
				Ui = Ux = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset,
						Ui_offset, Ux_offset, k, len) ;
				ukk = Udiag [Udiag_offset + k] ;
				//DIV (x [0], X [2*k], ukk) ;
				x [0] = X [X_offset + 2*k    ] / ukk ;
				//DIV (x [1], X [2*k + 1], ukk) ;
				x [1] = X [X_offset + 2*k + 1] / ukk ;

				X [X_offset + 2*k    ] = x [0] ;
				X [X_offset + 2*k + 1] = x [1] ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = Ui [Ui_offset[0] + p].toInt() ;
					uik = Ux [Ux_offset[0] + p] ;
					//MULT_SUB (X [2*i], uik, x [0]) ;
					X [X_offset + 2*i    ] -= uik * x [0] ;
					//MULT_SUB (X [2*i + 1], uik, x [1]) ;
					X [X_offset + 2*i + 1] -= uik * x [1] ;
				}
			}

			break ;

		case 3:

			for (k = n-1 ; k >= 0 ; k--)
			{
				Ui = Ux = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset,
						Ui_offset, Ux_offset, k, len) ;
				ukk = Udiag [Udiag_offset + k] ;

				//DIV (x [0], X [3*k], ukk) ;
				x [0] = X [X_offset + 3*k] / ukk ;
				//DIV (x [1], X [3*k + 1], ukk) ;
				x [1] = X [X_offset + 3*k + 1] / ukk ;
				//DIV (x [2], X [3*k + 2], ukk) ;
				x [2] = X [X_offset + 3*k + 2] / ukk ;

				X [3*k    ] = x [0] ;
				X [3*k + 1] = x [1] ;
				X [3*k + 2] = x [2] ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = Ui [Ui_offset[0] + p].toInt() ;
					uik = Ux [Ux_offset[0] + p] ;
					//MULT_SUB (X [3*i], uik, x [0]) ;
					X [X_offset + 3*i] -= uik * x [0] ;
					//MULT_SUB (X [3*i + 1], uik, x [1]) ;
					X [X_offset + 3*i + 1] -= uik * x [1] ;
					//MULT_SUB (X [3*i + 2], uik, x [2]) ;
					X [X_offset + 3*i + 2] -= uik * x [2] ;
				}
			}

			break ;

		case 4:

			for (k = n-1 ; k >= 0 ; k--)
			{
				Ui = Ux = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset,
						Ui_offset, Ux_offset, k, len) ;
				ukk = Udiag [Udiag_offset + k] ;

				//DIV (x [0], X [4*k], ukk) ;
				x [0] = X [X_offset + 4*k] / ukk ;
				//DIV (x [1], X [4*k + 1], ukk) ;
				x [1] = X [X_offset + 4*k + 1] / ukk ;
				//DIV (x [2], X [4*k + 2], ukk) ;
				x [2] = X [X_offset + 4*k + 2] / ukk ;
				//DIV (x [3], X [4*k + 3], ukk) ;
				x [3] = X [X_offset + 4*k + 3] / ukk ;

				X [X_offset + 4*k    ] = x [0] ;
				X [X_offset + 4*k + 1] = x [1] ;
				X [X_offset + 4*k + 2] = x [2] ;
				X [X_offset + 4*k + 3] = x [3] ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = Ui [Ui_offset[0] + p].toInt() ;
					uik = Ux [Ux_offset[0] + p] ;

					//MULT_SUB (X [4*i], uik, x [0]) ;
					X [X_offset + 4*i] -= uik * x [0] ;
					//MULT_SUB (X [4*i + 1], uik, x [1]) ;
					X [X_offset + 4*i + 1] -= uik * x [1] ;
					//MULT_SUB (X [4*i + 2], uik, x [2]) ;
					X [X_offset + 4*i + 2] -= uik * x [2] ;
					//MULT_SUB (X [4*i + 3], uik, x [3]) ;
					X [X_offset + 4*i + 3] -= uik * x [3] ;
				}
			}

			break ;

	}
}

/**
 * Solve L'x=b.  Assumes L is unit lower triangular and where the unit diagonal
 * entry is NOT stored.  Overwrites B with the solution X.  B is n-by-nrhs
 * and is stored in ROW form with row dimension nrhs.  nrhs must in the
 * range 1 to 4.
 *
 * @param n
 * @param Lip
 * @param Llen
 * @param LU
 * @param nrhs
 * @param X right-hand-side on input, solution to L'x=b on output
 */
void ltsolve(int n, List<int> Lip, int Lip_offset, List<int> Llen, int Llen_offset,
		List<double> LU, int nrhs, List<double> X, int X_offset)
{
	List<double> x = new List<double>(4) ;
	double lik ;
	/*List<int>*/List<double> Li ;
	List<double> Lx ;
	int k, p, i ;
	List<int> len = new List<int>(1) ;
	List<int> Li_offset = new List<int>(1) ;
	List<int> Lx_offset = new List<int>(1) ;

	switch (nrhs)
	{

		case 1:

			for (k = n-1 ; k >= 0 ; k--)
			{
				Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len) ;
				x [0] = X [X_offset + k] ;
				for (p = 0 ; p < len[0] ; p++)
				{
					{
						//MULT_SUB (x [0], Lx [p], X [Li [p]]) ;
						x [0] -= Lx [Lx_offset[0] + p] * X [X_offset + Li [Li_offset[0] + p].toInt()] ;
					}
				}
				X [X_offset + k] = x [0] ;
			}
			break ;

		case 2:

			for (k = n-1 ; k >= 0 ; k--)
			{
				x [0] = X [X_offset + 2*k    ] ;
				x [1] = X [X_offset + 2*k + 1] ;
				Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len) ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = Li [Li_offset[0] + p].toInt() ;
					{
						lik = Lx [Lx_offset[0] + p] ;
					}
					//MULT_SUB (x [0], lik, X [2*i]) ;
					x [0] -= lik * X [X_offset + 2*i] ;
					//MULT_SUB (x [1], lik, X [2*i + 1]) ;
					x [1] -= lik * X [X_offset + 2*i + 1] ;
				}
				X [X_offset + 2*k    ] = x [0] ;
				X [X_offset + 2*k + 1] = x [1] ;
			}
			break ;

		case 3:

			for (k = n-1 ; k >= 0 ; k--)
			{
				x [0] = X [X_offset + 3*k    ] ;
				x [1] = X [X_offset + 3*k + 1] ;
				x [2] = X [X_offset + 3*k + 2] ;
				Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len) ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = Li [Li_offset[0] + p].toInt() ;
					{
						lik = Lx [Lx_offset[0] + p] ;
					}
					//MULT_SUB (x [0], lik, X [3*i]) ;
					x [0] -= lik * X [X_offset + 3*i] ;
					//MULT_SUB (x [1], lik, X [3*i + 1]) ;
					x [1] -= lik * X [X_offset + 3*i + 1] ;
					//MULT_SUB (x [2], lik, X [3*i + 2]) ;
					x [2] -= lik * X [X_offset + 3*i + 2] ;
				}
				X [X_offset + 3*k    ] = x [0] ;
				X [X_offset + 3*k + 1] = x [1] ;
				X [X_offset + 3*k + 2] = x [2] ;
			}
			break ;

		case 4:

			for (k = n-1 ; k >= 0 ; k--)
			{
				x [0] = X [X_offset + 4*k    ] ;
				x [1] = X [X_offset + 4*k + 1] ;
				x [2] = X [X_offset + 4*k + 2] ;
				x [3] = X [X_offset + 4*k + 3] ;
				Li = Lx = GET_POINTER (LU, Lip, Lip_offset, Llen, Llen_offset, Li_offset, Lx_offset, k, len) ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = Li [Li_offset[0] + p].toInt() ;
					{
						lik = Lx [Lx_offset[0] + p] ;
					}
					//MULT_SUB (x [0], lik, X [4*i]) ;
					x [0] -= lik * X [X_offset + 4*i] ;
					//MULT_SUB (x [1], lik, X [4*i + 1]) ;
					x [1] -= lik * X [X_offset + 4*i + 1] ;
					//MULT_SUB (x [2], lik, X [4*i + 2]) ;
					x [2] -= lik * X [X_offset + 4*i + 2] ;
					//MULT_SUB (x [3], lik, X [4*i + 3]) ;
					x [3] -= lik * X [X_offset + 4*i + 3] ;
				}
				X [X_offset + 4*k    ] = x [0] ;
				X [X_offset + 4*k + 1] = x [1] ;
				X [X_offset + 4*k + 2] = x [2] ;
				X [X_offset + 4*k + 3] = x [3] ;
			}
			break ;
	}
}

/**
 * Solve U'x=b.  Assumes U is non-unit upper triangular and where the diagonal
 * entry is stored (and appears last in each column of U).  Overwrites B
 * with the solution X.  B is n-by-nrhs and is stored in ROW form with row
 * dimension nrhs.  nrhs must be in the range 1 to 4.
 *
 * @param n
 * @param Uip
 * @param Ulen
 * @param LU
 * @param Udiag
 * @param nrhs
 * @param X right-hand-side on input, solution to Ux=b on output
 */
void utsolve(int n, List<int> Uip, int Uip_offset,
		List<int> Ulen, int Ulen_offset, List<double> LU,
		List<double> Udiag, int Udiag_offset, int nrhs,
		List<double> X, int X_offset)
{
	List<double> x = new List<double>(4) ;
	double uik, ukk ;
	int k, p, i ;
	/*List<int>*/List<double> Ui ;
	List<double> Ux ;
	List<int> len = new List<int>(1) ;
	List<int> Ui_offset = new List<int>(1) ;
	List<int> Ux_offset = new List<int>(1) ;

	switch (nrhs)
	{

		case 1:

			for (k = 0 ; k < n ; k++)
			{
				Ui = Ux = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len) ;
				x [0] = X [X_offset + k] ;
				for (p = 0 ; p < len[0] ; p++)
				{
					{
						//MULT_SUB (x [0], Ux [p], X [Ui [p]]) ;
						x [0] -= Ux [Ux_offset[0] + p] * X [X_offset + Ui [Ui_offset[0] + p].toInt()] ;
					}
				}
				{
					ukk = Udiag [k] ;
				}
				//DIV (X [k], x [0], ukk) ;
				X [X_offset + k] = x [0] / ukk ;
			}
			break ;

		case 2:

			for (k = 0 ; k < n ; k++)
			{
				Ui = Ux = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len) ;
				x [0] = X [X_offset + 2*k    ] ;
				x [1] = X [X_offset + 2*k + 1] ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = Ui [Ui_offset[0] + p].toInt() ;
					{
						uik = Ux [Ux_offset[0] + p] ;
					}
					//MULT_SUB (x [0], uik, X [2*i]) ;
					x [0] -= uik * X [X_offset + 2*i] ;
					//MULT_SUB (x [1], uik, X [2*i + 1]) ;
					x [1] -= uik * X [X_offset + 2*i + 1] ;
				}
				{
					ukk = Udiag [k] ;
				}
				//DIV (X [2*k], x [0], ukk) ;
				X [X_offset + 2*k] = x [0] / ukk ;
				//DIV (X [2*k + 1], x [1], ukk) ;
				X [X_offset + 2*k + 1] = x [1] / ukk ;
			}
			break ;

		case 3:

			for (k = 0 ; k < n ; k++)
			{
				Ui = Ux = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len) ;
				x [0] = X [X_offset + 3*k    ] ;
				x [1] = X [X_offset + 3*k + 1] ;
				x [2] = X [X_offset + 3*k + 2] ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = Ui [Ui_offset[0] + p].toInt() ;
					{
						uik = Ux [Ux_offset[0] + p] ;
					}
					//MULT_SUB (x [0], uik, X [3*i]) ;
					x [0] -= uik * X [X_offset + 3*i] ;
					//MULT_SUB (x [1], uik, X [3*i + 1]) ;
					x [1] -= uik * X [X_offset + 3*i + 1] ;
					//MULT_SUB (x [2], uik, X [3*i + 2]) ;
					x [2] -= uik * X [X_offset + 3*i + 2] ;
				}
				{
					ukk = Udiag [k] ;
				}
				//DIV (X [3*k], x [0], ukk) ;
				X [X_offset + 3*k] = x [0] / ukk ;
				//DIV (X [3*k + 1], x [1], ukk) ;
				X [X_offset + 3*k + 1] = x [1] / ukk ;
				//DIV (X [3*k + 2], x [2], ukk) ;
				X [X_offset + 3*k + 2] = x [2] / ukk ;
			}
			break ;

		case 4:

			for (k = 0 ; k < n ; k++)
			{
				Ui = Ux = GET_POINTER (LU, Uip, Uip_offset, Ulen, Ulen_offset, Ui_offset, Ux_offset, k, len) ;
				x [0] = X [X_offset + 4*k    ] ;
				x [1] = X [X_offset + 4*k + 1] ;
				x [2] = X [X_offset + 4*k + 2] ;
				x [3] = X [X_offset + 4*k + 3] ;
				for (p = 0 ; p < len[0] ; p++)
				{
					i = Ui [Ui_offset[0] + p].toInt() ;
					{
						uik = Ux [Ux_offset[0] + p] ;
					}
					//MULT_SUB (x [0], uik, X [4*i]) ;
					x [0] -= uik * X [X_offset + 4*i] ;
					//MULT_SUB (x [1], uik, X [4*i + 1]) ;
					x [1] -= uik * X [X_offset + 4*i + 1] ;
					//MULT_SUB (x [2], uik, X [4*i + 2]) ;
					x [2] -= uik * X [X_offset + 4*i + 2] ;
					//MULT_SUB (x [3], uik, X [4*i + 3]) ;
					x [3] -= uik * X [X_offset + 4*i + 3] ;
				}
				{
					ukk = Udiag [k] ;
				}
				//DIV (X [4*k], x [0], ukk) ;
				X [X_offset + 4*k] = x [0] / ukk ;
				//DIV (X [4*k + 1], x [1], ukk) ;
				X [X_offset + 4*k + 1] = x [1] / ukk ;
				//DIV (X [4*k + 2], x [2], ukk) ;
				X [X_offset + 4*k + 2] = x [2] / ukk ;
				//DIV (X [4*k + 3], x [3], ukk) ;
				X [X_offset + 4*k + 3] = x [3] / ukk ;
			}
			break ;
	}
}
