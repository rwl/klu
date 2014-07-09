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

//public abstract class Dklu_version {

/*public static */final int KLU_OK = 0;
/** status > 0 is a warning, not an error */
/*public static */final int KLU_SINGULAR = 1;
/*public static */final int KLU_OUT_OF_MEMORY = -2;
/*public static */final int KLU_INVALID = -3;
/** integer overflow has occured */
/*public static */final int KLU_TOO_LARGE = -4;

/** enable diagnostic printing */
/*public static */boolean NPRINT = true ;

final int INT_MAX = 0x7fffffff ;

final String INT_ID = "%d" ;

//	int BYTES (Object type, double n)
//	{
//		return sizeof (type * n) ;
//	}
//
//	double CEILING (double b, double u)
//	{
//		return (b+u-1) / u ;
//	}
//
//	double UNITS (Object type, double n)
//	{
//		return CEILING (BYTES (type, n), sizeof (double)) ;
//	}
//
//	double DUNITS (Object type, int n)
//	{
//		return Math.ceil(BYTES (type, (double) n) / sizeof (double)) ;
//	}

List<double> GET_I_POINTER(List<double> LU, List<int> Xip,
    int Xip_offset, List<int> Xi_offset, int k)
{
  Xi_offset[0] = Xip [Xip_offset + k] ;
  return LU ;
}

//	void GET_X_POINTER(List<double> LU, List<int> Xip, int Xlen,
//			List<double> Xx, int k)
//	{
//		Xx = (List<double>) (LU + Xip [k] + UNITS (Int, Xlen [k])) ;
//	}

List<double> GET_POINTER(List<double> LU,
    List<int> Xip, int Xip_offset,
    List<int> Xlen, int Xlen_offset,
    List<int> Xi_offset,
    List<int> Xx_offset,
    int k, List<int> xlen)
{
  int xp = Xip [Xip_offset + k] ;
  xlen[0] = Xlen [Xlen_offset + k] ;
  Xi_offset[0] = xp ;
  Xx_offset[0] = xp + xlen[0] ;
  return LU ;
}

boolean SCALAR_IS_NAN (double x)
{
  return x != x ;
}

boolean SCALAR_IS_ZERO (double x)
{
  return x == 0.0 ;
}

boolean SCALAR_IS_NONZERO (double x)
{
  return x != 0.0 ;
}

boolean SCALAR_IS_LTZERO (double x)
{
  return x < 0.0 ;
}

/* scalar absolute value macro. If x is NaN, the result is NaN */
double SCALAR_ABS (double x)
{
  return SCALAR_IS_LTZERO (x) ? -x : x ;
}

void PRINTF(String format, Object args)
{
  if (!NPRINT) {
    System.out.printf(format, args) ;
  }
}

void PRINT_SCALAR (double a)
{
  if (!NPRINT) {

    if (SCALAR_IS_NONZERO (a))
    {
      PRINTF (" (%g)", a) ;
    }
    else
    {
      PRINTF (" (0)") ;
    }

  }
}

/* ---------------------------------------------------------------------- */
/* Real floating-point arithmetic */
/* ---------------------------------------------------------------------- */

/**
 * @return TRUE if a complex number is in split form, FALSE if in packed
 * form.
 */
int SPLIT (double s)
{
  return 1 ;
}

/**
 * @return real part of c
 */
//	double REAL (double c)
//	{
//		return c ;
//	}

/**
 * @return imag part of c
 */
//	double IMAG (double c)
//	{
//		return 0.0 ;
//	}

/**
 * c = (s1) + (s2)*i
 */
//	void ASSIGN (Double c, List<double> s1, List<double> s2, int p,
//			boolean split)
//	{
//		c = s1[p] ;
//	}

//	void CLEAR (Double c)
//	{
//		c = 0.0 ;
//	}

//	void CLEAR_AND_INCREMENT (Double p)
//	{
//		p = 0.0 ;
//		p++ ;
//	}

/**
 * @return True if a is NaN
 */
boolean IS_NAN (double a)
{
  return SCALAR_IS_NAN (a) ;
}

/**
 * @return True if a == 0
 */
boolean IS_ZERO (double a)
{
  return SCALAR_IS_ZERO (a) ;
}

/**
 * @return True if a != 0
 */
boolean IS_NONZERO (double a)
{
  return SCALAR_IS_NONZERO (a) ;
}

/**
 * c /= s
 */
double SCALE_DIV (double c, double s)
{
  c /= s ;
  return c ;
}

/**
 * a = c/s
 */
//	double SCALE_DIV_ASSIGN (double c, double s)
//	{
//		return c / s ;
//	}

/**
 * c *= s
 */
//	void SCALE (Double c, double s)
//	{
//		c *= s ;
//	}

/**
 * c += a
 */
//	void ASSEMBLE (Double c, double a)
//	{
//		c += a ;
//	}

/**
 * c += *p++
 */
//	void ASSEMBLE_AND_INCREMENT (Double c, double p)
//	{
//		c += p++ ;
//	}

/**
 * c -= a
 */
//	void DECREMENT (Double c, double a)
//	{
//		c -= a ;
//	}

/**
 * c = a*b
 */
//	void MULT (Double c, double a,  double b)
//	{
//		c = a * b ;
//	}

/**
 * c = a*conjugate(b)
 */
//	void MULT_CONJ (Double c, double a, double b)
//	{
//		c = a * b ;
//	}

/**
 * c -= a*b
 */
double MULT_SUB (double c, double a, double b)
{
  c -= a * b ;
  return c ;
}

//	/**
//	 * c -= a*conjugate(b)
//	 */
//	void MULT_SUB_CONJ (Double c, double a, double b)
//	{
//		c -= a * b ;
//	}
//
//	/**
//	 * c = a/b
//	 */
//	void DIV (Double c, double a, double b)
//	{
//		c = a / b ;
//	}
//
//	/**
//	 * c = 1/c
//	 */
//	void RECIPROCAL (Double c)
//	{
//		c = 1.0 / c ;
//	}
//
//	/**
//	 * c = a/conjugate(b)
//	 */
//	void DIV_CONJ (Double c, double a, double b)
//	{
//		c = a / b ;
//	}
//
//	/**
//	 * approximate absolute value, s = |r|+|i|
//	 */
//	void APPROX_ABS (Double s, double a)
//	{
//		s = SCALAR_ABS (a) ;
//	}
//
//	/**
//	 * exact absolute value, s = sqrt (a.real^2 + amag^2)
//	 */
//	void ABS (Double s, double a)
//	{
//		s = SCALAR_ABS (a) ;
//	}

void PRINT_ENTRY (double a)
{
  PRINT_SCALAR (a) ;
}

//	void CONJ (Double a, double x)
//	{
//		a = x ;
//	}

double ABS (double a)
{
  return Math.abs (a) ;
}

void CLEAR(List<double> A, int i)
{
  A [i] = 0.0 ;
}

/* for flop counts */
final double MULTSUB_FLOPS   = 2.0 ;      /* c -= a*b */
final double DIV_FLOPS       = 1.0 ;      /* c = a/b */
final double ABS_FLOPS       = 0.0 ;      /* c = abs(a) */
final double ASSEMBLE_FLOPS  = 1.0 ;      /* c += a */
final double DECREMENT_FLOPS = 1.0 ;      /* c -= a */
final double MULT_FLOPS      = 1.0 ;      /* c = a*b */
final double SCALE_FLOPS     = 1.0 ;      /* c = a/s */
